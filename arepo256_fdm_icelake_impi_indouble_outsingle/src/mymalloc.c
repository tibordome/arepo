/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/mymalloc.c
 * \date        MM/YYYY
 * \author
 * \brief       Manager for dynamic memory allocation
 * \details     This module handles the dynamic memory allocation for AREPO.
 *              To avoid memory allocation/deallocation overhead, a big chunk of memory
 *              (which will be the maximum amount of dynamically allocatable memory)
 *              is allocated upon initialization. This chunk is then filled by the memory
 *              blocks as in a stack structure. The blocks are automatically aligned to a
 *              64-bit boundary.
 *              Memory blocks come in two flavours: movable and non-movable. In non-movable
 *              blocks, the starting address is fixed once the block is allocated and cannot be
 *              changed. Due to the stack structure of the dynamic memory, this implies that the
 *              last (non-movable) block allocated must be the first block to be deallocated.
 *              If this condition is not met, an abort condition is triggered.
 *              If more flexibility is needed, movable memory blocks can be used. In this case,
 *              the starting address of the block is again fixed upon allocation but the block
 *              can be shifted (therefore its initial address changes) according to needs.
 *              For a movable block to be successfully shifted, it is required that all the
 *              subsequent allocated blocks are movable. Again, an abort condition is triggered
 *              if this condition is not met.
 *              Movable blocks can be deallocated in any order provided that the condition just
 *              described holds. The gap resulting from the deallocation of a block that is not in
 *              the last position will be automatically filled by shifting all the blocks coming
 *              after the deallocated block.
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

/* configuration options must be defined in order for the #ifdef check to make
 * sense */
#include "arepoconfig.h"
#if defined(MEMORY_MANAGER_USE_MPROTECT) && !defined(_GNU_SOURCE)
/* Need to define _GNU_SOURCE so that MAP_ANONYMOUS will be defined in the
 * headers. This has to be done before loading any system headers. */
#define _GNU_SOURCE
#endif

#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"

#ifdef HUGEPAGES
#include <hugetlbfs.h>
#endif

#ifdef MEMORY_MANAGER_USE_MPROTECT
#include <sys/mman.h>
#endif

#define CACHELINESIZE 64

#ifdef CHIMES_PTHREADS
#define MAXBLOCKS 250000
#else /* CHIMES_PTHREADS */
#ifdef GFM
#ifdef GFM_AGN_RADIATION
#define MAXBLOCKS 25000
#else
#ifdef RADCOOL_HOTHALO
#define MAXBLOCKS 500000
#else /* RADCOOL_HOTHALO */
#define MAXBLOCKS 50000
#endif /* RADCOOL_HOTHALO */
#endif /* GFM_AGN_RADIATION */
#else  /* GFM */
#define MAXBLOCKS 5000
#endif /* GFM */
#endif /* CHIMES_PTHREADS */

#define MAXCHARS 40
#define TAB_BUF_SIZE ((100 + 4 * MAXCHARS) * (MAXBLOCKS + 10))

static size_t AllocatedBytesGeneric;

static size_t HighMarkBytes;
static size_t HighMarkBytesWithoutGeneric;

static double OldGlobHighMarkMB;
static double OldGlobHighMarkMBWithoutGeneric;

static size_t TotBytes; /**< The total dimension (in bytes) of dynamic memory available to the current task. */
static void *Base;      /**< Base pointer to the start of the stack. */

static unsigned long Nblocks; /**< The current number of allocated memory blocks. */

static void **Table;         /**< Array containing the start addresses of all allocated memory blocks. */
static size_t *BlockSize;    /**< Array containing the size (in bytes) of all the allocated memory blocks. */
static bool *GenericFlag;    /**< Identifies whether a block has been allocated in the generic allocation templates
                                (generic_comm_helpers2.h, generic_comm_helpers_async.h). */
static void ***BasePointers; /**< For each movable block, keeps track of a pointer to the variable that holds the memory address of the
                                block (as returned by this function). This allows the variable to be updated with a new address when
                                memory blocks are moved. */
static char (*VarName)[MAXCHARS];        /**< The name of the variable with which the block has been allocated. */
static char (*FunctionName)[MAXCHARS];   /**< The function name that has allocated the memory block. */
static char (*ParentFileName)[MAXCHARS]; /**< The location from which the generic communication templates (generic_comm_helpers2.h,
                                            generic_comm_helpers_async.h) were called. */
static char (*FileName)[MAXCHARS];       /**< The file name where the function that has allocated the block is called. */
static int *LineNumber;      /**< The line number in FileName where the function that allocated the block has been called. */
static char *HighMarkTabBuf; /**< This is a buffer that holds the log file output corresponding to the largest memory use that has
                                occurred on this task. */
static char *HighMarkTabBufWithoutGeneric; /**< This is a buffer that holds the log file output corresponding to the largest memory use
                                              that has occurred on this task, without including allocations from the generic
                                              communication templates (generic_comm_helpers2.h, generic_comm_helpers_async.h). */

size_t FreeBytes(void) { return TotBytes - AllocatedBytes; }

static size_t roundup_to_multiple_of(size_t n, const size_t factor)
{
  if(n % factor > 0)
    n = (n / factor + 1) * factor;
  return n;
}

/*! \brief Allocation function wrapper to allocate memory pages (if necessary).
 *
 *  \param[in] size Size of the allocated memory.
 *
 *  \return void pointer to address in memory.
 */
static void *pmalloc(const size_t size)
{
#ifdef MEMORY_MANAGER_USE_MPROTECT
  /* mprotect() should be called on memory obtained from mmap() */
  return mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
#else
  return malloc(size);
#endif
}

/*! \brief Allocation function wrapper for hugepages usage.
 *
 *  \param[in] size Size of the allocated memory.
 *
 *  \return void pointer to address in memory.
 */
static void *hmalloc(const size_t size)
{
#ifdef HUGEPAGES
  void *p = get_hugepage_region(size, GHR_STRICT);
  if(!p)
    {
      warn("Failed to get_hugepage_region of size %g", size / pow(1024, 2));
      p = pmalloc(size);
      if(!p)
        terminate("Failed to allocate memory of size %g", size / pow(1024, 2));
    }
  /* mark pages as used */
  memset(p, 255, size);
  memset(p, 0, size);
  return p;
#else
  return pmalloc(size);
#endif
}

/*! \brief Updates memory protection using mprotect() to prevent access to
 *         unused memory blocks when a block is allocated or deallocated.
 *
 *  \param[in] allocated_bytes Number of bytes in currently allocated memory
 *                             blocks.
 *                             If offset > 0: allocated_bytes is before the
 *                             allocation.
 *                             If offset < 0: allocated_bytes is after the
 *                             deallocation.
 *  \param[in] offset Change in the number of allocated bytes due to the
 *                    allocation operation.
 *
 *  \return Return value of mprotect(), or 0 if mprotect() was not called.
 */
static int mprotect_on_allocation(const size_t allocated_bytes, const ptrdiff_t offset)
{
  if(offset == 0)
    return 0;
  int error = 0;
#ifdef MEMORY_MANAGER_USE_MPROTECT
  const size_t first_unused_page = roundup_to_multiple_of(allocated_bytes, sysconf(_SC_PAGESIZE));
  const ptrdiff_t extra_space    = first_unused_page - allocated_bytes;
  ptrdiff_t prot_len             = (ptrdiff_t)imaxabs(offset) - extra_space;
  if(first_unused_page + prot_len >= TotBytes)
    prot_len = TotBytes - first_unused_page;
  if(prot_len > 0)
    {
      const int prot_flag = offset > 0 ? PROT_READ | PROT_WRITE : PROT_NONE;
      error               = mprotect((char *)Base + first_unused_page, prot_len, prot_flag);
      if(error)
        warn("mprotect(%p + %zu, %td, %d) failed (%d): %s", Base, first_unused_page, prot_len, prot_flag, errno, strerror(errno));
    }
#endif
  return error;
}

/*! \brief Rounds up size to cache line size.
 *
 *  \param[in] n Size.
 *
 *  \return Rounded-up size.
 */
size_t roundup_to_multiple_of_cacheline_size(size_t n)
{
#ifdef _SC_LEVEL1_DCACHE_LINESIZE
  static bool check_cache_line_size = true;
  if(check_cache_line_size)
    {
      myassert(sysconf(_SC_LEVEL1_DCACHE_LINESIZE) == CACHELINESIZE);
      check_cache_line_size = false;
    }
#endif
  n = roundup_to_multiple_of(n, CACHELINESIZE);
  if(n < CACHELINESIZE)
    n = CACHELINESIZE;
  return n;
}

/*! \brief Initializes the memory manager.
 *
 *  This function initializes the memory manager. In particular, it sets
 *  the global variables of the module to their initial values and allocates
 *  the memory for the stack.
 *
 *  \return void
 */
void mymalloc_init(void)
{
  Table                        = (void **)hmalloc(MAXBLOCKS * sizeof(*Table));
  BlockSize                    = (size_t *)hmalloc(MAXBLOCKS * sizeof(*BlockSize));
  GenericFlag                  = (bool *)hmalloc(MAXBLOCKS * sizeof(*GenericFlag));
  BasePointers                 = (void ***)hmalloc(MAXBLOCKS * sizeof(*BasePointers));
  VarName                      = (char(*)[MAXCHARS])hmalloc(MAXBLOCKS * sizeof(*VarName));
  FunctionName                 = (char(*)[MAXCHARS])hmalloc(MAXBLOCKS * sizeof(*FunctionName));
  ParentFileName               = (char(*)[MAXCHARS])hmalloc(MAXBLOCKS * sizeof(*ParentFileName));
  FileName                     = (char(*)[MAXCHARS])hmalloc(MAXBLOCKS * sizeof(*FileName));
  LineNumber                   = (int *)hmalloc(MAXBLOCKS * sizeof(*LineNumber));
  HighMarkTabBuf               = (char *)hmalloc(TAB_BUF_SIZE);
  HighMarkTabBufWithoutGeneric = (char *)hmalloc(TAB_BUF_SIZE);

  memset(VarName, 0, MAXBLOCKS * sizeof(*VarName));
  memset(FunctionName, 0, MAXBLOCKS * sizeof(*FunctionName));
  memset(ParentFileName, 0, MAXBLOCKS * sizeof(*ParentFileName));
  memset(FileName, 0, MAXBLOCKS * sizeof(*FileName));

  size_t n = All.MaxMemSize * ((size_t)1024 * 1024);
  n        = roundup_to_multiple_of_cacheline_size(n);

  if(!(Base = hmalloc(n)))
    terminate("Failed to allocate memory for `Base' (%d Mbytes).", All.MaxMemSize);
#ifdef MEMORY_MANAGER_USE_MPROTECT
  /* make sure that the base pointer is page-aligned */
  myassert((uintptr_t)Base % sysconf(_SC_PAGESIZE) == 0);
  /* mark unused memory blocks as inaccessible */
  const int error = mprotect(Base, n, PROT_NONE);
  if(error)
    warn("mprotect(%p, %td, PROT_NONE) failed (%d): %s", Base, n, errno, strerror(errno));
#endif

  TotBytes = n;

  AllocatedBytes                  = 0;
  Nblocks                         = 0;
  HighMarkBytes                   = 0;
  HighMarkBytesWithoutGeneric     = 0;
  OldGlobHighMarkMB               = 0;
  OldGlobHighMarkMBWithoutGeneric = 0;
  HighMarkTabBuf[0]               = '\0';
  HighMarkTabBufWithoutGeneric[0] = '\0';
}

/*! \brief Writes memory usage in FdMemory.
 *
 *  \param[in] rank Number of tasks involved.
 *  \param[in] tabbuf Header message written in FdMemory.
 *
 *  \return void
 */
static void report_memory_usage(const int rank, const char *const tabbuf)
{
  if(ThisTask == rank)
    {
      const size_t buf_size = (100 + 4 * MAXCHARS) * (Nblocks + 10);
      char *buf             = (char *)mymalloc("buf", buf_size);
      int cc                = 0;
      cc += snprintf(buf + cc, buf_size - cc,
                     "\nMEMORY:  Largest Allocation = %g Mbyte  |  Largest Allocation Without Generic = %g Mbyte\n\n",
                     OldGlobHighMarkMB, OldGlobHighMarkMBWithoutGeneric);

      myassert(buf_size - cc > strlen(tabbuf));
      cc += snprintf(buf + cc, buf_size - cc, "%s", tabbuf);
      if(ThisTask == 0)
        {
          if(RestartFlag <= RESTART_SNAPSHOT)
            {
              fprintf(FdMemory, "%s", buf);
              fflush(FdMemory);
            }
        }
      else
        {
          MPI_Send(&cc, 1, MPI_INT, 0, TAG_N, MPI_COMM_WORLD);
          MPI_Send(buf, cc + 1, MPI_BYTE, 0, TAG_PDATA, MPI_COMM_WORLD);
        }
      myfree(buf);
    }

  if(ThisTask == 0 && rank > 0)
    {
      int cc;
      MPI_Recv(&cc, 1, MPI_INT, rank, TAG_N, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      char *buf = (char *)mymalloc("buf", cc + 1);
      MPI_Recv(buf, cc + 1, MPI_BYTE, rank, TAG_PDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if(RestartFlag <= RESTART_SNAPSHOT)
        {
          fprintf(FdMemory, "%s", buf);
          fflush(FdMemory);
        }
      myfree(buf);
    }
}

/*! \brief Writes memory useage of largest task in FdMemory.
 *
 *  \return void
 */
void report_detailed_memory_usage_of_largest_task(void)
{
  int flag = 0;

  struct
  {
    double mem;
    int rank;
  } local, global;

  local.mem  = HighMarkBytes / pow(1024, 2);
  local.rank = ThisTask;

  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  if(global.mem >= 1.05 * OldGlobHighMarkMB)
    {
      OldGlobHighMarkMB = global.mem;
      flag |= 1;
    }

  local.mem  = HighMarkBytesWithoutGeneric / pow(1024, 2);
  local.rank = ThisTask;

  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  if(global.mem >= 1.05 * OldGlobHighMarkMBWithoutGeneric)
    {
      OldGlobHighMarkMBWithoutGeneric = global.mem;
      flag |= 2;
    }

  if(flag & 2)
    report_memory_usage(global.rank, HighMarkTabBufWithoutGeneric);

  if(flag & 1)
    report_memory_usage(global.rank, HighMarkTabBuf);
}

/*! \brief Checks if a string is equal to any element in a set of other strings.
 *
 *  \param[in] str String to compare the other strings to.
 *  \param[in] size Maximum number of characters to compare (see strncmp()).
 *  \param[in] ... Set of strings that will be compared to str. The final
 *                 argument must be NULL;
 *
 *  \return true if any of the strings is equal to str (according to
 *          strncmp()), false otherwise.
 */
static bool strneq_any(const char *const str, const size_t size, ...)
{
  va_list args;
  va_start(args, size);
  bool result = false;
  char *next_str;
  while((next_str = va_arg(args, char *)))
    if(strncmp(str, next_str, size) == 0)
      {
        result = true;
        break;
      }
  va_end(args);
  return result;
}

/*! \brief Fills the output buffer with the memory log.
 *
 *  \param[out] p Output buffer.
 *  \param[in] buf_size Size of output buffer.
 *
 *  \return The number of characters written to p.
 */
static size_t dump_memory_table_buffer(char *const p, const size_t buf_size)
{
  size_t cc           = 0;
  size_t totBlocksize = 0;

  cc += snprintf(p + cc, buf_size - cc,
                 "-------------------------- Allocated Memory Blocks---- ( Step %8d )------------------\n"
                 "Task    Nr F                  Variable      MBytes   Cumulative  Function|File|Linenumber\n"
                 "------------------------------------------------------------------------------------------\n",
                 All.NumCurrentTiStep);
  myassert(cc < buf_size);
#ifdef GFM
  cc += snprintf(p + cc, buf_size - cc,
                 "-- Skipping: yieldsSNIa, yieldsSNII, yieldsAGB, netCoolingRate, netHeatingRate, netCoolingRateLowZ, "
                 "netHeatingRateLowZ, stellarLuminosityTable (if present) --\n");
  myassert(cc < buf_size);
#endif
  for(unsigned long i = 0; i < Nblocks; i++)
    {
      totBlocksize += BlockSize[i];

      if(strneq_any(VarName[i], sizeof(VarName[i]), "yieldsSNIa", "yieldsSNII", "yieldsAGB", "netCoolingRate", "netHeatingRate",
                    "netCoolingRateLowZ", "netHeatingRateLowZ", "stellarLuminosityTable", NULL))
        {
          continue;
        }
      cc += snprintf(p + cc, buf_size - cc, "%4d %5lu %d %40s  %10.4f   %10.4f  %s%s()|%s|%d\n", ThisTask, i, BasePointers[i] != NULL,
                     VarName[i], BlockSize[i] / pow(1024, 2), totBlocksize / pow(1024, 2), ParentFileName[i], FunctionName[i],
                     FileName[i], LineNumber[i]);
      myassert(cc < buf_size);
    }
  cc +=
      snprintf(p + cc, buf_size - cc, "------------------------------------------------------------------------------------------\n");
  myassert(cc < buf_size);

  return cc;
}

/*! \brief Dumps the buffer where the memory information is stored to the
 *         standard output.
 *
 *  \return void
 */
void dump_memory_table(void)
{
  const size_t buf_size = 200 * (Nblocks + 10);
  char *buf             = (char *)malloc(buf_size);
  const size_t nwritten = dump_memory_table_buffer(buf, buf_size);
  myassert(nwritten < buf_size);
  printf("%s", buf);
  free(buf);
}

/*! \brief Allocates a non-movable or movable memory block and stores the
 *         relevant information.
 *
 *  \param[in] ptr If NULL: allocate a non-movable memory block;
 *                 if not NULL: allocate a movable memory block, where ptr has
 *                 the same meaning as in #mymalloc_movable_fullinfo().
 *  \param[in] varname See #mymalloc_fullinfo() and
 *                     #mymalloc_movable_fullinfo().
 *  \param[in] n See #mymalloc_fullinfo() and #mymalloc_movable_fullinfo().
 *  \param[in] func See #mymalloc_fullinfo() and #mymalloc_movable_fullinfo().
 *  \param[in] file See #mymalloc_fullinfo() and #mymalloc_movable_fullinfo().
 *  \param[in] line See #mymalloc_fullinfo() and #mymalloc_movable_fullinfo().
 *  \param[in] clear_flag See #mymalloc_fullinfo().
 *  \param[in] callorigin See #mymalloc_fullinfo() and
 *                        #mymalloc_movable_fullinfo().
 *  \param[in] mymalloc_func Name of the calling function (#mymalloc_fullinfo()
 *                           or #mymalloc_movable_fullinfo()).
 *
 *  \return A pointer to the beginning of the allocated memory block.
 */
static void *mymalloc_impl(void **ptr, const char *varname, size_t n, const char *func, const char *file, const int line,
                           const bool clear_flag, const char *const callorigin, const char *const malloc_func)
{
  if(Nblocks >= MAXBLOCKS)
    terminate("No blocks left in %s() at %s()/%s/line %d. MAXBLOCKS=%d", malloc_func, func, file, line, MAXBLOCKS);

  n = roundup_to_multiple_of_cacheline_size(n);

  if(n > FreeBytes())
    {
      dump_memory_table();
      terminate("\nNot enough memory in %s() to allocate %g MB for variable '%s' at %s()/%s/line %d (FreeBytes=%g MB).", malloc_func,
                n / pow(1024, 2), varname, func, file, line, FreeBytes() / pow(1024, 2));
    }

  Table[Nblocks] = (char *)Base + AllocatedBytes;
  myassert(Nblocks == 0 || Table[Nblocks] == (char *)Table[Nblocks - 1] + BlockSize[Nblocks - 1]);
  mprotect_on_allocation(AllocatedBytes, n);

  strncpy(VarName[Nblocks], varname, sizeof(VarName[Nblocks]) - 1);
  strncpy(FunctionName[Nblocks], func, sizeof(FunctionName[Nblocks]) - 1);
  strncpy(FileName[Nblocks], file, sizeof(FileName[Nblocks]) - 1);
  LineNumber[Nblocks] = line;
  if(callorigin)
    {
      strncpy(ParentFileName[Nblocks], callorigin, sizeof(ParentFileName[Nblocks]) - 1);
      GenericFlag[Nblocks] = true;
      AllocatedBytesGeneric += n;
    }
  else
    {
      memset(ParentFileName[Nblocks], 0, sizeof(ParentFileName[Nblocks]));
      GenericFlag[Nblocks] = false;
    }

  AllocatedBytes += n;
  BlockSize[Nblocks]    = n;
  BasePointers[Nblocks] = ptr;

  Nblocks++;

  if(AllocatedBytes - AllocatedBytesGeneric > HighMarkBytesWithoutGeneric)
    {
      HighMarkBytesWithoutGeneric = AllocatedBytes - AllocatedBytesGeneric;
      const size_t nwritten       = dump_memory_table_buffer(HighMarkTabBufWithoutGeneric, TAB_BUF_SIZE);
      myassert(nwritten < TAB_BUF_SIZE);
    }
  if(AllocatedBytes > HighMarkBytes)
    {
      HighMarkBytes         = AllocatedBytes;
      const size_t nwritten = dump_memory_table_buffer(HighMarkTabBuf, TAB_BUF_SIZE);
      myassert(nwritten < TAB_BUF_SIZE);
    }

  if(clear_flag)
    memset(Table[Nblocks - 1], 0, n);

  return Table[Nblocks - 1];
}

/*! \brief Allocates a non-movable memory block and stores the relevant
 *         information.
 *
 *  \param[in] varname Name of the variable to be stored in the allocated
 *             block.
 *  \param[in] n Size of the memory block in bytes.
 *  \param[in] func Name of function that has called the allocation routine
 *                  (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the allocation
 *                  routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of file where the allocation routine was
 *                  called (usually given by the __LINE__ macro).
 *  \param[in] clear_flag If non-zero, all bytes in the returned memory block
 *                        will be initialized with 0 (analogous to #calloc()).
 *  \param[in] callorigin Used to record the location of the original call for
 *                        the generic communication templates
 *                        (generic_comm_helpers2.h, generic_comm_helpers_async.h)
 *                        with #mymalloc_g(). Should be NULL for normal use.
 *
 *  \return A pointer to the beginning of the allocated memory block.
 */
void *mymalloc_fullinfo(const char *const varname, const size_t n, const char *const func, const char *const file, const int line,
                        const int clear_flag, const char *const callorigin)
{
  return mymalloc_impl(NULL, varname, n, func, file, line, clear_flag, callorigin, __func__);
}

/*! \brief Allocates a movable memory block and stores the relevant
 *         information.
 *
 *  \param[in] ptr Pointer to the variable that will hold the memory address of
 *                 the block (as returned by this function). This pointer will
 *                 be used to update the variable with a new pointer in the
 *                 case that memory blocks are moved (i.e. when
 *                 #myfree_movable_fullinfo() or #myrealloc_movable_fullinfo()
 *                 are called).
 *                 Note that this can be problematic if the pointer to the
 *                 block is stored in variables other than the one indicated by
 *                 ptr (e.g. when passing the pointer to a function). Only the
 *                 pointer stored in the variable corresponding to ptr will be
 *                 updated when blocks are moved.
 *  \param[in] varname Name of the variable to be stored in the allocated block.
 *  \param[in] n Size of the memory block in bytes.
 *  \param[in] func Name of function that has called the allocation routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the allocation
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of file where the allocation routine was
 *             called (usually given by the __LINE__ macro).
 *  \param[in] clear_flag If non-zero, all bytes in the returned memory block
 *                        will be initialized with 0 (analogous to #calloc()).
 *  \param[in] callorigin Used to record the location of the original call for
 *                        the generic communication templates
 *                        (generic_comm_helpers2.h, generic_comm_helpers_async.h)
 *                        with #mymalloc_movable_g(). Should be NULL for normal
 *                        use.
 *
 *  \return A pointer to the beginning of the allocated memory block.
 */
void *mymalloc_movable_fullinfo(void *const ptr, const char *const varname, const size_t n, const char *const func,
                                const char *const file, const int line, const int clear_flag, const char *const callorigin)
{
  myassert(ptr);
  return mymalloc_impl((void **)ptr, varname, n, func, file, line, clear_flag, callorigin, __func__);
}

static void check_block_nonmovable(void *const p, const char *const func, const char *const file, const int line,
                                   const char *const mymalloc_func)
{
  if(p != Table[Nblocks - 1])
    {
      dump_memory_table();
      terminate("Wrong call of %s() at %s()/%s/line %d: not the last allocated block!", mymalloc_func, func, file, line);
    }
}

static unsigned long find_block(void *const p, const char *const func, const char *const file, const int line,
                                const char *const mymalloc_func)
{
  long nr;
  for(nr = Nblocks - 1; nr >= 0; nr--)
    if(p == Table[nr])
      break;
  if(nr < 0)
    {
      dump_memory_table();
      terminate("Wrong call of %s() from %s()/%s/line %d: this block has not been allocated!", mymalloc_func, func, file, line);
    }
  if((unsigned long)nr < Nblocks - 1)
    {
      /* the block is not the last allocated block;
       * check that all subsequent blocks are actually movable */
      for(unsigned long i = nr + 1; i < Nblocks; i++)
        if(!BasePointers[i])
          {
            dump_memory_table();
            myflush(stdout);
            terminate("Wrong call of %s() from %s()/%s/line %d: behind block=%lu, there are subsequent non-movable allocated blocks!",
                      mymalloc_func, func, file, line, nr);
          }
    }
  return nr;
}

/*! \brief Finds last allocated block.
 *
 *  \return void pointer to last allocated block.
 */
void *myfree_query_last_block(void)
{
  myassert(Nblocks > 0);
  return Table[Nblocks - 1];
}

/*! \brief Deallocates a non-movable or movable memory block.
 *
 *  \param[in] p See #myfree_fullinfo() and #myfree_movable_fullinfo().
 *  \param[in] func See #myfree_fullinfo() and #myfree_movable_fullinfo().
 *  \param[in] file See #myfree_fullinfo() and #myfree_movable_fullinfo().
 *  \param[in] line See #myfree_fullinfo() and #myfree_movable_fullinfo().
 *  \param[in] movable Whether the deallocation is supposed to be movable
 *                     (true) or non-movable (false).
 *  \param[in] free_func Name of the calling function (#myfree_fullinfo()
 *                       or #myfree_movable_fullinfo()).
 */
static void myfree_impl(void *p, const char *const func, const char *const file, const int line, const bool movable,
                        const char *const free_func)
{
  myassert(Nblocks > 0);
  if(!movable)
    check_block_nonmovable(p, func, file, line, free_func);
  /* first, let's find the block */
  const unsigned long nr = find_block(p, func, file, line, free_func);
  /* collect the block's relevant properties */
  const bool generic_flag = GenericFlag[nr];
  const size_t block_size = BlockSize[nr];
  const ptrdiff_t offset  = -block_size;
  size_t length           = 0;
  for(unsigned long i = nr + 1; i < Nblocks; i++)
    length += BlockSize[i];
  /* the following operations will move the memory block if necessary */
  if(nr < Nblocks - 1)
    {
      /* move the memory block data */
      myassert((char *)Table[nr + 1] + offset == Table[nr]);
      memmove(Table[nr], Table[nr + 1], length);
    }
  for(unsigned long i = nr + 1; i < Nblocks; i++)
    {
      /* update pointers to the new locations of the moved memory blocks */
      Table[i]         = (char *)Table[i] + offset;
      *BasePointers[i] = (char *)*BasePointers[i] + offset;
    }
  for(unsigned long i = nr + 1; i < Nblocks; i++)
    {
      /* shift array elements to the left to get rid of the gap created at
       * element nr */
      Table[i - 1]        = Table[i];
      BasePointers[i - 1] = BasePointers[i];
      BlockSize[i - 1]    = BlockSize[i];
      GenericFlag[i - 1]  = GenericFlag[i];
      memcpy(VarName[i - 1], VarName[i], sizeof(VarName[i]));
      memcpy(FunctionName[i - 1], FunctionName[i], sizeof(FunctionName[i]));
      memcpy(ParentFileName[i - 1], ParentFileName[i], sizeof(ParentFileName[i]));
      memcpy(FileName[i - 1], FileName[i], sizeof(FileName[i]));
      LineNumber[i - 1] = LineNumber[i];
    }
  /* update book-keeping variables */
  if(generic_flag)
    AllocatedBytesGeneric -= block_size;
  AllocatedBytes -= block_size;
  Nblocks--;
  mprotect_on_allocation(AllocatedBytes, offset);
}

/*! \brief Deallocates a non-movable memory block.
 *
 *  For this operation to be successful the block that has to be deallocated
 *  must be the last allocated one.
 *
 *  \param[in] p Pointer to the memory block to be deallocated.
 *  \param[in] func Name of function that has called the deallocation routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the deallocation
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of file where the deallocation routine was
 *             called (usually given by the __LINE__ macro).
 */
void myfree_fullinfo(void *const p, const char *const func, const char *const file, const int line)
{
  myfree_impl(p, func, file, line, false, __func__);
}

/*! \brief Deallocates a movable memory block.
 *
 *  For this operation to be successful, all the blocks allocated after the
 *  block that has to be freed must be of movable type.
 *
 *  \param[in] p pointer to the memory block to be deallocated.
 *  \param[in] func name of function that has called the deallocation routine
 *             (usually given by the __func__ macro).
 *  \param[in] file file where the function that has called the deallocation
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line line number of file where the deallocation routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return void
 */
void myfree_movable_fullinfo(void *const p, const char *const func, const char *const file, const int line)
{
  myfree_impl(p, func, file, line, true, __func__);
}

/*! \brief Reallocates an existing non-movable or movable memory block.
 *
 *  \param[in] p See #myrealloc_fullinfo() and #myrealloc_movable_fullinfo().
 *  \param[in] n See #myrealloc_fullinfo() and #myrealloc_movable_fullinfo().
 *  \param[in] func See #myrealloc_fullinfo() and
 *                  #myrealloc_movable_fullinfo().
 *  \param[in] file See #myrealloc_fullinfo() and
 *                  #myrealloc_movable_fullinfo().
 *  \param[in] line See #myrealloc_fullinfo() and
 *                  #myrealloc_movable_fullinfo().
 *  \param[in] movable Whether the reallocation is supposed to be movable
 *                     (true) or non-movable (false).
 *  \param[in] realloc_func Name of the calling function (#myrealloc_fullinfo()
 *                          or #myrealloc_movable_fullinfo()).
 *
 *  \return A pointer to the beginning of the newly allocated memory block.
 */
static void *myrealloc_impl(void *p, size_t n, const char *const func, const char *const file, const int line, const bool movable,
                            const char *const realloc_func)
{
  myassert(Nblocks > 0);
  if(!movable)
    check_block_nonmovable(p, func, file, line, realloc_func);

  n = roundup_to_multiple_of_cacheline_size(n);

  /* first, let's find the block */
  const unsigned long nr = find_block(p, func, file, line, realloc_func);

  /* movable reallocation is not used by the generic communication templates */
  myassert(!movable || !GenericFlag[nr]);

  if(GenericFlag[nr])
    AllocatedBytesGeneric -= BlockSize[nr];
  AllocatedBytes -= BlockSize[nr];

  if(n > FreeBytes())
    {
      dump_memory_table();
      terminate("Not enough memory in %s(n=%g MB) at %s()/%s/line %d. previous=%g FreeBytes=%g MB", realloc_func, n / pow(1024, 2),
                func, file, line, BlockSize[nr] / pow(1024, 2), FreeBytes() / pow(1024, 2));
    }

  const ptrdiff_t offset = n - BlockSize[nr];
  size_t length          = 0;
  for(unsigned long i = nr + 1; i < Nblocks; i++)
    length += BlockSize[i];

  if(offset > 0)
    mprotect_on_allocation(AllocatedBytes + BlockSize[nr], offset);

  if(nr < Nblocks - 1)
    memmove((char *)Table[nr + 1] + offset, Table[nr + 1], length);

  for(unsigned long i = nr + 1; i < Nblocks; i++)
    {
      Table[i]         = (char *)Table[i] + offset;
      *BasePointers[i] = (char *)*BasePointers[i] + offset;
    }

  myassert(movable || Table[nr] == (char *)Base + AllocatedBytes);

  AllocatedBytes += n;
  BlockSize[nr] = n;

  /* no need to check for #HighMarkBytesWithoutGeneric because movable
   * reallocation is not used by the generic communication templates
   * (see above) */
  if(AllocatedBytes > HighMarkBytes)
    {
      HighMarkBytes         = AllocatedBytes;
      const size_t nwritten = dump_memory_table_buffer(HighMarkTabBuf, TAB_BUF_SIZE);
      myassert(nwritten < TAB_BUF_SIZE);
    }

  if(offset < 0)
    mprotect_on_allocation(AllocatedBytes, offset);

  return Table[nr];
}

/*! \brief Reallocates an existing non-movable memory block.
 *
 *  For this operation to be successful this must be the last allocated block.
 *
 *  \param[in] p Pointer to the existing memory block to be reallocated.
 *  \param[in] n The new size of the memory block in bytes.
 *  \param[in] func Name of function that has called the reallocation routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the reallocation
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of file where the reallocation routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return A pointer to the beginning of the newly allocated memory block.
 */
void *myrealloc_fullinfo(void *const p, const size_t n, const char *const func, const char *const file, const int line)
{
  return myrealloc_impl(p, n, func, file, line, false, __func__);
}

/*! \brief Reallocates an existing movable memory block.
 *
 *  For this operation to be successful, all the blocks allocated after the
 *  block that has to be reallocated must be of movable type.
 *
 *  \param[in] p Pointer to the existing memory block to be reallocated.
 *  \param[in] n The new size of the memory block in bytes.
 *  \param[in] func Name of function that has called the reallocation routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the reallocation
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of file where the reallocation routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return A pointer to the beginning of the newly allocated memory block.
 */
void *myrealloc_movable_fullinfo(void *const p, size_t n, const char *const func, const char *const file, const int line)
{
  return myrealloc_impl(p, n, func, file, line, true, __func__);
}
