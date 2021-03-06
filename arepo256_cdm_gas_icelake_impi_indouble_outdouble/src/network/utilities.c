/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/network/utilities.c
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "../proto.h"

#include "utilities.h"

void myprintf(const char *format, ...)
{
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
}

char *util_fgets(char *str, int num, FILE *stream, char *file, int line)
{
  char *ret = fgets(str, num, stream);
  if(ret == NULL)
    {
      printf("error: fgets in file %s at line %d\n", file, line);
      terminate("stop");
    }

  return ret;
}

size_t util_fread(void *ptr, size_t size, size_t count, FILE *stream, char *file, int line)
{
  size_t result = fread(ptr, size, count, stream);
  if(result != count)
    {
      printf("error: fread in file %s at line %d\n", file, line);
      terminate("stop");
    }

  return result;
}
