/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/mpi_utils/sizelimited_sendrecv.c
 * \date        MM/YYYY
 * \author
 * \brief       implements a wrapper around MPI_Sendrecv that, if needed,
 *              transmits the data in smaller pieces than a prescribed maximum size
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

/*! \brief Self-made sendrecv function with limiter to the number of elements
 *         that can be sent in one go.
 *
 *  If the total message is longer, multiple MPI_Sendrecv calls are executed
 *  until the entire message has been communicated.
 *
 *  \param[in] sendb Initial address of send buffer.
 *  \param[in] sendcount Number of elements in send buffer.
 *  \param[in] sendtype Type of elements in send buffer (handle).
 *  \param[in] dest Rank of destination.
 *  \param[in] sendtag Send tag.
 *  \param[out] recvb Initial address of receive buffer.
 *  \param[in] recvcount Number of elements in receive buffer.
 *  \param[in] recvtype Type of elements in receive buffer (handle).
 *  \param[in] source Rank of source.
 *  \param[in] recvtag Receive tag.
 *  \param[in] comm MPI communicator.
 *  \param[out] status Status, referring to receive operation.
 *
 *  \return 0
 */
int myMPI_Sendrecv(void *sendb, size_t sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvb, size_t recvcount,
                   MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status)
{
  int iter      = 0, size_sendtype, size_recvtype, send_now, recv_now;
  char *sendbuf = (char *)sendb;
  char *recvbuf = (char *)recvb;

  if(dest != source)
    terminate("dest != source");

  MPI_Type_size(sendtype, &size_sendtype);
  MPI_Type_size(recvtype, &size_recvtype);

  int thistask;
  MPI_Comm_rank(comm, &thistask);

  if(dest == thistask)
    {
      memcpy(recvbuf, sendbuf, recvcount * size_recvtype);
      return 0;
    }

  size_t count_limit = MPI_MESSAGE_SIZELIMIT_IN_BYTES / size_sendtype;

  while(sendcount > 0 || recvcount > 0)
    {
      if(sendcount > count_limit)
        {
          send_now = count_limit;
          /*
             if(iter == 0)
             {
             printf("Imposing size limit on MPI_Sendrecv() on task=%d (send of size=%lld)\n", ThisTask, (long long) sendcount *
             size_sendtype); myflush(stdout);
             }
           */
          iter++;
        }
      else
        send_now = sendcount;

      if(recvcount > count_limit)
        recv_now = count_limit;
      else
        recv_now = recvcount;

      MPI_Sendrecv(sendbuf, send_now, sendtype, dest, sendtag, recvbuf, recv_now, recvtype, source, recvtag, comm, status);

      sendcount -= send_now;
      recvcount -= recv_now;

      sendbuf += send_now * size_sendtype;
      recvbuf += recv_now * size_recvtype;
    }

  return 0;
}
