/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/mpi_utils/hypercube_allgatherv.c
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

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef MPI_HYPERCUBE_ALLGATHERV

#define TAG 100

/*! \brief Allgatherv routine based on MPI_Sendrecv calls.
 *
 *  \param[in] sendbuf Starting address of send buffer.
 *  \param[in] sendcount Number of elements in send buffer.
 *  \param[in] sendtype Data type of send buffer elements.
 *  \param[out] recvbuf Address of receive buffer.
 *  \param[in] recvcount Integer array (of length group size) containing the
 *             number of elements that are to be received from each process.
 *  \param[in] displs Integer array (of length group size). Entry i specifies
 *             the displacement (relative to recvbuf ) at which to place the
 *             incoming data from process.
 *  \param[in] recvtype Data type of receive buffer elements.
 *  \param[in] comm Communicator.
 *
 *  \return 0
 */
int MPI_hypercube_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int *recvcount, int *displs,
                             MPI_Datatype recvtype, MPI_Comm comm)
{
  int ntask, thistask, ptask, ngrp, size_sendtype, size_recvtype;
  MPI_Status status;

  MPI_Comm_rank(comm, &thistask);
  MPI_Comm_size(comm, &ntask);

  MPI_Type_size(sendtype, &size_sendtype);
  MPI_Type_size(recvtype, &size_recvtype);

  for(ptask = 0; ntask > (1 << ptask); ptask++)
    ;

  for(ngrp = 1; ngrp < (1 << ptask); ngrp++)
    {
      int recvtask = thistask ^ ngrp;

      if(recvtask < ntask)
        MPI_Sendrecv(sendbuf, sendcount, sendtype, recvtask, TAG, recvbuf + displs[recvtask] * size_recvtype, recvcount[recvtask],
                     recvtype, recvtask, TAG, comm, &status);
    }

  if(sendbuf != recvbuf + displs[thistask] * size_recvtype)
    memcpy(recvbuf + displs[thistask] * size_recvtype, sendbuf, sendcount * size_sendtype);

  return 0;
}

#endif
