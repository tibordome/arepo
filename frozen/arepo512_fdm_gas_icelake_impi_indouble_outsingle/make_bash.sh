#!/bin/bash

module purge
module load rhel8/default-icl
module load intel/compilers/2018.4
#module load intel/impi/2018.4/intel # causes PMPI_Init 14 error
#module load python
