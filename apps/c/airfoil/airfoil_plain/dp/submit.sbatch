#!/bin/bash

#SBATCH --nodes=6
#SBATCH --job-name PGAS-OP2
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4096
#SBATCH --time=4:00

module purge

module load GCC/10.2.0 OpenMPI/4.0.5

mpirun ./airfoil_mpi_seq
