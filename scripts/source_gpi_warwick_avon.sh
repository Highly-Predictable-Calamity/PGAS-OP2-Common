module purge

module load GCC/10.2.0 OpenMPI/4.0.5 HDF5

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.local/lib64
export OP2_COMPILER=gnu
export GPI_INSTALL_PATH=$HOME/.local
export PARMETIS_INSTALL_PATH=$HOME/.local

