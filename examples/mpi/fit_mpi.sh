#!/bin/bash

CHAIN=${SGE_TASK_ID:-1}
if [ "$CHAIN" = undefined ];
then
    CHAIN=1
fi


# this is where my boost lives
export DYLD_LIBRARY_PATH=$HOME/local/lib

#time mpirun -np 4 ./oral_2cmt_mpi sample num_samples=20 num_warmup=150 save_warmup=0 algorithm=hmc stepsize=0.01 adapt delta=0.8 id=$CHAIN  data file=mpi_test_parallel4.R init=mpi_init.R random seed=12 output refresh=1 file=samples-$CHAIN-mpi.csv

time mpirun -np 1 ./oral_2cmt_mpi sample num_samples=20 num_warmup=150 save_warmup=0 algorithm=hmc stepsize=0.01 adapt delta=0.8 id=$CHAIN  data file=mpi_test_serial.R init=mpi_init.R random seed=12 output refresh=1 file=samples-$CHAIN-serial.csv

