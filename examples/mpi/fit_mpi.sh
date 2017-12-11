#!/bin/bash

CHAIN=${SGE_TASK_ID:-1}
if [ "$CHAIN" = undefined ];
then
    CHAIN=1
fi


# this is where my boost lives
export DYLD_LIBRARY_PATH=$HOME/local/lib


#mpirun -np 1 ./oral_2cmt_mpi3 diagnose data file=mpi_testgrad-0.R init=0.5 random seed=12
#mpirun -np 1 ./oral_2cmt_mpi3 diagnose data file=mpi_testgrad-1.R init=0.5 random seed=12
#mpirun -np 1 ./oral_2cmt_mpi3 diagnose data file=mpi_testgrad-2.R init=0.5 random seed=12

##time mpirun -np 4 ./oral_2cmt_mpi3 sample num_samples=20 num_warmup=150 save_warmup=0 algorithm=hmc stepsize=0.01 adapt delta=0.8 id=$CHAIN  data file=mpi_test_parallel.R init=mpi_init.R random seed=12 output refresh=1 file=samples-$CHAIN-mpi.csv

##exit 0

#time mpirun -np 1 ./oral_2cmt_mpi3 sample num_samples=20 num_warmup=150 save_warmup=1 algorithm=hmc stepsize=0.01 adapt delta=0.8 id=$CHAIN  data file=mpi_testgrad-0.R init=0.5 random seed=12 output refresh=10 file=samples-$CHAIN-rect-0.csv

#time mpirun -np 4 ./oral_2cmt_mpi3 sample num_samples=20 num_warmup=150 save_warmup=1 algorithm=hmc stepsize=0.01 adapt delta=0.8 id=$CHAIN  data file=mpi_testgrad-0.R init=0.5 random seed=12 output refresh=10 file=samples-$CHAIN-rect-0-2.csv

time mpirun -np 1 ./oral_2cmt_mpi3 sample num_samples=20 num_warmup=150 save_warmup=1 algorithm=hmc stepsize=0.01 adapt delta=0.8 id=$CHAIN  data file=mpi_testgrad-1.R init=0.5 random seed=12 output refresh=10 file=samples-$CHAIN-rect-1.csv

time mpirun -np 1 ./oral_2cmt_mpi3 sample num_samples=20 num_warmup=150 save_warmup=1 algorithm=hmc stepsize=0.01 adapt delta=0.8 id=$CHAIN  data file=mpi_testgrad-2.R init=0.5 random seed=12 output refresh=10 file=samples-$CHAIN-rect-2.csv

##mpi_init.R

#time mpirun -np 1 ./oral_2cmt_mpi3 sample num_samples=20 num_warmup=150 save_warmup=0 algorithm=hmc stepsize=0.01 adapt delta=0.8 id=$CHAIN  data file=mpi_test_serial.R init=mpi_init.R random seed=12 output refresh=1 file=samples-$CHAIN-serial.csv

