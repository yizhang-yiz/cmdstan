#!/bin/bash

CHAIN=${SGE_TASK_ID:-1}
if [ "$CHAIN" = undefined ];
then
    CHAIN=1
fi


# this is where my boost lives
export DYLD_LIBRARY_PATH=$HOME/local/lib


J=50
NITER=500

cat <(cat oral4_stan-base.R) <(printf "J <- $J\nuse_map_rect <- 0\n") > oral4_stan-$J-0.R
cat <(cat oral4_stan-base.R) <(printf "J <- $J\nuse_map_rect <- 1\n") > oral4_stan-$J-1.R
cat <(cat oral4_stan-base.R) <(printf "J <- $J\nuse_map_rect <- 2\n") > oral4_stan-$J-2.R

#exit 0
## MPI map_rect version
time mpirun -np 12 ./oral_2cmt_mpi4 sample num_samples=$NITER num_warmup=$NITER save_warmup=1 algorithm=hmc stepsize=0.01 adapt delta=0.9 id=$CHAIN  data file=oral4_stan-$J-0.R init=0.1 random seed=10 output refresh=10 file=samples5-$CHAIN-rect-$J-0.csv

exit 0

## we send the two serial programs into the background and wait for
## both of them to finish

## serial map_rect version
time ./oral_2cmt_mpi4 sample num_samples=$NITER num_warmup=$NITER save_warmup=1 algorithm=hmc stepsize=0.01 adapt delta=0.9 id=$CHAIN  data file=oral4_stan-$J-1.R init=0.1 random seed=10 output refresh=10 file=samples5-$CHAIN-rect-$J-1.csv &

## vanilla Stan version
time ./oral_2cmt_mpi4 sample num_samples=$NITER num_warmup=$NITER save_warmup=1 algorithm=hmc stepsize=0.01 adapt delta=0.9 id=$CHAIN  data file=oral4_stan-$J-2.R init=0.1 random seed=10 output refresh=10 file=samples5-$CHAIN-rect-$J-2.csv &

wait
