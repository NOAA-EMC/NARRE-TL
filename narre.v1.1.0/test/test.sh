#!/bin/ksh

export RUN_ENVIR=dev
export PDY=20191028
cyc=12
export model_ver=v1.1.0

export NDATE=/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.0/exec/ips/ndate
export RAPprod=/gpfs/hps/nco/ops/com/rap/prod
export COMINnam=/gpfs/dell1/nco/ops/com/nam/prod
./JNARRE_GETNAM.non_poe $cyc
./JNARRE_GETRAP.non_poe $cyc
./JNARRE_ENSPROD.non_poe $cyc

