#!/bin/sh
set -aux

source /etc/profile.d/lmod.sh

SORCnarre=`pwd`

module purge
module use $SORCnarre/../modulefiles
module load v1.1.0
#module unload g2
#module list
set -x


#export G2_INC4=/gpfs/dell1/nco/ops/nwprod/lib/g2/v3.1.0/ips/18.0.1/include/g2_v3.1.0_4
#export G2_LIB4=/gpfs/dell1/nco/ops/nwprod/lib/g2/v3.1.0/ips/18.0.1/libg2_v3.1.0_4.a


EXEC=`pwd`/../exec
mkdir -p $EXEC

#compile product generator sorc
cd $SORCnarre/narre_ensprod.fd
make
make copy
make clean 

#compile Precip code
cd $SORCnarre/narre_calPrecip.fd
make narre_precip
mv  narre_precip ${EXEC}/.
make clean

make narre_precip3hr
mv narre_precip3hr ${EXEC}/.
make clean

echo "All compilatios done"


