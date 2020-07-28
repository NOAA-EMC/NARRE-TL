
set -x
# Loading Intel Compiler Suite
module load ics/12.1

# Loding nceplibs modules
module load bacio/v2.0.2
module load w3emc/v2.2.0
module load w3nco/v2.0.6
module load g2/v2.5.0
module load jasper/v1.900.1
module load z/v1.2.6
module load png/v1.2.44

make
