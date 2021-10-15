#!/bin/sh
#
#  This makefile is to compile ensembel production generator
#  Author: Binbin Zhou
#  OS: Intel WCOSS
#

FC=     ftn -convert big_endian 
INC= ${G2_INC4}
LIBS= ${IP_LIB4} ${W3NCO_LIB4} ${W3EMC_LIB4} ${BACIO_LIB4}  ${SP_LIB4} ${G2_LIB4} ${JASPER_LIB} ${PNG_LIB} ${Z_LIB}
FFLAGS = -O3 -auto -I $(INC)
narre_precip:
	$(FC) $(FFLAGS) -o $(@) narre_1hr_precip.f $(LIBS)
narre_precip3hr:
	$(FC) $(FFLAGS) -o $(@) narre_3hr_precip.f $(LIBS)
clean:
	-rm -f narre_precip narre_precip3hr
#

