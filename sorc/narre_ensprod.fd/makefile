#!/bin/sh
#
#  This makefile is to compile ensembel production generator
#  Author: Binbin Zhou
#          Aug 1, 2005
#  OS: Linux (on WCOSS)
#  Modified for Cray: Binbin Zhou Dec. 16, 2016
#

FC=     ftn -convert big_endian

INC= ${G2_INC4}

LIBS= ${IP_LIB4} ${W3NCO_LIB4} ${W3EMC_LIB4} ${BACIO_LIB4}  ${SP_LIB4} ${G2_LIB4} ${JASPER_LIB} ${PNG_LIB} ${Z_LIB}

FFLAGS = -O3 -auto -I $(INC)

SRCS=	char_process.f util.f readGB2.f print_gfld.f Zputgb2.f \
	packGB2_mean.f packGB2_prob.f packGB2_mean_derv.f packGB2_prob_derv.f \
	preciptype.f wind.f get_icing.f get_cat.f getflight.f Hains_index.f getceiling.f \
	new_fog.f get_new_fog.f thickness.f llws.f getconv.f \
	get_mxp.f mxp.f packGB2_mxp.f \
	get_cptp_severe.f calibrate2_hrly_rgn3.f calibrate2_dryt.f  calibrate2_svr.f calibrate2.f \
	meanwind.f get_time_packGB2.f


OBJ=${SRCS:.f=.o}

narre: $(OBJ)
	$(FC) $(FFLAGS) -o narre_ensprod sref_ens_gen.f $(OBJ)  $(LIBS)
clean:
	/bin/rm  *.o
copy:
	/bin/cp narre_ensprod ../../exec/narre_ensprod
#

