###
###################################################################################
#
# CPMlib - Cartesian Partition Manager Library
# 
# Copyright (C) 2012-2014 Institute of Industrial Science, The University of Tokyo.
# All right reserved.
#
# Copyright (c) 2014-2015 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
###################################################################################
###
noinst_PROGRAMS=mconvp mconvp_LMR

DISTCLEANFILES=*~
CLEANFILES=*.mod data/*.sph

mconvp_SOURCES = \
    mconvp.f90

mconvp_LMR_SOURCES = \
    mconvp_LMR_main.C \
    mconvp_LMR_main.h \
    mconvp_LMR.f90 \
    README \
    ./data/mconv_ftcs.tp \
    ./data/mconv_jacobi.tp \
    ./data/leaf15.oct \
    ./data/run_ftcs.sh \
    ./data/run_jacobi.sh


mconvp_LMR_CXXFLAGS = \
    -I$(top_srcdir)/include \
    -I$(top_srcdir)/include/LMR \
    -I$(top_srcdir)/include/LMR/BCM \
    @MPICH_CFLAGS@ \
    @TP_CFLAGS@ \
    @REALOPT@

mconvp_LMR_FCFLAGS = \
    -I$(top_srcdir)/include \
    @FCFLAGS@ \
    @MPICH_CFLAGS@ \
    @FREALOPT@

#mconvp_LMR_LINK = @CPM_FC_LD@ $(FCFLAGS) \
#  $(AM_LDFLAGS) $(LDFLAGS) -o $@


mconvp_LMR_LDADD = \
    -L$(top_builddir)/src \
    -L$(top_builddir)/src/LMR \
    @CPMLMR_LIBS@ \
    @MPICH_LDFLAGS@ \
    @MPICH_FCLIBS@ \
    @TP_LDFLAGS@ \
    @TP_LIBS@

dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)

