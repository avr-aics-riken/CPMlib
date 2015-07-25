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
    ./data/mconv.tp \
    ./data/leaf15.oct \
    ./data/run.sh


mconvp_LMR_CXXFLAGS = \
    -I$(top_srcdir)/include \
    @MPICH_CFLAGS@ \
    @TP_CFLAGS@ \
    @REALOPT@

mconvp_LMR_FCFLAGS = \
    -I$(top_srcdir)/include \
    @MPICH_CFLAGS@ \
    @FREALOPT@

#mconvp_LMR_LINK = @CPM_FC_LD@ $(FCFLAGS) \
#  $(AM_LDFLAGS) $(LDFLAGS) -o $@


mconvp_LMR_LDADD = \
    -L$(top_builddir)/src \
    @CPM_LIBS@ \
    @MPICH_LDFLAGS@ \
    @MPICH_FCLIBS@ \
    @TP_LDFLAGS@ \
    @TP_LIBS@

dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)

