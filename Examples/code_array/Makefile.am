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
noinst_PROGRAMS=kernel_CPM

DISTCLEANFILES=*~
CLEANFILES=data/*.log

kernel_CPM_SOURCES = \
    main.cpp \
    voxinfo.cpp \
    voxinfo.h \
    LS.cpp \
    LS.h \
    FBdefine.h \
    FortFunc.h \
    kernel_def.h \
    linear_solver.f90 \
    utility.f90 \
    ffv_blas.f90 \
    cbc_f_params.h

kernel_CPM_CXXFLAGS = \
    -I$(top_srcdir)/include \
    @MPICH_CFLAGS@ \
    @TP_CFLAGS@ \
    @REALOPT@ \
    @PM_CFLAGS@

kernel_CPM_FCFLAGS = \
    -I$(top_srcdir)/include \
    @MPICH_CFLAGS@ \
    @FREALOPT@

kernel_CPM_LDADD = \
    -L$(top_builddir)/src \
    @CPM_LIBS@ \
    @MPICH_LDFLAGS@ \
    @MPICH_LIBS@ \
    @TP_LDFLAGS@ \
    @PM_LDFLAGS@

dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)

