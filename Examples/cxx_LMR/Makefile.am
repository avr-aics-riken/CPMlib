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
noinst_PROGRAMS=exampleCXX_LMR

DISTCLEANFILES=*~
CLEANFILES=data/*.log

exampleCXX_LMR_SOURCES = \
    main.C \
    ./include/_mpi_test.h \
    README \
    ./data/leaf15.tp \
    ./data/leaf15.oct \
    ./data/run.sh

exampleCXX_LMR_CXXFLAGS = \
    -I$(top_builddir)/include \
    -I$(top_builddir)/include/BCM \
    @MPICH_CFLAGS@ \
    @TP_CFLAGS@ \
    @REALOPT@

exampleCXX_LMR_LDADD = \
    -L$(top_builddir)/src \
    @CPM_LIBS@ \
    @MPICH_LDFLAGS@ \
    @MPICH_LIBS@ \
    @TP_LDFLAGS@

dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)

