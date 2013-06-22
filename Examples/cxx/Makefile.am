###
###################################################################################
#
# CPMlib - Cartesian Partition Manager Library
# 
# Copyright (C) 2012-2013 Institute of Industrial Science, The University of Tokyo.
# All right reserved.
#
###################################################################################
###
noinst_PROGRAMS=exampleCXX

DISTCLEANFILES=*~
CLEANFILES=data/*.log

exampleCXX_SOURCES = \
    main.C \
    ./include/_debug.h \
    ./include/_mpi_test.h \
    README \
    ./data/input1.txt \
    ./data/input2.txt \
    ./data/run.sh \
    ./data_subdomain/input.txt \
    ./data_subdomain/subdomain.dat \
    ./data_subdomain/run.sh

exampleCXX_CXXFLAGS = \
    -I$(top_builddir)/include \
    @MPICH_CFLAGS@ \
    @TP_CFLAGS@ \
    @REALOPT@

exampleCXX_LDADD = \
    -L$(top_builddir)/src \
    @CPM_LIBS@ \
    @MPICH_LDFLAGS@ \
    @MPICH_LIBS@ \
    @TP_LDFLAGS@ \
    @TP_LIBS@

dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)

