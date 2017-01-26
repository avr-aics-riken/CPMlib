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

AUTOMAKE_OPTIONS = subdir-objects

lib_LIBRARIES = libCPM.a

libCPM_a_CXXFLAGS = @MPICH_CFLAGS@ @TP_CFLAGS@ @PM_CFLAGS@ @REALOPT@ \
                    -I$(top_srcdir)/include 

libCPM_a_SOURCES = \
    cpm_DomainInfo.cpp \
    cpm_VoxelInfo.cpp \
    cpm_VoxelInfoCART.cpp \
    cpm_BaseParaManager.cpp \
    cpm_BaseParaManager_Alloc.cpp \
    cpm_BaseParaManager_MPI.cpp \
    cpm_ParaManager.cpp \
    cpm_ParaManager_Alloc.cpp \
    cpm_ParaManager_MPI.cpp \
    cpm_ParaManager_frtIF.cpp \
    cpm_TextParser.cpp \
    cpm_TextParserDomain.cpp 

CPM_includedir = $(includedir)
CPM_include_HEADERS = \
    $(top_srcdir)/include/cpm_Base.h \
    $(top_srcdir)/include/cpm_Define.h \
    $(top_srcdir)/include/cpm_DomainInfo.h \
    $(top_srcdir)/include/cpm_ObjList.h \
    $(top_srcdir)/include/cpm_BaseParaManager.h \
    $(top_srcdir)/include/cpm_VoxelInfo.h \
    $(top_srcdir)/include/cpm_TextParser.h \
    $(top_srcdir)/include/cpm_TextParserDomain.h \
    $(top_srcdir)/include/cpm_Version.h \
    $(top_srcdir)/include/cpm_fparam.fi \
    $(top_srcdir)/include/cpm_EndianUtil.h \
    $(top_srcdir)/include/cpm_PathUtil.h \
    $(top_srcdir)/include/cpm_ParaManager.h \
    $(top_srcdir)/include/cpm_VoxelInfoCART.h

CPM_inline_includedir = $(includedir)/inline

CPM_inline_include_HEADERS = \
    $(top_srcdir)/include/inline/cpm_BaseParaManager_inline.h \
    $(top_srcdir)/include/inline/cpm_ParaManager_BndComm.h \
    $(top_srcdir)/include/inline/cpm_ParaManager_BndCommEx.h 
