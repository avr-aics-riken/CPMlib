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

libCPM_a_CXXFLAGS = @MPICH_CFLAGS@ @TP_CFLAGS@ @PM_CFLAGS@ @REALOPT@ -I$(top_srcdir)/include -I$(top_srcdir)/include/BCM

libCPM_a_SOURCES = \
    cpm_DomainInfo.cpp \
    cpm_VoxelInfo.cpp \
    cpm_VoxelInfoCART.cpp \
    cpm_VoxelInfoLMR.cpp \
    cpm_ParaManager.cpp \
    cpm_ParaManagerCART.cpp \
    cpm_ParaManagerLMR.cpp \
    cpm_ParaManager_MPI.cpp \
    cpm_ParaManager_Alloc.cpp \
    cpm_ParaManager_frtIF.cpp \
    cpm_TextParser.cpp \
    cpm_TextParserDomain.cpp \
    cpm_TextParserDomainLMR.cpp \
    BCM/BCMOctree.cpp

CPM_includedir = $(includedir)
CPM_include_HEADERS = \
    $(top_srcdir)/include/cpm_Base.h \
    $(top_srcdir)/include/cpm_Define.h \
    $(top_srcdir)/include/cpm_DomainInfo.h \
    $(top_srcdir)/include/cpm_ObjList.h \
    $(top_srcdir)/include/cpm_ParaManager.h \
    $(top_srcdir)/include/cpm_ParaManagerCART.h \
    $(top_srcdir)/include/cpm_ParaManagerLMR.h \
    $(top_srcdir)/include/cpm_TextParser.h \
    $(top_srcdir)/include/cpm_TextParserDomain.h \
    $(top_srcdir)/include/cpm_TextParserDomainLMR.h \
    $(top_srcdir)/include/cpm_Version.h \
    $(top_srcdir)/include/cpm_VoxelInfo.h \
    $(top_srcdir)/include/cpm_VoxelInfoCART.h \
    $(top_srcdir)/include/cpm_VoxelInfoLMR.h \
    $(top_srcdir)/include/cpm_fparam.fi \
    $(top_srcdir)/include/cpm_EndianUtil.h \
    $(top_srcdir)/include/cpm_PathUtil.h \
    $(top_srcdir)/include/BCM/BCMFileCommon.h \
    $(top_srcdir)/include/BCM/BCMOctree.h \
    $(top_srcdir)/include/BCM/BCMTools.h \
    $(top_srcdir)/include/BCM/BitVoxel.h \
    $(top_srcdir)/include/BCM/Divider.h \
    $(top_srcdir)/include/BCM/NeighborInfo.h \
    $(top_srcdir)/include/BCM/Node.h \
    $(top_srcdir)/include/BCM/Partition.h \
    $(top_srcdir)/include/BCM/Pedigree.h \
    $(top_srcdir)/include/BCM/RootGrid.h \
    $(top_srcdir)/include/BCM/Vec3.h

CPM_inline_includedir = $(includedir)/inline

CPM_inline_include_HEADERS = \
    $(top_srcdir)/include/inline/cpm_ParaManager_inline.h \
    $(top_srcdir)/include/inline/cpm_ParaManager_BndComm.h \
    $(top_srcdir)/include/inline/cpm_ParaManager_BndCommEx.h \
    $(top_srcdir)/include/inline/cpm_ParaManager_BndComm_CART.h \
    $(top_srcdir)/include/inline/cpm_ParaManager_BndCommEx_CART.h \
    $(top_srcdir)/include/inline/cpm_ParaManager_BndComm_LMR.h \
    $(top_srcdir)/include/inline/cpm_ParaManager_BndCommEx_LMR.h
