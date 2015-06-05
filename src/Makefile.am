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

lib_LIBRARIES = libCPM.a

libCPM_a_CXXFLAGS = @MPICH_CFLAGS@ @TP_CFLAGS@ @PM_CFLAGS@ @REALOPT@ -I$(top_builddir)/include -I$(top_builddir)/include/BCM

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
    $(top_builddir)/include/cpm_Base.h \
    $(top_builddir)/include/cpm_Define.h \
    $(top_builddir)/include/cpm_DomainInfo.h \
    $(top_builddir)/include/cpm_ObjList.h \
    $(top_builddir)/include/cpm_ParaManager.h \
    $(top_builddir)/include/cpm_ParaManagerCART.h \
    $(top_builddir)/include/cpm_ParaManagerLMR.h \
    $(top_builddir)/include/cpm_TextParser.h \
    $(top_builddir)/include/cpm_TextParserDomain.h \
    $(top_builddir)/include/cpm_TextParserDomainLMR.h \
    $(top_builddir)/include/cpm_Version.h \
    $(top_builddir)/include/cpm_VoxelInfo.h \
    $(top_builddir)/include/cpm_VoxelInfoCART.h \
    $(top_builddir)/include/cpm_VoxelInfoLMR.h \
    $(top_builddir)/include/cpm_fparam.fi \
    $(top_builddir)/include/cpm_EndianUtil.h \
    $(top_builddir)/include/cpm_PathUtil.h \
    $(top_builddir)/include/BCM/BCMFileCommon.h \
    $(top_builddir)/include/BCM/BCMOctree.h \
    $(top_builddir)/include/BCM/BCMTools.h \
    $(top_builddir)/include/BCM/BitVoxel.h \
    $(top_builddir)/include/BCM/Divider.h \
    $(top_builddir)/include/BCM/NeighborInfo.h \
    $(top_builddir)/include/BCM/Node.h \
    $(top_builddir)/include/BCM/Partition.h \
    $(top_builddir)/include/BCM/Pedigree.h \
    $(top_builddir)/include/BCM/RootGrid.h \
    $(top_builddir)/include/BCM/Vec3.h

CPM_inline_includedir = $(includedir)/inline

CPM_inline_include_HEADERS = \
    $(top_builddir)/include/inline/cpm_ParaManager_inline.h \
    $(top_builddir)/include/inline/cpm_ParaManager_BndComm.h \
    $(top_builddir)/include/inline/cpm_ParaManager_BndCommEx.h \
    $(top_builddir)/include/inline/cpm_ParaManager_BndComm_CART.h \
    $(top_builddir)/include/inline/cpm_ParaManager_BndCommEx_CART.h \
    $(top_builddir)/include/inline/cpm_ParaManager_BndComm_LMR.h \
    $(top_builddir)/include/inline/cpm_ParaManager_BndCommEx_LMR.h
