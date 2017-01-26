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

lib_LIBRARIES = libCPMLMR.a

libCPMLMR_a_CXXFLAGS = @MPICH_CFLAGS@ @TP_CFLAGS@ @PM_CFLAGS@ @REALOPT@ \
                    -I$(top_srcdir)/include \
                    -I$(top_srcdir)/include/LMR \
                    -I$(top_srcdir)/include/LMR/BCM

libCPMLMR_a_SOURCES = \
    cpm_VoxelInfoLMR.cpp \
    cpm_ParaManagerLMR.cpp \
    cpm_ParaManagerLMR_MPI.cpp \
    cpm_ParaManagerLMR_frtIF.cpp \
    cpm_ParaManagerLMR_Alloc.cpp \
    cpm_LeafCommInfo.cpp \
    cpm_TextParserDomainLMR.cpp \
    BCM/BCMOctree.cpp

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
    $(top_srcdir)/include/LMR/cpm_DefineLMR.h \
    $(top_srcdir)/include/LMR/cpm_ParaManagerLMR.h \
    $(top_srcdir)/include/LMR/cpm_LeafCommInfo.h \
    $(top_srcdir)/include/LMR/cpm_VoxelInfoLMR.h \
    $(top_srcdir)/include/LMR/cpm_TextParserDomainLMR.h \
    $(top_srcdir)/include/LMR/BCM/BCMFileCommon.h \
    $(top_srcdir)/include/LMR/BCM/BCMOctree.h \
    $(top_srcdir)/include/LMR/BCM/BCMTools.h \
    $(top_srcdir)/include/LMR/BCM/BitVoxel.h \
    $(top_srcdir)/include/LMR/BCM/Divider.h \
    $(top_srcdir)/include/LMR/BCM/NeighborInfo.h \
    $(top_srcdir)/include/LMR/BCM/Node.h \
    $(top_srcdir)/include/LMR/BCM/Partition.h \
    $(top_srcdir)/include/LMR/BCM/Pedigree.h \
    $(top_srcdir)/include/LMR/BCM/RootGrid.h \
    $(top_srcdir)/include/LMR/BCM/Vec3.h

CPM_inline_includedir = $(includedir)/inline

CPM_inline_include_HEADERS = \
    $(top_srcdir)/include/inline/cpm_BaseParaManager_inline.h \
    $(top_srcdir)/include/LMR/inline/cpm_ParaManagerLMR_BndComm.h \
    $(top_srcdir)/include/LMR/inline/cpm_ParaManagerLMR_BndCommEx.h
