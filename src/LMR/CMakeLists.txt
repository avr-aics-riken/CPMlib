###################################################################################
#
# CPMlib - Computational space Partitioning Management library
#
# Copyright (c) 2012-2014 Institute of Industrial Science (IIS), The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2014-2016 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2016-2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################

# cpm_Version.hはPROJECT_BINARY_DIR
include_directories(
       ${PROJECT_BINARY_DIR}/include
       ${PROJECT_SOURCE_DIR}/include
       ${PROJECT_SOURCE_DIR}/include/LMR
       ${PROJECT_SOURCE_DIR}/include/inline
       ${PROJECT_SOURCE_DIR}/include/LMR/BCM
       ${TP_INC}
)


link_directories(
      ${TP_LIB}
)

set(cpm_files
    cpm_LeafCommInfo.cpp
    cpm_ParaManagerLMR_Alloc.cpp
    cpm_ParaManagerLMR_frtIF.cpp
    cpm_ParaManagerLMR_MPI.cpp
    cpm_ParaManagerLMR.cpp
    cpm_TextParserDomainLMR.cpp
    cpm_VoxelInfoLMR.cpp
    BCM/BCMOctree.cpp
)


if(NOT with_MPI)
  add_library(CPM_LMR STATIC ${cpm_files})
  target_link_libraries(CPM_LMR -lTP -lPM)
  install(TARGETS CPM_LMR DESTINATION lib)

else()
  add_library(CPM_LMRmpi STATIC ${cpm_files})
  target_link_libraries(CPM_LMRmpi -lTPmpi -lPMmpi)
  install(TARGETS CPM_LMRmpi DESTINATION lib)

endif()

install(FILES
        ${PROJECT_SOURCE_DIR}/include/LMR/cpm_DefineLMR.h
        ${PROJECT_SOURCE_DIR}/include/LMR/cpm_ParaManagerLMR.h
        ${PROJECT_SOURCE_DIR}/include/LMR/cpm_LeafCommInfo.h
        ${PROJECT_SOURCE_DIR}/include/LMR/cpm_VoxelInfoLMR.h
        ${PROJECT_SOURCE_DIR}/include/LMR/cpm_TextParserDomainLMR.h
        DESTINATION include/LMR
)

install(FILES
        ${PROJECT_SOURCE_DIR}/include/LMR/BCM/BCMFileCommon.h
        ${PROJECT_SOURCE_DIR}/include/LMR/BCM/BCMOctree.h
        ${PROJECT_SOURCE_DIR}/include/LMR/BCM/BCMTools.h
        ${PROJECT_SOURCE_DIR}/include/LMR/BCM/BitVoxel.h
        ${PROJECT_SOURCE_DIR}/include/LMR/BCM/Divider.h
        ${PROJECT_SOURCE_DIR}/include/LMR/BCM/NeighborInfo.h
        ${PROJECT_SOURCE_DIR}/include/LMR/BCM/Node.h
        ${PROJECT_SOURCE_DIR}/include/LMR/BCM/Partition.h
        ${PROJECT_SOURCE_DIR}/include/LMR/BCM/Pedigree.h
        ${PROJECT_SOURCE_DIR}/include/LMR/BCM/RootGrid.h
        ${PROJECT_SOURCE_DIR}/include/LMR/BCM/Vec3.h
        DESTINATION include/LMR/BCM
)

install(FILES
        ${PROJECT_SOURCE_DIR}/include/LMR/inline/cpm_ParaManagerLMR_BndComm.h
        ${PROJECT_SOURCE_DIR}/include/LMR/inline/cpm_ParaManagerLMR_BndCommEx.h
        DESTINATION include/LMR/inline
)
