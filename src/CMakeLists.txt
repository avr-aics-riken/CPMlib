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
       ${PROJECT_SOURCE_DIR}/include/inline
       ${TP_INC}
)


link_directories(
      ${TP_LIB}
)


set(cpm_files
    cpm_BaseParaManager_Alloc.cpp
    cpm_BaseParaManager_MPI.cpp
    cpm_BaseParaManager.cpp
    cpm_DomainInfo.cpp
    cpm_ParaManager_Alloc.cpp
    cpm_ParaManager_frtIF.cpp
    cpm_ParaManager_MPI.cpp
    cpm_ParaManager.cpp
    cpm_TextParser.cpp
    cpm_TextParserDomain.cpp
    cpm_VoxelInfo.cpp
    cpm_VoxelInfoCART.cpp
)


if(NOT with_MPI)
  add_library(CPM STATIC ${cpm_files})
  #target_link_libraries(CPM -lTP -lPM)
  target_link_libraries(CPM -lTP)
  install(TARGETS CPM DESTINATION lib)

else()
  add_library(CPMmpi STATIC ${cpm_files})
  #target_link_libraries(CPMmpi -lTPmpi -lPMmpi)
  target_link_libraries(CPMmpi -lTPmpi)
  install(TARGETS CPMmpi DESTINATION lib)

endif()

install(FILES
        ${PROJECT_SOURCE_DIR}/include/cpm_Base.h
        ${PROJECT_SOURCE_DIR}/include/cpm_BaseParaManager.h
        ${PROJECT_SOURCE_DIR}/include/cpm_Define.h
        ${PROJECT_SOURCE_DIR}/include/cpm_DomainInfo.h
        ${PROJECT_SOURCE_DIR}/include/cpm_EndianUtil.h
        ${PROJECT_SOURCE_DIR}/include/cpm_fparam.fi
        ${PROJECT_SOURCE_DIR}/include/cpm_ObjList.h
        ${PROJECT_SOURCE_DIR}/include/cpm_ParaManager.h
        ${PROJECT_SOURCE_DIR}/include/cpm_PathUtil.h
        ${PROJECT_SOURCE_DIR}/include/cpm_TextParser.h
        ${PROJECT_SOURCE_DIR}/include/cpm_TextParserDomain.h
        ${PROJECT_SOURCE_DIR}/include/cpm_VoxelInfo.h
        ${PROJECT_SOURCE_DIR}/include/cpm_VoxelInfoCART.h
        ${PROJECT_BINARY_DIR}/include/cpm_Version.h
        DESTINATION include
)

install(FILES
        ${PROJECT_SOURCE_DIR}/include/inline/cpm_ParaManager_BndComm.h
        ${PROJECT_SOURCE_DIR}/include/inline/cpm_ParaManager_BndCommEx.h
        ${PROJECT_SOURCE_DIR}/include/inline/cpm_BaseParaManager_inline.h
        DESTINATION include/inline
)
