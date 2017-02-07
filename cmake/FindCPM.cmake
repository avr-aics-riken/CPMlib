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

# - Try to find CPMlib
# Once done, this will define
#
#  CPM_FOUND - system has CPMlib
#  CPM_INCLUDE_DIRS - CPMlib include directories
#  CPM_LIBRARIES - link these to use CPMlib

include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(CPM_PKGCONF CPM)

if(CMAKE_PREFIX_PATH)
  set(CPM_CANDIDATE_PATH ${CMAKE_PREFIX_PATH})
  file(GLOB tmp "${CMAKE_PREFIX_PATH}/[Jj][Hh][Pp][Cc][Nn][Dd][Ff]*/")
  list(APPEND CPM_CANDIDATE_PATH ${tmp})
endif()

# Include dir
find_path(CPM_INCLUDE_DIR
  NAMES cpm_Base.h
  PATHS ${CPM_ROOT} ${CPM_PKGCONF_INCLUDE_DIRS} ${CPM_CANDIDATE_PATH}
  PATH_SUFFIXES include
)

# Finally the library itself
find_library(CPM_LIBRARY
  NAMES CPM
  PATHS ${CPM_ROOT} ${CPM_PKGCONF_LIBRARY_DIRS} ${CPM_CANDIDATE_PATH}
  PATH_SUFFIXES lib
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(CPM_PROCESS_INCLUDES CPM_INCLUDE_DIR)
set(CPM_PROCESS_LIBS CPM_LIBRARY)
libfind_process(CPM)
