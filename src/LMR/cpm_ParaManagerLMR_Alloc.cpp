/*
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
 */

/**
 * @file   cpm_ParaManagerLMR_Alloc.cpp
 * LMRパラレルマネージャクラスのソースファイル
 * @date   2012/05/31
 */
#include <stdlib.h>
#include "cpm_ParaManagerLMR.h"

////////////////////////////////////////////////////////////////////////////////
// 配列確保
double*
cpm_ParaManagerLMR::AllocDouble( int nmax, int sz[3], int vc, int procGrpNo )
{
  size_t nw = size_t(sz[0]+2*vc) * size_t(sz[1]+2*vc) * size_t(sz[2]+2*vc) * size_t(nmax);
  nw *= size_t(GetLocalNumLeaf(procGrpNo));
  if( nw == 0 ) return NULL;
  return new double[nw];
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
float*
cpm_ParaManagerLMR::AllocFloat( int nmax, int sz[3], int vc, int procGrpNo )
{
  size_t nw = size_t(sz[0]+2*vc) * size_t(sz[1]+2*vc) * size_t(sz[2]+2*vc) * size_t(nmax);
  nw *= size_t(GetLocalNumLeaf(procGrpNo));
  if( nw == 0 ) return NULL;
  return new float[nw];
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
int*
cpm_ParaManagerLMR::AllocInt( int nmax, int sz[3], int vc, int procGrpNo )
{
  size_t nw = size_t(sz[0]+2*vc) * size_t(sz[1]+2*vc) * size_t(sz[2]+2*vc) * size_t(nmax);
  nw *= size_t(GetLocalNumLeaf(procGrpNo));
  if( nw == 0 ) return NULL;
  return new int[nw];
}
