/*
 * CPMlib - Cartesian Partition Manager Library
 *
 * Copyright (C) 2012-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2014-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   cpm_ParaManager_Alloc.cpp
 * LMRパラレルマネージャクラスのソースファイル
 * @date   2012/05/31
 */
#include <stdlib.h>
#include "cpm_ParaManager.h"

////////////////////////////////////////////////////////////////////////////////
// 配列確保
double*
cpm_ParaManager::AllocDouble( int nmax, int sz[3], int vc, int procGrpNo )
{
  size_t nw = size_t(sz[0]+2*vc) * size_t(sz[1]+2*vc) * size_t(sz[2]+2*vc) * size_t(nmax);
  if( nw == 0 ) return NULL;
  return new double[nw];
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
float*
cpm_ParaManager::AllocFloat( int nmax, int sz[3], int vc, int procGrpNo )
{
  size_t nw = size_t(sz[0]+2*vc) * size_t(sz[1]+2*vc) * size_t(sz[2]+2*vc) * size_t(nmax);
  if( nw == 0 ) return NULL;
  return new float[nw];
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
int*
cpm_ParaManager::AllocInt( int nmax, int sz[3], int vc, int procGrpNo )
{
  size_t nw = size_t(sz[0]+2*vc) * size_t(sz[1]+2*vc) * size_t(sz[2]+2*vc) * size_t(nmax);
  if( nw == 0 ) return NULL;
  return new int[nw];
}

