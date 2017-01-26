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
 * @file   cpm_BaseParaManager_Alloc.cpp
 * パラレルマネージャクラスのソースファイル
 * @date   2012/05/31
 */
#include <stdlib.h>
#include "cpm_BaseParaManager.h"

////////////////////////////////////////////////////////////////////////////////
// 配列確保
double*
cpm_BaseParaManager::AllocDoubleS4D( int nmax, int vc, bool padding, int *pad_size, int procGrpNo )
{
  const int *sz = GetLocalArraySize( procGrpNo );
  if( !sz ) return NULL;

  int sz2[3] = {sz[0], sz[1], sz[2]};

  int psize[4];
  if( !pad_size ) pad_size = psize;
  for( int i=0;i<4;i++ ) pad_size[i] = 0;
  if( padding )
  {
    GetPaddingSize( CPM_ARRAY_S4D, sz2, vc, pad_size, nmax );
  }
  sz2[0] += pad_size[0];
  sz2[1] += pad_size[1];
  sz2[2] += pad_size[2];
  nmax   += pad_size[3];

  return AllocDouble( nmax, sz2, vc, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
float*
cpm_BaseParaManager::AllocFloatS4D( int nmax, int vc, bool padding, int *pad_size, int procGrpNo )
{
  const int *sz = GetLocalArraySize( procGrpNo );
  if( !sz ) return NULL;

  int sz2[3] = {sz[0], sz[1], sz[2]};

  int psize[4];
  if( !pad_size ) pad_size = psize;
  for( int i=0;i<4;i++ ) pad_size[i] = 0;
  if( padding )
  {
    GetPaddingSize( CPM_ARRAY_S4D, sz2, vc, pad_size, nmax );
  }
  sz2[0] += pad_size[0];
  sz2[1] += pad_size[1];
  sz2[2] += pad_size[2];
  nmax   += pad_size[3];

  return AllocFloat( nmax, sz2, vc, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
int*
cpm_BaseParaManager::AllocIntS4D( int nmax, int vc, bool padding, int *pad_size, int procGrpNo )
{
  const int *sz = GetLocalArraySize( procGrpNo );
  if( !sz ) return NULL;

  int sz2[3] = {sz[0], sz[1], sz[2]};

  int psize[4];
  if( !pad_size ) pad_size = psize;
  for( int i=0;i<4;i++ ) pad_size[i] = 0;
  if( padding )
  {
    GetPaddingSize( CPM_ARRAY_S4D, sz2, vc, pad_size, nmax );
  }
  sz2[0] += pad_size[0];
  sz2[1] += pad_size[1];
  sz2[2] += pad_size[2];
  nmax   += pad_size[3];

  return AllocInt( nmax, sz2, vc, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
double*
cpm_BaseParaManager::AllocDoubleS3D( int vc, bool padding, int *pad_size, int procGrpNo )
{
  int psz[4];
  double *ptr = AllocDoubleS4D(1, vc, padding, psz, procGrpNo);
  if( pad_size )
  {
    pad_size[0] = psz[0];
    pad_size[1] = psz[1];
    pad_size[2] = psz[2];
  }
  return ptr;
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
float*
cpm_BaseParaManager::AllocFloatS3D( int vc, bool padding, int *pad_size, int procGrpNo )
{
  int psz[4];
  float *ptr = AllocFloatS4D(1, vc, padding, psz, procGrpNo);
  if( pad_size )
  {
    pad_size[0] = psz[0];
    pad_size[1] = psz[1];
    pad_size[2] = psz[2];
  }
  return ptr;
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
int*
cpm_BaseParaManager::AllocIntS3D( int vc, bool padding, int *pad_size, int procGrpNo )
{
  int psz[4];
  int *ptr = AllocIntS4D(1, vc, padding, psz, procGrpNo);
  if( pad_size )
  {
    pad_size[0] = psz[0];
    pad_size[1] = psz[1];
    pad_size[2] = psz[2];
  }
  return ptr;
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
double*
cpm_BaseParaManager::AllocDoubleV3D( int vc, bool padding, int *pad_size, int procGrpNo )
{
  return AllocDoubleS4D(3, vc, padding, pad_size, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
float*
cpm_BaseParaManager::AllocFloatV3D( int vc, bool padding, int *pad_size, int procGrpNo )
{
  return AllocFloatS4D(3, vc, padding, pad_size, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
int*
cpm_BaseParaManager::AllocIntV3D( int vc, bool padding, int *pad_size, int procGrpNo )
{
  return AllocIntS4D(3, vc, padding, pad_size, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
double*
cpm_BaseParaManager::AllocDoubleS4DEx( int nmax, int vc, bool padding, int *pad_size, int procGrpNo )
{
  const int *sz = GetLocalArraySize( procGrpNo );
  if( !sz ) return NULL;

  int sz2[3] = {sz[0], sz[1], sz[2]};

  int psize[4];
  if( !pad_size ) pad_size = psize;
  for( int i=0;i<4;i++ ) pad_size[i] = 0;
  if( padding )
  {
    GetPaddingSize( CPM_ARRAY_S4DEX, sz2, vc, pad_size, nmax );
  }
  nmax   += pad_size[0];
  sz2[0] += pad_size[1];
  sz2[1] += pad_size[2];
  sz2[2] += pad_size[3];

  return AllocDouble( nmax, sz2, vc, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
float*
cpm_BaseParaManager::AllocFloatS4DEx( int nmax, int vc, bool padding, int *pad_size, int procGrpNo )
{
  const int *sz = GetLocalArraySize( procGrpNo );
  if( !sz ) return NULL;

  int sz2[3] = {sz[0], sz[1], sz[2]};

  int psize[4];
  if( !pad_size ) pad_size = psize;
  for( int i=0;i<4;i++ ) pad_size[i] = 0;
  if( padding )
  {
    GetPaddingSize( CPM_ARRAY_S4DEX, sz2, vc, pad_size, nmax );
  }
  nmax   += pad_size[0];
  sz2[0] += pad_size[1];
  sz2[1] += pad_size[2];
  sz2[2] += pad_size[3];

  return AllocFloat( nmax, sz2, vc, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
int*
cpm_BaseParaManager::AllocIntS4DEx( int nmax, int vc, bool padding, int *pad_size, int procGrpNo )
{
  const int *sz = GetLocalArraySize( procGrpNo );
  if( !sz ) return NULL;

  int sz2[3] = {sz[0], sz[1], sz[2]};

  int psize[4];
  if( !pad_size ) pad_size = psize;
  for( int i=0;i<4;i++ ) pad_size[i] = 0;
  if( padding )
  {
    GetPaddingSize( CPM_ARRAY_S4DEX, sz2, vc, pad_size, nmax );
  }
  nmax   += pad_size[0];
  sz2[0] += pad_size[1];
  sz2[1] += pad_size[2];
  sz2[2] += pad_size[3];

  return AllocInt( nmax, sz2, vc, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
double*
cpm_BaseParaManager::AllocDoubleV3DEx( int vc, bool padding, int *pad_size, int procGrpNo )
{
  return AllocDoubleS4DEx( 3, vc, padding, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
float*
cpm_BaseParaManager::AllocFloatV3DEx( int vc, bool padding, int *pad_size, int procGrpNo )
{
  return AllocFloatS4DEx( 3, vc, padding, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
int*
cpm_BaseParaManager::AllocIntV3DEx( int vc, bool padding, int *pad_size, int procGrpNo )
{
  return AllocIntS4DEx( 3, vc, padding, pad_size, procGrpNo );
}

