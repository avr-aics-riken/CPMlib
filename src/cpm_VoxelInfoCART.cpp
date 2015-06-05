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
 * @file   cpm_VoxelInfoCART.cpp
 * カーテシアン用のVOXEL空間情報クラスのソースファイル
 * @author University of Tokyo
 * @date   2015/03/27
 */
#include "cpm_VoxelInfoCART.h"

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_VoxelInfoCART::cpm_VoxelInfoCART()
  : cpm_VoxelInfo()
{
  m_rankMap = NULL;
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_VoxelInfoCART::~cpm_VoxelInfoCART()
{
  if( m_rankMap ) delete [] m_rankMap;
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割情報の生成
cpm_ErrorCode
cpm_VoxelInfoCART::Init( MPI_Comm comm, cpm_GlobalDomainInfo* dInfo )
{
  // 入力チェック
  if( IsCommNull(comm) )
  {
    return CPM_ERROR_MPI_INVALID_COMM;
  }
  if( !dInfo )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 入力をコピー
  m_comm = comm;
  m_globalDomainInfo = *dInfo;
   
  // ランク数、ランク番号をセット
  MPI_Comm_size(m_comm, &m_nRank);
  MPI_Comm_rank(m_comm, &m_rankNo);

  // ランクマップを生成
  if( !CreateRankMap() )
  {
    return CPM_ERROR_CREATE_RANKMAP;
  }

  // ローカル領域情報を生成
  if( !CreateLocalDomainInfo() )
  {
    return CPM_ERROR_CREATE_LOCALDOMAIN;
  }

  // 隣接ランク情報を生成
  if( !CreateNeighborRankInfo() )
  {
    return CPM_ERROR_CREATE_NEIGHBOR;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// ランクマップを生成
bool
cpm_VoxelInfoCART::CreateRankMap()
{
  // 領域分割数を取得
  const int* div = m_globalDomainInfo.GetDivNum();
  if( !div )
  {
    return false;
  }

  // マップ領域を確保(初期値NULL)
  size_t ndiv = size_t(div[0]) * size_t(div[1]) * size_t(div[2]);
  int *rankMap = new int[ndiv];
  if( !rankMap )
  {
    return false;
  }
  for( size_t i=0;i<ndiv;i++ ) rankMap[i] = getRankNull();

  // 活性サブドメイン情報配置位置に0をセット
  for( int i=0;i<m_globalDomainInfo.GetSubdomainNum();i++ )
  {
    //サブドメイン情報
    const cpm_ActiveSubdomainInfo* dom = m_globalDomainInfo.GetSubdomainInfo(i);
    if( !dom )
    {
      delete [] rankMap;
      return false;
    }

    // 位置を取得
    const int *pos = dom->GetPos();
    if( !pos )
    {
      delete [] rankMap;
      return false;
    }

    // 0をセット
    rankMap[_IDX_S3D(pos[0],pos[1],pos[2],div[0],div[1],div[2],0)] = 0;
  }

  // i->j->kの優先順で活性サブドメインにランク番号をセット
  int rankCount = 0;
  for( int k=0;k<div[2];k++ ){
  for( int j=0;j<div[1];j++ ){
  for( int i=0;i<div[0];i++ ){
    if( rankMap[_IDX_S3D(i,j,k,div[0],div[1],div[2],0)] == 0 )
    {
      rankMap[_IDX_S3D(i,j,k,div[0],div[1],div[2],0)] = rankCount;
      rankCount++;
    }
  }}}

  // ランクマップをセット
  if( m_rankMap ) delete [] m_rankMap;
  m_rankMap = rankMap;
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// 隣接ランク情報を生成
bool
cpm_VoxelInfoCART::CreateNeighborRankInfo()
{
  //自ランクの位置を取得
  const int *pos = m_localDomainInfo.GetPos();
  if( !pos )
  {
    return false;
  }

  //領域分割数を取得
  const int *div = m_globalDomainInfo.GetDivNum();
  if( !div )
  {
     return false;
  }

  //ランクマップと整合性が取れているかをチェック
  if( !m_rankMap )
  {
    return false;
  }
  if( m_rankMap[_IDX_S3D(pos[0],pos[1],pos[2],div[0],div[1],div[2],0)] != m_rankNo )
  {
    return false;
  }

  //隣接ランクを取得

  // -X face
  if( pos[0] != 0 )
  {
    m_neighborRankID[X_MINUS] = m_rankMap[_IDX_S3D(pos[0]-1,pos[1],pos[2],div[0],div[1],div[2],0)];
    m_periodicRankID[X_MINUS] = getRankNull();
  }
  else
  {
    m_neighborRankID[X_MINUS] = getRankNull();
    m_periodicRankID[X_MINUS] = m_rankMap[_IDX_S3D(div[0]-1,pos[1],pos[2],div[0],div[1],div[2],0)];
  }

  // -Y face
  if( pos[1] != 0 )
  {
    m_neighborRankID[Y_MINUS] = m_rankMap[_IDX_S3D(pos[0],pos[1]-1,pos[2],div[0],div[1],div[2],0)];
    m_periodicRankID[Y_MINUS] = getRankNull();
  }
  else
  {
    m_neighborRankID[Y_MINUS] = getRankNull();
    m_periodicRankID[Y_MINUS] = m_rankMap[_IDX_S3D(pos[0],div[1]-1,pos[2],div[0],div[1],div[2],0)];
  }

  // -Z face
  if( pos[2] != 0 )
  {
    m_neighborRankID[Z_MINUS] = m_rankMap[_IDX_S3D(pos[0],pos[1],pos[2]-1,div[0],div[1],div[2],0)];
    m_periodicRankID[Z_MINUS] = getRankNull();
  }
  else
  {
    m_neighborRankID[Z_MINUS] = getRankNull();
    m_periodicRankID[Z_MINUS] = m_rankMap[_IDX_S3D(pos[0],pos[1],div[2]-1,div[0],div[1],div[2],0)];
  }

  // +X face
  if( pos[0] != div[0]-1 )
  {
    m_neighborRankID[X_PLUS] = m_rankMap[_IDX_S3D(pos[0]+1,pos[1],pos[2],div[0],div[1],div[2],0)];
    m_periodicRankID[X_PLUS] = getRankNull();
  }
  else
  {
    m_neighborRankID[X_PLUS] = getRankNull();
    m_periodicRankID[X_PLUS] = m_rankMap[_IDX_S3D(0       ,pos[1],pos[2],div[0],div[1],div[2],0)];
  }

  // +Y face
  if( pos[1] != div[1]-1 )
  {
    m_neighborRankID[Y_PLUS] = m_rankMap[_IDX_S3D(pos[0],pos[1]+1,pos[2],div[0],div[1],div[2],0)];
    m_periodicRankID[Y_PLUS] = getRankNull();
  }
  else
  {
    m_neighborRankID[Y_PLUS] = getRankNull();
    m_periodicRankID[Y_PLUS] = m_rankMap[_IDX_S3D(pos[0],0       ,pos[2],div[0],div[1],div[2],0)];
  }

  // +Z face
  if( pos[2] != div[2]-1 )
  {
    m_neighborRankID[Z_PLUS] = m_rankMap[_IDX_S3D(pos[0],pos[1],pos[2]+1,div[0],div[1],div[2],0)];
    m_periodicRankID[Z_PLUS] = getRankNull();
  }
  else
  {
    m_neighborRankID[Z_PLUS] = getRankNull();
    m_periodicRankID[Z_PLUS] = m_rankMap[_IDX_S3D(pos[0],pos[1],0       ,div[0],div[1],div[2],0)];
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////
// ローカル領域情報を生成
bool
cpm_VoxelInfoCART::CreateLocalDomainInfo()
{
  if( !m_rankMap )
  {
    return false;
  }

  // 領域分割数
  const int *div = m_globalDomainInfo.GetDivNum();
  if( !div )
  {
    return false;
  }

  // 全体ボクセル数
  const int *gvox = m_globalDomainInfo.GetVoxNum();
  if( !gvox )
  {
    return false;
  }

  // ランクマップから、自ランクの位置を取得、セット
  int pos[3];
  for( int k=0;k<div[2];k++ ){
  for( int j=0;j<div[1];j++ ){
  for( int i=0;i<div[0];i++ ){
    int rankNo = m_rankMap[_IDX_S3D(i,j,k,div[0],div[1],div[2],0)];
    if( rankNo == m_rankNo )
    {
      pos[0] = i;
      pos[1] = j;
      pos[2] = k;
      break;
    }
  }}}
  m_localDomainInfo.SetPos(pos);

  // ローカルのVOXEL数
  int *nvX = new int[div[0]];
  int *nvY = new int[div[1]];
  int *nvZ = new int[div[2]];
  int *nv[3] = {nvX,nvY,nvZ};
  for( int n=0;n<3;n++ )
  {
    int *nvd = nv[n];
    //基準のボクセル数
    int nbase = gvox[n] / div[n];

    //余り
    int amari = gvox[n] % div[n];

    // ボクセル数をセット
    for( int i=0;i<div[n];i++ )
    {
      nvd[i] = nbase;
      if( i<amari ) nvd[i]++;
    }
  }

  // VOXEL数と始点インデクスをセット
  int lvox[3] = {0,0,0};
  int head[3] = {0,0,0};
  for( int n=0;n<3;n++ )
  {
    int *nvd = nv[n];
    int hd = 0;

    for( int i=0;i<div[n];i++ )
    {
      if( i == pos[n] )
      {
        lvox[n] = nvd[i];
        head[n] = hd;
        break;
      }
      hd += nvd[i];
    }
  }
  m_localDomainInfo.SetVoxNum(lvox);
  for( int n=0;n<3;n++ )
  {
    m_voxelHeadIndex[n] = head[n];
    m_voxelTailIndex[n] = head[n] + lvox[n] - 1;
  }

  // ローカルの原点、ピッチ、空間サイズをセット
  const double *gorg = m_globalDomainInfo.GetOrigin();
  const double *gpch = m_globalDomainInfo.GetPitch();
  double org[3], pch[3], rgn[3];
  for( int n=0;n<3;n++ )
  {
    pch[n] = gpch[n];
    org[n] = gorg[n] + double(head[n]) * pch[n];
    rgn[n] = double(lvox[n]) * pch[n];
  }
  m_localDomainInfo.SetOrigin(org);
  m_localDomainInfo.SetPitch (pch);
  m_localDomainInfo.SetRegion(rgn);

  delete [] nvX;
  delete [] nvY;
  delete [] nvZ;
  return true;
}

