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
 * @file   cpm_VoxelInfo.cpp
 * VOXEL空間情報クラスのソースファイル
 * @date   2012/05/31
 */
#include "cpm_VoxelInfo.h"

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_VoxelInfo::cpm_VoxelInfo()
  : cpm_Base()
{
  for( int i=0;i<3;i++ )
  {
    m_voxelHeadIndex[i] = 0;
    m_voxelTailIndex[i] = 0;
  }

// 2016/01/22 FEAST add.s
  for( int i=0;i<3;i++ )
  {
    m_nodeHeadIndex[i] = 0;
    m_nodeTailIndex[i] = 0;
  }
// 2016/01/22 FEAST add.e
 
  m_comm   = MPI_COMM_NULL;
  m_nRank  = 1;
  m_rankNo = 0;
  for( int i=0;i<6;i++ )
  {
    m_neighborRankID[i] = getRankNull();
    m_periodicRankID[i] = getRankNull();
  }
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_VoxelInfo::~cpm_VoxelInfo()
{
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割数を取得
const int*
cpm_VoxelInfo::GetDivNum() const
{
  return m_globalDomainInfo.GetDivNum();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの領域分割位置を取得
const int*
cpm_VoxelInfo::GetDivPos() const
{
  return m_localDomainInfo.GetPos();
}

////////////////////////////////////////////////////////////////////////////////
// ピッチを取得
const double*
cpm_VoxelInfo::GetPitch() const
{
  return m_localDomainInfo.GetPitch();
}

////////////////////////////////////////////////////////////////////////////////
// グローバルピッチを取得
const double*
cpm_VoxelInfo::GetGlobalPitch() const
{
  return m_globalDomainInfo.GetPitch();
}

////////////////////////////////////////////////////////////////////////////////
// 全体ボクセル数を取得
const int*
cpm_VoxelInfo::GetGlobalVoxelSize() const
{
  return m_globalDomainInfo.GetVoxNum();
}

// 2016/01/22 FEAST add.s

////////////////////////////////////////////////////////////////////////////////
// 全体頂点数を取得
const int*
cpm_VoxelInfo::GetGlobalNodeSize() const
{
  return m_globalDomainInfo.GetNodNum();
}

////////////////////////////////////////////////////////////////////////////////
// 全体ボクセル数または頂点数を取得
const int*
cpm_VoxelInfo::GetGlobalArraySize( cpm_DefPointType dtype ) const
{
  //定義点タイプがボクセルのとき
  if( dtype == CPM_DEFPOINTTYPE_FVM ) {
    return m_globalDomainInfo.GetVoxNum();
  }
  //定義点タイプが頂点のとき
  if( dtype == CPM_DEFPOINTTYPE_FDM ) {
    return m_globalDomainInfo.GetNodNum();
  }
  return NULL;
}

// 2016/01/22 FEAST add.e

////////////////////////////////////////////////////////////////////////////////
// 全体空間の原点を取得
const double*
cpm_VoxelInfo::GetGlobalOrigin() const
{
  return m_globalDomainInfo.GetOrigin();
}

////////////////////////////////////////////////////////////////////////////////
// 全体空間サイズを取得
const double*
cpm_VoxelInfo::GetGlobalRegion() const
{
  return m_globalDomainInfo.GetRegion();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクのボクセル数を取得
const int*
cpm_VoxelInfo::GetLocalVoxelSize() const
{
  return m_localDomainInfo.GetVoxNum();
}

// 2016/01/22 FEAST add.s

////////////////////////////////////////////////////////////////////////////////
// 自ランクの頂点数を取得
const int*
cpm_VoxelInfo::GetLocalNodeSize() const
{
  return m_localDomainInfo.GetNodNum();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクのボクセル数または頂点数を取得
const int*
cpm_VoxelInfo::GetLocalArraySize( cpm_DefPointType dtype ) const
{
  //定義点タイプがボクセルのとき
  if( dtype == CPM_DEFPOINTTYPE_FVM ) {
    return m_localDomainInfo.GetVoxNum();
  }
  //定義点タイプが頂点のとき
  if( dtype == CPM_DEFPOINTTYPE_FDM ) {
    return m_localDomainInfo.GetNodNum();
  }
  return NULL;
}

// 2016/01/22 FEAST add.e

////////////////////////////////////////////////////////////////////////////////
// 自ランクの空間原点を取得
const double*
cpm_VoxelInfo::GetLocalOrigin() const
{
  return m_localDomainInfo.GetOrigin();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの空間サイズを取得
const double*
cpm_VoxelInfo::GetLocalRegion() const
{
  return m_localDomainInfo.GetRegion();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの始点VOXELの全体空間でのインデクスを取得
const int*
cpm_VoxelInfo::GetVoxelHeadIndex() const
{
  return m_voxelHeadIndex;
}

// 2016/01/22 FEAST add.s

////////////////////////////////////////////////////////////////////////////////
// 自ランクの始点VOXELの全体空間でのインデクスを取得
const int*
cpm_VoxelInfo::GetNodeHeadIndex() const
{
  return m_nodeHeadIndex;
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの始点VOXELまたは頂点の全体空間でのインデクスを取得
const int*
cpm_VoxelInfo::GetArrayHeadIndex( cpm_DefPointType dtype ) const
{
  //定義点タイプがボクセルのとき
  if( dtype == CPM_DEFPOINTTYPE_FVM ) {
    return m_voxelHeadIndex;
  }
  //定義点タイプが頂点のとき
  if( dtype == CPM_DEFPOINTTYPE_FDM ) {
    return m_nodeHeadIndex;
  }
  return NULL;
}

// 2016/01/22 FEAST add.e

////////////////////////////////////////////////////////////////////////////////
// 自ランクの終点VOXELの全体空間でのインデクスを取得
const int*
cpm_VoxelInfo::GetVoxelTailIndex() const
{
  return m_voxelTailIndex;
}

// 2016/01/22 FEAST add.s

////////////////////////////////////////////////////////////////////////////////
// 自ランクの終点頂点の全体空間でのインデクスを取得
const int*
cpm_VoxelInfo::GetNodeTailIndex() const
{
  return m_nodeTailIndex;
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの終点頂点の全体空間でのインデクスを取得
const int*
cpm_VoxelInfo::GetArrayTailIndex( cpm_DefPointType dtype ) const
{
  //定義点タイプがボクセルのとき
  if( dtype == CPM_DEFPOINTTYPE_FVM ) {
    return m_voxelTailIndex;
  }
  //定義点タイプが頂点のとき
  if( dtype == CPM_DEFPOINTTYPE_FDM ) {
    return m_nodeTailIndex;
  }
  return NULL;
}

// 2016/01/22 FEAST add.e

////////////////////////////////////////////////////////////////////////////////
// 自ランクの隣接ランク番号を取得
const int*
cpm_VoxelInfo::GetNeighborRankID() const
{
  return m_neighborRankID;
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの周期境界の隣接ランク番号を取得
const int*
cpm_VoxelInfo::GetPeriodicRankID() const
{
  return m_periodicRankID;
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの境界が外部境界かどうかを判定
bool
cpm_VoxelInfo::IsOuterBoundary( cpm_FaceFlag face ) const
{
  // 隣のドメインが存在する
  if( !IsRankNull(m_neighborRankID[face]) )
  {
    return false;
  }

  // 分割数とポジション
  const int *div = GetDivNum();
  const int *pos = GetDivPos();

  // 分割数とポジションから外部境界かどうか判定
  if( face == X_MINUS && pos[0] == 0        ) return true;
  if( face == X_PLUS  && pos[0] == div[0]-1 ) return true;
  if( face == Y_MINUS && pos[1] == 0        ) return true;
  if( face == Y_PLUS  && pos[1] == div[1]-1 ) return true;
  if( face == Z_MINUS && pos[2] == 0        ) return true;
  if( face == Z_PLUS  && pos[2] == div[2]-1 ) return true;

  // ここに到達した場合は内部境界
  return false; 
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの境界が内部境界(隣が不活性ドメイン)かどうかを判定
bool
cpm_VoxelInfo::IsInnerBoundary( cpm_FaceFlag face ) const
{
  // 隣のドメインが存在する
  if( !IsRankNull(m_neighborRankID[face]) )
  {
    return false;
  }

  // 分割数とポジション
  const int *div = GetDivNum();
  const int *pos = GetDivPos();

  // 分割数とポジションから外部境界かどうか判定
  if( face == X_MINUS && pos[0] == 0        ) return false;
  if( face == X_PLUS  && pos[0] == div[0]-1 ) return false;
  if( face == Y_MINUS && pos[1] == 0        ) return false;
  if( face == Y_PLUS  && pos[1] == div[1]-1 ) return false;
  if( face == Z_MINUS && pos[2] == 0        ) return false;
  if( face == Z_PLUS  && pos[2] == div[2]-1 ) return false;

  // ここに到達した場合は内部境界
  return true;
}

