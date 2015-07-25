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

////////////////////////////////////////////////////////////////////////////////
// 自ランクの終点VOXELの全体空間でのインデクスを取得
const int*
cpm_VoxelInfo::GetVoxelTailIndex() const
{
  return m_voxelTailIndex;
}

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
// 指定面における自ランクの隣接ランク番号を取得
const int*
cpm_VoxelInfo::GetNeighborRankList( cpm_FaceFlag face, int &num ) const
{
  num = 1;
  return &(m_neighborRankID[face]);
}

////////////////////////////////////////////////////////////////////////////////
// 指定面における自ランクの周期境界の隣接ランク番号を取得
const int*
cpm_VoxelInfo::GetPeriodicRankList( cpm_FaceFlag face, int &num ) const
{
  num = 1;
  return &(m_periodicRankID[face]);
}

////////////////////////////////////////////////////////////////////////////////
// 指定面におけるレベル差を取得
int
cpm_VoxelInfo::GetNeighborLevelDiff( cpm_FaceFlag face ) const
{
  return 0;
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

