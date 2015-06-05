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
 * @file   cpm_VoxelInfo.h
 * VOXEL空間情報クラスのヘッダーファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */

#ifndef _CPM_VOXELINFO_H_
#define _CPM_VOXELINFO_H_

#include "cpm_Base.h"
#include "cpm_DomainInfo.h"

/** CPMのVOXEL空間情報管理クラス
 */
class cpm_VoxelInfo : public cpm_Base
{
friend class cpm_ParaManager;
friend class cpm_ParaManagerCART;
////////////////////////////////////////////////////////////////////////////////
// メンバー関数 
////////////////////////////////////////////////////////////////////////////////
public:


protected:
  /** コンストラクタ */
  cpm_VoxelInfo();

  /** デストラクタ */
  virtual ~cpm_VoxelInfo();

  /** 領域分割数を取得
   *  LMRのときは最大レベルにおける分割数を返す
   *  @return 領域分割数整数配列のポインタ
   */
  const int* GetDivNum() const;

  /** 自ランクの領域分割位置を取得
   *  LMRのときは最大レベルにおける分割位置を返す
   *  @return 自ランクの領域分割位置整数配列のポインタ
   */
  const int* GetDivPos() const;

  /** ローカルピッチを取得
   *  @return ピッチ実数配列のポインタ
   */
  const double* GetPitch() const;

  /** グローバルピッチを取得
   *  カーテシアンのときはGetPitchと同じ
   *  LMRのときは最大レベルにおけるピッチ
   *  @return ピッチ実数配列のポインタ
   */
  const double* GetGlobalPitch() const;

  /** 全体ボクセル数を取得
   *  @return 全体ボクセル数整数配列のポインタ
   */
  const int* GetGlobalVoxelSize() const;

  /** 全体空間の原点を取得
   *  @return 全体空間の原点実数配列のポインタ
   */
  const double* GetGlobalOrigin() const;

  /** 全体空間サイズを取得
   *  @return 全体空間サイズ実数配列のポインタ
   */
  const double* GetGlobalRegion() const;

  /** 自ランクのボクセル数を取得
   *  @return 自ランクのボクセル数整数配列のポインタ
   */
  const int* GetLocalVoxelSize() const;

  /** 自ランクの空間原点を取得
   *  @return 自ランクの空間原点実数配列のポインタ
   */
  const double* GetLocalOrigin() const;

  /** 自ランクの空間サイズを取得
   *  @return 自ランクの空間サイズ実数配列のポインタ
   */
  const double* GetLocalRegion() const;

  /** 自ランクの始点VOXELの全体空間でのインデクスを取得
   *  LMRのときは最大レベルにおける始点インデクスを返す
   *  @return 自ランクの始点インデクス整数配列のポインタ
   */
  const int* GetVoxelHeadIndex() const;

  /** 自ランクの終点VOXELの全体空間でのインデクスを取得
   *  LMRのときは最大レベルにおける終点インデクスを返す
   *  @return 自ランクの終点インデクス整数配列のポインタ
   */
  const int* GetVoxelTailIndex() const;

  /** 自ランクの隣接ランク番号を取得
   *  LMRで隣接ランクが4つの場合は、1番目のランクを返す
   *  @return 自ランクの隣接ランク番号整数配列のポインタ
   */
  const int* GetNeighborRankID() const;

  /** 自ランクの周期境界の隣接ランク番号を取得
   *  LMRで隣接ランクが4つの場合は、1番目のランクを返す
   *  @return 自ランクの周期境界の隣接ランク番号整数配列のポインタ
   */
  const int* GetPeriodicRankID() const;

  /** 指定面における自ランクの隣接ランク番号を取得
   *  @param[in]  face 面方向
   *  @param[out] num  面の数(CARTのとき1)
   *  @return 指定面における自ランクの隣接ランク番号整数配列のポインタ
   */
  virtual
  const int* GetNeighborRankList( cpm_FaceFlag face, int &num ) const;

  /** 指定面における自ランクの周期境界の隣接ランク番号を取得
   *  @param[in]  face 面方向
   *  @param[out] num  面の数(CARTのとき1)
   *  @return 指定面における自ランクの周期境界の隣接ランク番号整数配列のポインタ
   */
  virtual
  const int* GetPeriodicRankList( cpm_FaceFlag face, int &num ) const;

  /** 指定面におけるレベル差を取得
   *  @param[in]  face 面方向
   *  @return     レベル差(0:同じレベル, 1:fine, -1:coarse)
   */
  virtual
  int GetNeighborLevelDiff( cpm_FaceFlag face ) const;

  /** 自ランクの境界が外部境界かどうかを判定
   *  @param[in] face  面方向
   *  @retval    true  外部境界
   *  @retval    false 外部境界でない
   */
  virtual
  bool IsOuterBoundary( cpm_FaceFlag face ) const;

  /** 自ランクの境界が内部境界(隣が不活性ドメイン)かどうかを判定
   *  @param[in] face  面方向
   *  @retval    true  内部境界
   *  @retval    false 内部境界でない
   */
  virtual
  bool IsInnerBoundary( cpm_FaceFlag face ) const;



////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


protected:
  /**** 全体空間の情報 ****/
  cpm_GlobalDomainInfo m_globalDomainInfo; ///< 空間全体の領域情報
 
  /**** 自ランクの情報 ****/
  cpm_LocalDomainInfo m_localDomainInfo;   ///< 自ランクの領域情報
  int m_voxelHeadIndex[3];                 ///< 自ランクの始点ボクセルインデックス
  int m_voxelTailIndex[3];                 ///< 自ランクの終点ボクセルインデックス

  /**** 並列情報 ****/
  MPI_Comm  m_comm;        ///< MPIコミュニケータ
  int m_nRank;             ///< コミュニケータ内のランク数(=プロセス並列数)
  int m_rankNo;            ///< コミュニケータ内でのランク番号
  int m_neighborRankID[6]; ///< 隣接ランク番号(外部境界は負の値)
  int m_periodicRankID[6]; ///< 周期境界の隣接ランク番号
};

#endif /* _CPM_VOXELINFO_H_ */

