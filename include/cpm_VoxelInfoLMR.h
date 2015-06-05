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
 * @file   cpm_VoxelInfoLMR.h
 * LMR用のVOXEL空間情報クラスのヘッダーファイル
 * @author University of Tokyo
 * @date   2015/03/27
 */

#ifndef _CPM_VOXELINFO_LMR_H_
#define _CPM_VOXELINFO_LMR_H_

#include "cpm_VoxelInfo.h"
#include "BCMOctree.h"
#include "BCMFileCommon.h"
#include "cpm_TextParserDomainLMR.h"

/** LMR用のVOXEL空間情報管理クラス
 */
class cpm_VoxelInfoLMR : public cpm_VoxelInfo
{
friend class cpm_ParaManager;
friend class cpm_ParaManagerLMR;
////////////////////////////////////////////////////////////////////////////////
// メンバー関数 
////////////////////////////////////////////////////////////////////////////////
public:


protected:
  /** コンストラクタ */
  cpm_VoxelInfoLMR();

  /** デストラクタ */
  virtual ~cpm_VoxelInfoLMR();


////// 領域分割情報生成関係 //////

  /** CPM領域分割情報の生成
   *  - MPI_COMM_WORLDを使用した領域を生成する。
   *  @param[in]  comm   MPIコミュニケータ
   *  @param[in]  tpFile 領域情報ファイル
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Init( MPI_Comm comm, std::string treeFile );

  /** 木情報ファイルの読み込み
   *  @param[in]  octFile   木情報ファイル
   *  @param[out] header    ヘッダー情報
   *  @param[out] pedigrees ぺディグリー情報
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  static
  cpm_ErrorCode LoadOctreeFile( std::string octFile, BCMFileIO::OctHeader &header, std::vector<Pedigree> &pedigrees );

  /** 木情報ファイルのヘッダー読み込み
   *  ヘッダーのみを読み込み、ファイルをクローズする
   *  @param[in]  octFile 木情報ファイル
   *  @param[out] header  ヘッダー情報
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  static
  cpm_ErrorCode LoadOctreeHeader( std::string octFile, BCMFileIO::OctHeader &header );

  /** 木情報ファイルのヘッダー読み込み
   *  @param[in]  fp         木情報ファイルポインタ
   *  @param[out] header     ヘッダー情報
   *  @param[out] isNeedSwap エンディアン変換フラグ(true:要変換)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  static
  cpm_ErrorCode LoadOctreeHeader( FILE *fp, BCMFileIO::OctHeader &header, bool &isNeedSwap );

  /** グローバルの領域情報をセット
   *  @param[in] dInfo 領域情報
   */
  void SetGlobaliDomainInfo( S_OCT_DOMAIN_INFO &dInfo );

  /** ローカルの領域情報をセット
   *  @param[in] dInfo 領域情報
   */
  void SetLocalDomainInfo( S_OCT_DOMAIN_INFO &dInfo );

  /** 隣接情報の取得
   */
  void SetNeighborInfo();

  /** 木情報ファイルからリーフ数を取得する
   *  @param[in] treefile  木情報ファイル
   *  @return    リーフ数
   */
  static
  int GetNumLeaf( std::string treeFile );


////// 領域取得関係 //////

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
  /**** 木情報 ****/
  BCMFileIO::OctHeader m_octHeader;  ///< 木情報ファイルのヘッダー情報
  BCMOctree *m_octree;               ///< 生成された木情報
  Node      *m_node;                 ///< 自ランクが担当するリーフノード

  /**** 隣接情報 ****/
  const NeighborInfo *m_neighborInfo; ///< BCMOctreeから生成した隣接情報
  int m_neighborRankID_LMR[6][4]; ///< 隣接ランク番号(外部境界は負の値)
  int m_periodicRankID_LMR[6][4]; ///< 周期境界の隣接ランク番号
  int m_neighborLevelDiff[6]; ///< 隣接ランクとのレベル差(-1/0/1)

};

#endif /* _CPM_VOXELINFO_LMR_H_ */

