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
 * @file   cpm_VoxelInfoCART.h
 * カーテシアン用のVOXEL空間情報クラスのヘッダーファイルo
 * @date   2015/03/27
 */

#ifndef _CPM_VOXELINFO_CART_H_
#define _CPM_VOXELINFO_CART_H_

#include "cpm_VoxelInfo.h"

/** カーテシアン用のVOXEL空間情報管理クラス
 */
class cpm_VoxelInfoCART : public cpm_VoxelInfo
{
friend class cpm_ParaManager;
friend class cpm_ParaManagerCART;
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:


protected:
  /** コンストラクタ */
  cpm_VoxelInfoCART();

  /** デストラクタ */
  virtual ~cpm_VoxelInfoCART();

  /** CPM領域分割情報の生成
   *  - MPI_COMM_WORLDを使用した領域を生成する。
   *
   *  @param[in]  comm  MPIコミュニケータ
   *  @param[in]  dInfo 領域分割情報
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Init( MPI_Comm comm, cpm_GlobalDomainInfo* dInfo );

  /** ランクマップを生成
   *  @retval true  正常終了
   *  @retval false エラー
   */
  bool CreateRankMap();

  /** 隣接ランク情報を生成
   *  @retval true  正常終了
   *  @retval false エラー
   */
  bool CreateNeighborRankInfo();

  /** ローカル領域情報を生成
   *  @retval true  正常終了
   *  @retval false エラー
   */
  bool CreateLocalDomainInfo();







////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


protected:
  /**** 並列情報 ****/
  int *m_rankMap; ///< ランクマップ
};

#endif /* _CPM_VOXELINFO_CART_H_ */
