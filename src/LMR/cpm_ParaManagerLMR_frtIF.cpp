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
 * @file   cpm_ParaManager_frtIF.cpp
 * パラレルマネージャクラスのFortranインターフェイスのソースファイル
 * @date   2012/05/31
 */
#include "cpm_ParaManagerLMR.h"

/** extern宣言 */
#define CPM_EXTERN extern "C"

#if 0
/** S3D,V3D通信でS4D版を使う */
#define _USE_S4D_
#endif

#ifndef CPM_WINDOWS
  #define cpm_Initialize_LMR_          cpm_initialize_lmr_
  #define cpm_IsParallel_LMR_          cpm_isparallel_lmr_
  #define cpm_GetNumLeaf_LMR_          cpm_getnumleaf_lmr_
  #define cpm_GetLocalNumLeaf_LMR_     cpm_getlocalnumleaf_lmr_
  #define cpm_GetLeafID_LMR_           cpm_getleafid_lmr_
  #define cpm_GetDivNum_LMR_           cpm_getdivnum_lmr_
  #define cpm_GetPitch_LMR_            cpm_getpitch_lmr_
  #define cpm_GetGlobalVoxelSize_LMR_  cpm_getglobalvoxelsize_lmr_
  #define cpm_GetGlobalNodeSize_LMR_   cpm_getglobalnodesize_lmr_
  #define cpm_GetGlobalArraySize_LMR_  cpm_getglobalarraysize_lmr_
  #define cpm_GetGlobalOrigin_LMR_     cpm_getglobalorigin_lmr_
  #define cpm_GetGlobalRegion_LMR_     cpm_getglobalregion_lmr_
  #define cpm_GetLocalVoxelSize_LMR_   cpm_getlocalvoxelsize_lmr_
  #define cpm_GetLocalNodeSize_LMR_    cpm_getlocalnodesize_lmr_
  #define cpm_GetLocalArraySize_LMR_   cpm_getlocalarraysize_lmr_
  #define cpm_GetLocalOrigin_LMR_      cpm_getlocalorigin_lmr_
  #define cpm_GetLocalRegion_LMR_      cpm_getlocalregion_lmr_
  #define cpm_GetDivPos_LMR_           cpm_getdivpos_lmr_
  #define cpm_GetVoxelHeadIndex_LMR_   cpm_getvoxelheadindex_lmr_
  #define cpm_GetVoxelTailIndex_LMR_   cpm_getvoxeltailindex_lmr_
  #define cpm_GetNodeHeadIndex_LMR_    cpm_getnodeheadindex_lmr_
  #define cpm_GetNodeTailIndex_LMR_    cpm_getnodetailindex_lmr_
  #define cpm_GetArrayHeadIndex_LMR_   cpm_getarrayheadindex_lmr_
  #define cpm_GetArrayTailIndex_LMR_   cpm_getarraytailindex_lmr_
  #define cpm_GetDefPointType_LMR_     cpm_getdefpointtype_lmr_
  #define cpm_GetNeighborRankList_LMR_ cpm_getneighborranklist_lmr_
  #define cpm_GetPeriodicRankList_LMR_ cpm_getperiodicranklist_lmr_
  #define cpm_GetNeighborLeafList_LMR_ cpm_getneighborleaflist_lmr_
  #define cpm_GetPeriodicLeafList_LMR_ cpm_getperiodicleaflist_lmr_
  #define cpm_GetMyRankID_LMR_         cpm_getmyrankid_lmr_
  #define cpm_GetNumRank_LMR_          cpm_getnumrank_lmr_
  #define cpm_Abort_LMR_               cpm_abort_lmr_
  #define cpm_Barrier_LMR_             cpm_barrier_lmr_
  #define cpm_Wait_LMR_                cpm_wait_lmr_
  #define cpm_Waitall_LMR_             cpm_waitall_lmr_
  #define cpm_Bcast_LMR_               cpm_bcast_lmr_
  #define cpm_Send_LMR_                cpm_send_lmr_
  #define cpm_Recv_LMR_                cpm_recv_lmr_
  #define cpm_Isend_LMR_               cpm_isend_lmr_
  #define cpm_Irecv_LMR_               cpm_irecv_lmr_
  #define cpm_Allreduce_LMR_           cpm_allreduce_lmr_
  #define cpm_Gather_LMR_              cpm_gather_lmr_
  #define cpm_Allgather_LMR_           cpm_allgather_lmr_
  #define cpm_Gatherv_LMR_             cpm_gatherv_lmr_
  #define cpm_Allgatherv_LMR_          cpm_allgatherv_lmr_
  #define cpm_BndCommS3D_LMR_          cpm_bndcomms3d_lmr_
  #define cpm_BndCommV3D_LMR_          cpm_bndcommv3d_lmr_
  #define cpm_BndCommS4D_LMR_          cpm_bndcomms4d_lmr_
//  #define cpm_wait_BndCommS3D_LMR_   cpm_wait_bndcomms3d_lmr_
//  #define cpm_wait_BndCommV3D_LMR_   cpm_wait_bndcommv3d_lmr_
//  #define cpm_wait_BndCommS4D_LMR_   cpm_wait_bndcomms4d_lmr_
  #define cpm_BndCommV3DEx_LMR_        cpm_bndcommv3dex_lmr_
  #define cpm_BndCommS4DEx_LMR_        cpm_bndcomms4dex_lmr_
//  #define cpm_wait_BndCommV3DEx_LMR_ cpm_wait_bndcommv3dex_lmr_
//  #define cpm_wait_BndCommS4DEx_LMR_ cpm_wait_bndcomms4dex_lmr_
  #define cpm_PeriodicCommS3D_LMR_     cpm_periodiccomms3d_lmr_
  #define cpm_PeriodicCommV3D_LMR_     cpm_periodiccommv3d_lmr_
  #define cpm_PeriodicCommS4D_LMR_     cpm_periodiccomms4d_lmr_
  #define cpm_PeriodicCommV3DEx_LMR_   cpm_periodiccommv3dex_lmr_
  #define cpm_PeriodicCommS4DEx_LMR_   cpm_periodiccomms4dex_lmr_
#else
  #define cpm_Initialize_LMR_          CPM_INITIALIZE_LMR
  #define cpm_IsParallel_LMR_          CPM_ISPARALLEL_LMR
  #define cpm_GetNumLeaf_LMR_          CPM_GETNUMLEAF_LMR
  #define cpm_GetLocalNumLeaf_LMR_     CPM_GETLOCALNUMLEAF_LMR
  #define cpm_GetLeafID_LMR_           CPM_GETLEAFID_LMR
  #define cpm_GetDivNum_LMR_           CPM_GETDIVNUM_LMR
  #define cpm_GetPitch_LMR_            CPM_GETPITCH_LMR
  #define cpm_GetGlobalVoxelSize_LMR_  CPM_GETGLOBALVOXELSIZE_LMR
  #define cpm_GetGlobalNodeSize_LMR_   CPM_GETGLOBALNODESIZE_LMR
  #define cpm_GetGlobalArraySize_LMR_  CPM_GETGLOBALARRAYSIZE_LMR
  #define cpm_GetGlobalOrigin_LMR_     CPM_GETGLOBALORIGIN_LMR
  #define cpm_GetGlobalRegion_LMR_     CPM_GETGLOBALREGION_LMR
  #define cpm_GetLocalVoxelSize_LMR_   CPM_GETLOCALVOXELSIZE_LMR
  #define cpm_GetLocalNodeSize_LMR_    CPM_GETLOCALNODESIZE_LMR
  #define cpm_GetLocalArraySize_LMR_   CPM_GETLOCALARRAYSIZE_LMR
  #define cpm_GetLocalOrigin_LMR_      CPM_GETLOCALORIGIN_LMR
  #define cpm_GetLocalRegion_LMR_      CPM_GETLOCALREGION_LMR
  #define cpm_GetDivPos_LMR_           CPM_GETDIVPOS_LMR
  #define cpm_GetVoxelHeadIndex_LMR_   CPM_GETVOXELHEADINDEX_LMR
  #define cpm_GetVoxelTailIndex_LMR_   CPM_GETVOXELTAILINDEX_LMR
  #define cpm_GetNodeHeadIndex_LMR_    CPM_GETNODEHEADINDEX_LMR
  #define cpm_GetNodeTailIndex_LMR_    CPM_GETNODETAILINDEX_LMR
  #define cpm_GetArrayHeadIndex_LMR_   CPM_GETARRAYHEADINDEX_LMR
  #define cpm_GetArrayTailIndex_LMR_   CPM_GETARRAYTAILINDEX_LMR
  #define cpm_GetDefPointType_LMR_     CPM_GETDEFPOINTTYPE_LMR
  #define cpm_GetNeighborRankList_LMR_ CPM_GETNEIGHBORRANKLIST_LMR
  #define cpm_GetPeriodicRankList_LMR_ CPM_GETPERIODICRANKLIST_LMR
  #define cpm_GetNeighborLeafList_LMR_ CPM_GETNEIGHBORLEAFLIST_LMR
  #define cpm_GetPeriodicLeafList_LMR_ CPM_GETPERIODICLEAFLIST_LMR
  #define cpm_GetMyRankID_LMR_         CPM_GETMYRANKID_LMR
  #define cpm_GetNumRank_LMR_          CPM_GETNUMRANK_LMR
  #define cpm_Abort_LMR_               CPM_ABORT_LMR
  #define cpm_Barrier_LMR_             CPM_BARRIER_LMR
  #define cpm_Wait_LMR_                CPM_WAIT_LMR
  #define cpm_Waitall_LMR_             CPM_WAITALL_LMR
  #define cpm_Bcast_LMR_               CPM_BCAST_LMR
  #define cpm_Send_LMR_                CPM_SEND_LMR
  #define cpm_Recv_LMR_                CPM_RECV_LMR
  #define cpm_Isend_LMR_               CPM_ISEND_LMR
  #define cpm_Irecv_LMR_               CPM_IRECV_LMR
  #define cpm_Allreduce_LMR_           CPM_ALLREDUCE_LMR
  #define cpm_Gather_LMR_              CPM_GATHER_LMR
  #define cpm_Allgather_LMR_           CPM_ALLGATHER_LMR
  #define cpm_Gatherv_LMR_             CPM_GATHERV_LMR
  #define cpm_Allgatherv_LMR_          CPM_ALLGATHERV_LMR
  #define cpm_BndCommS3D_LMR_          CPM_BNDCOMMS3D_LMR
  #define cpm_BndCommV3D_LMR_          CPM_BNDCOMMV3D_LMR
  #define cpm_BndCommS4D_LMR_          CPM_BNDCOMMS4D_LMR
//  #define cpm_wait_BndCommS3D_LMR_     CPM_WAIT_BNDCOMMS3D_LMR
//  #define cpm_wait_BndCommV3D_LMR_     CPM_WAIT_BNDCOMMV3D_LMR
//  #define cpm_wait_BndCommS4D_LMR_     CPM_WAIT_BNDCOMMS4D_LMR
  #define cpm_BndCommV3DEx_LMR_        CPM_BNDCOMMV3DEX_LMR
  #define cpm_BndCommS4DEx_LMR_        CPM_BNDCOMMS4DEX_LMR
//  #define cpm_wait_BndCommV3DEx_LMR_   CPM_WAIT_BNDCOMMV3DEX_LMR
//  #define cpm_wait_BndCommS4DEx_LMR_   CPM_WAIT_BNDCOMMS4DEX_LMR
  #define cpm_PeriodicCommS3D_LMR_     CPM_PERIODICCOMMS3D_LMR
  #define cpm_PeriodicCommV3D_LMR_     CPM_PERIODICCOMMV3D_LMR
  #define cpm_PeriodicCommS4D_LMR_     CPM_PERIODICCOMMS4D_LMR
  #define cpm_PeriodicCommV3DEx_LMR_   CPM_PERIODICCOMMV3DEX_LMR
  #define cpm_PeriodicCommS4DEx_LMR_   CPM_PERIODICCOMMS4DEX_LMR
#endif

////////////////////////////////////////////////////////////////////////////////
/** 初期化処理(MPI_Initは実行済みの場合)
 *  - InitializeのFortranインターフェイス関数
 *  - FortranでMPI_Initがコールされている必要がある
 *  @param[out] ierr       終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Initialize_LMR_( int *ierr )
{
  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // 戻り値
  *ierr = CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
/** 並列実行であるかチェックする
 *  - IsParallelのFortranインターフェイス関数
 *  @param[out] ipara 並列実行フラグ(1=並列実行、1以外=逐次実行)
 *  @param[out] ierr  終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_IsParallel_LMR_( int *ipara, int *ierr )
{
  if( !ipara || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  *ipara = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // IsParallel
  if( paraMngr->IsParallel() )
    *ipara = 1;
  else
    *ipara = 0;
}

////////////////////////////////////////////////////////////////////////////////
/** 全リーフ数を取得する
 *  @param[out] numLeaf   全リーフ数
 *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 *  @return    全リーフ数
 */
CPM_EXTERN
void
cpm_GetNumLeaf_LMR_( int *numLeaf, int *procGrpNo, int *ierr )
{
  if( !numLeaf || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  *numLeaf = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  *numLeaf = paraMngr->GetNumLeaf( *procGrpNo );

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 自ランクが担当するリーフ数を取得する
 *  @param[out] numLeaf   リーフ数
 *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 *  @return    自ランクが担当するリーフ数
 */
CPM_EXTERN
void
cpm_GetLocalNumLeaf_LMR_( int *numLeaf, int *procGrpNo, int *ierr )
{
  if( !numLeaf || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  *numLeaf = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  *numLeaf = paraMngr->GetLocalNumLeaf( *procGrpNo );

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフのリーフIDを取得する
 *  - GetLeafIDのFortranインターフェイス関数
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[out] leafID    リーフID
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetLeafID_LMR_( int *leafIndex, int *leafID, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !procGrpNo || !leafID || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  *leafID = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  *leafID = paraMngr->GetLeafID( (*leafIndex)-1, *procGrpNo );

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 領域分割数を取得
 *  - GetDivNumのFortranインターフェイス関数
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[out] div       領域分割数(3wordの整数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetDivNum_LMR_( int *leafIndex, int *div, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !procGrpNo || !div || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  div[0] = div[1] = div[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *dd = paraMngr->GetDivNum( (*leafIndex)-1, *procGrpNo );
  if( !dd )
  {
    *ierr = CPM_ERROR_GET_DIVNUM;
    return;
  }

  div[0] = dd[0];
  div[1] = dd[1];
  div[2] = dd[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** ピッチを取得
 *  - GetPitchのFortranインターフェイス関数
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] pch       ピッチ(3wordの実数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetPitch_LMR_( int *leafIndex, double *pch, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !procGrpNo || !pch || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  pch[0] = pch[1] = pch[2] = 0.0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const double *pp = paraMngr->GetPitch( (*leafIndex)-1, *procGrpNo );
  if( !pp )
  {
    *ierr = CPM_ERROR_GET_PITCH;
    return;
  }

  pch[0] = pp[0];
  pch[1] = pp[1];
  pch[2] = pp[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 全体ボクセル数を取得
 *  - GetGlobalVoxelSizeのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] wsz       全体ボクセル数(3wordの整数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetGlobalVoxelSize_LMR_( int *wsz, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !wsz || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  wsz[0] = wsz[1] = wsz[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *sz = paraMngr->GetGlobalVoxelSize( *procGrpNo );
  if( !sz )
  {
    *ierr = CPM_ERROR_GET_GLOBALVOXELSIZE;
    return;
  }

  wsz[0] = sz[0];
  wsz[1] = sz[1];
  wsz[2] = sz[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 全体頂点数を取得
 *  - GetGlobalNodeSizeのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] wsz       全体頂点数(3wordの整数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetGlobalNodeSize_LMR_( int *wsz, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !wsz || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  wsz[0] = wsz[1] = wsz[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *sz = paraMngr->GetGlobalNodeSize( *procGrpNo );
  if( !sz )
  {
    *ierr = CPM_ERROR_GET_GLOBALNODESIZE;
    return;
  }

  wsz[0] = sz[0];
  wsz[1] = sz[1];
  wsz[2] = sz[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 全体ボクセル数または頂点数を取得
 *  - FVMのときはボクセル数、FDMのときは頂点数
 *  - GetGlobalArraySizeのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] wsz       全体ボクセル数または頂点数(3wordの整数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetGlobalArraySize_LMR_( int *wsz, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !wsz || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  wsz[0] = wsz[1] = wsz[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *sz = paraMngr->GetGlobalArraySize( *procGrpNo );
  if( !sz )
  {
    *ierr = CPM_ERROR_GET_GLOBALARRAYSIZE;
    return;
  }

  wsz[0] = sz[0];
  wsz[1] = sz[1];
  wsz[2] = sz[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 全体空間の原点を取得
 *  - GetGlobalOriginのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] worg      全体空間の原点(3wordの実数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetGlobalOrigin_LMR_( double *worg, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !worg || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  worg[0] = worg[1] = worg[2] = 0.0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const double *org = paraMngr->GetGlobalOrigin( *procGrpNo );
  if( !org )
  {
    *ierr = CPM_ERROR_GET_GLOBALORIGIN;
    return;
  }

  worg[0] = org[0];
  worg[1] = org[1];
  worg[2] = org[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 全体空間サイズを取得
 *  - GetGlobalRegionのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] wrgn      全体空間サイズ(3wordの実数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetGlobalRegion_LMR_( double *wrgn, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !wrgn || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  wrgn[0] = wrgn[1] = wrgn[2] = 0.0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const double *rgn = paraMngr->GetGlobalRegion( *procGrpNo );
  if( !rgn )
  {
    *ierr = CPM_ERROR_GET_GLOBALREGION;
    return;
  }

  wrgn[0] = rgn[0];
  wrgn[1] = rgn[1];
  wrgn[2] = rgn[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** １リーフのボクセル数を取得
 *  - GetLocalVoxelSizeのFortranインターフェイス関数
 *  @param[out] lsz       ボクセル数(3wordの整数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetLocalVoxelSize_LMR_( int *lsz, int *procGrpNo, int *ierr )
{
  if( !lsz || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  lsz[0] = lsz[1] = lsz[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *sz = paraMngr->GetLocalVoxelSize( *procGrpNo );
  if( !sz )
  {
    *ierr = CPM_ERROR_GET_LOCALVOXELSIZE;
    return;
  }

  lsz[0] = sz[0];
  lsz[1] = sz[1];
  lsz[2] = sz[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** １リーフの頂点数を取得
 *  - GetLocalNodeSizeのFortranインターフェイス関数
 *  @param[out] lsz       頂点数(3wordの整数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetLocalNodeSize_LMR_( int *lsz, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !lsz || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  lsz[0] = lsz[1] = lsz[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *sz = paraMngr->GetLocalNodeSize( *procGrpNo );
  if( !sz )
  {
    *ierr = CPM_ERROR_GET_LOCALNODESIZE;
    return;
  }

  lsz[0] = sz[0];
  lsz[1] = sz[1];
  lsz[2] = sz[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** １リーフのボクセル数または頂点数を取得
 *  - FVMのときはボクセル数、FDMのときは頂点数
 *  - GetLocalArraySizeのFortranインターフェイス関数
 *  @param[out] lsz       頂点数(3wordの整数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetLocalArraySize_LMR_( int *lsz, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !lsz || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  lsz[0] = lsz[1] = lsz[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *sz = paraMngr->GetLocalArraySize( *procGrpNo );
  if( !sz )
  {
    *ierr = CPM_ERROR_GET_LOCALARRAYSIZE;
    return;
  }

  lsz[0] = sz[0];
  lsz[1] = sz[1];
  lsz[2] = sz[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフの空間原点を取得
 *  - GetLocalOriginのFortranインターフェイス関数
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[out] lorg      空間原点(3wordの実数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetLocalOrigin_LMR_( int *leafIndex, double *lorg, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !procGrpNo || !lorg || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  lorg[0] = lorg[1] = lorg[2] = 0.0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const double *org = paraMngr->GetLocalOrigin( (*leafIndex)-1, *procGrpNo );
  if( !org )
  {
    *ierr = CPM_ERROR_GET_LOCALORIGIN;
    return;
  }

  lorg[0] = org[0];
  lorg[1] = org[1];
  lorg[2] = org[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフの空間サイズを取得
 *  - GetLocalRegionのFortranインターフェイス関数
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[out] lrgn      空間サイズ(3wordの実数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetLocalRegion_LMR_( int *leafIndex, double *lrgn, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !procGrpNo || !lrgn || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  lrgn[0] = lrgn[1] = lrgn[2] = 0.0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const double *rgn = paraMngr->GetLocalRegion( (*leafIndex)-1, *procGrpNo );
  if( !rgn )
  {
    *ierr = CPM_ERROR_GET_LOCALREGION;
    return;
  }

  lrgn[0] = rgn[0];
  lrgn[1] = rgn[1];
  lrgn[2] = rgn[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフの領域分割位置を取得
 *  - GetDivPosのFortranインターフェイス関数
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[out] pos       領域分割位置(3wordの整数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetDivPos_LMR_( int *leafIndex, int *pos, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !procGrpNo || !pos || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  pos[0] = pos[1] = pos[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *pp = paraMngr->GetDivPos( (*leafIndex)-1, *procGrpNo );
  if( !pp )
  {
    *ierr = CPM_ERROR_GET_DIVPOS;
    return;
  }

  pos[0] = pp[0];
  pos[1] = pp[1];
  pos[2] = pp[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフの始点VOXELの全体空間でのインデクスを取得
 *  - GetVoxelHeadIndexのFortranインターフェイス関数
 *  - 全体空間の先頭インデクスを0としたC型のインデクス
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[out] idx       始点VOXELインデクス(3wordの整数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetVoxelHeadIndex_LMR_( int *leafIndex, int *idx, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !procGrpNo || !idx || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

 idx[0] = idx[1] = idx[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *id = paraMngr->GetVoxelHeadIndex( (*leafIndex)-1, *procGrpNo );
  if( !id )
  {
    *ierr = CPM_ERROR_GET_HEADINDEX;
    return;
  }

  idx[0] = id[0];
  idx[1] = id[1];
  idx[2] = id[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフの終点VOXELの全体空間でのインデクスを取得
 *  - GetVoxelTailIndexのFortranインターフェイス関数
 *  - 全体空間の先頭インデクスを0としたC型のインデクス
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[out] idx       終点VOXELインデクス(3wordの整数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetVoxelTailIndex_LMR( int *leafIndex, int *idx, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !procGrpNo || !idx || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

 idx[0] = idx[1] = idx[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *id = paraMngr->GetVoxelTailIndex( (*leafIndex)-1, *procGrpNo );
  if( !id )
  {
    *ierr = CPM_ERROR_GET_TAILINDEX;
    return;
  }

  idx[0] = id[0];
  idx[1] = id[1];
  idx[2] = id[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフの始点頂点の全体空間でのインデクスを取得
 *  - GetNodeHeadIndexのFortranインターフェイス関数
 *  - 全体空間の先頭インデクスを0としたC型のインデクス
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[out] idx       始点頂点インデクス(3wordの整数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetNodeHeadIndex_LMR_( int *leafIndex,int *idx, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !procGrpNo || !idx || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  idx[0] = idx[1] = idx[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *id = paraMngr->GetNodeHeadIndex( (*leafIndex)-1, *procGrpNo );
  if( !id )
  {
    *ierr = CPM_ERROR_GET_HEADINDEX;
    return;
  }

  idx[0] = id[0];
  idx[1] = id[1];
  idx[2] = id[2];

  *ierr = CPM_SUCCESS;

  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフの終点頂点の全体空間でのインデクスを取得
 *  - GetNodeTailIndexのFortranインターフェイス関数
 *  - 全体空間の先頭インデクスを0としたC型のインデクス
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[out] idx       始点頂点インデクス(3wordの整数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetNodeTailIndex_LMR_( int *leafIndex, int *idx, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !procGrpNo || !idx || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  idx[0] = idx[1] = idx[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *id = paraMngr->GetNodeTailIndex( (*leafIndex)-1, *procGrpNo );
  if( !id )
  {
    *ierr = CPM_ERROR_GET_HEADINDEX;
    return;
  }

  idx[0] = id[0];
  idx[1] = id[1];
  idx[2] = id[2];

  *ierr = CPM_SUCCESS;

  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフの始点ボクセルまたは頂点の全体空間でのインデクスを取得
 *  - FVMのときはボクセル、FDMのときは頂点始点インデックスを取得
 *  - GetArrayHeadIndexのFortranインターフェイス関数
 *  - 全体空間の先頭インデクスを0としたC型のインデクス
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[out] idx       始点ボクセルまたは頂点インデクス(3wordの整数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetArrayHeadIndex_LMR_( int *leafIndex, int *idx, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !procGrpNo || !idx || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  idx[0] = idx[1] = idx[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *id = paraMngr->GetArrayHeadIndex( (*leafIndex)-1, *procGrpNo );
  if( !id )
  {
    *ierr = CPM_ERROR_GET_HEADINDEX;
    return;
  }

  idx[0] = id[0];
  idx[1] = id[1];
  idx[2] = id[2];

  *ierr = CPM_SUCCESS;

  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフの終点ボクセルまたは頂点の全体空間でのインデクスを取得
 *  - FVMのときはボクセル、FDMのときは頂点終点インデックスを取得
 *  - GetArrayTailIndexのFortranインターフェイス関数
 *  - 全体空間の先頭インデクスを0としたC型のインデクス
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[out] idx       終点ボクセルまたは頂点インデクス(3wordの整数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetArrayTailIndex_LMR_( int *leafIndex, int *idx, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !procGrpNo || !idx || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  idx[0] = idx[1] = idx[2] = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *id = paraMngr->GetArrayTailIndex( (*leafIndex)-1, *procGrpNo );
  if( !id )
  {
    *ierr = CPM_ERROR_GET_HEADINDEX;
    return;
  }

  idx[0] = id[0];
  idx[1] = id[1];
  idx[2] = id[2];

  *ierr = CPM_SUCCESS;

  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 定義点タイプを取得
 *  - GetDefPointTypeのFortranインターフェイス関数
 *  - 全体空間の先頭インデクスを0としたC型のインデクス
 *  @param[out] ideftyp   定義点タイプ(-1=未定義、0=FVM、1=FDM )
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetDefPointType_LMR_( int *ideftyp, int *procGrpNo, int *ierr )
{

  if( !procGrpNo || !ideftyp || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  *ideftyp = paraMngr->GetDefPointType( *procGrpNo );

  return;
}

////////////////////////////////////////////////////////////////////////////////
/** ランク番号の取得
 *  - GetMyRankIDのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] id        ランク番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetMyRankID_LMR_( int *id, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !id || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  *id = 0;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  *id = paraMngr->GetMyRankID( *procGrpNo );
  if( !id )
  {
    *ierr = CPM_ERROR_GET_MYRANK;
    return;
  }

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** ランク数の取得
 *  - GetNumRankのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] nrank     ランク数
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetNumRank_LMR_( int *nrank, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !nrank || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  *nrank = 1;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  *nrank = paraMngr->GetNumRank( *procGrpNo );
  if( !nrank )
  {
    *ierr = CPM_ERROR_GET_NUMRANK;
    return;
  }

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフの指定面における自リーフの隣接ランク番号を取得
 *  - GetNeighborRankListのFortranインターフェイス関数
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[in]  face      面方向
 *  @param[out] rankList  指定リーフの指定面における自リーフの隣接ランク番号整数配列
 *  @param[out] numRank   面の数(0 or 1 or 4)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetNeighborRankList_LMR_( int *leafIndex, int *face, int *rankList, int *numRank, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !face || !rankList || !numRank || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  *numRank = 0;
  rankList[0] = -1;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *rtmp = paraMngr->GetNeighborRankList( (*leafIndex)-1, (cpm_FaceFlag)(*face), *numRank, *procGrpNo);
  if( !rtmp )
  {
    *ierr = CPM_ERROR_GET_NUMRANK;
    return;
  }

  for( int i=0;i<*numRank;i++ )
  {
    rankList[i] = rtmp[i];
  }

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフの指定面における自リーフの周期境界の隣接ランク番号を取得
 *  - GetPeriodicRankListのFortranインターフェイス関数
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[in]  face      面方向
 *  @param[out] rankList  指定リーフの指定面における自リーフの周期境界の隣接ランク番号整数配列
 *  @param[out] numRank   面の数(0 or 1 or 4)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetPeriodicRankList_LMR_( int *leafIndex, int *face, int *rankList, int *numRank, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !face || !rankList || !numRank || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  *numRank = 0;
  rankList[0] = -1;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *rtmp = paraMngr->GetPeriodicRankList( (*leafIndex)-1, (cpm_FaceFlag)(*face), *numRank, *procGrpNo);
  if( !rtmp )
  {
    *ierr = CPM_ERROR_GET_NUMRANK;
    return;
  }

  for( int i=0;i<*numRank;i++ )
  {
    rankList[i] = rtmp[i];
  }

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフの指定面における隣接リーフ番号を取得
 *  - GetNeighborLeafListのFortranインターフェイス関数
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[in]  face      面方向
 *  @param[out] leafList  指定リーフの指定面における隣接リーフ番号配列
 *  @param[out] numLeaf   面の数(0 or 1 or 4)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetNeighborLeafList_LMR_( int *leafIndex, int *face, int *leafList, int *numLeaf, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !face || !leafList || !numLeaf || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  *numLeaf = 0;
  leafList[0] = -1;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *rtmp = paraMngr->GetNeighborLeafList( (*leafIndex)-1, (cpm_FaceFlag)(*face), *numLeaf, *procGrpNo);
  if( !rtmp )
  {
    *ierr = CPM_ERROR_GET_NUMRANK;
    return;
  }

  for( int i=0;i<*numLeaf;i++ )
  {
    leafList[i] = rtmp[i];
  }

  *ierr = CPM_SUCCESS;
  return;
}

  #define cpm_GetPeriodicLeafList_LMR_ CPM_GETPERIODICLEAFLIST_LMR

////////////////////////////////////////////////////////////////////////////////
/** 指定リーフの指定面における周期境界の隣接リーフ番号を取得
 *  - GetPeriodicLeafListのFortranインターフェイス関数
 *  @param[in]  leafIndex リーフ順番号(1~)
 *  @param[in]  face      面方向
 *  @param[out] leafList  指定リーフの指定面における周期境界の隣接リーフ番号配列
 *  @param[out] numLeaf   面の数(0 or 1 or 4)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetPeriodicLeafList_LMR_( int *leafIndex, int *face, int *leafList, int *numLeaf, int *procGrpNo, int *ierr )
{
  if( !leafIndex || !face || !leafList || !numLeaf || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  *numLeaf = 0;
  leafList[0] = -1;

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *rtmp = paraMngr->GetPeriodicLeafList( (*leafIndex)-1, (cpm_FaceFlag)(*face), *numLeaf, *procGrpNo);
  if( !rtmp )
  {
    *ierr = CPM_ERROR_GET_NUMRANK;
    return;
  }

  for( int i=0;i<*numLeaf;i++ )
  {
    leafList[i] = rtmp[i];
  }

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** Abort
 *  - AbortのFortranインターフェイス関数
 *  @param[in]  errorcode MPI_Abortに渡すエラーコード
 */
CPM_EXTERN
void
cpm_Abort_LMR_( int *errorcode )
{
  int err = 0;
  if( errorcode )
  {
    err = *errorcode;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    MPI_Abort( MPI_COMM_WORLD, err );
    exit(err);
    return;
  }

  paraMngr->Abort( err );
}

////////////////////////////////////////////////////////////////////////////////
/** Barrier
 *  - BarrierのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Barrier_LMR_( int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // Barrier
  *ierr = paraMngr->Barrier( *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Wait
 *  - WaitのFortranインターフェイス関数
 *  @param[in]  reqNo     リクエスト番号(0以上の整数)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Wait_LMR_( int *reqNo, int *ierr )
{
  if( !reqNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_Wait
  *ierr = paraMngr->cpm_Wait( *reqNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Waitall
 *  - WaitallのFortranインターフェイス関数
 *  @param[in]  count     リクエストの数
 *  @param[in]  reqlist   リクエスト番号のリスト(0以上の整数)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Waitall_LMR_( int *count, int *reqlist, int *ierr )
{
  if( !count || !reqlist || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_Waitall
  *ierr = paraMngr->cpm_Waitall( *count, reqlist );
}

////////////////////////////////////////////////////////////////////////////////
/** Bcast
 *  - BcastのFortranインターフェイス関数
 *  @param[inout] buf       送受信バッファ
 *  @param[in]    count     送信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    root      送信元のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Bcast_LMR_( void *buf, int *count, int *datatype, int *root, int *procGrpNo, int *ierr )
{
  if( !buf || !count || !datatype || !root || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Bcast
  *ierr = paraMngr->Bcast( dtype, buf, *count, *root, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Send
 *  - SendのFortranインターフェイス関数
 *  @param[inout] buf       送信バッファ
 *  @param[in]    count     送信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    dest      送信先のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Send_LMR_( void *buf, int *count, int *datatype, int *dest, int *procGrpNo, int *ierr )
{
  if( !buf || !count || !datatype || !dest || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Send
  *ierr = paraMngr->Send( dtype, buf, *count, *dest, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Recv
 *  - RecvのFortranインターフェイス関数
 *  @param[inout] buf       受信バッファ
 *  @param[in]    count     受信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    source    送信元のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Recv_LMR_( void *buf, int *count, int *datatype, int *source, int *procGrpNo, int *ierr )
{
  if( !buf || !count || !datatype || !source || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Recv
  *ierr = paraMngr->Recv( dtype, buf, *count, *source, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Isend
 *  - IsendのFortranインターフェイス関数
 *  @param[inout] buf       送信バッファ
 *  @param[in]    count     送信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    dest      送信先のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   reqNo     リクエスト番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Isend_LMR_( void *buf, int *count, int *datatype, int *dest, int *procGrpNo, int *reqNo, int *ierr )
{
  if( !buf || !count || !datatype || !dest || !procGrpNo || !reqNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_Isend
  if( (*ierr = paraMngr->cpm_Isend( buf, *count, *datatype, *dest, reqNo, *procGrpNo ) ) != CPM_SUCCESS )
  {
    return;
  }

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** Irecv
 *  - IrecvのFortranインターフェイス関数
 *  @param[inout] buf       受信バッファ
 *  @param[in]    count     受信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    source    送信元先のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[in]    reqNo     リクエスト番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Irecv_LMR_( void *buf, int *count, int *datatype, int *source, int *procGrpNo, int *reqNo, int *ierr )
{
  if( !buf || !count || !datatype || !source || !procGrpNo || !reqNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_Irecv
  if( (*ierr = paraMngr->cpm_Irecv( buf, *count, *datatype, *source, reqNo, *procGrpNo ) ) != CPM_SUCCESS )
  {
    return;
  }

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** MPI_AllreduceのFortranインターフェイス
 *  - MPI_AllreduceのFortranインターフェイス関数
 *  @param[in]  sendbuf   送信データ
 *  @param[out] recvbuf   受信データ
 *  @param[in]  count     送受信データのサイズ
 *  @param[in]  datatype  送受信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]  op        オペレータ
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Allreduce_LMR_( void *sendbuf, void *recvbuf, int *count, int *datatype, int *op, int *procGrpNo, int *ierr )
{
  if( !sendbuf || !recvbuf || !count || !datatype || !op || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // MPI_Op
  MPI_Op ope = cpm_BaseParaManager::GetMPI_Op( *op );
  if( ope == MPI_OP_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_OPERATOR;
    return;
  }

  // Allreduce
  *ierr = paraMngr->Allreduce( dtype, sendbuf, recvbuf, *count, ope, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** MPI_GatherのFortranインターフェイス
 *  - MPI_GatherのFortranインターフェイス関数
 *  @param[in]  sendbuf   送信データ
 *  @param[in]  sendcnt   送信データのサイズ
 *  @param[in]  sendtype  送信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out] recvbuf   受信データ
 *  @param[in]  recvcnt   受信データのサイズ
 *  @param[in]  recvtype  受信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]  root      受信するランク番号(procGrpNo内でのランク番号)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Gather_LMR_( void *sendbuf, int *sendcnt, int *sendtype, void *recvbuf, int *recvcnt, int *recvtype
           , int *root, int *procGrpNo, int *ierr )
{
  if( !sendbuf || !sendcnt || !sendtype || !recvbuf || !recvcnt || !recvtype ||
      !root || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype stype = cpm_BaseParaManager::GetMPI_Datatype( *sendtype );
  if( stype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }
  MPI_Datatype rtype = cpm_BaseParaManager::GetMPI_Datatype( *recvtype );
  if( rtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Gather
  *ierr = paraMngr->Gather( stype, sendbuf, *sendcnt, rtype, recvbuf, *recvcnt, *root, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** MPI_AllgatherのFortranインターフェイス
 *  - MPI_AllgatherのFortranインターフェイス関数
 *  @param[in]  sendbuf   送信データ
 *  @param[in]  sendcnt   送信データのサイズ
 *  @param[in]  sendtype  送信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out] recvbuf   受信データ
 *  @param[in]  recvcnt   受信データのサイズ
 *  @param[in]  recvtype  受信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Allgather_LMR_( void *sendbuf, int *sendcnt, int *sendtype, void *recvbuf, int *recvcnt, int *recvtype
              , int *procGrpNo, int *ierr )
{
  if( !sendbuf || !sendcnt || !sendtype || !recvbuf || !recvcnt || !recvtype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype stype = cpm_BaseParaManager::GetMPI_Datatype( *sendtype );
  if( stype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }
  MPI_Datatype rtype = cpm_BaseParaManager::GetMPI_Datatype( *recvtype );
  if( rtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Allgather
  *ierr = paraMngr->Allgather( stype, sendbuf, *sendcnt, rtype, recvbuf, *recvcnt, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** MPI_GathervのFortranインターフェイス
 *  - MPI_GathervのFortranインターフェイス関数
 *  @param[in]  sendbuf   送信データ
 *  @param[in]  sendcnt   送信データのサイズ
 *  @param[in]  sendtype  送信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out] recvbuf   受信データ
 *  @param[in]  recvcnts  各ランクからの受信データサイズ
 *  @param[in]  displs    各ランクからの受信データ配置位置
 *  @param[in]  recvtype  受信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]  root      受信するランク番号(procGrpNo内でのランク番号)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Gatherv_LMR_( void *sendbuf, int *sendcnt, int *sendtype, void *recvbuf, int *recvcnts, int *displs, int *recvtype
            , int *root, int *procGrpNo, int *ierr )
{
  if( !sendbuf || !sendcnt || !sendtype || !recvbuf || !recvcnts || !displs || !recvtype ||
      !root || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype stype = cpm_BaseParaManager::GetMPI_Datatype( *sendtype );
  if( stype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }
  MPI_Datatype rtype = cpm_BaseParaManager::GetMPI_Datatype( *recvtype );
  if( rtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Gatherv
  *ierr = paraMngr->Gatherv( stype, sendbuf, *sendcnt, rtype, recvbuf, recvcnts, displs, *root, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** MPI_AllgathervのFortranインターフェイス
 *  - MPI_AllgathervのFortranインターフェイス関数
 *  @param[in]  sendbuf   送信データ
 *  @param[in]  sendcnt   送信データのサイズ
 *  @param[in]  sendtype  送信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out] recvbuf   受信データ
 *  @param[in]  recvcnts  各ランクからの受信データサイズ
 *  @param[in]  displs    各ランクからの受信データ配置位置
 *  @param[in]  recvtype  受信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Allgatherv_LMR_( void *sendbuf, int *sendcnt, int *sendtype, void *recvbuf, int *recvcnts, int *displs, int *recvtype
               , int *procGrpNo, int *ierr )
{
  if( !sendbuf || !sendcnt || !sendtype || !recvbuf || !recvcnts || !displs || !recvtype ||
      !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype stype = cpm_BaseParaManager::GetMPI_Datatype( *sendtype );
  if( stype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }
  MPI_Datatype rtype = cpm_BaseParaManager::GetMPI_Datatype( *recvtype );
  if( rtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Gatherv
  *ierr = paraMngr->Allgatherv( stype, sendbuf, *sendcnt, rtype, recvbuf, recvcnts, displs, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 袖通信(Scalar4D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,nmax,nLeaf)の形式の配列の袖通信を行う
 *  - BndCommS4DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    nmax      配列サイズ(成分数)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommS4D_LMR_( void *array, int *imax, int *jmax, int *kmax, int *nmax, int *vc, int *vc_comm
               , int *datatype, int *procGrpNo, int *ierr )
{
  if( !array || !imax || !jmax || !kmax || !nmax || !vc || !vc_comm || !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // BndCommS4D
   *ierr = paraMngr->BndCommS4D( dtype, array, *imax, *jmax, *kmax, *nmax, *vc, *vc_comm, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 袖通信(Scalar3D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,nLeaf)の形式の配列の袖通信を行う
 *  - BndCommS3DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommS3D_LMR_( void *array, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm, int *datatype
               , int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 1;
  cpm_BndCommS4D_LMR_( array, imax, jmax, kmax, &nmax, vc, vc_comm, datatype, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // BndCommS3D
  *ierr = paraMngr->BndCommS3D( dtype, array, *imax, *jmax, *kmax, *vc, *vc_comm, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 袖通信(Vector3D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,3,nLeaf)の形式の配列の袖通信を行う
 *  - BndCommV3DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommV3D_LMR_( void *array, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm, int *datatype
               , int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 3;
  cpm_BndCommS4D_LMR_( array, imax, jmax, kmax, &nmax, vc, vc_comm, datatype, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // BndCommV3D
  *ierr = paraMngr->BndCommV3D( dtype, array, *imax, *jmax, *kmax, *vc, *vc_comm, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 袖通信(Scalar4DEx版)のFortranインターフェイス
 *  - (nmax,imax,jmax,kmax,nLeaf)の形式の配列の袖通信を行う
 *  - BndCommS4DExのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    nmax      配列サイズ(成分数)
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommS4DEx_LMR_( void *array, int *nmax, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm
                 , int *datatype, int *procGrpNo, int *ierr )
{
  if( !array || !imax || !jmax || !kmax || !nmax || !vc || !vc_comm || !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // BndCommS4DEx
   *ierr = paraMngr->BndCommS4DEx( dtype, array, *nmax, *imax, *jmax, *kmax, *vc, *vc_comm, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 袖通信(Vector3DEx版)のFortranインターフェイス
 *  - (3,imax,jmax,kmax,nLeaf)の形式の配列の袖通信を行う
 *  - BndCommV3DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommV3DEx_LMR_( void *array, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm, int *datatype
                 , int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 3;
  cpm_BndCommS4DEx_LMR_( array, &nmax, imax, jmax, kmax, vc, vc_comm, datatype, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // BndCommV3DEx
  *ierr = paraMngr->BndCommV3DEx( dtype, array, *imax, *jmax, *kmax, *vc, *vc_comm, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 周期境界袖通信(Scalar4D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,nmax,nLeaf)の形式の配列の周期境界方向の袖通信を行う
 *  - PeriodicCommS4DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    nmax      配列サイズ(成分数)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
 *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_PeriodicCommS4D_LMR_( void *array, int *imax, int *jmax, int *kmax, int *nmax, int *vc, int *vc_comm
                    , int *dir, int *pm, int *datatype, int *procGrpNo, int *ierr )
{
  if( !array || !imax || !jmax || !kmax || !nmax || !vc || !vc_comm || !dir || !pm ||
      !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // cpm_DirFlag
  if( *dir < 0 || *dir > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_DIR;
    return;
  }

  // cpm_PMFlag
  if( *pm < 0 || *pm > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_PM;
    return;
  }

  // PeriodicCommS4D
  *ierr = paraMngr->PeriodicCommS4D( dtype, array, *imax, *jmax, *kmax, *nmax, *vc, *vc_comm
                                   , (cpm_DirFlag)*dir, (cpm_PMFlag)*pm, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 周期境界袖通信(Scalar3D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,nLeaf)の形式の配列の周期境界方向の袖通信を行う
 *  - PeriodicCommS3DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
 *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_PeriodicCommS3D_LMR_( void *array, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm
                    , int *dir, int *pm, int *datatype, int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 1;
  cpm_PeriodicCommS4D_LMR_( array, imax, jmax, kmax, &nmax, vc, vc_comm, dir, pm
                      , datatype, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !dir || !pm ||
      !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // cpm_DirFlag
  if( *dir < 0 || *dir > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_DIR;
    return;
  }

  // cpm_PMFlag
  if( *pm < 0 || *pm > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_PM;
    return;
  }

  // PeriodicCommS3D
  *ierr = paraMngr->PeriodicCommS3D( dtype, array, *imax, *jmax, *kmax, *vc, *vc_comm
                                   , (cpm_DirFlag)*dir, (cpm_PMFlag)*pm, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 周期境界袖通信(Vector3D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,3,nLeaf)の形式の配列の周期境界方向の袖通信を行う
 *  - PeriodicCommV3DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
 *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_PeriodicCommV3D_LMR_( void *array, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm
                    , int *dir, int *pm, int *datatype, int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 3;
  cpm_PeriodicCommS4D_LMR_( array, imax, jmax, kmax, &nmax, vc, vc_comm, dir, pm
                      , datatype, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !dir || !pm ||
      !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // cpm_DirFlag
  if( *dir < 0 || *dir > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_DIR;
    return;
  }

  // cpm_PMFlag
  if( *pm < 0 || *pm > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_PM;
    return;
  }

  // PeriodicCommV3D
  *ierr = paraMngr->PeriodicCommV3D( dtype, array, *imax, *jmax, *kmax, *vc, *vc_comm
                                   , (cpm_DirFlag)*dir, (cpm_PMFlag)*pm, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 周期境界袖通信(Scalar4DEx版)のFortranインターフェイス
 *  - (nmax,imax,jmax,kmax,nLeaf)の形式の配列の周期境界方向の袖通信を行う
 *  - PeriodicCommS4DExのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    nmax      配列サイズ(成分数)
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
 *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_PeriodicCommS4DEx_LMR_( void *array, int *nmax, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm
                      , int *dir, int *pm, int *datatype, int *procGrpNo, int *ierr )
{
  if( !array || !imax || !jmax || !kmax || !nmax || !vc || !vc_comm || !dir || !pm ||
      !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // cpm_DirFlag
  if( *dir < 0 || *dir > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_DIR;
    return;
  }

  // cpm_PMFlag
  if( *pm < 0 || *pm > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_PM;
    return;
  }

  // PeriodicCommS4DEx
  *ierr = paraMngr->PeriodicCommS4DEx( dtype, array, *nmax, *imax, *jmax, *kmax, *vc, *vc_comm
                                     , (cpm_DirFlag)*dir, (cpm_PMFlag)*pm, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 周期境界袖通信(Vector3DEx版)のFortranインターフェイス
 *  - (3,imax,jmax,kmax,nLeaf)の形式の配列の周期境界方向の袖通信を行う
 *  - PeriodicCommV3DExのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
 *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_PeriodicCommV3DEx_LMR_( void *array, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm
                      , int *dir, int *pm, int *datatype, int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 3;
  cpm_PeriodicCommS4DEx_LMR_( array, &nmax, imax, jmax, kmax, vc, vc_comm, dir, pm
                        , datatype, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !dir || !pm ||
      !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManagerLMR *paraMngr = cpm_ParaManagerLMR::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_BaseParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // cpm_DirFlag
  if( *dir < 0 || *dir > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_DIR;
    return;
  }

  // cpm_PMFlag
  if( *pm < 0 || *pm > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_PM;
    return;
  }

  // PeriodicCommV3DEx
  *ierr = paraMngr->PeriodicCommV3DEx( dtype, array, *imax, *jmax, *kmax, *vc, *vc_comm
                                     , (cpm_DirFlag)*dir, (cpm_PMFlag)*pm, *procGrpNo );
#endif
}

