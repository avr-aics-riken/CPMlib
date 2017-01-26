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
 * @file   cpm_ParaManagerLMR.h
 * LMR用のパラレルマネージャクラスのヘッダーファイル
 * @date   2015/03/27
 */

#ifndef _CPM_PARAMANAGER_LMR_H_
#define _CPM_PARAMANAGER_LMR_H_

#include "cpm_DefineLMR.h"
#include "cpm_BaseParaManager.h"
#include "cpm_VoxelInfoLMR.h"
#include "cpm_LeafCommInfo.h"

/** プロセスグループ毎のVOXEL空間情報管理マップ */
//typedef std::map<int, cpm_VoxelInfoLMR*> LeafMap; //map<leafID,VoxelInfo*> -> cpm_VoxelInfoLMR.h
typedef std::map<int, LeafMap> VoxelInfoMapLMR; //map<procGrpID, LeafMap>

/** プロセスグループ内の袖通信情報マップ */
typedef std::map<int, cpm_LeafCommInfo*> LeafCommInfoMap; //map<distRankNo,cpm_LeafCommInfo*>

/** 全プロセスグループの袖通信情報マップ */
typedef std::map<int, LeafCommInfoMap> BndCommInfoMap;  //map<procGrpID,LeafCommInfoMap>


/** LMR用の並列管理クラス
 *  - 現時点ではユーザがインスタンスすることを許していない
 *  - get_instance静的関数を用いて唯一のインスタンスを取得する
 */
class cpm_ParaManagerLMR : public cpm_BaseParaManager
{
friend class cpm_BaseParaManager;

////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:

  /** インスタンスの取得
   * @return インスタンスのポインタ
   */
  static cpm_ParaManagerLMR* get_instance();

  /** インスタンスの取得(initialize処理も実行)
   *  @param[in] argc       プログラム実行時引数の数
   *  @param[in] argv       プログラム実行時引数
   *  @return インスタンスのポインタ
   */
  static cpm_ParaManagerLMR* get_instance(int &argc, char**& argv);

  /** LMR用の領域分割
   *  - FXgen出力の領域情報ファイル、木情報ファイルを渡して領域分割情報を生成する
   *
   *  @param[in] treeFile  木情報ファイル
   *  @param[in] maxVC     最大の袖数(袖通信用)
   *  @param[in] maxN      最大の成分数(袖通信用)
   *  @param[in] procGrpNo 領域分割を行うプロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  virtual
  cpm_ErrorCode VoxelInit_LMR( std::string treeFile
                             , size_t maxVC=1, size_t maxN=3, int procGrpNo=0 );

  /** 木情報ファイルからリーフ数を取得する
   *  @param[in] treeFile  木情報ファイル
   *  @return    リーフ数
   */
  static
  int GetNumLeaf( std::string treeFile );




////// 領域情報の取得関数 //////

  /** VOXEL空間マップを検索(0番目のリーフ情報を取得)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return VOXEL空間情報ポインタ
   */
  virtual const cpm_VoxelInfo* FindVoxelInfo( int procGrpNo=0 );

  /** VOXEL空間マップを検索(リーフID指定)
   *  @param[in] leafID リーフID
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return VOXEL空間情報ポインタ
   */
  const cpm_VoxelInfoLMR* FindLeafVoxelInfo_byID( int leafID, int procGrpNo=0 );

  /** VOXEL空間マップを検索(ランク内のリーフ順番号指定)
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return VOXEL空間情報ポインタ
   */
  const cpm_VoxelInfoLMR* FindLeafVoxelInfo( int leafIndex, int procGrpNo=0 );

  /** 全リーフ数を取得する
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return    リーフ数
   */
  int GetNumLeaf( int procGrpNo=0 );

  /** 自ランクが担当するリーフ数を取得する
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return    自ランクが担当するリーフ数
   */
  int GetLocalNumLeaf( int procGrpNo=0 );

  /** 自ランクが担当するリーフIDリストを取得する
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return    自ランクが担当するリーフ数
   */
  std::vector<int> GetLocalLeafIDs( int procGrpNo=0 );

  /** 指定リーフのリーフIDを取得する
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return リーフID
   */
  int GetLeafID( int leafIndex, int procGrpNo=0 );

  /** 自ランクが担当するリーフIDのインデクスを取得する
   *  @param[in] leafID    リーフID
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return    自ランク内での順番号(0～、IDが存在しない場合-1)
   */
  int GetLocalLeafIndex_byID( int leafID, int procGrpNo=0 );

  /** 領域分割数を取得(指定リーフのレベルにおける値)
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフのリーフレベルにおける領域分割数整数配列のポインタ
   */
  const int* GetDivNum( int leafIndex, int procGrpNo=0 );

  /** 指定リーフのピッチを取得
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフのピッチ実数配列のポインタ
   */
  const double* GetPitch( int leafIndex, int procGrpNo=0 );

  /** 指定リーフの空間原点を取得
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフの空間原点実数配列のポインタ
   */
  const double* GetLocalOrigin( int leafIndex, int procGrpNo=0 );

  /** 指定リーフの空間サイズを取得
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフの空間サイズ実数配列のポインタ
   */
  const double* GetLocalRegion( int leafIndex, int procGrpNo=0 );

  /** 指定リーフの領域分割位置を取得
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフの領域分割位置整数配列のポインタ
   */
  const int* GetDivPos( int leafIndex, int procGrpNo=0 );

  /** 指定リーフの始点VOXELの全体空間でのインデクスを取得
   *  - 全体空間の先頭インデクスを0としたC型のインデクス
   *
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフの始点インデクス整数配列のポインタ
   */
  const int* GetVoxelHeadIndex( int leafIndex, int procGrpNo=0 );

  /** 指定リーフの始点頂点の全体空間でのインデクスを取得
   *  - 全体空間の先頭インデクスを0としたC型のインデクス
   *
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフの始点インデクス整数配列のポインタ
   */
  const int* GetNodeHeadIndex( int leafIndex, int procGrpNo=0 );

  /** 指定リーフの始点VOXELまたは頂点の全体空間でのインデクスを取得
   *  - FVMのときはボクセル数、FDMのときは頂点数を取得
   *  - 全体空間の先頭インデクスを0としたC型のインデクス
   *
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフの始点インデクス整数配列のポインタ
   */
  const int* GetArrayHeadIndex( int leafIndex, int procGrpNo=0 );

  /** 指定リーフの終点VOXELの全体空間でのインデクスを取得
   *  - 全体空間の先頭インデクスを0としたC型のインデクス
   *
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフの終点インデクス整数配列のポインタ
   */
  const int* GetVoxelTailIndex( int leafIndex, int procGrpNo=0 );

  /** 指定リーフの終点頂点の全体空間でのインデクスを取得
   *  - 全体空間の先頭インデクスを0としたC型のインデクス
   *
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフの終点インデクス整数配列のポインタ
   */
  const int* GetNodeTailIndex( int leafIndex, int procGrpNo=0 );

  /** 指定リーフの終点VOXELまたは頂点の全体空間でのインデクスを取得
   *  - FVMのときはボクセル数、FDMのときは頂点数を取得
   *  - 全体空間の先頭インデクスを0としたC型のインデクス
   *
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフの終点インデクス整数配列のポインタ
   */
  const int* GetArrayTailIndex( int leafIndex, int procGrpNo=0 );

  /** 指定面における自リーフの隣接リーフ番号を取得
   *  @param[in]  leafIndex リーフ順番号(0~)
   *  @param[in]  face 面方向
   *  @param[out] num  面の数(0 or 1 or 4)
   *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定面における自リーフの隣接リーフ番号整数配列のポインタ
   */
  const int* GetNeighborLeafList( int leafIndex, cpm_FaceFlag face, int &num, int procGrpNo=0 );

  /** 指定リーフの指定面における自リーフの周期境界の隣接リーフ番号を取得
   *  @param[in]  leafIndex リーフ順番号(0~)
   *  @param[in]  face 面方向
   *  @param[out] num  面の数(0 or 1 or 4)
   *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフの指定面における自リーフの周期境界の隣接リーフ番号整数配列のポインタ
   */
  const int* GetPeriodicLeafList( int leafIndex, cpm_FaceFlag face, int &num, int procGrpNo=0 );

  /** 指定リーフの指定面における自リーフの隣接ランク番号を取得
   *  @param[in]  leafIndex リーフ順番号(0~)
   *  @param[in]  face 面方向
   *  @param[out] num  面の数(0 or 1 or 4)
   *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフの指定面における自リーフの隣接ランク番号整数配列のポインタ
   */
  const int* GetNeighborRankList( int leafIndex, cpm_FaceFlag face, int &num, int procGrpNo=0 );

  /** 指定リーフの指定面における自リーフの周期境界の隣接ランク番号を取得
   *  @param[in]  leafIndex リーフ順番号(0~)
   *  @param[in]  face 面方向
   *  @param[out] num  面の数(0 or 1 or 4)
   *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定リーフの指定面における自リーフの周期境界の隣接ランク番号整数配列のポインタ
   */
  const int* GetPeriodicRankList( int leafIndex, cpm_FaceFlag face, int &num, int procGrpNo=0 );

  /** 指定リーフの指定面におけるレベル差を取得
   *  @param[in]  leafIndex リーフ順番号(0~)
   *  @param[in]  face 面方向
   *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
   *  @return     レベル差(0:同じレベル, 1:fine, -1:coarse)
   */
  int GetNeighborLevelDiff( int leafIndex, cpm_FaceFlag face, int procGrpNo=0 );

  /** 指定リーフの境界が外部境界かどうかを判定
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] face      面方向
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @retval    true      外部境界
   *  @retval    false     外部境界でない
   */
  bool IsOuterBoundary( int leafIndex, cpm_FaceFlag face, int procGrpNo=0 );

  /** 指定リーフの境界が内部境界(隣が不活性ドメイン)かどうかを判定
   *  @param[in] leafIndex リーフ順番号(0~)
   *  @param[in] face      面方向
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @retval    true      内部境界
   *  @retval    false     内部境界でない
   */
  bool IsInnerBoundary( int leafIndex, cpm_FaceFlag face, int procGrpNo=0 );






////// 袖通信関数 //////

  /** 袖通信バッファサイズの取得
   *  - 袖通信バッファとして確保されている配列サイズ(byte)を返す
   *
   *  @param[in] procGrpNo プロセスグループ番号(負の場合、全プロセスグループでのトータルを返す)
   *  @return バッファサイズ(byte)
   */
  virtual
  size_t GetBndCommBufferSize( int procGrpNo=0 ); 

  /** 袖通信バッファのセット(LMR用の袖通信情報も生成する)
   *  @param[in] maxVC     送受信バッファの最大袖数
   *  @param[in] maxN      送受信バッファの最大成分数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode
  SetBndCommBuffer( size_t maxVC, size_t maxN, int procGrpNo=0 );

  /** 袖通信(Scalar3D版)
   *  - (imax,jmax,kmax,nLeaf)の形式の配列の袖通信を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                          , int procGrpNo=0 ); 

  /** 袖通信(Scalar3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,nLeaf)の形式の配列の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                          , int vc, int vc_comm, int procGrpNo=0 ); 

  /** 袖通信(Vector3D版)
   *  - (imax,jmax,kmax,3,nLeaf)の形式の配列の袖通信を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                          , int procGrpNo=0 ); 

  /** 袖通信(Vector3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,3,nLeaf)の形式の配列の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                          , int vc, int vc_comm, int procGrpNo=0 );
 
  /** 袖通信(Vector3DEx版)
   *  - (3,imax,jmax,kmax,nLeaf)の形式の配列の袖通信を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm, int procGrpNo=0 ); 

  /** 袖通信(Vector3DEx版, MPI_Datatype指定)
   *  - (3,imax,jmax,kmax,nLeaf)の形式の配列の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *   @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                            , int vc, int vc_comm, int procGrpNo=0 );

  /** 非同期版袖通信(Scalar3D版)
   *  - (imax,jmax,kmax,nLeaf)の形式の配列の非同期袖通信を行う
   *  - waitと展開は行わない
   *  - wait、展開はwait_BndCommS3Dをコールする
   *
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS3D_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , int procGrpNo=0 );

  /** 非同期版袖通信(Scalar3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,nLeaf)の形式の配列の非同期袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  - waitと展開は行わない
   *  - wait、展開はwait_BndCommS3Dをコールする
   *
   *  @param[in]   dtype     袖通信データのMPI_Datatype
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS3D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, int procGrpNo=0 ); 


  /** 非同期版袖通信(Vector3D版)
   *  - (imax,jmax,kmax,3,nLeaf)の形式の配列の非同期袖通信を行う
   *  - waitと展開は行わない
   *  - wait、展開はwait_BndCommV3Dをコールする
   *
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommV3D_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , int procGrpNo=0 ); 

  /** 非同期版袖通信(Vector3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,3,nLeaf)の形式の配列の非同期袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  - waitと展開は行わない
   *  - wait、展開はwait_BndCommV3Dをコールする
   *
   *  @param[in]   dtype     袖通信データのMPI_Datatype
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommV3D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, int procGrpNo=0 ); 

  /** 非同期版袖通信(Vector3DEx版)
   *  - (3,imax,jmax,kmax,nLeaf)の形式の配列の非同期袖通信を行う
   *  - waitと展開は行わない
   *  - wait、展開はwait_BndCommV3DExをコールする
   *
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommV3DEx_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                   , int procGrpNo=0 );

  /** 非同期版袖通信(Vector3DEx版, MPI_Datatype指定)
   *  - (3,imax,jmax,kmax,nLeaf)の形式の配列の非同期袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  - waitと展開は行わない
   *  - wait、展開はwait_BndCommV3DExをコールする
   *
   *  @param[in]   dtype     袖通信データのMPI_Datatype
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommV3DEx_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                   , int vc, int vc_comm, int procGrpNo=0 );


  /** 非同期版袖通信のwait、展開(Scalar3D版)
   *  - (imax,jmax,kmax,nLeaf)の形式の配列の非同期版袖通信のwaitと展開を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                               , int procGrpNo=0 );

  /** 非同期版袖通信のwait、展開(Vector3D版)
   *  - (imax,jmax,kmax,3,nLeaf)の形式の配列の非同期版袖通信のwaitと展開を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                               , int procGrpNo=0 ); 

 /** 非同期版袖通信のwait、展開(Vector3D版, MPI_Datatype指定)
  *  - (imax,jmax,kmax,3,nLeaf)の形式の配列の非同期版袖通信のwaitと展開を行う
  *  - MPI_Datatypeを指定するバージョン
  *
  *  @param[in]    dtype     袖通信データのMPI_Datatype
  *  @param[inout] array     袖通信をする配列の先頭ポインタ
  *  @param[in]    imax      配列サイズ(I方向)
  *  @param[in]    jmax      配列サイズ(J方向)
  *  @param[in]    kmax      配列サイズ(K方向)
  *  @param[in]    vc        仮想セル数
  *  @param[in]    vc_comm   通信する仮想セル数
  *  @param[in]    procGrpNo プロセスグループ番号
  *  @return 終了コード(CPM_SUCCESS=正常終了)
  */
  cpm_ErrorCode wait_BndCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                               , int vc, int vc_comm, int procGrpNo=0 ); 

  /** 非同期版袖通信のwait、展開(Vector3DEx版)
   *  - (3,imax,jmax,kmax,nLeaf)の形式の配列の非同期版袖通信のwaitと展開を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , int procGrpNo=0 );

  /** 周期境界袖通信(Scalar3D版)
   *  - (imax,jmax,kmax,nLeaf)の形式の配列の周期境界方向の袖通信を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                               , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 ); 

  /** 周期境界袖通信(Scalar3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,nLeaf)の形式の配列の周期境界方向の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 ); 

  /** 周期境界袖通信(Vector3D版)
   *  - (imax,jmax,kmax,3,nLeaf)の形式の配列の周期境界方向の袖通信を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                               , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 ); 

  /** 周期境界袖通信(Vector3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,3,nLeaf)の形式の配列の周期境界方向の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 ); 

  /** 周期境界袖通信(Vector3DEx版)
   *  - (3,imax,jmax,kmax,nLeaf)の形式の配列の周期境界方向の袖通信を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 ); 

  /** 周期境界袖通信(Vector3DEx版, MPI_Datatype指定)
   *  - (3,imax,jmax,kmax,nLeaf)の形式の配列の周期境界方向の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );
 

// 2016/01/22 FEAST add.e 
 

  /** 袖通信(Scalar4D版)
   *  - (imax,jmax,kmax,nmax,nLeaf)の形式の配列の袖通信を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                          , int procGrpNo=0 );

// 2016/01/22 FEAST add.s

  /** 袖通信(Scalar4D版)
   *  - (imax,jmax,kmax,nmax,nLeaf)の形式の配列の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                          , int vc, int vc_comm, int procGrpNo=0 );

  /** 非同期版袖通信(Scalar4D版)
   *  - (imax,jmax,kmax,nmax,nLeaf)の形式の配列の非同期袖通信を行う
   *  - waitと展開は行わない
   *  - wait、展開はwait_BndCommS4Dをコールする
   *
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4D_nowait( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                 , int procGrpNo=0 );

  /** 非同期版袖通信(Scalar4D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,nmax,nLeaf)の形式の配列の非同期袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  - waitと展開は行わない
   *  - wait、展開はwait_BndCommS4Dをコールする
   *
   *  @param[in]   dtype     袖通信データのMPI_Datatype
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                                 , int vc, int vc_comm, int procGrpNo=0 ); 


  /** 非同期版袖通信のwait、展開(Scalar4D版)
   *  - (imax,jmax,kmax,nmax,nLeaf)の形式の配列の非同期版袖通信のwaitと展開を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                               , int procGrpNo=0 );

  /** 周期境界袖通信(Scalar4D版)
   *  - (imax,jmax,kmax,nmax,nLeaf)の形式の配列の周期境界方向の袖通信を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                               , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );

  /** 周期境界袖通信(Scalar4D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,nmax,nLeaf)の形式の配列の周期境界方向の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 ); 

  /** 袖通信(Scalar4DEx版)
   *  - (nmax,imax,jmax,kmax,nLeaf)の形式の配列の袖通信を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                            , int procGrpNo=0 );

  /** 袖通信(Scalar4DEx版, MPI_Datatype指定)
   *  - (nmax,imax,jmax,kmax,nLeaf)の形式の配列の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                            , int vc, int vc_comm, int procGrpNo=0 );

  /** 非同期版袖通信(Scalar4DEx版)
   *  - (nmax,imax,jmax,kmax,nLeaf)の形式の配列の非同期袖通信を行う
   *  - waitと展開は行わない
   *  - wait、展開はwait_BndCommS4DExをコールする
   *
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4DEx_nowait( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                   , int procGrpNo=0 );

  /** 非同期版袖通信(Scalar4DEx版, MPI_Datatype指定)
   *  - (nmax,imax,jmax,kmax,nLeaf)の形式の配列の非同期袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  - waitと展開は行わない
   *  - wait、展開はwait_BndCommS4DExをコールする
   *
   *  @param[in]   dtype     袖通信データのMPI_Datatype
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4DEx_nowait( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                   , int vc, int vc_comm, int procGrpNo=0 ); 

  /** 非同期版袖通信のwait、展開(Scalar4DEx版)
   *  - (nmax,imax,jmax,kmax,nLeaf)の形式の配列の非同期版袖通信のwaitと展開を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , int procGrpNo=0 );

  /** 周期境界袖通信(Scalar4DEx版)
   *  - (nmax,imax,jmax,kmax,nLeaf)の形式の配列の周期境界方向の袖通信を行う
   *
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );

  /** 周期境界袖通信(Scalar4DEx版, MPI_Datatype指定)
   *  - (nmax,imax,jmax,kmax,nLeaf)の形式の配列の周期境界方向の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );




protected:

  /** コンストラクタ */
  cpm_ParaManagerLMR();

  /** デストラクタ */
  virtual ~cpm_ParaManagerLMR();

  /** 配列確保(double)
   *  @param[in] nmax 成分数
   *  @param[in] sz   配列サイズ
   *  @param[in] vc   仮想セル数
   *  @param[in] procGrpNo  領域分割を行うプロセスグループ番号
   *  @return         確保した配列のポインタ
   */
  virtual
  double* AllocDouble( int nmax, int sz[3], int vc, int procGrpNo );

  /** 配列確保(float)
   *  @param[in] nmax 成分数
   *  @param[in] sz   配列サイズ
   *  @param[in] vc   仮想セル数
   *  @param[in] procGrpNo  領域分割を行うプロセスグループ番号
   *  @return         確保した配列のポインタ
   */
  virtual
  float* AllocFloat( int nmax, int sz[3], int vc, int procGrpNo );

  /** 配列確保(int)
   *  @param[in] nmax 成分数
   *  @param[in] sz   配列サイズ
   *  @param[in] vc   仮想セル数
   *  @param[in] procGrpNo  領域分割を行うプロセスグループ番号
   *  @return         確保した配列のポインタ
   */
  virtual
  int* AllocInt( int nmax, int sz[3], int vc, int procGrpNo );

  /** １方向の非同期受信処理
   *  @param[in]  commInfoMap 通信情報マップ
   *  @param[in]  sz_face     面内の格子数
   *  @param[in]  nmax        成分数
   *  @param[in]  vc_comm     通信層数
   *  @param[in]  bPeriodic   周期境界フラグ(true:周期境界通信、false:内部袖通信のみ)
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode recv_LMR( LeafCommInfoMap &commInfoMap
                        , size_t sz_face[2], int nmax, int vc_comm
                        , bool bPeriodic, int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の１面の送信データのパックと送信
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfoMap 通信情報マップ
   *  @param[in]  bPeriodic   周期境界フラグ(true:周期境界通信、false:内部袖通信のみ)
   *  @param[in]  face        送信方向
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode send_LMR( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , LeafCommInfoMap &commInfoMap, bool bPeriodic, cpm_FaceFlag face
                        , int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)のランク内コピー処理
   *  @param[inout] array        袖通信をする配列の先頭ポインタ
   *  @param[in]    imax         配列サイズ(I方向)
   *  @param[in]    jmax         配列サイズ(J方向)
   *  @param[in]    kmax         配列サイズ(K方向)
   *  @param[in]    nmax         配列サイズ(成分数)
   *  @param[in]    vc           仮想セル数
   *  @param[in]    vc_comm      通信する仮想セル数
   *  @param[in]    commInfoMapM 通信情報マップ(マイナス側)
   *  @param[in]    commInfoMapP 通信情報マップ(プラス側)
   *  @param[in]    bPeriodic    周期境界フラグ(true:周期境界通信、false:内部袖通信のみ)
   *  @param[in]    dir          通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm           通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo    プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode copy_LMR( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , LeafCommInfoMap &commInfoMapM, LeafCommInfoMap &commInfoMapP
                        , bool bPeriodic, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の１面の受信待機とデータの展開
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfoMap 通信情報マップ
   *  @param[in]  bPeriodic   周期境界フラグ(true:周期境界通信、false:内部袖通信のみ)
   *  @param[in]  face        受信方向
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode recv_LMR_wait( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                             , LeafCommInfoMap &commInfoMap, bool bPeriodic, cpm_FaceFlag face
                             , int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の１面の送信待機
   *  @param[in]  commInfoMap 通信情報マップ
   */
  template<class T>
  cpm_ErrorCode send_LMR_wait( LeafCommInfoMap &commInfoMap );

  /** 袖通信(Scalar3D,4D,Vector3D版)の-X面への送信データのパック(通信面毎)
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfo    リーフ間の通信情報
   *  @param[out] sendbuf     送信バッファ 
   *  @param[in]  nw          送信バッファサイズ
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode packMX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                       , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw
                       , int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の+X面への送信データのパック(通信面毎)
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfo    リーフ間の通信情報
   *  @param[out] sendbuf     送信バッファ 
   *  @param[in]  nw          送信バッファサイズ
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode packPX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                       , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw
                       , int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の-Y面への送信データのパック(通信面毎)
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfo    リーフ間の通信情報
   *  @param[out] sendbuf     送信バッファ 
   *  @param[in]  nw          送信バッファサイズ
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode packMY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                       , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw
                       , int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の+Y面への送信データのパック(通信面毎)
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfo    リーフ間の通信情報
   *  @param[out] sendbuf     送信バッファ 
   *  @param[in]  nw          送信バッファサイズ
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode packPY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                       , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw
                       , int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の-Z面への送信データのパック(通信面毎)
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfo    リーフ間の通信情報
   *  @param[out] sendbuf     送信バッファ 
   *  @param[in]  nw          送信バッファサイズ
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode packMZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                       , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw
                       , int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の+Z面への送信データのパック(通信面毎)
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfo    リーフ間の通信情報
   *  @param[out] sendbuf     送信バッファ 
   *  @param[in]  nw          送信バッファサイズ
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode packPZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                       , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw
                       , int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の-X面からの受信データの展開(通信面毎)
   *  @param[inout] array       袖通信をする配列の先頭ポインタ
   *  @param[in]    imax        配列サイズ(I方向)
   *  @param[in]    jmax        配列サイズ(J方向)
   *  @param[in]    kmax        配列サイズ(K方向)
   *  @param[in]    nmax        配列サイズ(成分数)
   *  @param[in]    vc          仮想セル数
   *  @param[in]    vc_comm     通信する仮想セル数
   *  @param[in]    commInfo    リーフ間の通信情報
   *  @param[in]    recvbuf     受信バッファ 
   *  @param[in]    procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode unpackMX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の+X面からの受信データの展開(通信面毎)
   *  @param[inout] array       袖通信をする配列の先頭ポインタ
   *  @param[in]    imax        配列サイズ(I方向)
   *  @param[in]    jmax        配列サイズ(J方向)
   *  @param[in]    kmax        配列サイズ(K方向)
   *  @param[in]    nmax        配列サイズ(成分数)
   *  @param[in]    vc          仮想セル数
   *  @param[in]    vc_comm     通信する仮想セル数
   *  @param[in]    commInfo    リーフ間の通信情報
   *  @param[in]    recvbuf     受信バッファ 
   *  @param[in]    procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode unpackPX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の-Y面からの受信データの展開(通信面毎)
   *  @param[inout] array       袖通信をする配列の先頭ポインタ
   *  @param[in]    imax        配列サイズ(I方向)
   *  @param[in]    jmax        配列サイズ(J方向)
   *  @param[in]    kmax        配列サイズ(K方向)
   *  @param[in]    nmax        配列サイズ(成分数)
   *  @param[in]    vc          仮想セル数
   *  @param[in]    vc_comm     通信する仮想セル数
   *  @param[in]    commInfo    リーフ間の通信情報
   *  @param[in]    recvbuf     受信バッファ 
   *  @param[in]    procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode unpackMY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の+Y面からの受信データの展開(通信面毎)
   *  @param[inout] array       袖通信をする配列の先頭ポインタ
   *  @param[in]    imax        配列サイズ(I方向)
   *  @param[in]    jmax        配列サイズ(J方向)
   *  @param[in]    kmax        配列サイズ(K方向)
   *  @param[in]    nmax        配列サイズ(成分数)
   *  @param[in]    vc          仮想セル数
   *  @param[in]    vc_comm     通信する仮想セル数
   *  @param[in]    commInfo    リーフ間の通信情報
   *  @param[in]    recvbuf     受信バッファ 
   *  @param[in]    procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode unpackPY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の-Z面からの受信データの展開(通信面毎)
   *  @param[inout] array       袖通信をする配列の先頭ポインタ
   *  @param[in]    imax        配列サイズ(I方向)
   *  @param[in]    jmax        配列サイズ(J方向)
   *  @param[in]    kmax        配列サイズ(K方向)
   *  @param[in]    nmax        配列サイズ(成分数)
   *  @param[in]    vc          仮想セル数
   *  @param[in]    vc_comm     通信する仮想セル数
   *  @param[in]    commInfo    リーフ間の通信情報
   *  @param[in]    recvbuf     受信バッファ 
   *  @param[in]    procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode unpackMZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の+Z面からの受信データの展開(通信面毎)
   *  @param[inout] array       袖通信をする配列の先頭ポインタ
   *  @param[in]    imax        配列サイズ(I方向)
   *  @param[in]    jmax        配列サイズ(J方向)
   *  @param[in]    kmax        配列サイズ(K方向)
   *  @param[in]    nmax        配列サイズ(成分数)
   *  @param[in]    vc          仮想セル数
   *  @param[in]    vc_comm     通信する仮想セル数
   *  @param[in]    commInfo    リーフ間の通信情報
   *  @param[in]    recvbuf     受信バッファ 
   *  @param[in]    procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode unpackPZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の１面の送信データのパックと送信
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfoMap 通信情報マップ
   *  @param[in]  bPeriodic   周期境界フラグ(true:周期境界通信、false:内部袖通信のみ)
   *  @param[in]  face        送信方向
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode send_LMR_Ex( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                           , LeafCommInfoMap &commInfoMap, bool bPeriodic, cpm_FaceFlag face
                           , int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)のランク内コピー処理
   *  @param[inout] array        袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax         配列サイズ(成分数)
   *  @param[in]    imax         配列サイズ(I方向)
   *  @param[in]    jmax         配列サイズ(J方向)
   *  @param[in]    kmax         配列サイズ(K方向)
   *  @param[in]    vc           仮想セル数
   *  @param[in]    vc_comm      通信する仮想セル数
   *  @param[in]    commInfoMapM 通信情報マップ(マイナス側)
   *  @param[in]    commInfoMapP 通信情報マップ(プラス側)
   *  @param[in]    bPeriodic    周期境界フラグ(true:周期境界通信、false:内部袖通信のみ)
   *  @param[in]    dir          通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm           通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo    プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode copy_LMR_Ex( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                           , LeafCommInfoMap &commInfoMapM, LeafCommInfoMap &commInfoMapP
                           , bool bPeriodic, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の１面の受信待機とデータの展開
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfoMap 通信情報マップ
   *  @param[in]  bPeriodic   周期境界フラグ(true:周期境界通信、false:内部袖通信のみ)
   *  @param[in]  face        受信方向
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode recv_LMR_Ex_wait( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                , LeafCommInfoMap &commInfoMap, bool bPeriodic, cpm_FaceFlag face
                                , int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の-X面への送信データのパック(通信面毎)
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfo    リーフ間の通信情報
   *  @param[out] sendbuf     送信バッファ 
   *  @param[in]  nw          送信バッファサイズ
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode packMXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw
                        , int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の+X面への送信データのパック(通信面毎)
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfo    リーフ間の通信情報
   *  @param[out] sendbuf     送信バッファ 
   *  @param[in]  nw          送信バッファサイズ
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode packPXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw
                        , int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の-Y面への送信データのパック(通信面毎)
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfo    リーフ間の通信情報
   *  @param[out] sendbuf     送信バッファ 
   *  @param[in]  nw          送信バッファサイズ
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode packMYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw
                        , int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の+Y面への送信データのパック(通信面毎)
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfo    リーフ間の通信情報
   *  @param[out] sendbuf     送信バッファ 
   *  @param[in]  nw          送信バッファサイズ
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode packPYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw
                        , int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の-Z面への送信データのパック(通信面毎)
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfo    リーフ間の通信情報
   *  @param[out] sendbuf     送信バッファ 
   *  @param[in]  nw          送信バッファサイズ
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode packMZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw
                        , int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の+Z面への送信データのパック(通信面毎)
   *  @param[in]  array       袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax        配列サイズ(成分数)
   *  @param[in]  imax        配列サイズ(I方向)
   *  @param[in]  jmax        配列サイズ(J方向)
   *  @param[in]  kmax        配列サイズ(K方向)
   *  @param[in]  vc          仮想セル数
   *  @param[in]  vc_comm     通信する仮想セル数
   *  @param[in]  commInfo    リーフ間の通信情報
   *  @param[out] sendbuf     送信バッファ 
   *  @param[in]  nw          送信バッファサイズ
   *  @param[in]  procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode packPZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw
                        , int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の-X面からの受信データの展開(通信面毎)
   *  @param[inout] array       袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax        配列サイズ(成分数)
   *  @param[in]    imax        配列サイズ(I方向)
   *  @param[in]    jmax        配列サイズ(J方向)
   *  @param[in]    kmax        配列サイズ(K方向)
   *  @param[in]    vc          仮想セル数
   *  @param[in]    vc_comm     通信する仮想セル数
   *  @param[in]    commInfo    リーフ間の通信情報
   *  @param[in]    recvbuf     受信バッファ 
   *  @param[in]    procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode unpackMXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の+X面からの受信データの展開(通信面毎)
   *  @param[inout] array       袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax        配列サイズ(成分数)
   *  @param[in]    imax        配列サイズ(I方向)
   *  @param[in]    jmax        配列サイズ(J方向)
   *  @param[in]    kmax        配列サイズ(K方向)
   *  @param[in]    vc          仮想セル数
   *  @param[in]    vc_comm     通信する仮想セル数
   *  @param[in]    commInfo    リーフ間の通信情報
   *  @param[in]    recvbuf     受信バッファ 
   *  @param[in]    procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode unpackPXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の-Y面からの受信データの展開(通信面毎)
   *  @param[inout] array       袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax        配列サイズ(成分数)
   *  @param[in]    imax        配列サイズ(I方向)
   *  @param[in]    jmax        配列サイズ(J方向)
   *  @param[in]    kmax        配列サイズ(K方向)
   *  @param[in]    vc          仮想セル数
   *  @param[in]    vc_comm     通信する仮想セル数
   *  @param[in]    commInfo    リーフ間の通信情報
   *  @param[in]    recvbuf     受信バッファ 
   *  @param[in]    procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode unpackMYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の+Y面からの受信データの展開(通信面毎)
   *  @param[inout] array       袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax        配列サイズ(成分数)
   *  @param[in]    imax        配列サイズ(I方向)
   *  @param[in]    jmax        配列サイズ(J方向)
   *  @param[in]    kmax        配列サイズ(K方向)
   *  @param[in]    vc          仮想セル数
   *  @param[in]    vc_comm     通信する仮想セル数
   *  @param[in]    commInfo    リーフ間の通信情報
   *  @param[in]    recvbuf     受信バッファ 
   *  @param[in]    procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode unpackPYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の-Z面からの受信データの展開(通信面毎)
   *  @param[inout] array       袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax        配列サイズ(成分数)
   *  @param[in]    imax        配列サイズ(I方向)
   *  @param[in]    jmax        配列サイズ(J方向)
   *  @param[in]    kmax        配列サイズ(K方向)
   *  @param[in]    vc          仮想セル数
   *  @param[in]    vc_comm     通信する仮想セル数
   *  @param[in]    commInfo    リーフ間の通信情報
   *  @param[in]    recvbuf     受信バッファ 
   *  @param[in]    procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode unpackMZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo=0 );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の+Z面からの受信データの展開(通信面毎)
   *  @param[inout] array       袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax        配列サイズ(成分数)
   *  @param[in]    imax        配列サイズ(I方向)
   *  @param[in]    jmax        配列サイズ(J方向)
   *  @param[in]    kmax        配列サイズ(K方向)
   *  @param[in]    vc          仮想セル数
   *  @param[in]    vc_comm     通信する仮想セル数
   *  @param[in]    commInfo    リーフ間の通信情報
   *  @param[in]    recvbuf     受信バッファ 
   *  @param[in]    procGrpNo   プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode unpackPZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo=0 );




////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:

protected:

  /** プロセスグループ毎のVOXEL空間情報マップ
   *  - VOXEL空間番号をキーとしたVOXEL空間情報マップ
   *  - 自ランクが含まれるVOXEL空間のみを管理する
   *  - １プロセス複数リーフに対応
   */
  VoxelInfoMapLMR m_voxelInfoMap;

  /** -X方向袖通信情報 */
  BndCommInfoMap m_bndCommInfoMapMX;

  /** -Y方向袖通信情報 */
  BndCommInfoMap m_bndCommInfoMapMY;

  /** -Z方向袖通信情報 */
  BndCommInfoMap m_bndCommInfoMapMZ;

  /** +X方向袖通信情報 */
  BndCommInfoMap m_bndCommInfoMapPX;

  /** +Y方向袖通信情報 */
  BndCommInfoMap m_bndCommInfoMapPY;

  /** +Z方向袖通信情報 */
  BndCommInfoMap m_bndCommInfoMapPZ;

};

//インライン関数
#include "inline/cpm_ParaManagerLMR_BndComm.h"
#include "inline/cpm_ParaManagerLMR_BndCommEx.h"

#endif /* _CPM_PARAMANAGER_LMR_H_ */
