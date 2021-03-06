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
 * @file   cpm_BaseParaManager.h
 * パラレルマネージャ基底クラスのヘッダーファイル
 * @date   2012/05/31
 */

#ifndef _CPM_BASEPARAMANAGER_H_
#define _CPM_BASEPARAMANAGER_H_

#include "cpm_Base.h"
#include <map>
#include <vector>
#include <typeinfo>
#include "cpm_DomainInfo.h"
#include "cpm_VoxelInfo.h"
#include "cpm_ObjList.h"
#include <string.h> // for memset()

/** プロセスグループ毎の定義点タイプ管理マップ */
typedef std::map<int, cpm_DefPointType> DefPointMap;

/** CPMの並列管理クラス
 *  - 現時点ではユーザがインスタンスすることを許していない
 *  - get_instance静的関数を用いて唯一のインスタンスを取得する
 */
class cpm_BaseParaManager : public cpm_Base
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:

  /** 初期化処理(MPI_Initは実行済みの場合)
   *  - MPI_Initは既に実行済みである必要がある
   *  - 並列数、自ランク番号を取得
   *
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Initialize();

  /** 初期化処理(MPI_Initも実行する)
   *  - MPI_Initが実行されていない場合、実行する
   *  - 並列数、自ランク番号を取得
   *
   *  @param[in] argc プログラム実行時引数の数
   *  @param[in] argv プログラム実行時引数
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Initialize(int &argc, char**& argv);

  /** 並列実行であるかチェックする
   *  - 並列実行であっても、並列数が1のときはfalseとなる
   *
   *  @retval true  並列実行
   *  @retval false 逐次実行
   */
  bool IsParallel();

  /** 並列実行であるかチェックする(const)
   *  - 並列実行であっても、並列数が1のときはfalseとなる
   *
   *  @retval true  並列実行
   *  @retval false 逐次実行
   */
  bool IsParallel() const;

  /** プロセスグループの作成
   *  - 指定されたプロセスリストを使用してプロセスグループを生成する
   *
   *  @param[in]  nproc           使用するプロセスの数
   *  @param[in]  proclist        使用するプロセスのリスト(親プロセスグループでのランク番号)
   *  @param[in]  parentProcGrpNo 親とするプロセスグループ番号(省略時0)
   *  @retval 0以上 生成されたプロセスグループ番号
   *  @retval -1    エラー
   */
  int CreateProcessGroup( int nproc, int *proclist, int parentProcGrpNo=0 );





////// 領域情報の取得関数 //////

  /** 領域分割タイプを取得
   *  @return 領域分割タイプ
   */
  cpm_DomainType GetDomainType();

  /** 定義点タイプを取得
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 定義点情報(enum)
   */
  cpm_DefPointType GetDefPointType( int procGrpNo=0 );

  /** VOXEL空間マップを検索 \n
   *  (LMRのときは自プロセス内先頭リーフのVoxelInfoを返す)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return VOXEL空間情報ポインタ
   */
  virtual
  const cpm_VoxelInfo* FindVoxelInfo( int procGrpNo=0 ) = 0;

  /** 全体ボクセル数を取得 \n
   *  (LMRのときは最大レベルでの全体ボクセル数を返す)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 全体ボクセル数の整数配列ポインタ(3word)
   */
  const int* GetGlobalVoxelSize( int procGrpNo=0 );

  /** 全体頂点数を取得 \n
   *  (LMRのときは最大レベルでの全体頂点数を返す)
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 全体頂点数の整数配列ポインタ(3word)
   */
  const int* GetGlobalNodeSize( int procGrpNo=0 );

  /** 全体ボクセル数または頂点数を取得
   *  - FVMのときはボクセル数、FDMのときは頂点数を取得
   *  (LMRのときは最大レベルでの数を返す)
   *
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 全体ボクセル数または頂点数の整数配列ポインタ(3word)
   */
  const int* GetGlobalArraySize( int procGrpNo=0 );

  /** 全体空間の原点を取得
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 全体空間の原点実数配列のポインタ
   */
  const double* GetGlobalOrigin( int procGrpNo=0 );

  /** 全体空間サイズを取得
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 全体空間サイズ実数配列のポインタ
   */
  const double* GetGlobalRegion( int procGrpNo=0 );

  /** 自ランクのボクセル数を取得
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return ローカルボクセル数の整数配列ポインタ(3word)
   */
  const int* GetLocalVoxelSize( int procGrpNo=0 );

  /** 自ランクの頂点数を取得
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return ローカル頂点数の整数配列ポインタ(3word)
   */
  const int* GetLocalNodeSize( int procGrpNo=0 );

  /** 自ランクのボクセル数または頂点数を取得
   *  - FVMのときはボクセル数、FDMのときは頂点数を取得
   *
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return ローカル頂点数の整数配列ポインタ(3word)
   */
  const int* GetLocalArraySize( int procGrpNo=0 );




////// MPI処理のインターフェイス関数 //////

  /** MPI_Datatypeを取得
   *  @param[in] ptr 取得したいデータのポインタ
   *  @return MPI_Datatype
   */
  template<class T> CPM_INLINE
  static MPI_Datatype GetMPI_Datatype(T *ptr);

  /** MPI_Datatypeを取得
   *  - FortranデータタイプからMPI_Datatypeを取得
   *
   *  @param[in] datatype 取得したいデータのポインタ
   *  @return MPI_Datatype
   */
  static MPI_Datatype GetMPI_Datatype(int datatype);

  /** MPI_Opを取得
   *  - FortranオペレータタイプからMPI_Opを取得
   *
   *  @param[in] op 取得したいデータのポインタ
   *  @return MPI_Op
   */
  static MPI_Op GetMPI_Op(int op);

  /** ランク番号の取得
   *  - MPI_PROC_NULLが返ってきた場合は、\n
   *    1.プロセスグループが存在しない、 \n
   *    2.プロセスグループに自ランクが含まれていない、\n
   *    のいずれか
   *
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return ランク番号
   */
  int GetMyRankID( int procGrpNo=0 );

  /** ランク数の取得
   *  - プロセスグループのランク数を取得する
   *
   *  @param[in] procGrpNo プロセスグループ番号(省略時0)
   *  @return ランク数
   */
  int GetNumRank( int procGrpNo=0 );

  /** ホスト名の取得
   *  - 自ランクのホスト名を取得
   *
   *  @return ホスト名
   */
  std::string GetHostName();

  /** MPIコミュニケータの取得
   *  - MPI_COMM_NULLが返ってきた場合は、\n
   *    1.プロセスグループが存在しない、 \n
   *    2.プロセスグループに自ランクが含まれていない、\n
   *    のいずれか
   *
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return MPIコミュニケータ
   */
  MPI_Comm GetMPI_Comm( int procGrpNo=0 );

  /** Abort
   *  - MPI_Abortのインターフェイス
   *
   *  @param[in] errorcode MPI_Abortに渡すエラーコード
   */
  void Abort( int errorcode );

  /** Barrier
   *  - MPI_Barrierのインターフェイス
   *
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Barrier( int procGrpNo=0 );

  /** Wait
   *  - MPI_Waitのインターフェイス
   *
   *  @param[in] request リクエストハンドル
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Wait( MPI_Request *request );

  /** Waitall
   *  - MPI_Waitallのインターフェイス
   *
   *  @param[in] count    リクエストの数
   *  @param[in] requests リクエストハンドル配列
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Waitall( int count, MPI_Request requests[] );

  /** Bcast
   *  - MPI_Bcastのインターフェイス
   *
   *  @param[inout] buf       送受信バッファ
   *  @param[in]    count     送信バッファのサイズ(ワード数)
   *  @param[in]    root      送信元のランク番号(procGrpNo内でのランク番号)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode Bcast( T *buf, int count, int root, int procGrpNo=0 );

  /** Bcast
   *  - MPI_Bcastのインターフェイス
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]    dtype     送信バッファのMPI_Datatype
   *  @param[inout] buf       送受信バッファ
   *  @param[in]    count     送信バッファのサイズ(ワード数)
   *  @param[in]    root      送信元のランク番号(procGrpNo内でのランク番号)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Bcast( MPI_Datatype dtype, void *buf, int count, int root
                     , int procGrpNo=0 );

  /** Send
   *  - MPI_Sendのインターフェイス
   *
   *  @param[in] buf       送信データ
   *  @param[in] count     送信データのサイズ
   *  @param[in] dest      送信先のランク番号(procGrpNo内でのランク番号)
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode Send( T *buf, int count, int dest, int procGrpNo=0 );

  /** Send
   *  - MPI_Sendのインターフェイス
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in] dtype     送信データのMPI_Datatype
   *  @param[in] buf       送信データ
   *  @param[in] count     送信データのサイズ
   *  @param[in] dest      送信先のランク番号(procGrpNo内でのランク番号)
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Send( MPI_Datatype dtype, void *buf, int count, int dest
                    , int procGrpNo=0 );

  /** Recv
   *  - MPI_Recvのインターフェイス
   *
   *  @param[out] buf      受信データ
   *  @param[in] count     受信データのサイズ
   *  @param[in] source    送信元のランク番号(procGrpNo内でのランク番号)
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode Recv( T *buf, int count, int source, int procGrpNo=0 );

  /** Recv
   *  - MPI_Recvのインターフェイス
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]  dtype     送信データのMPI_Datatype
   *  @param[out] buf       受信データ
   *  @param[in]  count     受信データのサイズ
   *  @param[in]  source    送信元のランク番号(procGrpNo内でのランク番号)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Recv( MPI_Datatype dtype, void *buf, int count, int source
                    , int procGrpNo=0 );

  /** Isend
   *  - MPI_Isendのインターフェイス
   *
   *  @param[in]  buf       送信データ
   *  @param[in]  count     送信データのサイズ
   *  @param[in]  dest      送信先のランク番号(procGrpNo内でのランク番号)
   *  @param[out] request   リクエストハンドル
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode Isend( T *buf, int count, int dest, MPI_Request *request
                     , int procGrpNo=0 );

  /** Isend
   *  - MPI_Isendのインターフェイス
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]  dtype     送信データのMPI_Datatype
   *  @param[in]  buf       送信データ
   *  @param[in]  count     送信データのサイズ
   *  @param[in]  dest      送信先のランク番号(procGrpNo内でのランク番号)
   *  @param[out] request   リクエストハンドル
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Isend( MPI_Datatype dtype, void *buf, int count, int dest
                     , MPI_Request *request, int procGrpNo=0 );

  /** Irecv
   *  - MPI_Irecvのインターフェイス
   *
   *  @param[out] buf       受信データ
   *  @param[in]  count     受信データのサイズ
   *  @param[in]  source    送信元のランク番号(procGrpNo内でのランク番号)
   *  @param[out] request   リクエストハンドル
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode Irecv( T *buf, int count, int source, MPI_Request *request
                     , int procGrpNo=0 );

  /** Irecv
   *  - MPI_Irecvのインターフェイス
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]  dtype     送信データのMPI_Datatype
   *  @param[out] buf       受信データ
   *  @param[in]  count     受信データのサイズ
   *  @param[in]  source    送信元のランク番号(procGrpNo内でのランク番号)
   *  @param[out] request   リクエストハンドル
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Irecv( MPI_Datatype dtype, void *buf, int count, int source
                     , MPI_Request *request, int procGrpNo=0 );

  /** Allreduce
   *  - MPI_Allreduceのインターフェイス
   *
   *  @param[in]  sendbuf   送信データ
   *  @param[out] recvbuf   受信データ
   *  @param[in]  count     送受信データのサイズ
   *  @param[in]  op        オペレータ
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode Allreduce( T *sendbuf, T *recvbuf, int count, MPI_Op op
                         , int procGrpNo=0 );

  /** Allreduce
   *  - MPI_Allreduceのインターフェイス
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]  dtype     送信データのMPI_Datatype
   *  @param[in]  sendbuf   送信データ
   *  @param[out] recvbuf   受信データ
   *  @param[in]  count     送受信データのサイズ
   *  @param[in]  op        オペレータ
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Allreduce( MPI_Datatype dtype, void *sendbuf, void *recvbuf
                         , int count, MPI_Op op, int procGrpNo=0 );

  /** Gather
   *  - MPI_Gatherのインターフェイス
   *
   *  @param[in]  sendbuf   送信データ
   *  @param[in]  sendcnt   送信データのサイズ
   *  @param[out] recvbuf   受信データ
   *  @param[in]  recvcnt   送信データのサイズ
   *  @param[in]  root      受信するランク番号(procGrpNo内でのランク番号)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class Ts, class Tr> CPM_INLINE
  cpm_ErrorCode Gather( Ts *sendbuf, int sendcnt, Tr *recvbuf, int recvcnt
                      , int root, int procGrpNo=0 );

  /** Gather
   *  - MPI_Gatherのインターフェイス
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]  stype     送信データのMPI_Datatype
   *  @param[in]  sendbuf   送信データ
   *  @param[in]  sendcnt   送信データのサイズ
   *  @param[in]  rtype     受信データのMPI_Datatype
   *  @param[out] recvbuf   受信データ
   *  @param[in]  recvcnt   送信データのサイズ
   *  @param[in]  root      受信するランク番号(procGrpNo内でのランク番号)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Gather( MPI_Datatype stype, void *sendbuf, int sendcnt
                      , MPI_Datatype rtype, void *recvbuf, int recvcnt
                      , int root, int procGrpNo=0 );

  /** Allgather
   *  - MPI_Allgatherのインターフェイス
   *
   *  @param[in]  sendbuf   送信データ
   *  @param[in]  sendcnt   送信データのサイズ
   *  @param[out] recvbuf   受信データ
   *  @param[in]  recvcnt   送信データのサイズ
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class Ts, class Tr> CPM_INLINE
  cpm_ErrorCode Allgather( Ts *sendbuf, int sendcnt, Tr *recvbuf, int recvcnt
                         , int procGrpNo=0 );

  /** Allgather
   *  - MPI_Allgatherのインターフェイス
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]  stype     送信データのMPI_Datatype
   *  @param[in]  sendbuf   送信データ
   *  @param[in]  sendcnt   送信データのサイズ
   *  @param[in]  rtype     受信データのMPI_Datatype
   *  @param[out] recvbuf   受信データ
   *  @param[in]  recvcnt   送信データのサイズ
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Allgather( MPI_Datatype stype, void *sendbuf, int sendcnt
                         , MPI_Datatype rtype, void *recvbuf, int recvcnt
                         , int procGrpNo=0 );

  /** Gatherv
   *  - MPI_Gathervのインターフェイス
   *
   *  @param[in]  sendbuf   送信データ
   *  @param[in]  sendcnt   送信データのサイズ
   *  @param[out] recvbuf   受信データ
   *  @param[in]  recvcnts  各ランクからの受信データサイズ
   *  @param[in]  displs    各ランクからの受信データ配置位置
   *  @param[in]  root      受信するランク番号(procGrpNo内でのランク番号)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class Ts, class Tr> CPM_INLINE
  cpm_ErrorCode Gatherv( Ts *sendbuf, int sendcnt, Tr *recvbuf, int *recvcnts, int *displs, int root, int procGrpNo=0 );

  /** Gatherv
   *  - MPI_Gathervのインターフェイス
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]  stype     送信データのMPI_Datatype
   *  @param[in]  sendbuf   送信データ
   *  @param[in]  sendcnt   送信データのサイズ
   *  @param[in]  rtype     受信データのMPI_Datatype
   *  @param[out] recvbuf   受信データ
   *  @param[in]  recvcnts  各ランクからの受信データサイズ
   *  @param[in]  displs    各ランクからの受信データ配置位置
   *  @param[in]  root      受信するランク番号(procGrpNo内でのランク番号)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Gatherv( MPI_Datatype stype, void *sendbuf, int sendcnt
                       , MPI_Datatype rtype, void *recvbuf, int *recvcnts
                       , int *displs, int root, int procGrpNo=0 );

  /** Allgatherv
   *  - MPI_Allgathervのインターフェイス
   *
   *  @param[in]  sendbuf   送信データ
   *  @param[in]  sendcnt   送信データのサイズ
   *  @param[out] recvbuf   受信データ
   *  @param[in]  recvcnts  各ランクからの受信データサイズ
   *  @param[in]  displs    各ランクからの受信データ配置位置
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class Ts, class Tr> CPM_INLINE
  cpm_ErrorCode Allgatherv( Ts *sendbuf, int sendcnt, Tr *recvbuf, int *recvcnts, int *displs, int procGrpNo=0 );

  /** Allgatherv
   *  - MPI_Allgathervのインターフェイス
   *  - MPI_Datatypeを指定するバージョン
   *
   *  @param[in]  stype     送信データのMPI_Datatype
   *  @param[in]  sendbuf   送信データ
   *  @param[in]  sendcnt   送信データのサイズ
   *  @param[in]  rtype     受信データのMPI_Datatype
   *  @param[out] recvbuf   受信データ
   *  @param[in]  recvcnts  各ランクからの受信データサイズ
   *  @param[in]  displs    各ランクからの受信データ配置位置
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Allgatherv( MPI_Datatype stype, void *sendbuf, int sendcnt
                          , MPI_Datatype rtype, void *recvbuf, int *recvcnts
                          , int *displs, int procGrpNo=0 );





////// MPI処理のFortran用インターフェイス関数 //////

  /** cpm_Wait
   *  - MPI_Waitのインターフェイス
   *  - Fortranインターフェイス用
   *
   *  @param[in] reqNo リクエスト番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_Wait( int reqNo );

  /** cpm_Waitall
   *  - MPI_Waitallのインターフェイス
   *
   *  @param[in] count     リクエストの数
   *  @param[in] reqNoList リクエスト番号のリスト
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_Waitall( int count, int reqNoList[] );

  /** cpm_Isend
   *  - MPI_Isendのインターフェイス
   *  - Fortranインターフェイス用
   *
   *  @param[in]  buf       送信データ
   *  @param[in]  count     送信データのサイズ
   *  @param[in]  datatype  受信データのデータタイプ(cpm_fparam.fi参照)
   *  @param[in]  dest      送信先のランク番号(procGrpNo内でのランク番号)
   *  @param[out] reqNo     リクエスト番号(Fortran用)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_Isend( void *buf, int count, int datatype, int dest, int *reqNo, int procGrpNo=0 );

  /** cpm_Irecv
   *  - MPI_Irecvのインターフェイス
   *  - Fortranインターフェイス用
   *
   *  @param[out] buf       受信データ
   *  @param[in]  count     受信データのサイズ
   *  @param[in]  datatype  受信データのデータタイプ(cpm_fparam.fi参照)
   *  @param[in]  source    送信元のランク番号(procGrpNo内でのランク番号)
   *  @param[out] reqNo     リクエスト番号(Fortran用)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_Irecv( void *buf, int count, int datatype, int source, int *reqNo, int procGrpNo=0 );

#if 0
  /** cpm_BndCommS3D_nowait
   *  - BndCommS3D_nowaitのインターフェイス
   *  - Fortranインターフェイス用
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   datatype  袖通信データのデータタイプ(cpm_fparam.fi参照)
   *  @param[out]  reqNo     リクエスト番号配列(サイズ48,CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_BndCommS3D_nowait( void *array, int imax, int jmax, int kmax
                                     , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo=0 );

  /** cpm_BndCommV3D_nowait
   *  - BndCommV3D_nowaitのインターフェイス
   *  - Fortranインターフェイス用
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   datatype  袖通信データのデータタイプ(cpm_fparam.fi参照)
   *  @param[out]  reqNo     リクエスト番号配列(サイズ48,CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_BndCommV3D_nowait( void *array, int imax, int jmax, int kmax
                                     , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo=0 );

  /** cpm_BndCommS4D_nowait
   *  - BndCommS4D_nowaitのインターフェイス
   *  - Fortranインターフェイス用
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   datatype  袖通信データのデータタイプ(cpm_fparam.fi参照)
   *  @param[out]  reqNo     リクエスト番号配列(サイズ48,CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_BndCommS4D_nowait( void *array, int imax, int jmax, int kmax, int nmax
                                     , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo=0 );

  /** cpm_wait_BndCommS3D
   *  - wait_BndCommS3Dのインターフェイス
   *  - Fortranインターフェイス用
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   datatype  袖通信データのデータタイプ(cpm_fparam.fi参照)
   *  @param[in]   reqNo     リクエスト番号配列(サイズ48,CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_wait_BndCommS3D( void *array, int imax, int jmax, int kmax
                                   , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo=0 );

  /** cpm_wait_BndCommV3D
   *  - wait_BndCommV3Dのインターフェイス
   *  - Fortranインターフェイス用
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   datatype  袖通信データのデータタイプ(cpm_fparam.fi参照)
   *  @param[in]   reqNo     リクエスト番号配列(サイズ48,CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_wait_BndCommV3D( void *array, int imax, int jmax, int kmax
                                   , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo=0 );

  /** cpm_wait_BndCommS4D
   *  - wait_BndCommS4Dのインターフェイス
   *  - Fortranインターフェイス用
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   datatype  袖通信データのデータタイプ(cpm_fparam.fi参照)
   *  @param[in]   reqNo     リクエスト番号配列(サイズ48,CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_wait_BndCommS4D( void *array, int imax, int jmax, int kmax, int nmax
                                   , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo=0 );

  /** cpm_BndCommV3DEx_nowait
   *  - BndCommV3DEx_nowaitのインターフェイス
   *  - Fortranインターフェイス用
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   datatype  袖通信データのデータタイプ(cpm_fparam.fi参照)
   *  @param[out]  reqNo     リクエスト番号配列(サイズ48,CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_BndCommV3DEx_nowait( void *array, int imax, int jmax, int kmax
                                       , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo=0 );

  /** cpm_BndCommS4DEx_nowait
   *  - BndCommS4DEx_nowaitのインターフェイス
   *  - Fortranインターフェイス用
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   datatype  袖通信データのデータタイプ(cpm_fparam.fi参照)
   *  @param[out]  reqNo     リクエスト番号配列(サイズ48,CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_BndCommS4DEx_nowait( void *array, int nmax, int imax, int jmax, int kmax
                                       , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo=0 );

  /** cpm_wait_BndCommV3DEx
   *  - wait_BndCommV3DExのインターフェイス
   *  - Fortranインターフェイス用
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   datatype  袖通信データのデータタイプ(cpm_fparam.fi参照)
   *  @param[in]   reqNo     リクエスト番号配列(サイズ48,CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_wait_BndCommV3DEx( void *array, int imax, int jmax, int kmax
                                     , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo=0 );

  /** cpm_wait_BndCommS4DEx
   *  - wait_BndCommS4DExのインターフェイス
   *  - Fortranインターフェイス用
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[in]   datatype  袖通信データのデータタイプ(cpm_fparam.fi参照)
   *  @param[in]   reqNo     リクエスト番号配列(サイズ48,CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_wait_BndCommS4DEx( void *array, int nmax, int imax, int jmax, int kmax
                                     , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo=0 );
#endif





////// 袖通信関数 //////

  /** 袖通信バッファのセット
   *  - 6face分の送受信バッファを確保する
   *
   *  @param[in] maxVC     送受信バッファの最大袖数
   *  @param[in] maxN      送受信バッファの最大成分数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  virtual
  cpm_ErrorCode SetBndCommBuffer( size_t maxVC, size_t maxN, int procGrpNo=0 ) = 0;

  /** 袖通信バッファサイズの取得
   *  - 袖通信バッファとして確保されている配列サイズ(byte)を返す
   *
   *  @param[in] procGrpNo プロセスグループ番号(負の場合、全プロセスグループでのトータルを返す)
   *  @return バッファサイズ(byte)
   */
  virtual
  size_t GetBndCommBufferSize( int procGrpNo=0 ) = 0;

#if 0
  /** 袖通信(Scalar3D版)
   *  - (imax,jmax,kmax)の形式の配列の袖通信を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                          , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 袖通信(Scalar3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax)の形式の配列の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                          , int vc, int vc_comm, int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 袖通信(Vector3D版)
   *  - (imax,jmax,kmax,3)の形式の配列の袖通信を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                          , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 袖通信(Vector3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,3)の形式の配列の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                          , int vc, int vc_comm, int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 袖通信(Scalar4D版)
   *  - (imax,jmax,kmax,nmax)の形式の配列の袖通信を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                          , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 袖通信(Scalar4D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                          , int vc, int vc_comm, int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 袖通信(Scalar4D版, パディングサイズ指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の袖通信を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    pad_size  パディングサイズ(i,j,k,n)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                          , int pad_size[4], int procGrpNo );

  /** 袖通信(Scalar4D版, MPI_Datatype指定, パディングサイズ指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    pad_size  パディングサイズ(i,j,k,n)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                          , int vc, int vc_comm, int pad_size[4], int procGrpNo );

  /** 非同期版袖通信(Scalar3D版)
   *  - (imax,jmax,kmax)の形式の配列の非同期袖通信を行う
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommS3Dをコールする
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @param[in]   padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS3D_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , MPI_Request req[48], int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信(Scalar3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax)の形式の配列の非同期袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommS3Dをコールする
   *  @param[in]   dtype     袖通信データのMPI_Datatype
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @param[in]   padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS3D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, MPI_Request req[48]
                                 , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信(Vector3D版)
   *  - (imax,jmax,kmax,3)の形式の配列の非同期袖通信を行う
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommV3Dをコールする
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @param[in]   padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommV3D_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , MPI_Request req[48], int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信(Vector3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,3)の形式の配列の非同期袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommV3Dをコールする
   *  @param[in]   dtype     袖通信データのMPI_Datatype
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @param[in]   padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommV3D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, MPI_Request req[48]
                                 , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信(Scalar4D版)
   *  - (imax,jmax,kmax,nmax)の形式の配列の非同期袖通信を行う
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommS4Dをコールする
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @param[in]   padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4D_nowait( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                 , MPI_Request req[48], int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信(Scalar4D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の非同期袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommS4Dをコールする
   *  @param[in]   dtype     袖通信データのMPI_Datatype
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @param[in]   padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                                 , int vc, int vc_comm, MPI_Request req[48]
                                 , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信(Scalar4D版, パディングサイズ指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の非同期袖通信を行う
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommS4Dをコールする
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   pad_size  パディングサイズ(i,j,k,n)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4D_nowait( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                 , MPI_Request req[48], int pad_size[4], int procGrpNo );

  /** 非同期版袖通信(Scalar4D版, MPI_Datatype指定, パディングサイズ指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の非同期袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommS4Dをコールする
   *  @param[in]   dtype     袖通信データのMPI_Datatype
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   pad_size  パディングサイズ(i,j,k,n)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                                 , int vc, int vc_comm, MPI_Request req[48]
                                 , int pad_size[4], int procGrpNo );

  /** 非同期版袖通信のwait、展開(Scalar3D版)
   *  - (imax,jmax,kmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                               , MPI_Request req[48], int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信のwait、展開(Scalar3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode wait_BndCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                               , int vc, int vc_comm, MPI_Request req[48]
                               , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信のwait、展開(Vector3D版)
   *  - (imax,jmax,kmax,3)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                               , MPI_Request req[48], int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信のwait、展開(Vector3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,3)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode wait_BndCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                               , int vc, int vc_comm, MPI_Request req[48]
                               , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信のwait、展開(Scalar4D版)
   *  - (imax,jmax,kmax,nmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                               , MPI_Request req[48], int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信のwait、展開(Scalar4D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode wait_BndCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                               , int vc, int vc_comm, MPI_Request req[48]
                               , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信のwait、展開(Scalar4D版, パディングサイズ指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    pad_size  パディングサイズ(n,i,j,k)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                               , MPI_Request req[48], int pad_size[4], int procGrpNo );

  /** 非同期版袖通信のwait、展開(Scalar4D版, MPI_Datatype指定, パディングサイズ指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    pad_size  パディングサイズ(n,i,j,k)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode wait_BndCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                               , int vc, int vc_comm, MPI_Request req[48]
                               , int pad_size[4], int procGrpNo );

  /** 周期境界袖通信(Scalar3D版)
   *  - (imax,jmax,kmax)の形式の配列の周期境界方向の袖通信を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                               , cpm_DirFlag dir, cpm_PMFlag pm
                               , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Scalar3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax)の形式の配列の周期境界方向の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                               , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Vector3D版)
   *  - (imax,jmax,kmax,3)の形式の配列の周期境界方向の袖通信を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                               , cpm_DirFlag dir, cpm_PMFlag pm
                               , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Vector3D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,3)の形式の配列の周期境界方向の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                               , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Scalar4D版)
   *  - (imax,jmax,kmax,nmax)の形式の配列の周期境界方向の袖通信を行う
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                               , cpm_DirFlag dir, cpm_PMFlag pm
                               , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Scalar4D版, MPI_Datatype指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の周期境界方向の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                               , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Scalar4D版, パディングサイズ指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の周期境界方向の袖通信を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    pad_size  パディングサイズ(i,j,k,n)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                               , cpm_DirFlag dir, cpm_PMFlag pm
                               , int pad_size[4], int procGrpNo );

  /** 周期境界袖通信(Scalar4D版, MPI_Datatype指定, パディングサイズ指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の周期境界方向の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
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
   *  @param[in]    pad_size  パディングサイズ(i,j,k,n)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                               , int pad_size[4], int procGrpNo );

  /** 袖通信(Vector3DEx版)
   *  - (3,imax,jmax,kmax)の形式の配列の袖通信を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                            , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 袖通信(Vector3DEx版, MPI_Datatype指定)
   *  - (3,imax,jmax,kmax)の形式の配列の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                            , int vc, int vc_comm, int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 袖通信(Scalar4DEx版)
   *  - (nmax,imax,jmax,kmax)の形式の配列の袖通信を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                            , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 袖通信(Scalar4DEx版, MPI_Datatype指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                            , int vc, int vc_comm, int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 袖通信(Scalar4DEx版, パディングサイズ指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の袖通信を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    pad_size  パディングサイズ(n,i,j,k)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                            , int pad_size[4], int procGrpNo=0 );

  /** 袖通信(Scalar4DEx版, MPI_Datatype指定, パディングサイズ指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    pad_size  パディングサイズ(n,i,j,k)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                            , int vc, int vc_comm, int pad_size[4], int procGrpNo=0 );

  /** 非同期版袖通信(Vector3DEx版)
   *  - (3,imax,jmax,kmax)の形式の配列の非同期袖通信を行う
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommV3DExをコールする
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @param[in]   padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommV3DEx_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                   , MPI_Request req[48], int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信(Vector3DEx版, MPI_Datatype指定)
   *  - (3,imax,jmax,kmax)の形式の配列の非同期袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommV3DExをコールする
   *  @param[in]   dtype     袖通信データのMPI_Datatype
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @param[in]   padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommV3DEx_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                   , int vc, int vc_comm, MPI_Request req[48]
                                   , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信(Scalar4DEx版)
   *  - (nmax,imax,jmax,kmax)の形式の配列の非同期袖通信を行う
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommS4DExをコールする
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @param[in]   padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4DEx_nowait( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                   , MPI_Request req[48], int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信(Scalar4DEx版, MPI_Datatype指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の非同期袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommS4DExをコールする
   *  @param[in]   dtype     袖通信データのMPI_Datatype
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @param[in]   padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4DEx_nowait( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                   , int vc, int vc_comm, MPI_Request req[48]
                                   , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信(Scalar4DEx版, パディングサイズ指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の非同期袖通信を行う
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommS4DExをコールする
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   pad_size  パディングサイズ(n,i,j,k)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4DEx_nowait( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                   , MPI_Request req[48], int pad_size[4], int procGrpNo=0 );

  /** 非同期版袖通信(Scalar4DEx版, MPI_Datatype指定, パディングサイズ指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の非同期袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
   *  - waitと展開は行わず、requestを返す
   *  - wait、展開はwait_BndCommS4DExをコールする
   *  @param[in]   dtype     袖通信データのMPI_Datatype
   *  @param[in]   array     袖通信をする配列の先頭ポインタ
   *  @param[in]   nmax      配列サイズ(成分数)
   *  @param[in]   imax      配列サイズ(I方向)
   *  @param[in]   jmax      配列サイズ(J方向)
   *  @param[in]   kmax      配列サイズ(K方向)
   *  @param[in]   vc        仮想セル数
   *  @param[in]   vc_comm   通信する仮想セル数
   *  @param[out]  req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]   pad_size  パディングサイズ(n,i,j,k)
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4DEx_nowait( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                   , int vc, int vc_comm, MPI_Request req[48]
                                   , int pad_size[4], int procGrpNo=0 );

  /** 非同期版袖通信のwait、展開(Vector3DEx版)
   *  - (3,imax,jmax,kmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , MPI_Request req[48], int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信のwait、展開(Vector3DEx版, MPI_Datatype指定)
   *  - (3,imax,jmax,kmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode wait_BndCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, MPI_Request req[48]
                                 , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信のwait、展開(Scalar4DEx版)
   *  - (nmax,imax,jmax,kmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , MPI_Request req[48], int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信のwait、展開(Scalar4DEx版, MPI_Datatype指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode wait_BndCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, MPI_Request req[48]
                                 , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 非同期版袖通信のwait、展開(Scalar4DEx版, パディングサイズ指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    pad_size  パディングサイズ(n,i,j,k)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , MPI_Request req[48], int pad_size[4], int procGrpNo=0 );

  /** 非同期版袖通信のwait、展開(Scalar4DEx版, MPI_Datatype指定, パディングサイズ指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  - MPI_Datatypeを指定するバージョン
   *  @param[in]    dtype     袖通信データのMPI_Datatype
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト(サイズ48、CARTの場合12でも良い)
   *  @param[in]    pad_size  パディングサイズ(n,i,j,k)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode wait_BndCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, MPI_Request req[48]
                                 , int pad_size[4], int procGrpNo=0 );

  /** 周期境界袖通信(Vector3DEx版)
   *  - (3,imax,jmax,kmax)の形式の配列の周期境界方向の袖通信を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , cpm_DirFlag dir, cpm_PMFlag pm
                                 , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Vector3DEx版, MPI_Datatype指定)
   *  - (3,imax,jmax,kmax)の形式の配列の周期境界方向の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                                 , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Scalar4DEx版)
   *  - (nmax,imax,jmax,kmax)の形式の配列の周期境界方向の袖通信を行う
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , cpm_DirFlag dir, cpm_PMFlag pm
                                 , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Scalar4DEx版, MPI_Datatype指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の周期境界方向の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                                 , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Scalar4DEx版, パディングサイズ指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の周期境界方向の袖通信を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
   *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
   *  @param[in]    pad_size  パディングサイズ(n,i,j,k)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , cpm_DirFlag dir, cpm_PMFlag pm
                                 , int pad_size[4], int procGrpNo=0 );

  /** 周期境界袖通信(Scalar4DEx版, MPI_Datatype指定, パディングサイズ指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の周期境界方向の袖通信を行う
   *  - MPI_Datatypeを指定するバージョン
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
   *  @param[in]    pad_size  パディングサイズ(n,i,j,k)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                                 , int pad_size[4], int procGrpNo=0 );
#endif




////// ローカルボクセルサイズでの配列確保関数 //////

  /** 配列の初期化処理
   *  @param[out] array 初期化する配列のポインタ
   *  @param[in]  size  配列サイズ
   */
  template<class T>
  void InitArray( T *array, size_t size );

  /** 配列のコピー
   *  @param[in]  source コピー元の配列のポインタ
   *  @param[out] dist   コピー先の配列のポインタ
   *  @param[in]  size   配列サイズ
   */
  template<class T>
  void CopyArray( T *source, T *dist, size_t size );

  /** 配列確保 double(imax,jmax,kmax)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(3word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  double* AllocDoubleS3D( int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 float(imax,jmax,kmax)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(3word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  float* AllocFloatS3D( int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 int(imax,jmax,kmax)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(3word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  int* AllocIntS3D( int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 double(imax,jmax,kmax,3)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(4word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  double* AllocDoubleV3D( int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 float(imax,jmax,kmax,3)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(4word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  float* AllocFloatV3D( int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 int(imax,jmax,kmax,3)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(4word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  int* AllocIntV3D( int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 double(3,imax,jmax,kmax)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(4word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  double* AllocDoubleV3DEx( int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 float(3,imax,jmax,kmax)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(4word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  float* AllocFloatV3DEx( int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 int(3,imax,jmax,kmax)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(4word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  int* AllocIntV3DEx( int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 double(imax,jmax,kmax,nmax)
   *  @param[in]  nmax      成分数
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(4word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  double* AllocDoubleS4D( int nmax, int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 float(imax,jmax,kmax,nmax)
   *  @param[in]  nmax      成分数
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(4word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  float* AllocFloatS4D( int nmax, int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 int(imax,jmax,kmax,nmax)
   *  @param[in]  nmax      成分数
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(4word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  int* AllocIntS4D( int nmax, int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 double(nmax,imax,jmax,kmax)
   *  @param[in]  nmax      成分数
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(4word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  double* AllocDoubleS4DEx( int nmax, int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 float(nmax,imax,jmax,kmax)
   *  @param[in]  nmax      成分数
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(4word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  float* AllocFloatS4DEx( int nmax, int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** 配列確保 int(nmax,imax,jmax,kmax)
   *  @param[in]  nmax      成分数
   *  @param[in]  vc        仮想セル数
   *  @param[in]  padding   パディングフラグ(true:する、false:しない)
   *  @param[out] pad_size  各次元のパディング数(4word, NULLのとき格納しない)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  int* AllocIntS4DEx( int nmax, int vc, bool padding=false, int *pad_size=NULL, int procGrpNo=0 );

  /** パディングサイズ取得処理(静的関数)
   *  @param[in] size  配列サイズ
   *  @param[in] vc    仮想セル数
   *  @return          パディングサイズ
   */
  static int GetPaddingSize1D( const int size, const int vc );

  /** パディングサイズ取得処理(静的関数)
   *  @param[in]  atype    配列形状タイプ
   *  @param[in]  size     配列サイズ{imax,jmax,kmax}
   *  @param[in]  vc       仮想セル数
   *  @param[out] pad_size パディングサイズ(S3Dのとき3word、V3D,S4Dのとき4word{px,py,pz,pn}、V3DEx,S4DExのとき4word{pn,px,py,pz})
   *  @param[in]  nmax     成分数(S4D,S4DExのとき必須)
   */
  static void GetPaddingSize( CPM_ARRAY_SHAPE atype, const int *size, const int vc, int *pad_size, int nmax=0 );




  /** flush */
  void flush(std::ostream &out, int procGrpNo=0);

  /** flush */
  void flush(FILE *fp, int procGrpNo=0);

#ifdef _DEBUG
  void printVoxelInfo(int myrank=-1);
#endif

protected:

  /** コンストラクタ */
  cpm_BaseParaManager();

  /** デストラクタ */
  virtual ~cpm_BaseParaManager();

  /** 配列確保(double)
   *  @param[in] nmax 成分数
   *  @param[in] sz   配列サイズ
   *  @param[in] vc   仮想セル数
   *  @param[in] procGrpNo  領域分割を行うプロセスグループ番号
   *  @return         確保した配列のポインタ
   */
  virtual
  double* AllocDouble( int nmax, int sz[3], int vc, int procGrpNo ) = 0;

  /** 配列確保(float)
   *  @param[in] nmax 成分数
   *  @param[in] sz   配列サイズ
   *  @param[in] vc   仮想セル数
   *  @param[in] procGrpNo  領域分割を行うプロセスグループ番号
   *  @return         確保した配列のポインタ
   */
  virtual
  float* AllocFloat( int nmax, int sz[3], int vc, int procGrpNo ) = 0;

  /** 配列確保(int)
   *  @param[in] nmax 成分数
   *  @param[in] sz   配列サイズ
   *  @param[in] vc   仮想セル数
   *  @param[in] procGrpNo  領域分割を行うプロセスグループ番号
   *  @return         確保した配列のポインタ
   */
  virtual
  int* AllocInt( int nmax, int sz[3], int vc, int procGrpNo ) = 0;





////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:

protected:
  /** プロセス並列数 */
  int m_nRank;

  /** MPI_COMM_WORLDでの自ランク番号 */
  int m_rankNo;

  /** 領域分割タイプ */
  cpm_DomainType m_domainType;

  /** プロセスグループのリスト
   *  - VOXEL空間番号をインデクスとしたVOXEL空間のMPIコミュニケータを格納
   *  - vectorのインデクス=プロセスグループ番号とする
   *  - [0]には必ずMPI_COMM_WORLDを格納
   *  - 自ランクが含まれるプロセスグループのみを管理する
   *    (同じプロセスグループでもプロセス毎に異なるプロセスグループ番号になる場合もある)
   */
  std::vector<MPI_Comm> m_procGrpList;

  /** MPI_Requestの管理マップ
   *  - Fortranインターフェイス用
   */
  cpm_ObjList<MPI_Request> m_reqList;

  /** プロセスグループ毎の定義点タイプ管理マップ
   *  - プロセスグループ番号をキーとした定義点タイプマップ
   *  - 自ランクが含まれるプロセスグループのみを管理する
   */
  DefPointMap m_defPointMap;
};

//インライン関数
#include "inline/cpm_BaseParaManager_inline.h"

#endif /* _CPM_BASEPARAMANAGER_H_ */
