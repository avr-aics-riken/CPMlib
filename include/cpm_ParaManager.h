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
 * @file   cpm_ParaManager.h
 * パラレルマネージャクラスのヘッダーファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */

#ifndef _CPM_PARAMANAGER_H_
#define _CPM_PARAMANAGER_H_

#include "cpm_Base.h"
#include <map>
#include <vector>
#include <typeinfo>
#include "cpm_DomainInfo.h"
#include "cpm_VoxelInfo.h"
#include "cpm_ObjList.h"
#include <string.h> // for memset()

/** プロセスグループ毎のVOXEL空間情報管理マップ */
typedef std::map<int, cpm_VoxelInfo*> VoxelInfoMap;

/** CPMの並列管理クラス
 *  - 現時点ではユーザがインスタンスすることを許していない
 *  - get_instance静的関数を用いて唯一のインスタンスを取得する
 */
class cpm_ParaManager : public cpm_Base
{
friend class C_PARAMANAGER;
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:

  /** 唯一のインスタンスの取得
   *  @return インスタンスのポインタ
   */
  static cpm_ParaManager* get_instance();

  /** 唯一のインスタンスの取得(initialize処理も実行)
   *  @param[in] argc       プログラム実行時引数の数
   *  @param[in] argv       プログラム実行時引数
   *  @param[in] domainType インスタンスする領域タイプ(デフォルト=カーテシアン)
   *  @return インスタンスのポインタ
   */
  static cpm_ParaManager* get_instance(int &argc, char**& argv, cpm_DomainType domainType=CPM_DOMAIN_CARTESIAN);

  /** 唯一のインスタンスの取得(initialize処理も実行. fortanインターフェイス用)
   *  @param[in] domainType インスタンスする領域タイプ(デフォルト=カーテシアン)
   *  @return インスタンスのポインタ
   */
  static cpm_ParaManager* get_instance(cpm_DomainType domainType);

  /** 初期化処理(MPI_Initは実行済みの場合)
   *  - MPI_Initは既に実行済みである必要がある
   *  - 並列数、自ランク番号を取得
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Initialize();

  /** 初期化処理(MPI_Initも実行する)
   *  - MPI_Initが実行されていない場合、実行する
   *  - 並列数、自ランク番号を取得
   *  @param[in] argc プログラム実行時引数の数
   *  @param[in] argv プログラム実行時引数
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Initialize(int &argc, char**& argv);

  /** 並列実行であるかチェックする
   *  並列実行であっても、並列数が1のときはfalseとなる
   *  @retval true  並列実行
   *  @retval false 逐次実行
   */
  bool IsParallel();

  /** 並列実行であるかチェックする(const)
   *  - 並列実行であっても、並列数が1のときはfalseとなる
   *  @retval true  並列実行
   *  @retval false 逐次実行
   */
  bool IsParallel() const;

  /** カーテシアン用の領域分割
   *  - 既に作成済みの領域分割情報を用いた領域分割処理
   *  @param[in] domainInfo 領域分割情報
   *  @param[in] maxVC      最大の袖数(袖通信用)
   *  @param[in] maxN       最大の成分数(袖通信用)
   *  @param[in] procGrpNo  領域分割を行うプロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  virtual
  cpm_ErrorCode VoxelInit( cpm_GlobalDomainInfo* domainInfo
                         , size_t maxVC=1, size_t maxN=3
                         , int procGrpNo=0 );

  /** カーテシアン用の領域分割
   *  - 領域分割の各種情報を引数で渡して領域分割を行う
   *  - プロセスグループの全てのランクが活性ドメインになる
   *  - I,J,K方向の領域分割数を指定するバージョン
   *  @param[in] div        領域分割数
   *  @param[in] vox        空間全体のボクセル数
   *  @param[in] origin     空間全体の原点
   *  @param[in] region     空間全体のサイズ
   *  @param[in] maxVC      最大の袖数(袖通信用)
   *  @param[in] maxN       最大の成分数(袖通信用)
   *  @param[in] divPolicy  自動分割ポリシー
   *  @param[in] procGrpNo  領域分割を行うプロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  virtual
  cpm_ErrorCode VoxelInit( int div[3], int vox[3], double origin[3], double region[3]
                         , size_t maxVC=1, size_t maxN=3, cpm_DivPolicy divPolicy=DIV_COMM_SIZE
                         , int procGrpNo=0 );

  /** カーテシアン用の領域分割
   *  - 領域分割の各種情報を引数で渡して領域分割を行う
   *  - プロセスグループの全てのランクが活性ドメインになる
   *  - 並列数=プロセスグループの並列数とし、内部で自動的に領域分割をするバージョン
   *  @param[in] vox        空間全体のボクセル数
   *  @param[in] origin     空間全体の原点
   *  @param[in] region     空間全体のサイズ
   *  @param[in] maxVC      最大の袖数(袖通信用)
   *  @param[in] maxN       最大の成分数(袖通信用)
   *  @param[in] divPolicy  自動分割ポリシー
   *  @param[in] procGrpNo  領域分割を行うプロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  virtual
  cpm_ErrorCode VoxelInit( int vox[3], double origin[3], double region[3]
                         , size_t maxVC=1, size_t maxN=3, cpm_DivPolicy divPolicy=DIV_COMM_SIZE
                         , int procGrpNo=0 );

  /** カーテシアン用の領域分割(ActiveSubdomain指定)
   *  - 領域分割の各種情報を引数で渡して領域分割を行う
   *  - ActiveSubdomainファイルで指定される領域分割位置のランクが活性ドメインになる
   *  - I,J,K方向の領域分割数を指定するバージョン
   *  - 指定の領域分割数とActiveSubdomainファイルで指定されている領域分割数が一致している必要がある
   *  - ActiveSubdomain数と並列数が一致している必要がある
   *  @param[in] div           領域分割数
   *  @param[in] vox           空間全体のボクセル数
   *  @param[in] origin        空間全体の原点
   *  @param[in] region        空間全体のサイズ
   *  @param[in] subDomainFile ActiveSubdomainファイル名
   *  @param[in] maxVC         最大の袖数(袖通信用)
   *  @param[in] maxN          最大の成分数(袖通信用)
   *  @param[in] procGrpNo     領域分割を行うプロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  virtual
  cpm_ErrorCode VoxelInit_Subdomain( int div[3], int vox[3], double origin[3], double region[3]
                                   , std::string subDomainFile
                                   , size_t maxVC=1, size_t maxN=3, int procGrpNo=0 );

  /** カーテシアン用の領域分割(ActiveSubdomain指定)
   *  - 領域分割の各種情報を引数で渡して領域分割を行う
   *  - ActiveSubdomainファイルで指定される領域分割位置のランクが活性ドメインになる
   *  - ActiveSubdomainファイルで指定されている領域分割数で領域分割を行う
   *  - ActiveSubdomain数と並列数が一致している必要がある
   *  @param[in] vox           空間全体のボクセル数
   *  @param[in] origin        空間全体の原点
   *  @param[in] region        空間全体のサイズ
   *  @param[in] subDomainFile ActiveSubdomainファイル名
   *  @param[in] maxVC         最大の袖数(袖通信用)
   *  @param[in] maxN          最大の成分数(袖通信用)
   *  @param[in] procGrpNo     領域分割を行うプロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  virtual
  cpm_ErrorCode VoxelInit_Subdomain( int vox[3], double origin[3], double region[3]
                                   , std::string subDomainFile
                                   , size_t maxVC=1, size_t maxN=3, int procGrpNo=0 );

  /** LMR用の領域分割
   *  - FXgen出力の領域情報ファイル、木情報ファイルを渡して領域分割情報を生成する
   *  @param[in] treefile  木情報ファイル
   *  @param[in] maxVC     最大の袖数(袖通信用)
   *  @param[in] maxN      最大の成分数(袖通信用)
   *  @param[in] procGrpNo 領域分割を行うプロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  virtual
  cpm_ErrorCode VoxelInit_LMR( std::string treeFile
                             , size_t maxVC=1, size_t maxN=3, int procGrpNo=0 );

  /** プロセスグループの作成
   *  - 指定されたプロセスリストを使用してプロセスグループを生成する
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

  /** VOXEL空間マップを検索
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return VOXEL空間情報ポインタ
   */
  const cpm_VoxelInfo* FindVoxelInfo( int procGrpNo=0 );

  /** 領域分割数を取得
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 領域分割数整数配列のポインタ
   */
  const int* GetDivNum( int procGrpNo=0 );

  /** ピッチを取得
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return ピッチ実数配列のポインタ
   */
  const double* GetPitch( int procGrpNo=0 );

  /** 全体ボクセル数を取得
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 全体ボクセル数の整数配列ポインタ(3word)
   */
  const int* GetGlobalVoxelSize( int procGrpNo=0 );

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

  /** 自ランクの空間原点を取得
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 自ランクの空間原点実数配列のポインタ
   */
  const double* GetLocalOrigin( int procGrpNo=0 );

  /** 自ランクの空間サイズを取得
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 自ランクの空間サイズ実数配列のポインタ
   */
  const double* GetLocalRegion( int procGrpNo=0 );

  /** 自ランクの領域分割位置を取得
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 自ランクの領域分割位置整数配列のポインタ
   */
  const int* GetDivPos( int procGrpNo=0 );

  /** 自ランクの始点VOXELの全体空間でのインデクスを取得
   *  - 全体空間の先頭インデクスを0としたC型のインデクス
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 自ランクの始点インデクス整数配列のポインタ
   */
  const int* GetVoxelHeadIndex( int procGrpNo=0 );

  /** 自ランクの終点VOXELの全体空間でのインデクスを取得
   *  - 全体空間の先頭インデクスを0としたC型のインデクス
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 自ランクの終点インデクス整数配列のポインタ
   */
  const int* GetVoxelTailIndex( int procGrpNo=0 );

  /** 自ランクの隣接ランク番号を取得
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 自ランクの隣接ランク番号整数配列のポインタ
   */
  const int* GetNeighborRankID( int procGrpNo=0 );

  /** 自ランクの周期境界の隣接ランク番号を取得
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 自ランクの周期境界の隣接ランク番号整数配列のポインタ
   */
  const int* GetPeriodicRankID( int procGrpNo=0 );

  /** 指定面における自ランクの隣接ランク番号を取得
   *  @param[in]  face 面方向
   *  @param[out] num  面の数(CARTのとき1)
   *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定面における自ランクの隣接ランク番号整数配列のポインタ
   */
  const int* GetNeighborRankList( cpm_FaceFlag face, int &num, int procGrpNo=0 );

  /** 指定面における自ランクの周期境界の隣接ランク番号を取得
   *  @param[in]  face 面方向
   *  @param[out] num  面の数(CARTのとき1)
   *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
   *  @return 指定面における自ランクの周期境界の隣接ランク番号整数配列のポインタ
   */
  const int* GetPeriodicRankList( cpm_FaceFlag face, int &num, int procGrpNo=0 );

  /** 指定面におけるレベル差を取得
   *  @param[in]  face 面方向
   *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
   *  @return     レベル差(0:同じレベル, 1:fine, -1:coarse)
   */
  int GetNeighborLevelDiff( cpm_FaceFlag face, int procGrpNo=0 );

  /** 指定idを含む全体ボクセル空間のインデクス範囲を取得
   *  - 全体空間実セルのスタートインデクスを0としたときの，i,j,k各方向の
   *    スタートインデクスと長さを取得する．
   *  @param[in]  id        判定するid
   *  @param[in]  array     判定対象の配列ポインタ
   *  @param[in]  vc        仮想セル数
   *  @param[out] ista      I方向範囲のスタートインデクス
   *  @param[out] jsta      J方向範囲のスタートインデクス
   *  @param[out] ksta      K方向範囲のスタートインデクス
   *  @param[out] ilen      I方向範囲の長さ
   *  @param[out] jlen      J方向範囲の長さ
   *  @param[out] klen      K方向範囲の長さ
   *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
   *  @retval true  指定idを含むセルが存在した
   *  @retval false 指定idを含むセルが存在しない
   */
  virtual
  bool GetBndIndexExtGc( int id, int *array, int vc
                       , int &ista, int &jsta, int &ksta, int &ilen, int &jlen, int &klen
                       , int procGrpNo=0 );

  /** 指定idを含む全体ボクセル空間のインデクス範囲を取得
   *  - 全体空間実セルのスタートインデクスを0としたときの，i,j,k各方向の
   *    スタートインデクスと長さを取得する．
   *  @param[in]  id        判定するid
   *  @param[in]  array     判定対象の配列ポインタ
   *  @param[in]  imax      配列サイズ(I方向)
   *  @param[in]  jmax      配列サイズ(J方向)
   *  @param[in]  kmax      配列サイズ(K方向)
   *  @param[in]  vc        仮想セル数
   *  @param[out] ista      I方向範囲のスタートインデクス
   *  @param[out] jsta      J方向範囲のスタートインデクス
   *  @param[out] ksta      K方向範囲のスタートインデクス
   *  @param[out] ilen      I方向範囲の長さ
   *  @param[out] jlen      J方向範囲の長さ
   *  @param[out] klen      K方向範囲の長さ
   *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
   *  @retval true  指定idを含むセルが存在した
   *  @retval false 指定idを含むセルが存在しない
   */
  virtual
  bool GetBndIndexExtGc( int id, int *array, int imax, int jmax, int kmax, int vc
                       , int &ista, int &jsta, int &ksta, int &ilen, int &jlen, int &klen
                       , int procGrpNo=0 );

  /** 自ランクの境界が外部境界かどうかを判定
   *  @param[in] face      面方向
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @retval    true      外部境界
   *  @retval    false     外部境界でない
   */
  bool IsOuterBoundary( cpm_FaceFlag face, int procGrpNo=0 );

  /** 自ランクの境界が内部境界(隣が不活性ドメイン)かどうかを判定
   *  @param[in] face      面方向
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @retval    true      内部境界
   *  @retval    false     内部境界でない
   */
  bool IsInnerBoundary( cpm_FaceFlag face, int procGrpNo=0 );



////// MPI処理のインターフェイス関数 //////

  /** MPI_Datatypeを取得
   *  @param[in] ptr 取得したいデータのポインタ
   *  @return MPI_Datatype
   */
  template<class T> CPM_INLINE
  static MPI_Datatype GetMPI_Datatype(T *ptr);

  /** MPI_Datatypeを取得
   *  - FortranデータタイプからMPI_Datatypeを取得
   *  @param[in] datatype 取得したいデータのポインタ
   *  @return MPI_Datatype
   */
  static MPI_Datatype GetMPI_Datatype(int datatype);

  /** MPI_Opを取得
   *  - FortranオペレータタイプからMPI_Opを取得
   *  @param[in] op 取得したいデータのポインタ
   *  @return MPI_Op
   */
  static MPI_Op GetMPI_Op(int op);

  /** ランク番号の取得
   *  - MPI_PROC_NULLが返ってきた場合は、
   *    1.プロセスグループが存在しない、
   *    2.プロセスグループに自ランクが含まれていない、
   *    のいずれか
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return ランク番号
   */
  int GetMyRankID( int procGrpNo=0 );

  /** ランク数の取得
   *  - プロセスグループのランク数を取得する
   *  @param[in] procGrpNo プロセスグループ番号(省略時0)
   *  @return ランク数
   */
  int GetNumRank( int procGrpNo=0 );

  /** ホスト名の取得
   *  - 自ランクのホスト名を取得
   *  @return ホスト名
   */
  std::string GetHostName();

  /** MPIコミュニケータの取得
   *  - MPI_COMM_NULLが返ってきた場合は、
   *    1.プロセスグループが存在しない、
   *    2.プロセスグループに自ランクが含まれていない、
   *    のいずれか
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return MPIコミュニケータ 
   */
  MPI_Comm GetMPI_Comm( int procGrpNo=0 );

  /** Abort
   *  - MPI_Abortのインターフェイス
   *  @param[in] errorcode MPI_Abortに渡すエラーコード
   */
  void Abort( int errorcode );

  /** Barrier
   *  - MPI_Barrierのインターフェイス
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Barrier( int procGrpNo=0 );

  /** Wait
   *  - MPI_Waitのインターフェイス
   *  @param[in] request リクエストハンドル
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Wait( MPI_Request *request );

  /** Waitall
   *  - MPI_Waitallのインターフェイス
   *  @param[in] count    リクエストの数
   *  @param[in] requests リクエストハンドル配列
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Waitall( int count, MPI_Request requests[] );

  /** Bcast
   *  - MPI_Bcastのインターフェイス
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
   *  @param[in] reqNo リクエスト番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_Wait( int reqNo );

  /** cpm_Waitall
   *  - MPI_Waitallのインターフェイス
   *  @param[in] count     リクエストの数
   *  @param[in] reqNoList リクエスト番号のリスト
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_Waitall( int count, int reqNoList[] );

  /** cpm_Isend
   *  - MPI_Isendのインターフェイス
   *  - Fortranインターフェイス用
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
   *  @param[out] buf       受信データ
   *  @param[in]  count     受信データのサイズ
   *  @param[in]  datatype  受信データのデータタイプ(cpm_fparam.fi参照)
   *  @param[in]  source    送信元のランク番号(procGrpNo内でのランク番号)
   *  @param[out] reqNo     リクエスト番号(Fortran用)
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode cpm_Irecv( void *buf, int count, int datatype, int source, int *reqNo, int procGrpNo=0 );

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





////// 袖通信関数 //////

  /** 袖通信バッファのセット
   *  - 6face分の送受信バッファを確保する
   *  @param[in] maxVC     送受信バッファの最大袖数
   *  @param[in] maxN      送受信バッファの最大成分数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  virtual
  cpm_ErrorCode SetBndCommBuffer( size_t maxVC, size_t maxN, int procGrpNo=0 );

  /** 袖通信バッファサイズの取得
   *  - 袖通信バッファとして確保されている配列サイズ(byte)を返す
   *  @param[in] procGrpNo プロセスグループ番号(負の場合、全プロセスグループでのトータルを返す)
   *  @return バッファサイズ(byte)
   */
  virtual
  size_t GetBndCommBufferSize( int procGrpNo=0 ); 

  /** 袖通信(Scalar3D版)
   *  - (imax,jmax,kmax)の形式の配列の袖通信を行う
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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                          , int vc, int vc_comm, int procGrpNo=0 );

  /** 袖通信(Vector3D版)
   *  - (imax,jmax,kmax,3)の形式の配列の袖通信を行う
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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                          , int vc, int vc_comm, int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                          , int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                          , int vc, int vc_comm, int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS3D_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS3D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommV3D_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommV3D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4D_nowait( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                 , MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                                 , int vc, int vc_comm, MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                               , MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode wait_BndCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                               , int vc, int vc_comm, MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                               , MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode wait_BndCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                               , int vc, int vc_comm, MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                               , MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode wait_BndCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                               , int vc, int vc_comm, MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                               , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                               , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                               , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );

  /** 袖通信(Vector3DEx版)
   *  - (3,imax,jmax,kmax)の形式の配列の袖通信を行う
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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                            , int vc, int vc_comm, int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                            , int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                            , int vc, int vc_comm, int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommV3DEx_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                   , MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommV3DEx_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                   , int vc, int vc_comm, MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4DEx_nowait( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                   , MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4DEx_nowait( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                   , int vc, int vc_comm, MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode wait_BndCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode wait_BndCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, MPI_Request req[48], int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode PeriodicCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );

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
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo=0 );





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
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  double* AllocDoubleS3D( int vc, int procGrpNo=0 );

  /** 配列確保 float(imax,jmax,kmax)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  float* AllocFloatS3D( int vc, int procGrpNo=0 );

  /** 配列確保 int(imax,jmax,kmax)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  int* AllocIntS3D( int vc, int procGrpNo=0 );

  /** 配列確保 double(imax,jmax,kmax,3)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  double* AllocDoubleV3D( int vc, int procGrpNo=0 );

  /** 配列確保 float(imax,jmax,kmax,3)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  float* AllocFloatV3D( int vc, int procGrpNo=0 );

  /** 配列確保 int(imax,jmax,kmax,3)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  int* AllocIntV3D( int vc, int procGrpNo=0 );

  /** 配列確保 double(3,imax,jmax,kmax)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  double* AllocDoubleV3DEx( int vc, int procGrpNo=0 );

  /** 配列確保 float(3,imax,jmax,kmax)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  float* AllocFloatV3DEx( int vc, int procGrpNo=0 );

  /** 配列確保 int(3,imax,jmax,kmax)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  int* AllocIntV3DEx( int vc, int procGrpNo=0 );

  /** 配列確保 double(imax,jmax,kmax,nmax)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  double* AllocDoubleS4D( int nmax, int vc, int procGrpNo=0 );

  /** 配列確保 float(imax,jmax,kmax,nmax)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  float* AllocFloatS4D( int nmax, int vc, int procGrpNo=0 );

  /** 配列確保 int(imax,jmax,kmax,nmax)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  int* AllocIntS4D( int nmax, int vc, int procGrpNo=0 );

  /** 配列確保 double(nmax,imax,jmax,kmax)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  double* AllocDoubleS4DEx( int nmax, int vc, int procGrpNo=0 );

  /** 配列確保 float(nmax,imax,jmax,kmax)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  float* AllocFloatS4DEx( int nmax, int vc, int procGrpNo=0 );

  /** 配列確保 int(nmax,imax,jmax,kmax)
   *  @param[in] vc        仮想セル数
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 配列ポインタ
   */
  int* AllocIntS4DEx( int nmax, int vc, int procGrpNo=0 );




  /** flush */
  void flush(std::ostream &out, int procGrpNo=0);

  /** flush */
  void flush(FILE *fp, int procGrpNo=0);

#ifdef _DEBUG
  void printVoxelInfo(int myrank=-1);
#endif

protected:

  /** コンストラクタ */
  cpm_ParaManager();

  /** デストラクタ */
  virtual ~cpm_ParaManager();





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

  /** プロセスグループ毎のVOXEL空間情報マップ
   *  - VOXEL空間番号をキーとしたVOXEL空間情報マップ
   *  - 自ランクが含まれるVOXEL空間のみを管理する
   */
  VoxelInfoMap m_voxelInfoMap;

  /** MPI_Requestの管理マップ
   *  - Fortranインターフェイス用
   */
  cpm_ObjList<MPI_Request> m_reqList;
};

#include "cpm_ParaManagerCART.h"
#include "cpm_ParaManagerLMR.h"

//インライン関数
#include "inline/cpm_ParaManager_inline.h"
#include "inline/cpm_ParaManager_BndComm.h"
#include "inline/cpm_ParaManager_BndCommEx.h"

#endif /* _CPM_PARAMANAGER_H_ */
