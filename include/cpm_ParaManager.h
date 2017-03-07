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
 * @file   cpm_ParaManager.h
 * カーテシアン用のパラレルマネージャクラスのヘッダーファイル
 * @date   2015/03/27
 */

#ifndef _CPM_PARAMANAGER_H_
#define _CPM_PARAMANAGER_H_

#include "cpm_BaseParaManager.h"

/** プロセスグループ毎のVOXEL空間情報管理マップ */
typedef std::map<int, cpm_VoxelInfo*> VoxelInfoMap;

/** 袖通信バッファ情報 */
struct S_BNDCOMM_BUFFER
{
  size_t m_maxVC; ///< 最大袖数
  size_t m_maxN;  ///< 最大成分数
  size_t m_nwX;   ///< バッファサイズ
  size_t m_nwY;   ///< バッファサイズ
  size_t m_nwZ;   ///< バッファサイズ
  REAL_BUF_TYPE *m_bufX[4]; ///< バッファ
  REAL_BUF_TYPE *m_bufY[4]; ///< バッファ
  REAL_BUF_TYPE *m_bufZ[4]; ///< バッファ

  S_BNDCOMM_BUFFER()
  {
    m_maxVC = m_maxN = 0;
    m_nwX = m_nwY = m_nwZ = 0;
    for( int i=0;i<4;i++ )
    {
      m_bufX[i] = NULL;
      m_bufY[i] = NULL;
      m_bufZ[i] = NULL;
    }
  }

  ~S_BNDCOMM_BUFFER()
  {
    for( int i=0;i<4;i++ )
    {
      if( m_bufX[i] ) delete [] m_bufX[i];
      if( m_bufY[i] ) delete [] m_bufY[i];
      if( m_bufZ[i] ) delete [] m_bufZ[i];
    }
  }

  /** バッファサイズの計算
   *  @return バッファサイズ[Byte]
   */
  size_t CalcBufferSize()
  {
    return (m_nwX*4 + m_nwY*4 + m_nwZ*4) * sizeof(REAL_BUF_TYPE);
  }
};

/** カーテシアン用の並列管理クラス
 *  - cpm_BaseParaManagerクラスからの派生
 *  - get_instance関数の引数のdomainTypeがCPM_DOMAIN_CARTESIANのとき、
 *  このクラスがインスタンスされる
 */
class cpm_ParaManager : public cpm_BaseParaManager
{
friend class cpm_BaseParaManager;
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
   *  @return インスタンスのポインタ
   */
  static cpm_ParaManager* get_instance(int &argc, char**& argv);


public:
  /** 領域分割(FVM用)
   *  - 既に作成済みの領域分割情報を用いた領域分割処理
   *
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

  /** 領域分割(FVM用)
   *  - 領域分割の各種情報を引数で渡して領域分割を行う
   *  - プロセスグループの全てのランクが活性ドメインになる
   *  - I,J,K方向の領域分割数を指定するバージョン
   *
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

  /** 領域分割(FVM用)
   *  - 領域分割の各種情報を引数で渡して領域分割を行う
   *  - プロセスグループの全てのランクが活性ドメインになる
   *  - 並列数=プロセスグループの並列数とし、内部で自動的に領域分割をするバージョン
   *
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

// 2016/01/22 FEAST add.s
  /** 領域分割(FDM用)
   *  - 領域分割の各種情報を引数で渡して領域分割を行う
   *  - プロセスグループの全てのランクが活性ドメインになる
   *  - I,J,K方向の領域分割数を指定するバージョン
   *
   *  @param[in] div        領域分割数
   *  @param[in] nod        空間全体の頂点数
   *  @param[in] origin     空間全体の原点
   *  @param[in] region     空間全体のサイズ
   *  @param[in] maxVC      最大の袖数(袖通信用)
   *  @param[in] maxN       最大の成分数(袖通信用)
   *  @param[in] divPolicy  自動分割ポリシー
   *  @param[in] procGrpNo  領域分割を行うプロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  virtual
  cpm_ErrorCode NodeInit ( int div[3], int nod[3], double origin[3], double region[3]
                         , size_t maxVC=1, size_t maxN=3, cpm_DivPolicy divPolicy=DIV_COMM_SIZE
                         , int procGrpNo=0 );

  /** 領域分割(FDM用)
   *  - 領域分割の各種情報を引数で渡して領域分割を行う
   *  - プロセスグループの全てのランクが活性ドメインになる
   *  - 並列数=プロセスグループの並列数とし、内部で自動的に領域分割をするバージョン
   *
   *  @param[in] nod        空間全体の頂点数
   *  @param[in] origin     空間全体の原点
   *  @param[in] region     空間全体のサイズ
   *  @param[in] maxVC      最大の袖数(袖通信用)
   *  @param[in] maxN       最大の成分数(袖通信用)
   *  @param[in] divPolicy  自動分割ポリシー
   *  @param[in] procGrpNo  領域分割を行うプロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  virtual
  cpm_ErrorCode NodeInit ( int nod[3], double origin[3], double region[3]
                         , size_t maxVC=1, size_t maxN=3, cpm_DivPolicy divPolicy=DIV_COMM_SIZE
                         , int procGrpNo=0 );
// 2016/01/22 FEAST add.e

  /** 領域分割(ActiveSubdomain指定)(FVM用)
   *  - 領域分割の各種情報を引数で渡して領域分割を行う
   *  - ActiveSubdomainファイルで指定される領域分割位置のランクが活性ドメインになる
   *  - I,J,K方向の領域分割数を指定するバージョン
   *  - 指定の領域分割数とActiveSubdomainファイルで指定されている領域分割数が一致している必要がある
   *  - ActiveSubdomain数と並列数が一致している必要がある
   *
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

  /** 領域分割(ActiveSubdomain指定)(FVM用)
   *  - 領域分割の各種情報を引数で渡して領域分割を行う
   *  - ActiveSubdomainファイルで指定される領域分割位置のランクが活性ドメインになる
   *  - ActiveSubdomainファイルで指定されている領域分割数で領域分割を行う
   *  - ActiveSubdomain数と並列数が一致している必要がある
   *
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

// 2016/01/22 FEAST add.s

  /** 領域分割(ActiveSubdomain指定)(FDM用)
   *  - 領域分割の各種情報を引数で渡して領域分割を行う
   *  - ActiveSubdomainファイルで指定される領域分割位置のランクが活性ドメインになる
   *  - I,J,K方向の領域分割数を指定するバージョン
   *  - 指定の領域分割数とActiveSubdomainファイルで指定されている領域分割数が一致している必要がある
   *  - ActiveSubdomain数と並列数が一致している必要がある
   *
   *  @param[in] div           領域分割数
   *  @param[in] nod           空間全体の頂点数
   *  @param[in] origin        空間全体の原点
   *  @param[in] region        空間全体のサイズ
   *  @param[in] subDomainFile ActiveSubdomainファイル名
   *  @param[in] maxVC         最大の袖数(袖通信用)
   *  @param[in] maxN          最大の成分数(袖通信用)
   *  @param[in] procGrpNo     領域分割を行うプロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  virtual
  cpm_ErrorCode NodeInit_Subdomain( int div[3], int nod[3], double origin[3], double region[3]
                                   , std::string subDomainFile
                                   , size_t maxVC=1, size_t maxN=3, int procGrpNo=0 );

  /** 領域分割(ActiveSubdomain指定)(FDM用)
   *  - 領域分割の各種情報を引数で渡して領域分割を行う
   *  - ActiveSubdomainファイルで指定される領域分割位置のランクが活性ドメインになる
   *  - ActiveSubdomainファイルで指定されている領域分割数で領域分割を行う
   *  - ActiveSubdomain数と並列数が一致している必要がある
   *
   *  @param[in] nod           空間全体の頂点数
   *  @param[in] origin        空間全体の原点
   *  @param[in] region        空間全体のサイズ
   *  @param[in] subDomainFile ActiveSubdomainファイル名
   *  @param[in] maxVC         最大の袖数(袖通信用)
   *  @param[in] maxN          最大の成分数(袖通信用)
   *  @param[in] procGrpNo     領域分割を行うプロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  virtual
  cpm_ErrorCode NodeInit_Subdomain( int nod[3], double origin[3], double region[3]
                                   , std::string subDomainFile
                                   , size_t maxVC=1, size_t maxN=3, int procGrpNo=0 );

// 2016/01/22 FEAST add.e




////// 領域情報の取得関数 //////

  /** VOXEL空間マップを検索
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return VOXEL空間情報ポインタ
   */
  virtual
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
   *
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 自ランクの始点インデクス整数配列のポインタ
   */
  const int* GetVoxelHeadIndex( int procGrpNo=0 );

  /** 自ランクの始点頂点の全体空間でのインデクスを取得
   *  - 全体空間の先頭インデクスを0としたC型のインデクス
   *
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 自ランクの始点インデクス整数配列のポインタ
   */
  const int* GetNodeHeadIndex( int procGrpNo=0 );

  /** 自ランクの始点VOXELまたは頂点の全体空間でのインデクスを取得
   *  - FVMのときはボクセル数、FDMのときは頂点数を取得
   *  - 全体空間の先頭インデクスを0としたC型のインデクス
   *
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 自ランクの始点インデクス整数配列のポインタ
   */
  const int* GetArrayHeadIndex( int procGrpNo=0 );

  /** 自ランクの終点VOXELの全体空間でのインデクスを取得
   *  - 全体空間の先頭インデクスを0としたC型のインデクス
   *
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 自ランクの終点インデクス整数配列のポインタ
   */
  const int* GetVoxelTailIndex( int procGrpNo=0 );

  /** 自ランクの終点頂点の全体空間でのインデクスを取得
   *  - 全体空間の先頭インデクスを0としたC型のインデクス
   *
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 自ランクの終点インデクス整数配列のポインタ
   */
  const int* GetNodeTailIndex( int procGrpNo=0 );

  /** 自ランクの終点ボクセルまたは頂点の全体空間でのインデクスを取得
   *  - FVMのときはボクセル数、FDMのときは頂点数を取得
   *  - 全体空間の先頭インデクスを0としたC型のインデクス
   *
   *  @param[in] procGrpNo プロセスグループ番号(省略時=0)
   *  @return 自ランクの終点インデクス整数配列のポインタ
   */
  const int* GetArrayTailIndex( int procGrpNo=0 );

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

  /** 指定idを含む全体ボクセル空間のインデクス範囲を取得
   *  - 全体空間実セルのスタートインデクスを0としたときの，i,j,k各方向の
   *    スタートインデクスと長さを取得する．
   *
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
   *  @param[in]  padding   パディングフラグ
   *  @retval true  指定idを含むセルが存在した
   *  @retval false 指定idを含むセルが存在しない
   */
  bool GetBndIndexExtGc( int id, int *array, int vc
                       , int &ista, int &jsta, int &ksta, int &ilen, int &jlen, int &klen
                       , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 指定idを含む全体ボクセル空間のインデクス範囲を取得
   *  - 全体空間実セルのスタートインデクスを0としたときの，i,j,k各方向の
   *    スタートインデクスと長さを取得する．
   *
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
   *  @param[in]  pad_size  パディングサイズ(i,j,k)
   *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
   *  @retval true  指定idを含むセルが存在した
   *  @retval false 指定idを含むセルが存在しない
   */
  bool GetBndIndexExtGc( int id, int *array, int imax, int jmax, int kmax, int vc
                       , int &ista, int &jsta, int &ksta, int &ilen, int &jlen, int &klen
                       , int pad_size[3], int procGrpNo=0 );

  /** グローバルインデクスからローカルインデクスを計算
   *  - 全体空間実セルのスタートインデクスを0としたときの，i,j,k各方向の
   *    ローカルインデクスの計算と自領域に含まれるかの判定を行う．
   *
   *  @param[in]  iG        グローバルのiインデクス
   *  @param[in]  jG        グローバルのjインデクス
   *  @param[in]  kG        グローバルのkインデクス
   *  @param[out] iL        ローカルのiインデクス
   *  @param[out] jL        ローカルのiインデクス
   *  @param[out] kL        ローカルのiインデクス
   *  @param[in]  procGrpNo プロセスグループ番号(省略時=0)
   *  @retval     true      自領域の実計算セルに含まれる
   *  @retval     false     自領域の実計算セルに含まれない
   */
  bool Global2LocalIndex( int iG, int jG, int kG, int &iL, int &jL, int &kL, int procGrpNo=0 );

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
  cpm_ErrorCode SetBndCommBuffer( size_t maxVC, size_t maxN, int procGrpNo=0 );

  /** 袖通信バッファサイズの取得
   *  - 袖通信バッファとして確保されている配列サイズ(byte)を返す
   *
   *  @param[in] procGrpNo プロセスグループ番号(負の場合、全プロセスグループでのトータルを返す)
   *  @return バッファサイズ(byte)
   */
  virtual
  size_t GetBndCommBufferSize( int procGrpNo=0 );

  /** 袖通信(Scalar3D版)
   *  - (imax,jmax,kmax)の形式の配列の袖通信を行う
   *
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
   *
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
   *
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
   *
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
   *
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                          , int vc, int vc_comm, int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 袖通信(Scalar4D版, パディングサイズ指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の袖通信を行う
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                               , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Vector3D版)
   *  - (imax,jmax,kmax,3)の形式の配列の周期境界方向の袖通信を行う
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                               , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Scalar4D版)
   *  - (imax,jmax,kmax,nmax)の形式の配列の周期境界方向の袖通信を行う
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                               , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Scalar4D版, パディングサイズ指定)
   *  - (imax,jmax,kmax,nmax)の形式の配列の周期境界方向の袖通信を行う
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
   *  @param[in]    pad_size  パディングサイズ(i,j,k,n)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                               , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                               , int pad_size[4], int procGrpNo );

  /** 袖通信(Vector3DEx版)
   *  - (3,imax,jmax,kmax)の形式の配列の袖通信を行う
   *
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
   *
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
   *
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode BndCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                            , int vc, int vc_comm, int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 袖通信(Scalar4DEx版, パディングサイズ指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の袖通信を行う
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                                 , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Scalar4DEx版)
   *  - (nmax,imax,jmax,kmax)の形式の配列の周期境界方向の袖通信を行う
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
   *  @param[in]    padding   パディングフラグ(true:ON、false:OFF)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                                 , int procGrpNo=0, CPM_PADDING padding=CPM_PADDING_OFF );

  /** 周期境界袖通信(Scalar4DEx版, パディングサイズ指定)
   *  - (nmax,imax,jmax,kmax)の形式の配列の周期境界方向の袖通信を行う
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
   *  @param[in]    pad_size  パディングサイズ(n,i,j,k)
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode PeriodicCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                 , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm
                                 , int pad_size[4], int procGrpNo=0 );





////// MPI処理のFortran用インターフェイス関数 //////

  /** cpm_BndCommS3D_nowait
   *  - BndCommS3D_nowaitのインターフェイス
   *  - Fortranインターフェイス用
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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
   *
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





protected:

  /** プロセスグループ毎の袖通信バッファ情報マップの定義 */
  typedef std::map<int, S_BNDCOMM_BUFFER*> BndCommInfoMap;


  /** コンストラクタ */
  cpm_ParaManager();

  /** デストラクタ */
  virtual ~cpm_ParaManager();

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

  /** 並列プロセス数からI,J,K方向の分割数を取得する
   *  - 通信面のトータルサイズが小さい分割パターンを採用する
   *
   *  @param[in]  divNum  ランク数
   *  @param[in]  voxSize 空間全体のボクセル数
   *  @param[out] divPttn 領域分割数
   *  @return             終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode
  DecideDivPattern_CommSize( int divNum
                           , int voxSize[3]
                           , int divPttn[3] ) const;

  /** I,J,K分割を行った時の通信点数の総数を取得する
   *  @param[in] iDiv    i方向領域分割数
   *  @param[in] jDiv    j方向領域分割数
   *  @param[in] kDiv    k方向領域分割数
   *  @param[in] voxSize 空間全体のボクセル数
   *  @return            袖通信点数
   */
  unsigned long long
  CalcCommSize( unsigned long long iDiv
              , unsigned long long jDiv
              , unsigned long long kDiv
              , unsigned long long voxSize[3] ) const;

  /** 並列プロセス数からI,J,K方向の分割数を取得する
   *  - １つのサブドメインが立方体に一番近い分割パターンを採用する
   *
   *  @param[in]  divNum  ランク数
   *  @param[in]  voxSize 空間全体のボクセル数
   *  @param[out] divPttn 領域分割数
   *  @return             終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode
  DecideDivPattern_Cube( int divNum
                       , int voxSize[3]
                       , int divPttn[3] ) const;

  /** I,J,K分割を行った時のI,J,Kボクセル数の最大/最小の差を取得する
   *  @param[in] iDiv    i方向領域分割数
   *  @param[in] jDiv    j方向領域分割数
   *  @param[in] kDiv    k方向領域分割数
   *  @param[in] voxSize 空間全体のボクセル数
   *  @retval 0以上      I,J,Kボクセル数の最大/最小の差
   *  @retval 負値       領域分割不可のパターン
   */
  long long
  CheckCube( unsigned long long iDiv
           , unsigned long long jDiv
           , unsigned long long kDiv
           , unsigned long long voxSize[3] ) const;

  /** 袖通信バッファの取得
   *  - 袖通信バッファ情報の取得
   *
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 袖通信バッファ情報のポインタ
   */
  CPM_INLINE
  S_BNDCOMM_BUFFER* GetBndCommBuffer( int procGrpNo=0 )
  {
    BndCommInfoMap::iterator it = m_bndCommInfoMap.find(procGrpNo);
    if( it == m_bndCommInfoMap.end() ) return NULL;
    return it->second;
  }

  /** 袖通信(Scalar3D,4D,Vector3D版)のX方向送信バッファのセット
   *  @param[in]  array    袖通信をする配列の先頭ポインタ
   *  @param[in]  imax     配列サイズ(I方向)
   *  @param[in]  jmax     配列サイズ(J方向)
   *  @param[in]  kmax     配列サイズ(K方向)
   *  @param[in]  nmax     配列サイズ(成分数)
   *  @param[in]  vc       仮想セル数
   *  @param[in]  vc_comm  通信する仮想セル数
   *  @param[in]  pad_size パディングサイズ(i,j,k,n)
   *  @param[out] sendm    マイナス方向の送信バッファ
   *  @param[out] sendp    プラス方向の送信バッファ
   *  @param[in]  nIDm     マイナス方向の隣接ランク番号
   *  @param[in]  nIDp     プラス方向の隣接ランク番号
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode packX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int pad_size[4]
// 2016/01/22 FEAST mod.s
//                   , T *sendm, T *sendp, int nIDm, int nIDp );
                     , T *sendm, T *sendp, int nIDm, int nIDp , int procGrpNo);
// 2016/01/22 FEAST mod.e

  /** 袖通信(Scalar3D,4D,Vector3D版)のX方向受信バッファを元に戻す
   *  @param[inout] array    袖通信をした配列の先頭ポインタ
   *  @param[in]    imax     配列サイズ(I方向)
   *  @param[in]    jmax     配列サイズ(J方向)
   *  @param[in]    kmax     配列サイズ(K方向)
   *  @param[in]    nmax     配列サイズ(成分数)
   *  @param[in]    vc       仮想セル数
   *  @param[in]    vc_comm  通信する仮想セル数
   *  @param[in]    pad_size パディングサイズ(i,j,k,n)
   *  @param[in]    recvm    マイナス方向の受信バッファ
   *  @param[in]    recvp    プラス方向の受信バッファ
   *  @param[in]    nIDm     マイナス方向の隣接ランク番号
   *  @param[in]    nIDp     プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode unpackX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int pad_size[4]
                       , T *recvm, T *recvp, int nIDm, int nIDp );

  /** 袖通信(Scalar3D,4D,Vector3D版)のY方向送信バッファのセット
   *  @param[in]  array    袖通信をする配列の先頭ポインタ
   *  @param[in]  imax     配列サイズ(I方向)
   *  @param[in]  jmax     配列サイズ(J方向)
   *  @param[in]  kmax     配列サイズ(K方向)
   *  @param[in]  nmax     配列サイズ(成分数)
   *  @param[in]  vc       仮想セル数
   *  @param[in]  vc_comm  通信する仮想セル数
   *  @param[in]  pad_size パディングサイズ(i,j,k,n)
   *  @param[out] sendm    マイナス方向の送信バッファ
   *  @param[out] sendp    プラス方向の送信バッファ
   *  @param[in]  nIDm     マイナス方向の隣接ランク番号
   *  @param[in]  nIDp     プラス方向の隣接ランク番号
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode packY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int pad_size[4]
// 2016/01/22 FEAST mod.s
//                   , T *sendm, T *sendp, int nIDm, int nIDp );
                     , T *sendm, T *sendp, int nIDm, int nIDp, int procGrpNo );
// 2016/01/22 FEAST mod.e

  /** 袖通信(Scalar3D,4D,Vector3D版)のY方向受信バッファを元に戻す
   *  @param[inout] array    袖通信をした配列の先頭ポインタ
   *  @param[in]    imax     配列サイズ(I方向)
   *  @param[in]    jmax     配列サイズ(J方向)
   *  @param[in]    kmax     配列サイズ(K方向)
   *  @param[in]    nmax     配列サイズ(成分数)
   *  @param[in]    vc       仮想セル数
   *  @param[in]    vc_comm  通信する仮想セル数
   *  @param[in]    pad_size パディングサイズ(i,j,k,n)
   *  @param[in]    recvm    マイナス方向の受信バッファ
   *  @param[in]    recvp    プラス方向の受信バッファ
   *  @param[in]    nIDm     マイナス方向の隣接ランク番号
   *  @param[in]    nIDp     プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode unpackY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int pad_size[4]
                       , T *recvm, T *recvp, int nIDm, int nIDp );

  /** 袖通信(Scalar3D,4D,Vector3D版)のZ方向送信バッファのセット
   *  @param[in]  array    袖通信をする配列の先頭ポインタ
   *  @param[in]  imax     配列サイズ(I方向)
   *  @param[in]  jmax     配列サイズ(J方向)
   *  @param[in]  kmax     配列サイズ(K方向)
   *  @param[in]  nmax     配列サイズ(成分数)
   *  @param[in]  vc       仮想セル数
   *  @param[in]  vc_comm  通信する仮想セル数
   *  @param[in]  pad_size パディングサイズ(i,j,k,n)
   *  @param[out] sendm    マイナス方向の送信バッファ
   *  @param[out] sendp    プラス方向の送信バッファ
   *  @param[in]  nIDm     マイナス方向の隣接ランク番号
   *  @param[in]  nIDp     プラス方向の隣接ランク番号
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode packZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int pad_size[4]
// 2016/01/22 FEAST mod.s
//                   , T *sendm, T *sendp, int nIDm, int nIDp );
                     , T *sendm, T *sendp, int nIDm, int nIDp, int procGrpNo );
// 2016/01/22 FEAST mod.e

  /** 袖通信(Scalar3D,4D,Vector3D版)のZ方向受信バッファを元に戻す
   *  @param[inout] array    袖通信をした配列の先頭ポインタ
   *  @param[in]    imax     配列サイズ(I方向)
   *  @param[in]    jmax     配列サイズ(J方向)
   *  @param[in]    kmax     配列サイズ(K方向)
   *  @param[in]    nmax     配列サイズ(成分数)
   *  @param[in]    vc       仮想セル数
   *  @param[in]    vc_comm  通信する仮想セル数
   *  @param[in]    pad_size パディングサイズ(i,j,k,n)
   *  @param[in]    recvm    マイナス方向の受信バッファ
   *  @param[in]    recvp    プラス方向の受信バッファ
   *  @param[in]    nIDm     マイナス方向の隣接ランク番号
   *  @param[in]    nIDp     プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode unpackZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int pad_size[4]
                       , T *recvm, T *recvp, int nIDm, int nIDp );

  /** 袖通信(Scalar4DEx,Vector3DEx版)のX方向送信バッファのセット
   *  @param[in]  array   袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax    配列サイズ(成分数)
   *  @param[in]  imax    配列サイズ(I方向)
   *  @param[in]  jmax    配列サイズ(J方向)
   *  @param[in]  kmax    配列サイズ(K方向)
   *  @param[in]  vc      仮想セル数
   *  @param[in]  vc_comm 通信する仮想セル数
   *  @param[in]  pad_size パディングサイズ(n,i,j,k)
   *  @param[out] sendm   マイナス方向の送信バッファ
   *  @param[out] sendp   プラス方向の送信バッファ
   *  @param[in]  nIDm    マイナス方向の隣接ランク番号
   *  @param[in]  nIDp    プラス方向の隣接ランク番号
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode packXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm, int pad_size[4]
// 2016/01/22 FEAST mod.s
//                     , T *sendm, T *sendp, int nIDm, int nIDp );
                       , T *sendm, T *sendp, int nIDm, int nIDp, int procGrpNo );
// 2016/01/22 FEAST mod.e

  /** 袖通信(Scalar4DEx,Vector3DEx版)のX方向受信バッファを元に戻す
   *  @param[inout] array   袖通信をした配列の先頭ポインタ
   *  @param[in]    nmax    配列サイズ(成分数)
   *  @param[in]    imax    配列サイズ(I方向)
   *  @param[in]    jmax    配列サイズ(J方向)
   *  @param[in]    kmax    配列サイズ(K方向)
   *  @param[in]    vc      仮想セル数
   *  @param[in]    vc_comm 通信する仮想セル数
   *  @param[in]    pad_size パディングサイズ(n,i,j,k)
   *  @param[in]    recvm   マイナス方向の受信バッファ
   *  @param[in]    recvp   プラス方向の受信バッファ
   *  @param[in]    nIDm    マイナス方向の隣接ランク番号
   *  @param[in]    nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode unpackXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm, int pad_size[4]
                         , T *recvm, T *recvp, int nIDm, int nIDp );

  /** 袖通信(Scalar4DEx,Vector3DEx版)のY方向送信バッファのセット
   *  @param[in]  array   袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax    配列サイズ(成分数)
   *  @param[in]  imax    配列サイズ(I方向)
   *  @param[in]  jmax    配列サイズ(J方向)
   *  @param[in]  kmax    配列サイズ(K方向)
   *  @param[in]  vc      仮想セル数
   *  @param[in]  vc_comm 通信する仮想セル数
   *  @param[in]  pad_size パディングサイズ(n,i,j,k)
   *  @param[out] sendm   マイナス方向の送信バッファ
   *  @param[out] sendp   プラス方向の送信バッファ
   *  @param[in]  nIDm    マイナス方向の隣接ランク番号
   *  @param[in]  nIDp    プラス方向の隣接ランク番号
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode packYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm, int pad_size[4]
// 2016/01/22 FEAST mod.s
//                     , T *sendm, T *sendp, int nIDm, int nIDp );
                       , T *sendm, T *sendp, int nIDm, int nIDp, int procGrpNo );
// 2016/01/22 FEAST mod.e

  /** 袖通信(Scalar4DEx,Vector3DEx版)のY方向受信バッファを元に戻す
   *  @param[inout] array   袖通信をした配列の先頭ポインタ
   *  @param[in]    nmax    配列サイズ(成分数)
   *  @param[in]    imax    配列サイズ(I方向)
   *  @param[in]    jmax    配列サイズ(J方向)
   *  @param[in]    kmax    配列サイズ(K方向)
   *  @param[in]    vc      仮想セル数
   *  @param[in]    vc_comm 通信する仮想セル数
   *  @param[in]    pad_size パディングサイズ(n,i,j,k)
   *  @param[in]    recvm   マイナス方向の受信バッファ
   *  @param[in]    recvp   プラス方向の受信バッファ
   *  @param[in]    nIDm    マイナス方向の隣接ランク番号
   *  @param[in]    nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode unpackYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm, int pad_size[4]
                         , T *recvm, T *recvp, int nIDm, int nIDp );

  /** 袖通信(Scalar4DEx,Vector3DEx版)のZ方向送信バッファのセット
   *  @param[in]  array   袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax    配列サイズ(成分数)
   *  @param[in]  imax    配列サイズ(I方向)
   *  @param[in]  jmax    配列サイズ(J方向)
   *  @param[in]  kmax    配列サイズ(K方向)
   *  @param[in]  vc      仮想セル数
   *  @param[in]  vc_comm 通信する仮想セル数
   *  @param[in]  pad_size パディングサイズ(n,i,j,k)
   *  @param[out] sendm   マイナス方向の送信バッファ
   *  @param[out] sendp   プラス方向の送信バッファ
   *  @param[in]  nIDm    マイナス方向の隣接ランク番号
   *  @param[in]  nIDp    プラス方向の隣接ランク番号
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode packZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm, int pad_size[4]
// 2016/01/22 FEAST mod.s
//                     , T *sendm, T *sendp, int nIDm, int nIDp );
                       , T *sendm, T *sendp, int nIDm, int nIDp, int procGrpNo );
// 2016/01/22 FEAST mod.e

  /** 袖通信(Scalar4DEx,Vector3DEx版)のZ方向受信バッファを元に戻す
   *  @param[inout] array   袖通信をした配列の先頭ポインタ
   *  @param[in]    imax    配列サイズ(I方向)
   *  @param[in]    jmax    配列サイズ(J方向)
   *  @param[in]    kmax    配列サイズ(K方向)
   *  @param[in]    nmax    配列サイズ(成分数)
   *  @param[in]    vc      仮想セル数
   *  @param[in]    vc_comm 通信する仮想セル数
   *  @param[in]    pad_size パディングサイズ(n,i,j,k)
   *  @param[in]    recvm   マイナス方向の受信バッファ
   *  @param[in]    recvp   プラス方向の受信バッファ
   *  @param[in]    nIDm    マイナス方向の隣接ランク番号
   *  @param[in]    nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode unpackZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm, int pad_size[4]
                         , T *recvm, T *recvp, int nIDm, int nIDp );

  /** １方向(プラス、マイナス)の双方向袖通信処理
   *  @param[in]  sendm     マイナス方向の送信バッファ
   *  @param[in]  sendp     プラス方向の送信バッファ
   *  @param[in]  recvm     マイナス方向の受信バッファ
   *  @param[in]  recvp     プラス方向の受信バッファ
   *  @param[in]  nw        送受信サイズ
   *  @param[out] req       MPI_Request配列のポインタ(サイズ4)
   *  @param[in]  nIDsm     マイナス方向受信用の隣接ランク番号
   *  @param[in]  nIDrm     マイナス方向送信用の隣接ランク番号
   *  @param[in]  nIDsp     プラス方向受信用の隣接ランク番号
   *  @param[in]  nIDrp     プラス方向送信用の隣接ランク番号
   *  @param[in]  procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode sendrecv( T *sendm, T *recvm, T *sendp, T *recvp, size_t nw, MPI_Request *req
                        , int nIDsm, int nIDrm, int nIDsp, int nIDrp, int procGrpNo=0 );




////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:

protected:

  /** プロセスグループ毎のVOXEL空間情報マップ
   *  - VOXEL空間番号をキーとしたVOXEL空間情報マップ
   *  - 自ランクが含まれるVOXEL空間のみを管理する
   */
  VoxelInfoMap m_voxelInfoMap;

  /** プロセスグループ毎の袖通信バッファ情報
   */
  BndCommInfoMap m_bndCommInfoMap;
};

//インライン関数
#include "inline/cpm_ParaManager_BndComm.h"
#include "inline/cpm_ParaManager_BndCommEx.h"

#endif /* _CPM_PARAMANAGER_H_ */
