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
 * @file   cpm_ParaManagerCART.h
 * カーテシアン用のパラレルマネージャクラスのヘッダーファイル
 * @date   2015/03/27
 */

#ifndef _CPM_PARAMANAGER_CART_H_
#define _CPM_PARAMANAGER_CART_H_

#include "cpm_ParaManager.h"

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
 *  cpm_ParaManagerクラスからの派生
 *  get_instance関数の引数のdomainTypeがCPM_DOMAIN_CARTESIANのとき、
 *  このクラスがインスタンスされる
 */
class cpm_ParaManagerCART : public cpm_ParaManager
{
friend class cpm_ParaManager;
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:

  /** 領域分割
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

  /** 領域分割
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

  /** 領域分割
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

  /** 領域分割(ActiveSubdomain指定)
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

  /** 領域分割(ActiveSubdomain指定)
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





////// 領域情報の取得関数 //////

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
   *  @param[out]  req       MPIリクエスト
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4D_nowait( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                 , MPI_Request req[12], int procGrpNo=0 );

  /** 非同期版袖通信のwait、展開(Scalar4D版)
   *  - (imax,jmax,kmax,nmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                               , MPI_Request req[12], int procGrpNo=0 );

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
   *  @param[out]  req       MPIリクエスト
   *  @param[in]   procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode BndCommS4DEx_nowait( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                   , MPI_Request req[12], int procGrpNo=0 );

  /** 非同期版袖通信のwait、展開(Scalar4DEx版)
   *  - (nmax,imax,jmax,kmax)の形式の配列の非同期版袖通信のwaitと展開を行う
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    req       MPIリクエスト
   *  @param[in]    procGrpNo プロセスグループ番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T> CPM_INLINE
  cpm_ErrorCode wait_BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                 , MPI_Request req[12], int procGrpNo=0 );

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





protected:

  /** プロセスグループ毎の袖通信バッファ情報マップの定義 */
  typedef std::map<int, S_BNDCOMM_BUFFER*> BndCommInfoMap;


  /** コンストラクタ */
  cpm_ParaManagerCART();

  /** デストラクタ */
  virtual ~cpm_ParaManagerCART();

  /** 並列プロセス数からI,J,K方向の分割数を取得する
   *  通信面のトータルサイズが小さい分割パターンを採用する
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
              , unsigned long long voxsize[3] ) const;

  /** 並列プロセス数からI,J,K方向の分割数を取得する
   *  １つのサブドメインが立方体に一番近い分割パターンを採用する
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
           , unsigned long long voxsize[3] ) const;

  /** 袖通信バッファの取得
   *  - 袖通信バッファ情報の取得
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
   *  @param[in]  array   袖通信をする配列の先頭ポインタ
   *  @param[in]  imax    配列サイズ(I方向)
   *  @param[in]  jmax    配列サイズ(J方向)
   *  @param[in]  kmax    配列サイズ(K方向)
   *  @param[in]  nmax    配列サイズ(成分数)
   *  @param[in]  vc      仮想セル数
   *  @param[in]  vc_comm 通信する仮想セル数
   *  @param[out] sendm   マイナス方向の送信バッファ
   *  @param[out] sendp   プラス方向の送信バッファ
   *  @param[in]  nIDm    マイナス方向の隣接ランク番号
   *  @param[in]  nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode packX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                     , T *sendm, T *sendp, int nIDm, int nIDp );

  /** 袖通信(Scalar3D,4D,Vector3D版)のX方向受信バッファを元に戻す
   *  @param[inout] array   袖通信をした配列の先頭ポインタ
   *  @param[in]    imax    配列サイズ(I方向)
   *  @param[in]    jmax    配列サイズ(J方向)
   *  @param[in]    kmax    配列サイズ(K方向)
   *  @param[in]    nmax    配列サイズ(成分数)
   *  @param[in]    vc      仮想セル数
   *  @param[in]    vc_comm 通信する仮想セル数
   *  @param[in]    recvm   マイナス方向の受信バッファ
   *  @param[in]    recvp   プラス方向の受信バッファ
   *  @param[in]    nIDm    マイナス方向の隣接ランク番号
   *  @param[in]    nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode unpackX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                       , T *recvm, T *recvp, int nIDm, int nIDp );

  /** 袖通信(Scalar3D,4D,Vector3D版)のY方向送信バッファのセット 
   *  @param[in]  array   袖通信をする配列の先頭ポインタ
   *  @param[in]  imax    配列サイズ(I方向)
   *  @param[in]  jmax    配列サイズ(J方向)
   *  @param[in]  kmax    配列サイズ(K方向)
   *  @param[in]  nmax    配列サイズ(成分数)
   *  @param[in]  vc      仮想セル数
   *  @param[in]  vc_comm 通信する仮想セル数
   *  @param[out] sendm   マイナス方向の送信バッファ
   *  @param[out] sendp   プラス方向の送信バッファ
   *  @param[in]  nIDm    マイナス方向の隣接ランク番号
   *  @param[in]  nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode packY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                     , T *sendm, T *sendp, int nIDm, int nIDp );

  /** 袖通信(Scalar3D,4D,Vector3D版)のY方向受信バッファを元に戻す
   *  @param[inout] array   袖通信をした配列の先頭ポインタ
   *  @param[in]    imax    配列サイズ(I方向)
   *  @param[in]    jmax    配列サイズ(J方向)
   *  @param[in]    kmax    配列サイズ(K方向)
   *  @param[in]    nmax    配列サイズ(成分数)
   *  @param[in]    vc      仮想セル数
   *  @param[in]    vc_comm 通信する仮想セル数
   *  @param[in]    recvm   マイナス方向の受信バッファ
   *  @param[in]    recvp   プラス方向の受信バッファ
   *  @param[in]    nIDm    マイナス方向の隣接ランク番号
   *  @param[in]    nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode unpackY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                       , T *recvm, T *recvp, int nIDm, int nIDp );

  /** 袖通信(Scalar3D,4D,Vector3D版)のZ方向送信バッファのセット 
   *  @param[in]  array   袖通信をする配列の先頭ポインタ
   *  @param[in]  imax    配列サイズ(I方向)
   *  @param[in]  jmax    配列サイズ(J方向)
   *  @param[in]  kmax    配列サイズ(K方向)
   *  @param[in]  nmax    配列サイズ(成分数)
   *  @param[in]  vc      仮想セル数
   *  @param[in]  vc_comm 通信する仮想セル数
   *  @param[out] sendm   マイナス方向の送信バッファ
   *  @param[out] sendp   プラス方向の送信バッファ
   *  @param[in]  nIDm    マイナス方向の隣接ランク番号
   *  @param[in]  nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode packZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                     , T *sendm, T *sendp, int nIDm, int nIDp );

  /** 袖通信(Scalar3D,4D,Vector3D版)のZ方向受信バッファを元に戻す
   *  @param[inout] array   袖通信をした配列の先頭ポインタ
   *  @param[in]    imax    配列サイズ(I方向)
   *  @param[in]    jmax    配列サイズ(J方向)
   *  @param[in]    kmax    配列サイズ(K方向)
   *  @param[in]    nmax    配列サイズ(成分数)
   *  @param[in]    vc      仮想セル数
   *  @param[in]    vc_comm 通信する仮想セル数
   *  @param[in]    recvm   マイナス方向の受信バッファ
   *  @param[in]    recvp   プラス方向の受信バッファ
   *  @param[in]    nIDm    マイナス方向の隣接ランク番号
   *  @param[in]    nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode unpackZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                       , T *recvm, T *recvp, int nIDm, int nIDp );

  /** 袖通信(Scalar4DEx,Vector3DEx版)のX方向送信バッファのセット 
   *  @param[in]  array   袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax    配列サイズ(成分数)
   *  @param[in]  imax    配列サイズ(I方向)
   *  @param[in]  jmax    配列サイズ(J方向)
   *  @param[in]  kmax    配列サイズ(K方向)
   *  @param[in]  vc      仮想セル数
   *  @param[in]  vc_comm 通信する仮想セル数
   *  @param[out] sendm   マイナス方向の送信バッファ
   *  @param[out] sendp   プラス方向の送信バッファ
   *  @param[in]  nIDm    マイナス方向の隣接ランク番号
   *  @param[in]  nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode packXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                       , T *sendm, T *sendp, int nIDm, int nIDp );

  /** 袖通信(Scalar4DEx,Vector3DEx版)のX方向受信バッファを元に戻す
   *  @param[inout] array   袖通信をした配列の先頭ポインタ
   *  @param[in]    nmax    配列サイズ(成分数)
   *  @param[in]    imax    配列サイズ(I方向)
   *  @param[in]    jmax    配列サイズ(J方向)
   *  @param[in]    kmax    配列サイズ(K方向)
   *  @param[in]    vc      仮想セル数
   *  @param[in]    vc_comm 通信する仮想セル数
   *  @param[in]    recvm   マイナス方向の受信バッファ
   *  @param[in]    recvp   プラス方向の受信バッファ
   *  @param[in]    nIDm    マイナス方向の隣接ランク番号
   *  @param[in]    nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode unpackXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                         , T *recvm, T *recvp, int nIDm, int nIDp );

  /** 袖通信(Scalar4DEx,Vector3DEx版)のY方向送信バッファのセット 
   *  @param[in]  array   袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax    配列サイズ(成分数)
   *  @param[in]  imax    配列サイズ(I方向)
   *  @param[in]  jmax    配列サイズ(J方向)
   *  @param[in]  kmax    配列サイズ(K方向)
   *  @param[in]  vc      仮想セル数
   *  @param[in]  vc_comm 通信する仮想セル数
   *  @param[out] sendm   マイナス方向の送信バッファ
   *  @param[out] sendp   プラス方向の送信バッファ
   *  @param[in]  nIDm    マイナス方向の隣接ランク番号
   *  @param[in]  nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode packYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                       , T *sendm, T *sendp, int nIDm, int nIDp );

  /** 袖通信(Scalar4DEx,Vector3DEx版)のY方向受信バッファを元に戻す
   *  @param[inout] array   袖通信をした配列の先頭ポインタ
   *  @param[in]    nmax    配列サイズ(成分数)
   *  @param[in]    imax    配列サイズ(I方向)
   *  @param[in]    jmax    配列サイズ(J方向)
   *  @param[in]    kmax    配列サイズ(K方向)
   *  @param[in]    vc      仮想セル数
   *  @param[in]    vc_comm 通信する仮想セル数
   *  @param[in]    recvm   マイナス方向の受信バッファ
   *  @param[in]    recvp   プラス方向の受信バッファ
   *  @param[in]    nIDm    マイナス方向の隣接ランク番号
   *  @param[in]    nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode unpackYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                         , T *recvm, T *recvp, int nIDm, int nIDp );

  /** 袖通信(Scalar4DEx,Vector3DEx版)のZ方向送信バッファのセット 
   *  @param[in]  array   袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax    配列サイズ(成分数)
   *  @param[in]  imax    配列サイズ(I方向)
   *  @param[in]  jmax    配列サイズ(J方向)
   *  @param[in]  kmax    配列サイズ(K方向)
   *  @param[in]  vc      仮想セル数
   *  @param[in]  vc_comm 通信する仮想セル数
   *  @param[out] sendm   マイナス方向の送信バッファ
   *  @param[out] sendp   プラス方向の送信バッファ
   *  @param[in]  nIDm    マイナス方向の隣接ランク番号
   *  @param[in]  nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode packZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                       , T *sendm, T *sendp, int nIDm, int nIDp );

  /** 袖通信(Scalar4DEx,Vector3DEx版)のZ方向受信バッファを元に戻す
   *  @param[inout] array   袖通信をした配列の先頭ポインタ
   *  @param[in]    imax    配列サイズ(I方向)
   *  @param[in]    jmax    配列サイズ(J方向)
   *  @param[in]    kmax    配列サイズ(K方向)
   *  @param[in]    nmax    配列サイズ(成分数)
   *  @param[in]    vc      仮想セル数
   *  @param[in]    vc_comm 通信する仮想セル数
   *  @param[in]    recvm   マイナス方向の受信バッファ
   *  @param[in]    recvp   プラス方向の受信バッファ
   *  @param[in]    nIDm    マイナス方向の隣接ランク番号
   *  @param[in]    nIDp    プラス方向の隣接ランク番号
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  template<class T>
  cpm_ErrorCode unpackZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
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

  /** プロセスグループ毎の袖通信バッファ情報
   */
  BndCommInfoMap m_bndCommInfoMap;
};

//インライン関数
#include "inline/cpm_ParaManager_BndComm_CART.h"
#include "inline/cpm_ParaManager_BndCommEx_CART.h"

#endif /* _CPM_PARAMANAGER_CART_H_ */
