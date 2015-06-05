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
 * @author University of Tokyo
 * @date   2015/03/27
 */

#ifndef _CPM_PARAMANAGER_LMR_H_
#define _CPM_PARAMANAGER_LMR_H_

#include "cpm_ParaManager.h"

/** 袖通信バッファ情報 */
struct S_BNDCOMM_BUFFER_LMR
{
  size_t m_maxVC; ///< 最大袖数
  size_t m_maxN;  ///< 最大成分数
  size_t m_nsend[6];   ///< 送信バッファサイズ
  size_t m_nrecv[6];   ///< 受信バッファサイズ
  size_t m_nface[6];   ///< 通信相手の面数
  REAL_BUF_TYPE *m_bufSend[6][4]; ///< 送信バッファ
  REAL_BUF_TYPE *m_bufRecv[6][4]; ///< 受信バッファ

  S_BNDCOMM_BUFFER_LMR()
  {
    m_maxVC = m_maxN = 0;
    for( int i=0;i<6;i++ )
    {
      m_nsend[i] = 0;
      m_nrecv[i] = 0;
      m_nface[i] = 0;
      for( int j=0;j<4;j++ )
      {
        m_bufSend[i][j] = NULL;
        m_bufRecv[i][j] = NULL;
      }
    }
  }

  ~S_BNDCOMM_BUFFER_LMR()
  {
    for( int i=0;i<6;i++ )
    {
      for( int j=0;j<4;j++ )
      {
        delete [] m_bufSend[i][j];
        delete [] m_bufRecv[i][j];
      }
    }
  }

  /** バッファサイズの取得
   *  @return バッファサイズ[Byte]
   */
  size_t CalcBufferSize()
  {
    size_t sz = 0;
    for( int i=0;i<6;i++ )
    {
      sz += m_nsend[i];
      sz += m_nrecv[i];
    }
    sz *= sizeof(REAL_BUF_TYPE);
    return sz;
  }
};

/** LMR用の並列管理クラス
 *  - 現時点ではユーザがインスタンスすることを許していない
 *  - get_instance静的関数を用いて唯一のインスタンスを取得する
 */
class cpm_ParaManagerLMR : public cpm_ParaManager
{
friend class cpm_ParaManager;
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:

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

  /** 木情報ファイルからリーフ数を取得する
   *  @param[in] treefile  木情報ファイル
   *  @return    リーフ数
   */
  static
  int GetNumLeaf( std::string treeFile );





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
                                 , MPI_Request req[48], int procGrpNo=0 );

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
                               , MPI_Request req[48], int procGrpNo=0 );

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
                                   , MPI_Request req[48], int procGrpNo=0 );

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
                                 , MPI_Request req[48], int procGrpNo=0 );

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
  typedef std::map<int, S_BNDCOMM_BUFFER_LMR*> BndCommInfoMapLMR;

  /** コンストラクタ */
  cpm_ParaManagerLMR();

  /** デストラクタ */
  virtual ~cpm_ParaManagerLMR();

  /** 袖通信バッファの取得
   *  - 袖通信バッファ情報の取得
   *  @param[in] procGrpNo プロセスグループ番号
   *  @return 袖通信バッファ情報のポインタ
   */
  CPM_INLINE
  S_BNDCOMM_BUFFER_LMR* GetBndCommBuffer( int procGrpNo=0 )
  {
    BndCommInfoMapLMR::iterator it = m_bndCommInfoMap.find(procGrpNo);
    if( it == m_bndCommInfoMap.end() ) return NULL;
    return it->second;
  }

  /** １方向(プラス、マイナス)の非同期受信処理
   *  @param[in]  nID       隣接ランクリスト([0]:マイナス側、[1]:プラス側)
   *  @param[in]  nFace     隣接ランク数  ([0]:マイナス側、[1]:プラス側)
   *  @param[in]  levelDiff 隣接領域とのレベル差([0]:マイナス側、[1]:プラス側)
   *  @param[in]  nw        1面あたりの受信サイズ([0]:マイナス側、[1]:プラス側)
   *  @param[in]  recvm     マイナス方向の受信バッファ
   *  @param[out] reqm      マイナス方向のMPI_Request配列のポインタ(サイズ4)
   *  @param[in]  recvp     プラス方向の受信バッファ
   *  @param[out] reqp      プラス方向のMPI_Request配列のポインタ(サイズ4)
   *  @param[in]  procGrpNo プロセスグループ番号
   */
  template<class T>
  cpm_ErrorCode recv_LMR( const int nID[2][4], int nFace[2], int levelDiff[2], size_t nw[2]
                        , T* recvm[4], MPI_Request *reqm, T* recvp[4], MPI_Request *reqp
                        , int procGrpNo=0 );

  /** 袖通信(Scalar3D,4D,Vector3D版)の-X面への送信データのパック(通信面毎)
   *  @param[in]  array     袖通信をする配列の先頭ポインタ
   *  @param[in]  imax      配列サイズ(I方向)
   *  @param[in]  jmax      配列サイズ(J方向)
   *  @param[in]  kmax      配列サイズ(K方向)
   *  @param[in]  nmax      配列サイズ(成分数)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  vc_comm   通信する仮想セル数
   *  @param[in]  nID       隣接ランク番号
   *  @param[in]  faceNo    面番号(0/1/2/3)
   *  @param[in]  levelDiff 隣接領域とのレベル差
   *  @param[out] sendbuf   送信バッファ
   *  @param[in]  nw        送信バッファサイズ
   */
  template<class T>
  cpm_ErrorCode packMX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                      , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw );

  /** 袖通信(Scalar3D,4D,Vector3D版)の+X面への送信データのパック(通信面毎)
   *  @param[in]  array     袖通信をする配列の先頭ポインタ
   *  @param[in]  imax      配列サイズ(I方向)
   *  @param[in]  jmax      配列サイズ(J方向)
   *  @param[in]  kmax      配列サイズ(K方向)
   *  @param[in]  nmax      配列サイズ(成分数)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  vc_comm   通信する仮想セル数
   *  @param[in]  nID       隣接ランク番号
   *  @param[in]  faceNo    面番号(0/1/2/3)
   *  @param[in]  levelDiff 隣接領域とのレベル差
   *  @param[out] sendbuf   送信バッファ
   *  @param[in]  nw        送信バッファサイズ
   */
  template<class T>
  cpm_ErrorCode packPX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                      , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw );

  /** 袖通信(Scalar3D,4D,Vector3D版)の-X面からの受信データの展開(通信面毎)
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    nID       隣接ランク番号
   *  @param[in]    faceNo    面番号(0/1/2/3)
   *  @param[in]    levelDiff 隣接領域とのレベル差
   *  @param[in]    recvbuf   受信バッファ
   */
  template<class T>
  cpm_ErrorCode unpackMX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , const int nID, int faceNo, int levelDiff, T* recvbuf );

  /** 袖通信(Scalar3D,4D,Vector3D版)の+X面からの受信データの展開(通信面毎)
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    nID       隣接ランク番号
   *  @param[in]    faceNo    面番号(0/1/2/3)
   *  @param[in]    levelDiff 隣接領域とのレベル差
   *  @param[in]    recvbuf   受信バッファ
   */
  template<class T>
  cpm_ErrorCode unpackPX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , const int nID, int faceNo, int levelDiff, T* recvbuf );

  /** 袖通信(Scalar3D,4D,Vector3D版)の-Y面への送信データのパック(通信面毎)
   *  @param[in]  array     袖通信をする配列の先頭ポインタ
   *  @param[in]  imax      配列サイズ(I方向)
   *  @param[in]  jmax      配列サイズ(J方向)
   *  @param[in]  kmax      配列サイズ(K方向)
   *  @param[in]  nmax      配列サイズ(成分数)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  vc_comm   通信する仮想セル数
   *  @param[in]  nID       隣接ランク番号
   *  @param[in]  faceNo    面番号(0/1/2/3)
   *  @param[in]  levelDiff 隣接領域とのレベル差
   *  @param[out] sendbuf   送信バッファ
   *  @param[in]  nw        送信バッファサイズ
   */
  template<class T>
  cpm_ErrorCode packMY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                      , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw );

  /** 袖通信(Scalar3D,4D,Vector3D版)の+Y面への送信データのパック(通信面毎)
   *  @param[in]  array     袖通信をする配列の先頭ポインタ
   *  @param[in]  imax      配列サイズ(I方向)
   *  @param[in]  jmax      配列サイズ(J方向)
   *  @param[in]  kmax      配列サイズ(K方向)
   *  @param[in]  nmax      配列サイズ(成分数)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  vc_comm   通信する仮想セル数
   *  @param[in]  nID       隣接ランク番号
   *  @param[in]  faceNo    面番号(0/1/2/3)
   *  @param[in]  levelDiff 隣接領域とのレベル差
   *  @param[out] sendbuf   送信バッファ
   *  @param[in]  nw        送信バッファサイズ
   */
  template<class T>
  cpm_ErrorCode packPY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                      , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw );

  /** 袖通信(Scalar3D,4D,Vector3D版)の-Y面からの受信データの展開(通信面毎)
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    nID       隣接ランク番号
   *  @param[in]    faceNo    面番号(0/1/2/3)
   *  @param[in]    levelDiff 隣接領域とのレベル差
   *  @param[in]    recvbuf   受信バッファ
   */
  template<class T>
  cpm_ErrorCode unpackMY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , const int nID, int faceNo, int levelDiff, T* recvbuf );

  /** 袖通信(Scalar3D,4D,Vector3D版)の+Y面からの受信データの展開(通信面毎)
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    nID       隣接ランク番号
   *  @param[in]    faceNo    面番号(0/1/2/3)
   *  @param[in]    levelDiff 隣接領域とのレベル差
   *  @param[in]    recvbuf   受信バッファ
   */
  template<class T>
  cpm_ErrorCode unpackPY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , const int nID, int faceNo, int levelDiff, T* recvbuf );

  /** 袖通信(Scalar3D,4D,Vector3D版)の-Z面への送信データのパック(通信面毎)
   *  @param[in]  array     袖通信をする配列の先頭ポインタ
   *  @param[in]  imax      配列サイズ(I方向)
   *  @param[in]  jmax      配列サイズ(J方向)
   *  @param[in]  kmax      配列サイズ(K方向)
   *  @param[in]  nmax      配列サイズ(成分数)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  vc_comm   通信する仮想セル数
   *  @param[in]  nID       隣接ランク番号
   *  @param[in]  faceNo    面番号(0/1/2/3)
   *  @param[in]  levelDiff 隣接領域とのレベル差
   *  @param[out] sendbuf   送信バッファ
   *  @param[in]  nw        送信バッファサイズ
   */
  template<class T>
  cpm_ErrorCode packMZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                      , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw );

  /** 袖通信(Scalar3D,4D,Vector3D版)の+Z面への送信データのパック(通信面毎)
   *  @param[in]  array     袖通信をする配列の先頭ポインタ
   *  @param[in]  imax      配列サイズ(I方向)
   *  @param[in]  jmax      配列サイズ(J方向)
   *  @param[in]  kmax      配列サイズ(K方向)
   *  @param[in]  nmax      配列サイズ(成分数)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  vc_comm   通信する仮想セル数
   *  @param[in]  nID       隣接ランク番号
   *  @param[in]  faceNo    面番号(0/1/2/3)
   *  @param[in]  levelDiff 隣接領域とのレベル差
   *  @param[out] sendbuf   送信バッファ
   *  @param[in]  nw        送信バッファサイズ
   */
  template<class T>
  cpm_ErrorCode packPZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                      , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw );

  /** 袖通信(Scalar3D,4D,Vector3D版)の-Z面からの受信データの展開(通信面毎)
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    nID       隣接ランク番号
   *  @param[in]    faceNo    面番号(0/1/2/3)
   *  @param[in]    levelDiff 隣接領域とのレベル差
   *  @param[in]    recvbuf   受信バッファ
   */
  template<class T>
  cpm_ErrorCode unpackMZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , const int nID, int faceNo, int levelDiff, T* recvbuf );

  /** 袖通信(Scalar3D,4D,Vector3D版)の+Z面からの受信データの展開(通信面毎)
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    nID       隣接ランク番号
   *  @param[in]    faceNo    面番号(0/1/2/3)
   *  @param[in]    levelDiff 隣接領域とのレベル差
   *  @param[in]    recvbuf   受信バッファ
   */
  template<class T>
  cpm_ErrorCode unpackPZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                        , const int nID, int faceNo, int levelDiff, T* recvbuf );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の-X面への送信データのパック(通信面毎)
   *  @param[in]  array     袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax      配列サイズ(成分数)
   *  @param[in]  imax      配列サイズ(I方向)
   *  @param[in]  jmax      配列サイズ(J方向)
   *  @param[in]  kmax      配列サイズ(K方向)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  vc_comm   通信する仮想セル数
   *  @param[in]  nID       隣接ランク番号
   *  @param[in]  faceNo    面番号(0/1/2/3)
   *  @param[in]  levelDiff 隣接領域とのレベル差
   *  @param[out] sendbuf   送信バッファ
   *  @param[in]  nw        送信バッファサイズ
   */
  template<class T>
  cpm_ErrorCode packMXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の+X面への送信データのパック(通信面毎)
   *  @param[in]  array     袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax      配列サイズ(成分数)
   *  @param[in]  imax      配列サイズ(I方向)
   *  @param[in]  jmax      配列サイズ(J方向)
   *  @param[in]  kmax      配列サイズ(K方向)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  vc_comm   通信する仮想セル数
   *  @param[in]  nID       隣接ランク番号
   *  @param[in]  faceNo    面番号(0/1/2/3)
   *  @param[in]  levelDiff 隣接領域とのレベル差
   *  @param[out] sendbuf   送信バッファ
   *  @param[in]  nw        送信バッファサイズ
   */
  template<class T>
  cpm_ErrorCode packPXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の-X面からの受信データの展開(通信面毎)
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    nID       隣接ランク番号
   *  @param[in]    faceNo    面番号(0/1/2/3)
   *  @param[in]    levelDiff 隣接領域とのレベル差
   *  @param[in]    recvbuf   受信バッファ
   */
  template<class T>
  cpm_ErrorCode unpackMXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , const int nID, int faceNo, int levelDiff, T* recvbuf );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の+X面からの受信データの展開(通信面毎)
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    nID       隣接ランク番号
   *  @param[in]    faceNo    面番号(0/1/2/3)
   *  @param[in]    levelDiff 隣接領域とのレベル差
   *  @param[in]    recvbuf   受信バッファ
   */
  template<class T>
  cpm_ErrorCode unpackPXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , const int nID, int faceNo, int levelDiff, T* recvbuf );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の-Y面への送信データのパック(通信面毎)
   *  @param[in]  array     袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax      配列サイズ(成分数)
   *  @param[in]  imax      配列サイズ(I方向)
   *  @param[in]  jmax      配列サイズ(J方向)
   *  @param[in]  kmax      配列サイズ(K方向)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  vc_comm   通信する仮想セル数
   *  @param[in]  nID       隣接ランク番号
   *  @param[in]  faceNo    面番号(0/1/2/3)
   *  @param[in]  levelDiff 隣接領域とのレベル差
   *  @param[out] sendbuf   送信バッファ
   *  @param[in]  nw        送信バッファサイズ
   */
  template<class T>
  cpm_ErrorCode packMYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の+Y面への送信データのパック(通信面毎)
   *  @param[in]  array     袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax      配列サイズ(成分数)
   *  @param[in]  imax      配列サイズ(I方向)
   *  @param[in]  jmax      配列サイズ(J方向)
   *  @param[in]  kmax      配列サイズ(K方向)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  vc_comm   通信する仮想セル数
   *  @param[in]  nID       隣接ランク番号
   *  @param[in]  faceNo    面番号(0/1/2/3)
   *  @param[in]  levelDiff 隣接領域とのレベル差
   *  @param[out] sendbuf   送信バッファ
   *  @param[in]  nw        送信バッファサイズ
   */
  template<class T>
  cpm_ErrorCode packPYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の-Y面からの受信データの展開(通信面毎)
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    nID       隣接ランク番号
   *  @param[in]    faceNo    面番号(0/1/2/3)
   *  @param[in]    levelDiff 隣接領域とのレベル差
   *  @param[in]    recvbuf   受信バッファ
   */
  template<class T>
  cpm_ErrorCode unpackMYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , const int nID, int faceNo, int levelDiff, T* recvbuf );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の+Y面からの受信データの展開(通信面毎)
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    nID       隣接ランク番号
   *  @param[in]    faceNo    面番号(0/1/2/3)
   *  @param[in]    levelDiff 隣接領域とのレベル差
   *  @param[in]    recvbuf   受信バッファ
   */
  template<class T>
  cpm_ErrorCode unpackPYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , const int nID, int faceNo, int levelDiff, T* recvbuf );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の-Z面への送信データのパック(通信面毎)
   *  @param[in]  array     袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax      配列サイズ(成分数)
   *  @param[in]  imax      配列サイズ(I方向)
   *  @param[in]  jmax      配列サイズ(J方向)
   *  @param[in]  kmax      配列サイズ(K方向)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  vc_comm   通信する仮想セル数
   *  @param[in]  nID       隣接ランク番号
   *  @param[in]  faceNo    面番号(0/1/2/3)
   *  @param[in]  levelDiff 隣接領域とのレベル差
   *  @param[out] sendbuf   送信バッファ
   *  @param[in]  nw        送信バッファサイズ
   */
  template<class T>
  cpm_ErrorCode packMZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の+Z面への送信データのパック(通信面毎)
   *  @param[in]  array     袖通信をする配列の先頭ポインタ
   *  @param[in]  nmax      配列サイズ(成分数)
   *  @param[in]  imax      配列サイズ(I方向)
   *  @param[in]  jmax      配列サイズ(J方向)
   *  @param[in]  kmax      配列サイズ(K方向)
   *  @param[in]  vc        仮想セル数
   *  @param[in]  vc_comm   通信する仮想セル数
   *  @param[in]  nID       隣接ランク番号
   *  @param[in]  faceNo    面番号(0/1/2/3)
   *  @param[in]  levelDiff 隣接領域とのレベル差
   *  @param[out] sendbuf   送信バッファ
   *  @param[in]  nw        送信バッファサイズ
   */
  template<class T>
  cpm_ErrorCode packPZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の-Z面からの受信データの展開(通信面毎)
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    nID       隣接ランク番号
   *  @param[in]    faceNo    面番号(0/1/2/3)
   *  @param[in]    levelDiff 隣接領域とのレベル差
   *  @param[in]    recvbuf   受信バッファ
   */
  template<class T>
  cpm_ErrorCode unpackMZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , const int nID, int faceNo, int levelDiff, T* recvbuf );

  /** 袖通信(Scalar4DEx,Vector3DEx版)の+Z面からの受信データの展開(通信面毎)
   *  @param[inout] array     袖通信をする配列の先頭ポインタ
   *  @param[in]    nmax      配列サイズ(成分数)
   *  @param[in]    imax      配列サイズ(I方向)
   *  @param[in]    jmax      配列サイズ(J方向)
   *  @param[in]    kmax      配列サイズ(K方向)
   *  @param[in]    vc        仮想セル数
   *  @param[in]    vc_comm   通信する仮想セル数
   *  @param[in]    nID       隣接ランク番号
   *  @param[in]    faceNo    面番号(0/1/2/3)
   *  @param[in]    levelDiff 隣接領域とのレベル差
   *  @param[in]    recvbuf   受信バッファ
   */
  template<class T>
  cpm_ErrorCode unpackPZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , const int nID, int faceNo, int levelDiff, T* recvbuf );




////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:

protected:

  /** プロセスグループ毎の袖通信バッファ情報
   */
  BndCommInfoMapLMR m_bndCommInfoMap;

};

//インライン関数
#include "inline/cpm_ParaManager_BndComm_LMR.h"
#include "inline/cpm_ParaManager_BndCommEx_LMR.h"

#endif /* _CPM_PARAMANAGER_LMR_H_ */
