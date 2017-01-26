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
 * @file   cpm_LeafCommInfo.h
 * LMRの袖通信情報管理クラスヘッダー
 * @date   2016/07/29
 */

#ifndef _CPM_LEAFCOMMINFO_H_
#define _CPM_LEAFCOMMINFO_H_

#include "cpm_Base.h"
#include "cpm_DefineLMR.h"
#include <vector>

/** LMRの袖通信情報管理クラス
 *  - 現時点ではユーザがインスタンスすることを許していない
 *  - get_instance静的関数を用いて唯一のインスタンスを取得する
 */
class cpm_LeafCommInfo : public cpm_Base
{
public:
  /** １通信経路情報構造体 */
  struct stCommInfo
  {
    /// 自身のリーフID
    int iOwnLeafID;

    /// 通信相手のリーフID
    int iDistLeafID;

    /// 通信相手とのレベル差
    int iLevelDiff;

    /// 自身のfaceIdx
    int iFaceIdx;

    /// 周期境界フラグ
    bool bPeriodic;

    /// コンストラクタ
    stCommInfo()
    {
      iOwnLeafID  = -1;
      iDistLeafID = -1;
      iLevelDiff  = 0;
      iFaceIdx    = 0;
      bPeriodic   = false;
    }

    /** LeafIDを取得
     *  @param[in] type ソートタイプ(0:iOwnLeafID, 1:iDistLeafID)
     */
    int GetLeafID( int type )
    {
      if( type==0 )
      {
        return iOwnLeafID;
      }
      else if( type==1 )
      {
        return iDistLeafID;
      }
      return -1;
    }

    /** 必要な送信バッファサイズを計算
     *  @param[in] sz_face   1リーフの格子数(平面内２軸)
     *  @param[in] vc_comm   通信する仮想セル数
     *  @param[in] nmax      送信バッファの最大成分数
     */
    size_t CalcSendBufferSize(size_t sz_face[2], size_t vc_comm, size_t nmax)
    {
      size_t jmax = sz_face[0];
      size_t kmax = sz_face[1];

      size_t sz = 0;
      if( iLevelDiff == 0 )
      {
        sz = size_t(jmax  +2*vc_comm) * size_t(kmax  +2*vc_comm) *    vc_comm  * nmax;
      }
      else if( iLevelDiff == 1 )
      {
        sz = size_t(jmax/2+2*vc_comm) * size_t(kmax/2+2*vc_comm) *    vc_comm  * nmax;
      }
      else if( iLevelDiff == -1 )
      {
        sz = size_t(jmax  +4*vc_comm) * size_t(kmax  +4*vc_comm) * (2*vc_comm) * nmax;
      }

      return sz;
    }

    /** 必要な受信バッファサイズを計算
     *  @param[in] sz_face   1リーフの格子数(平面内２軸)
     *  @param[in] vc_comm   通信する仮想セル数 
     *  @param[in] nmax      受信バッファの最大成分数
     */
    size_t CalcRecvBufferSize(size_t sz_face[2], size_t vc_comm, size_t nmax)
    {
      size_t jmax = sz_face[0];
      size_t kmax = sz_face[1];

      size_t sz = 0;
      if( iLevelDiff == 0 )
      {
        sz = size_t(jmax  +2*vc_comm) * size_t(kmax  +2*vc_comm) *    vc_comm  * nmax;
      }
      else if( iLevelDiff == 1 )
      {
        sz = size_t(jmax  +4*vc_comm) * size_t(kmax  +4*vc_comm) * (2*vc_comm) * nmax;
      }
      else if( iLevelDiff == -1 )
      {
        sz = size_t(jmax/2+2*vc_comm) * size_t(kmax/2+2*vc_comm) *    vc_comm  * nmax;
      }

      return sz;
    }

  };

////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:

  /** コンストラクタ */
  cpm_LeafCommInfo(int distRankID);

  /** デストラクタ */
  virtual ~cpm_LeafCommInfo();

  /** CommInfoを追加 */
  void AddCommInfo(stCommInfo *commInfo);

  /** CommInfoリストのソート
   *  @param[in] type ソートタイプ(0:自身のリーフ番号でソート, 1:相手のリーフ番号でソート)
   */
  void Sort( int type );

  /** 袖通信バッファのセット
   *  @param[in] myRankNo  自身のランク番号
   *  @param[in] sz_face   1リーフの格子数(平面内２軸)
   *  @param[in] maxVC     送受信バッファの最大袖数
   *  @param[in] maxN      送受信バッファの最大成分数
   */
  bool SetBndCommBuffer( int myRankNo, size_t sz[2], size_t maxVC, size_t maxN );

  /** 袖通信の送信バッファの取得
   *  @return バッファのポインタ
   */
  void* GetBndCommSendBufferPtr() { return m_pCommSendBuf; };

  /** 袖通信の受信バッファの取得
   *  @return バッファのポインタ
   */
  void* GetBndCommRecvBufferPtr() { return m_pCommRecvBuf; };

  /** 袖通信バッファサイズの取得
   *  @return バッファサイズ(word)
   */
  size_t GetBndCommBufferSize() { return (m_CommSendBufSize + m_CommRecvBufSize); };

  /** 対となる通信情報を検索 \n
   *  (OwnLeafとDistLeafが入れ替わりでbPeriodicが同じ通信情報)
   *  @param[in] commInfo    検索したい元の通信情報
   *  @return 対となる通信情報(存在しない場合NULL)
   */
  stCommInfo* SearchDistCommInfo(stCommInfo *commInfo);



protected:

  
  /** CommInfoリストのクイックソート
   *  @param[in]    type        ソートタイプ(0:自身のリーフ番号でソート, 1:相手のリーフ番号でソート)
   *  @param[inout] vecCommInfo ソート対象の配列
   *  @param[in]    startIndex  開始インデクス
   *  @param[in]    endIndex    終了インデクス
   */
  void Qsort( int type, std::vector<stCommInfo*> &vecCommInfo, int startIndex, int endIndex );




////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:

  /// 通信相手のランク番号
  int m_iDistRankNo;

  /// 経路情報
  std::vector<stCommInfo*> m_vecCommInfo;

  // 送信リクエストID
  MPI_Request m_reqSend;

  // 受信リクエストID
  MPI_Request m_reqRecv;

protected:

  /// 送信バッファ
  REAL_BUF_TYPE *m_pCommSendBuf;

  /// 受信バッファ
  REAL_BUF_TYPE *m_pCommRecvBuf;

  /// 送信バッファサイズ(WORD数)
  size_t m_CommSendBufSize;

  /// 受信バッファサイズ(WORD数)
  size_t m_CommRecvBufSize;


};

#endif _CPM_LEAFCOMMINFO_H_

