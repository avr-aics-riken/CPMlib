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
 * @file   cpm_LeafCommInfo.cpp
 * LMRの袖通信情報管理クラスソースファイル
 * @date   2016/07/29
 */

#include "cpm_LeafCommInfo.h"

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_LeafCommInfo::cpm_LeafCommInfo(int distRankID)
{
  m_iDistRankNo = distRankID;
  m_vecCommInfo.clear();
  m_pCommSendBuf = NULL;
  m_pCommRecvBuf = NULL;
  m_reqSend = MPI_REQUEST_NULL;
  m_reqRecv = MPI_REQUEST_NULL;
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_LeafCommInfo::~cpm_LeafCommInfo()
{
  for( int i=0;i<m_vecCommInfo.size();i++ )
  {
    delete m_vecCommInfo[i];
  }
  m_vecCommInfo.clear();

  delete [] m_pCommSendBuf;
  delete [] m_pCommRecvBuf;
  m_pCommSendBuf = NULL;
  m_pCommRecvBuf = NULL;
}

////////////////////////////////////////////////////////////////////////////////
// CommInfoを追加
void
cpm_LeafCommInfo::AddCommInfo(stCommInfo *commInfo)
{
  m_vecCommInfo.push_back(commInfo);
}

////////////////////////////////////////////////////////////////////////////////
// CommInfoリストのソート
// type ソートタイプ(0:自身のリーフ番号でソート, 1:相手のリーフ番号でソート)
// 0のときは自身のリーフ番号を優先１、相手のリーフ番号を優先２としてソートする
// 1のときは相手のリーフ番号を優先１、自身のリーフ番号を優先２としてソートする
void
cpm_LeafCommInfo::Sort( int type )
{
  if( m_vecCommInfo.size() == 0 )
  {
    return;
  }

  // 優先１でソート
  Qsort(type, m_vecCommInfo, 0, (int)m_vecCommInfo.size()-1);

  // 優先２でソート
  int type2 = 1 - type;
  int iStart = 0;
  int iEnd   = iStart;
  while( iEnd < m_vecCommInfo.size()-1 )
  {
    // 優先１の同じリーフIDのインデクス範囲を検索
    int leafID = m_vecCommInfo[iStart]->GetLeafID(type);
    for( int i=iStart+1;i<(int)m_vecCommInfo.size();i++ )
    {
      if( m_vecCommInfo[i]->GetLeafID(type) == leafID )
      {
        iEnd = i;
      }
      else
      {
        break;
      }
    }

    // 優先２でソート
    if( iEnd > iStart )
    {
      Qsort(type2, m_vecCommInfo, iStart, iEnd);
    }

    // 次の範囲を検索
    iStart = iEnd + 1;
    iEnd = iStart;
  }
}

////////////////////////////////////////////////////////////////////////////////
// CommInfoリストのクイックソート
void
cpm_LeafCommInfo::Qsort( int type, std::vector<stCommInfo*> &vecCommInfo, int iStart, int iEnd )
{
  if( iStart >= iEnd ) return;                      //終了番号が開始番号以下の場合、関数を抜ける
  int iBaseNumber = (iStart + iEnd) / 2;            //中央のインデックスを求める
  stCommInfo *BaseValue = vecCommInfo[iBaseNumber]; //配列の真ん中の値を基準値にする
  vecCommInfo[iBaseNumber] = vecCommInfo[iStart];   //中央の要素に開始番号の値を格納
  int iCounter = iStart;                            //格納位置カウンタを開始番号と同じにする

  //大小比較と入れ替え
  for( int i=iStart+1;i<=iEnd;i++ ) //開始番号の次の要素から終了番号までループ
  {
    if( vecCommInfo[i]->GetLeafID(type) < BaseValue->GetLeafID(type) ) //大小で比較
    {
      iCounter++;                             //格納位置カウンタをインクリメント
      stCommInfo *Buffer = vecCommInfo[iCounter];         //[i] と [vntCounter] の値をスワップ
      vecCommInfo[iCounter] = vecCommInfo[i];
      vecCommInfo[i] = Buffer;
    }
  }

  vecCommInfo[iStart] = vecCommInfo[iCounter]; //[iCounter]を開始番号の値にする
  vecCommInfo[iCounter] = BaseValue;           //基準値を[iCounter]に格納

  Qsort(type, vecCommInfo, iStart, iCounter-1); //分割された配列をクイックソート(再帰)
  Qsort(type, vecCommInfo, iCounter+1, iEnd);   //分割された配列をクイックソート(再帰) 
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信バッファのセット
bool
cpm_LeafCommInfo::SetBndCommBuffer( int myRankNo, size_t sz_face[2], size_t maxVC, size_t maxN )
{
  m_CommSendBufSize = 0;
  m_CommRecvBufSize = 0;
  m_pCommSendBuf = NULL;
  m_pCommRecvBuf = NULL;

  // 各経路ごとに計算、加算
  for( int i=0;i<m_vecCommInfo.size();i++ )
  {
    // 送受信バッファサイズをそれぞれ加算
    m_CommSendBufSize += m_vecCommInfo[i]->CalcSendBufferSize(sz_face, maxVC, maxN);
    m_CommRecvBufSize += m_vecCommInfo[i]->CalcRecvBufferSize(sz_face, maxVC, maxN);
  }

  m_pCommSendBuf = new REAL_BUF_TYPE[m_CommSendBufSize];
  m_pCommRecvBuf = new REAL_BUF_TYPE[m_CommRecvBufSize];
  if( !m_pCommSendBuf || !m_pCommRecvBuf )
  {
    return false;
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////
// 対となる通信情報を検索
cpm_LeafCommInfo::stCommInfo* cpm_LeafCommInfo::SearchDistCommInfo(cpm_LeafCommInfo::stCommInfo *commInfo)
{
  for( size_t i=0;i<m_vecCommInfo.size();i++ )
  {
    stCommInfo *commInfoDist = m_vecCommInfo[i];

    // OwnLeafとDistLeafが入れ替わりでbPeriodicが同じ通信情報かどうか
    if( commInfoDist->iOwnLeafID  == commInfo->iDistLeafID &&
        commInfoDist->iDistLeafID == commInfo->iOwnLeafID  &&
        commInfoDist->bPeriodic   == commInfo->bPeriodic   )
    {
      return commInfoDist;
    }
  }

  return NULL;
}

