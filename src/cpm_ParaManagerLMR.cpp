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
 * @file   cpm_ParaManagerLMR.cpp
 * パラレルマネージャクラスのソースファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */
//#include "cpm_ParaManagerLMR.h"
#include "cpm_ParaManager.h"
#include "cpm_VoxelInfoLMR.h"


////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_ParaManagerLMR::cpm_ParaManagerLMR()
  : cpm_ParaManager()
{
  // 領域分割タイプ
  m_domainType = CPM_DOMAIN_LMR;

  // 袖通信バッファ情報のクリア
   m_bndCommInfoMap.clear();
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_ParaManagerLMR::~cpm_ParaManagerLMR()
{
  // 袖通信バッファ情報の削除、クリア
  {
    BndCommInfoMapLMR::iterator it  = m_bndCommInfoMap.begin();
    BndCommInfoMapLMR::iterator ite = m_bndCommInfoMap.end();
    for( ; it!=ite; it++ )
    {
      if( it->second ) delete it->second;
    }
    m_bndCommInfoMap.clear();
  }
}

////////////////////////////////////////////////////////////////////////////////
// LMR用の領域分割
cpm_ErrorCode
cpm_ParaManagerLMR::VoxelInit_LMR( std::string treeFile
                              , size_t maxVC, size_t maxN, int procGrpNo )
{
  cpm_ErrorCode ret = CPM_SUCCESS;
  cpm_VoxelInfoLMR *voxelInfo = NULL;

  // 既に領域分割済みか
  if( m_voxelInfoMap.find(procGrpNo) != m_voxelInfoMap.end() )
  {
    return CPM_ERROR_ALREADY_VOXELINIIT;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm( procGrpNo );
  if( IsCommNull( comm ) )
  {
    return CPM_ERROR_MPI_INVALID_COMM;
  }

  // インスタンス
  voxelInfo = new cpm_VoxelInfoLMR();
  if( !voxelInfo )
  {
    Abort(CPM_ERROR_INVALID_PTR);
    return CPM_ERROR_INVALID_PTR;
  }

  // 領域分割情報の生成
  if( (ret = voxelInfo->Init( comm, treeFile )) != CPM_SUCCESS )
  {
    delete voxelInfo;
    Abort(ret);
    return ret;
  }

  // VOXEL空間マップに登録
    if( !m_voxelInfoMap.insert(std::make_pair(procGrpNo, voxelInfo)).second )
  {
    m_procGrpList.pop_back();
    delete voxelInfo;
    return CPM_ERROR_INSERT_VOXELMAP;
  }

  // 袖通信バッファの設定
  SetBndCommBuffer( maxVC, maxN, procGrpNo );

  return ret;
}

////////////////////////////////////////////////////////////////////////////////
// 木情報ファイルからリーフ数を取得する
int
cpm_ParaManagerLMR::GetNumLeaf( std::string treeFile )
{
   return cpm_VoxelInfoLMR::GetNumLeaf( treeFile );
}





























////////////////////////////////////////////////////////////////////////////////
// 袖通信バッファのセット
cpm_ErrorCode
cpm_ParaManagerLMR::SetBndCommBuffer( size_t maxVC, size_t maxN, int procGrpNo )
{
  if( maxVC==0 || maxN==0 )
  {
    return CPM_ERROR_BNDCOMM;
  }

  // local voxel size
  const int *sz = GetLocalVoxelSize(procGrpNo);
  if( !sz ) return CPM_ERROR_BNDCOMM_VOXELSIZE;

  // 隣接面ごとにサイズを計算
  const cpm_FaceFlag face[6] = {X_MINUS, X_PLUS, Y_MINUS, Y_PLUS, Z_MINUS, Z_PLUS};
  size_t sz_face[6][2] = { {sz[1], sz[2]}, {sz[1], sz[2]}
                         , {sz[2], sz[0]}, {sz[2], sz[0]}
                         , {sz[0], sz[1]}, {sz[0], sz[1]} };
  size_t nsend[6], nrecv[6]; //send/recv
  int nface[6];
  for( int i=0;i<6;i++ )
  {
    // 隣接領域とのレベル差
    int levelDiff = GetNeighborLevelDiff(face[i], procGrpNo);

    // レベル差ごとにサイズ計算
    if( levelDiff == 0 ) // 同じレベル
    {
      // 送信バッファサイズ
      // 同じvoxel数を送信
      nsend[i] = size_t(sz_face[i][0]+2*maxVC) * size_t(sz_face[i][1]+2*maxVC) * maxVC * maxN;

      // 受信バッファサイズ
      // 同じvoxel数を受信
      nrecv[i] = size_t(sz_face[i][0]+2*maxVC) * size_t(sz_face[i][1]+2*maxVC) * maxVC * maxN;

      // 相手は1面
      nface[i] = 1;
    }
    else if( levelDiff == 1 ) // 隣がfine
    {
      // 送信バッファサイズ
      // 1/4面を送信(相手は4面)
      nsend[i] = size_t(sz_face[i][0]/2+2*maxVC) * size_t(sz_face[i][1]/2+2*maxVC) * maxVC * maxN;

      // 受信バッファサイズ
      // 同じvoxel数で層数2倍を受信(相手は4面)
      nrecv[i] = size_t(sz_face[i][0]+4*maxVC) * size_t(sz_face[i][1]+4*maxVC) * (2*maxVC) * maxN;

      // 相手は4面
      nface[i] = 4;
    }
    else if( levelDiff == -1 ) // 隣がcoarse
    {
      // 送信バッファサイズ
      // 層数2倍を送信
      nsend[i] = size_t(sz_face[i][0]+4*maxVC) * size_t(sz_face[i][1]+4*maxVC) * (2*maxVC) * maxN;

      // 受信バッファサイズ
      // 1/4面を受信
      nrecv[i] = size_t(sz_face[i][0]/2+2*maxVC) * size_t(sz_face[i][1]/2+2*maxVC) * maxVC * maxN;

      // 相手は1面
      nface[i] = 1;
    }
    else
    {
      return CPM_ERROR_BNDCOMM;
    }
  }

  // 送受信バッファ情報のインスタンス
  S_BNDCOMM_BUFFER_LMR *bufInfo = new S_BNDCOMM_BUFFER_LMR();
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_ALLOC_BUFFER;
  }

  bufInfo->m_maxVC = maxVC;
  bufInfo->m_maxN  = maxN;
  for( int i=0;i<6;i++ )
  {
    bufInfo->m_nsend[i] = nsend[i];
    bufInfo->m_nrecv[i] = nrecv[i];
    bufInfo->m_nface[i] = nface[i];
    for( int j=0;j<nface[i];j++ )
    {
      bufInfo->m_bufSend[i][j] = new REAL_BUF_TYPE[nsend[i]];
      if( !bufInfo->m_bufSend[i][j] )
      {
        delete bufInfo;
        return CPM_ERROR_BNDCOMM_ALLOC_BUFFER;
      }
      bufInfo->m_bufRecv[i][j] = new REAL_BUF_TYPE[nrecv[i]];
      if( !bufInfo->m_bufRecv[i][j] )
      {
        delete bufInfo;
        return CPM_ERROR_BNDCOMM_ALLOC_BUFFER;
      }
    }
  }

  //マップにセット
  BndCommInfoMapLMR::iterator it = m_bndCommInfoMap.find(procGrpNo);
  if( it != m_bndCommInfoMap.end() )
  {
    delete it->second;
    m_bndCommInfoMap.erase( it );
  }
  m_bndCommInfoMap.insert( std::make_pair(procGrpNo, bufInfo) );

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信バッファサイズの取得
size_t
cpm_ParaManagerLMR::GetBndCommBufferSize( int procGrpNo )
{
  size_t mem = 0;

  // マップを検索
  if( procGrpNo < 0 )
  {
    //全体
    BndCommInfoMapLMR::iterator it = m_bndCommInfoMap.begin();
    for( ; it!=m_bndCommInfoMap.end(); it++ )
    {
      if( !it->second ) continue;
      mem += it->second->CalcBufferSize();
    }
  }
  else
  {
    // 指定プロセスグループを検索
    S_BNDCOMM_BUFFER_LMR *bInfo = GetBndCommBuffer( procGrpNo );
    if( !bInfo ) return 0;
    mem = bInfo->CalcBufferSize();
  }

  return mem;
}

