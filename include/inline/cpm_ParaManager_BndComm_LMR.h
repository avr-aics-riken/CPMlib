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
 * @file   cpm_ParaManager_BndComm_LMR.h
 * LMR用パラレルマネージャクラスのインラインヘッダーファイル
 * @date   2015/03/27
 */

#ifndef _CPM_PARAMANAGER_BNDCOMM_LMR_H_
#define _CPM_PARAMANAGER_BNDCOMM_LMR_H_

#define _IDXFX(_I,_J,_K,_N,_IS,_NJ,_NK,_VC) \
( size_t(_N)     * size_t(_VC) * size_t(_NJ+2*_VC) * size_t(_NK+2*_VC) \
+ size_t(_K+_VC) * size_t(_VC) * size_t(_NJ+2*_VC) \
+ size_t(_J+_VC) * size_t(_VC) \
+ size_t(_I-(_IS)) \
)

#define _IDXFY(_I,_J,_K,_N,_NI,_JS,_NK,_VC) \
( size_t(_N)       * size_t(_NI+2*_VC) * size_t(_VC) * size_t(_NK+2*_VC) \
+ size_t(_K+_VC)   * size_t(_NI+2*_VC) * size_t(_VC) \
+ size_t(_J-(_JS)) * size_t(_NI+2*_VC) \
+ size_t(_I+_VC) \
)

#define _IDXFZ(_I,_J,_K,_N,_NI,_NJ,_KS,_VC) \
( size_t(_N)       * size_t(_NI+2*_VC) * size_t(_NJ+2*_VC) * size_t(_VC) \
+ size_t(_K-(_KS)) * size_t(_NI+2*_VC) * size_t(_NJ+2*_VC) \
+ size_t(_J+_VC)   * size_t(_NI+2*_VC) \
+ size_t(_I+_VC) \
)

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // リクエスト格納配列
  MPI_Request req[48];
  for( int i=0;i<48;i++ ) req[i] = MPI_REQUEST_NULL;


  // 隣接領域とのレベル差
  int levelDiff[6];
  const cpm_FaceFlag face[6] = {X_MINUS, X_PLUS, Y_MINUS, Y_PLUS, Z_MINUS, Z_PLUS};
  for( int i=0;i<6;i++ )
  {
    levelDiff[face[i]] = GetNeighborLevelDiff(face[i], procGrpNo);
  }

  // 通信バッファを取得
  S_BNDCOMM_BUFFER_LMR *bufInfo = GetBndCommBuffer(procGrpNo);
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }

  //// X方向の通信 ////
  {
    // リクエスト格納配列
    MPI_Request *reqx = &(req[0]);

    // ポインタ
    T **sendm = (T**)(bufInfo->m_bufSend[X_MINUS]);
    T **recvm = (T**)(bufInfo->m_bufRecv[X_MINUS]);
    T **sendp = (T**)(bufInfo->m_bufSend[X_PLUS]);
    T **recvp = (T**)(bufInfo->m_bufRecv[X_PLUS]);

    // 隣接面
    cpm_FaceFlag faceFlag[2] = {X_MINUS, X_PLUS};

    // 隣接ランクとの通信情報
    int nID[2][4];     //隣接ランク番号リスト
    int nFace[2];      //隣接ランク数
    int lDiff[2];      //レベル差
    size_t nw_send[2]; //送信サイズ
    size_t nw_recv[2]; //受信サイズ
    for( int i=0;i<2;i++ )
    {
      // 隣接ランク番号
      const int *nIDtmp = GetNeighborRankList(faceFlag[i], nFace[i], procGrpNo );
      for( int j=0;j<4;j++ )
      {
        nID[i][j] = nIDtmp[j];
      }

      // レベル差とバッファサイズ
      lDiff[i] = levelDiff[faceFlag[i]];
      if( lDiff[i] == 0 )
      {
        nw_send[i] = size_t(jmax  +2*vc_comm) * size_t(kmax  +2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(jmax  +2*vc_comm) * size_t(kmax  +2*vc_comm) *    vc_comm  * nmax;
      }
      else if( lDiff[i] == 1 )
      {
        nw_send[i] = size_t(jmax/2+2*vc_comm) * size_t(kmax/2+2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(jmax  +4*vc_comm) * size_t(kmax  +4*vc_comm) * (2*vc_comm) * nmax;
      }
      else if( lDiff[i] == -1 )
      {
        nw_send[i] = size_t(jmax  +4*vc_comm) * size_t(kmax  +4*vc_comm) * (2*vc_comm) * nmax;
        nw_recv[i] = size_t(jmax/2+2*vc_comm) * size_t(kmax/2+2*vc_comm) *    vc_comm  * nmax;
      }
      if( nw_send[i] > bufInfo->m_nsend[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
      if( nw_recv[i] > bufInfo->m_nrecv[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for recv error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
    }

    // irecv
    if( (ret = recv_LMR(nID, nFace, lDiff, nw_recv, recvm, &reqx[0], recvp, &reqx[4], procGrpNo)) != CPM_SUCCESS ) return ret;

    // +側通信面ループ(send)
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // pack
          if( (ret = packPX( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], sendp[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendp[j], nw_send[i], nID[i][j], &reqx[8+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // -側通信面ループ(send)
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // pack
          if( (ret = packMX( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], sendm[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendm[j], nw_send[i], nID[i][j], &reqx[12+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // wait(recv)
    if( (ret = Waitall( 8, &reqx[0] )) != CPM_SUCCESS ) return ret;

    // -側unpack
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // unpack
          if( (ret = unpackMX( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], recvm[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // +側unpack
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // unpack
          if( (ret = unpackPX( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], recvp[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // wait(send)
    if( (ret = Waitall( 8, &reqx[8] )) != CPM_SUCCESS ) return ret;

  } // X方向完了


  //// Y方向の通信 ////
  {
    // リクエスト格納配列
    MPI_Request *reqy = &(req[16]);

    // ポインタ
    T **sendm = (T**)(bufInfo->m_bufSend[Y_MINUS]);
    T **recvm = (T**)(bufInfo->m_bufRecv[Y_MINUS]);
    T **sendp = (T**)(bufInfo->m_bufSend[Y_PLUS]);
    T **recvp = (T**)(bufInfo->m_bufRecv[Y_PLUS]);

    // 隣接面
    cpm_FaceFlag faceFlag[2] = {Y_MINUS, Y_PLUS};

    // 隣接ランクとの通信情報
    int nID[2][4];     //隣接ランク番号リスト
    int nFace[2];      //隣接ランク数
    int lDiff[2];      //レベル差
    size_t nw_send[2]; //送信サイズ
    size_t nw_recv[2]; //受信サイズ
    for( int i=0;i<2;i++ )
    {
      // 隣接ランク番号
      const int *nIDtmp = GetNeighborRankList(faceFlag[i], nFace[i], procGrpNo );
      for( int j=0;j<4;j++ )
      {
        nID[i][j] = nIDtmp[j];
      }

      // レベル差とバッファサイズ
      lDiff[i] = levelDiff[faceFlag[i]];
      if( lDiff[i] == 0 )
      {
        nw_send[i] = size_t(kmax  +2*vc_comm) * size_t(imax  +2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(kmax  +2*vc_comm) * size_t(imax  +2*vc_comm) *    vc_comm  * nmax;
      }
      else if( lDiff[i] == 1 )
      {
        nw_send[i] = size_t(kmax/2+2*vc_comm) * size_t(imax/2+2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(kmax  +4*vc_comm) * size_t(imax  +4*vc_comm) * (2*vc_comm) * nmax;
      }
      else if( lDiff[i] == -1 )
      {
        nw_send[i] = size_t(kmax  +4*vc_comm) * size_t(imax  +4*vc_comm) * (2*vc_comm) * nmax;
        nw_recv[i] = size_t(kmax/2+2*vc_comm) * size_t(imax/2+2*vc_comm) *    vc_comm  * nmax;
      }
      if( nw_send[i] > bufInfo->m_nsend[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
      if( nw_recv[i] > bufInfo->m_nrecv[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for recv error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
    }

    // irecv
    if( (ret = recv_LMR(nID, nFace, lDiff, nw_recv, recvm, &reqy[0], recvp, &reqy[4], procGrpNo)) != CPM_SUCCESS ) return ret;

    // +側通信面ループ(send)
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // pack
          if( (ret = packPY( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], sendp[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendp[j], nw_send[i], nID[i][j], &reqy[8+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // -側通信面ループ(send)
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // pack
          if( (ret = packMY( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], sendm[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendm[j], nw_send[i], nID[i][j], &reqy[12+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // wait(recv)
    if( (ret = Waitall( 8, &reqy[0] )) != CPM_SUCCESS ) return ret;

    // -側unpack
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // unpack
          if( (ret = unpackMY( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], recvm[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // +側unpack
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // unpack
          if( (ret = unpackPY( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], recvp[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // wait(send)
    if( (ret = Waitall( 8, &reqy[8] )) != CPM_SUCCESS ) return ret;

  } // Y方向完了


  //// Z方向の通信 ////
  {
    // リクエスト格納配列
    MPI_Request *reqz = &(req[32]);

    // ポインタ
    T **sendm = (T**)(bufInfo->m_bufSend[Z_MINUS]);
    T **recvm = (T**)(bufInfo->m_bufRecv[Z_MINUS]);
    T **sendp = (T**)(bufInfo->m_bufSend[Z_PLUS]);
    T **recvp = (T**)(bufInfo->m_bufRecv[Z_PLUS]);

    // 隣接面
    cpm_FaceFlag faceFlag[2] = {Z_MINUS, Z_PLUS};

    // 隣接ランクとの通信情報
    int nID[2][4];     //隣接ランク番号リスト
    int nFace[2];      //隣接ランク数
    int lDiff[2];      //レベル差
    size_t nw_send[2]; //送信サイズ
    size_t nw_recv[2]; //受信サイズ
    for( int i=0;i<2;i++ )
    {
      // 隣接ランク番号
      const int *nIDtmp = GetNeighborRankList(faceFlag[i], nFace[i], procGrpNo );
      for( int j=0;j<4;j++ )
      {
        nID[i][j] = nIDtmp[j];
      }

      // レベル差とバッファサイズ
      lDiff[i] = levelDiff[faceFlag[i]];
      if( lDiff[i] == 0 )
      {
        nw_send[i] = size_t(imax  +2*vc_comm) * size_t(jmax  +2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(imax  +2*vc_comm) * size_t(jmax  +2*vc_comm) *    vc_comm  * nmax;
      }
      else if( lDiff[i] == 1 )
      {
        nw_send[i] = size_t(imax/2+2*vc_comm) * size_t(jmax/2+2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(imax  +4*vc_comm) * size_t(jmax  +4*vc_comm) * (2*vc_comm) * nmax;
      }
      else if( lDiff[i] == -1 )
      {
        nw_send[i] = size_t(imax  +4*vc_comm) * size_t(jmax  +4*vc_comm) * (2*vc_comm) * nmax;
        nw_recv[i] = size_t(imax/2+2*vc_comm) * size_t(jmax/2+2*vc_comm) *    vc_comm  * nmax;
      }
      if( nw_send[i] > bufInfo->m_nsend[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
      if( nw_recv[i] > bufInfo->m_nrecv[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for recv error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
    }

    // irecv
    if( (ret = recv_LMR(nID, nFace, lDiff, nw_recv, recvm, &reqz[0], recvp, &reqz[4], procGrpNo)) != CPM_SUCCESS ) return ret;

    // +側通信面ループ(send)
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // pack
          if( (ret = packPZ( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], sendp[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendp[j], nw_send[i], nID[i][j], &reqz[8+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // -側通信面ループ(send)
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // pack
          if( (ret = packMZ( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], sendm[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendm[j], nw_send[i], nID[i][j], &reqz[12+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // wait(recv)
    if( (ret = Waitall( 8, &reqz[0] )) != CPM_SUCCESS ) return ret;

    // -側unpack
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // unpack
          if( (ret = unpackMZ( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], recvm[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // +側unpack
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // unpack
          if( (ret = unpackPZ( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], recvp[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // wait(send)
    if( (ret = Waitall( 8, &reqz[8] )) != CPM_SUCCESS ) return ret;

  } // Z方向完了


  // 正常終了
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4D版、waitなし)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommS4D_nowait( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                       , MPI_Request req[48], int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // リクエスト格納配列
  for( int i=0;i<48;i++ ) req[i] = MPI_REQUEST_NULL;


  // 隣接領域とのレベル差
  int levelDiff[6];
  const cpm_FaceFlag face[6] = {X_MINUS, X_PLUS, Y_MINUS, Y_PLUS, Z_MINUS, Z_PLUS};
  for( int i=0;i<6;i++ )
  {
    levelDiff[face[i]] = GetNeighborLevelDiff(face[i], procGrpNo);
  }

  // 通信バッファを取得
  S_BNDCOMM_BUFFER_LMR *bufInfo = GetBndCommBuffer(procGrpNo);
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }

  //// X方向の通信 ////
  {
    // リクエスト格納配列
    MPI_Request *reqx = &(req[0]);

    // ポインタ
    T **sendm = (T**)(bufInfo->m_bufSend[X_MINUS]);
    T **recvm = (T**)(bufInfo->m_bufRecv[X_MINUS]);
    T **sendp = (T**)(bufInfo->m_bufSend[X_PLUS]);
    T **recvp = (T**)(bufInfo->m_bufRecv[X_PLUS]);

    // 隣接面
    cpm_FaceFlag faceFlag[2] = {X_MINUS, X_PLUS};

    // 隣接ランクとの通信情報
    int nID[2][4];     //隣接ランク番号リスト
    int nFace[2];      //隣接ランク数
    int lDiff[2];      //レベル差
    size_t nw_send[2]; //送信サイズ
    size_t nw_recv[2]; //受信サイズ
    for( int i=0;i<2;i++ )
    {
      // 隣接ランク番号
      const int *nIDtmp = GetNeighborRankList(faceFlag[i], nFace[i], procGrpNo );
      for( int j=0;j<4;j++ )
      {
        nID[i][j] = nIDtmp[j];
      }

      // レベル差とバッファサイズ
      lDiff[i] = levelDiff[faceFlag[i]];
      if( lDiff[i] == 0 )
      {
        nw_send[i] = size_t(jmax  +2*vc_comm) * size_t(kmax  +2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(jmax  +2*vc_comm) * size_t(kmax  +2*vc_comm) *    vc_comm  * nmax;
      }
      else if( lDiff[i] == 1 )
      {
        nw_send[i] = size_t(jmax/2+2*vc_comm) * size_t(kmax/2+2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(jmax  +4*vc_comm) * size_t(kmax  +4*vc_comm) * (2*vc_comm) * nmax;
      }
      else if( lDiff[i] == -1 )
      {
        nw_send[i] = size_t(jmax  +4*vc_comm) * size_t(kmax  +4*vc_comm) * (2*vc_comm) * nmax;
        nw_recv[i] = size_t(jmax/2+2*vc_comm) * size_t(kmax/2+2*vc_comm) *    vc_comm  * nmax;
      }
      if( nw_send[i] > bufInfo->m_nsend[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
      if( nw_recv[i] > bufInfo->m_nrecv[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for recv error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
    }

    // irecv
    if( (ret = recv_LMR(nID, nFace, lDiff, nw_recv, recvm, &reqx[0], recvp, &reqx[4], procGrpNo)) != CPM_SUCCESS ) return ret;

    // +側通信面ループ(send)
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // pack
          if( (ret = packPX( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], sendp[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendp[j], nw_send[i], nID[i][j], &reqx[8+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // -側通信面ループ(send)
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // pack
          if( (ret = packMX( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], sendm[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendm[j], nw_send[i], nID[i][j], &reqx[12+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }
  } // X方向完了


  //// Y方向の通信 ////
  {
    // リクエスト格納配列
    MPI_Request *reqy = &(req[16]);

    // ポインタ
    T **sendm = (T**)(bufInfo->m_bufSend[Y_MINUS]);
    T **recvm = (T**)(bufInfo->m_bufRecv[Y_MINUS]);
    T **sendp = (T**)(bufInfo->m_bufSend[Y_PLUS]);
    T **recvp = (T**)(bufInfo->m_bufRecv[Y_PLUS]);

    // 隣接面
    cpm_FaceFlag faceFlag[2] = {Y_MINUS, Y_PLUS};

    // 隣接ランクとの通信情報
    int nID[2][4];     //隣接ランク番号リスト
    int nFace[2];      //隣接ランク数
    int lDiff[2];      //レベル差
    size_t nw_send[2]; //送信サイズ
    size_t nw_recv[2]; //受信サイズ
    for( int i=0;i<2;i++ )
    {
      // 隣接ランク番号
      const int *nIDtmp = GetNeighborRankList(faceFlag[i], nFace[i], procGrpNo );
      for( int j=0;j<4;j++ )
      {
        nID[i][j] = nIDtmp[j];
      }

      // レベル差とバッファサイズ
      lDiff[i] = levelDiff[faceFlag[i]];
      if( lDiff[i] == 0 )
      {
        nw_send[i] = size_t(kmax  +2*vc_comm) * size_t(imax  +2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(kmax  +2*vc_comm) * size_t(imax  +2*vc_comm) *    vc_comm  * nmax;
      }
      else if( lDiff[i] == 1 )
      {
        nw_send[i] = size_t(kmax/2+2*vc_comm) * size_t(imax/2+2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(kmax  +4*vc_comm) * size_t(imax  +4*vc_comm) * (2*vc_comm) * nmax;
      }
      else if( lDiff[i] == -1 )
      {
        nw_send[i] = size_t(kmax  +4*vc_comm) * size_t(imax  +4*vc_comm) * (2*vc_comm) * nmax;
        nw_recv[i] = size_t(kmax/2+2*vc_comm) * size_t(imax/2+2*vc_comm) *    vc_comm  * nmax;
      }
      if( nw_send[i] > bufInfo->m_nsend[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
      if( nw_recv[i] > bufInfo->m_nrecv[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for recv error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
    }

    // irecv
    if( (ret = recv_LMR(nID, nFace, lDiff, nw_recv, recvm, &reqy[0], recvp, &reqy[4], procGrpNo)) != CPM_SUCCESS ) return ret;

    // +側通信面ループ(send)
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // pack
          if( (ret = packPY( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], sendp[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendp[j], nw_send[i], nID[i][j], &reqy[8+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // -側通信面ループ(send)
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // pack
          if( (ret = packMY( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], sendm[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendm[j], nw_send[i], nID[i][j], &reqy[12+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }
  } // Y方向完了


  //// Z方向の通信 ////
  {
    // リクエスト格納配列
    MPI_Request *reqz = &(req[32]);

    // ポインタ
    T **sendm = (T**)(bufInfo->m_bufSend[Z_MINUS]);
    T **recvm = (T**)(bufInfo->m_bufRecv[Z_MINUS]);
    T **sendp = (T**)(bufInfo->m_bufSend[Z_PLUS]);
    T **recvp = (T**)(bufInfo->m_bufRecv[Z_PLUS]);

    // 隣接面
    cpm_FaceFlag faceFlag[2] = {Z_MINUS, Z_PLUS};

    // 隣接ランクとの通信情報
    int nID[2][4];     //隣接ランク番号リスト
    int nFace[2];      //隣接ランク数
    int lDiff[2];      //レベル差
    size_t nw_send[2]; //送信サイズ
    size_t nw_recv[2]; //受信サイズ
    for( int i=0;i<2;i++ )
    {
      // 隣接ランク番号
      const int *nIDtmp = GetNeighborRankList(faceFlag[i], nFace[i], procGrpNo );
      for( int j=0;j<4;j++ )
      {
        nID[i][j] = nIDtmp[j];
      }

      // レベル差とバッファサイズ
      lDiff[i] = levelDiff[faceFlag[i]];
      if( lDiff[i] == 0 )
      {
        nw_send[i] = size_t(imax  +2*vc_comm) * size_t(jmax  +2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(imax  +2*vc_comm) * size_t(jmax  +2*vc_comm) *    vc_comm  * nmax;
      }
      else if( lDiff[i] == 1 )
      {
        nw_send[i] = size_t(imax/2+2*vc_comm) * size_t(jmax/2+2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(imax  +4*vc_comm) * size_t(jmax  +4*vc_comm) * (2*vc_comm) * nmax;
      }
      else if( lDiff[i] == -1 )
      {
        nw_send[i] = size_t(imax  +4*vc_comm) * size_t(jmax  +4*vc_comm) * (2*vc_comm) * nmax;
        nw_recv[i] = size_t(imax/2+2*vc_comm) * size_t(jmax/2+2*vc_comm) *    vc_comm  * nmax;
      }
      if( nw_send[i] > bufInfo->m_nsend[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
      if( nw_recv[i] > bufInfo->m_nrecv[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for recv error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
    }

    // irecv
    if( (ret = recv_LMR(nID, nFace, lDiff, nw_recv, recvm, &reqz[0], recvp, &reqz[4], procGrpNo)) != CPM_SUCCESS ) return ret;

    // +側通信面ループ(send)
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // pack
          if( (ret = packPZ( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], sendp[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendp[j], nw_send[i], nID[i][j], &reqz[8+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // -側通信面ループ(send)
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // pack
          if( (ret = packMZ( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], sendm[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendm[j], nw_send[i], nID[i][j], &reqz[12+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }
  } // Z方向完了


  // 正常終了
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信のwait、展開(Scalar4D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::wait_BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                    , MPI_Request req[48], int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 隣接領域とのレベル差
  int levelDiff[6];
  const cpm_FaceFlag face[6] = {X_MINUS, X_PLUS, Y_MINUS, Y_PLUS, Z_MINUS, Z_PLUS};
  for( int i=0;i<6;i++ )
  {
    levelDiff[face[i]] = GetNeighborLevelDiff(face[i], procGrpNo);
  }

  // 通信バッファを取得
  S_BNDCOMM_BUFFER_LMR *bufInfo = GetBndCommBuffer(procGrpNo);
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }

  //// X方向の通信 ////
  {
    // リクエスト格納配列
    MPI_Request *reqx = &(req[0]);

    // ポインタ
    T **recvm = (T**)(bufInfo->m_bufRecv[X_MINUS]);
    T **recvp = (T**)(bufInfo->m_bufRecv[X_PLUS]);

    // 隣接面
    cpm_FaceFlag faceFlag[2] = {X_MINUS, X_PLUS};

    // 隣接ランクとの通信情報
    int nID[2][4];     //隣接ランク番号リスト
    int nFace[2];      //隣接ランク数
    int lDiff[2];      //レベル差
    for( int i=0;i<2;i++ )
    {
      // 隣接ランク番号
      const int *nIDtmp = GetNeighborRankList(faceFlag[i], nFace[i], procGrpNo );
      for( int j=0;j<4;j++ )
      {
        nID[i][j] = nIDtmp[j];
      }

      // レベル差
      lDiff[i] = levelDiff[faceFlag[i]];
    }

    // wait(recv)
    if( (ret = Waitall( 8, &reqx[0] )) != CPM_SUCCESS ) return ret;

    // -側unpack
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // unpack
          if( (ret = unpackMX( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], recvm[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // +側unpack
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // unpack
          if( (ret = unpackPX( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], recvp[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }
  } // X方向完了


  //// Y方向の通信 ////
  {
    // リクエスト格納配列
    MPI_Request *reqy = &(req[16]);

    // ポインタ
    T **recvm = (T**)(bufInfo->m_bufRecv[Y_MINUS]);
    T **recvp = (T**)(bufInfo->m_bufRecv[Y_PLUS]);

    // 隣接面
    cpm_FaceFlag faceFlag[2] = {Y_MINUS, Y_PLUS};

    // 隣接ランクとの通信情報
    int nID[2][4];     //隣接ランク番号リスト
    int nFace[2];      //隣接ランク数
    int lDiff[2];      //レベル差
    for( int i=0;i<2;i++ )
    {
      // 隣接ランク番号
      const int *nIDtmp = GetNeighborRankList(faceFlag[i], nFace[i], procGrpNo );
      for( int j=0;j<4;j++ )
      {
        nID[i][j] = nIDtmp[j];
      }

      // レベル差
      lDiff[i] = levelDiff[faceFlag[i]];
    }

    // wait(recv)
    if( (ret = Waitall( 8, &reqy[0] )) != CPM_SUCCESS ) return ret;

    // -側unpack
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // unpack
          if( (ret = unpackMY( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], recvm[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // +側unpack
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // unpack
          if( (ret = unpackPY( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], recvp[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }
  } // Y方向完了


  //// Z方向の通信 ////
  {
    // リクエスト格納配列
    MPI_Request *reqz = &(req[32]);

    // ポインタ
    T **recvm = (T**)(bufInfo->m_bufRecv[Z_MINUS]);
    T **recvp = (T**)(bufInfo->m_bufRecv[Z_PLUS]);

    // 隣接面
    cpm_FaceFlag faceFlag[2] = {Z_MINUS, Z_PLUS};

    // 隣接ランクとの通信情報
    int nID[2][4];     //隣接ランク番号リスト
    int nFace[2];      //隣接ランク数
    int lDiff[2];      //レベル差
    for( int i=0;i<2;i++ )
    {
      // 隣接ランク番号
      const int *nIDtmp = GetNeighborRankList(faceFlag[i], nFace[i], procGrpNo );
      for( int j=0;j<4;j++ )
      {
        nID[i][j] = nIDtmp[j];
      }

      // レベル差
      lDiff[i] = levelDiff[faceFlag[i]];
    }

    // wait(recv)
    if( (ret = Waitall( 8, &reqz[0] )) != CPM_SUCCESS ) return ret;

    // -側unpack
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // unpack
          if( (ret = unpackMZ( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], recvm[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // +側unpack
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nID[i][j]) )
        {
          // unpack
          if( (ret = unpackPZ( array, imax, jmax, kmax, nmax, vc, vc_comm, nID[i][j], j, lDiff[i], recvp[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }
  } // Z方向完了

  // wait(send)
  if( (ret = Waitall( 8, &req[8]  )) != CPM_SUCCESS ) return ret;
  if( (ret = Waitall( 8, &req[24] )) != CPM_SUCCESS ) return ret;
  if( (ret = Waitall( 8, &req[40] )) != CPM_SUCCESS ) return ret;

  // 正常終了
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar4D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::PeriodicCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                    , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // リクエスト格納配列
  MPI_Request req[48];
  for( int i=0;i<48;i++ ) req[i] = MPI_REQUEST_NULL;


  // 隣接領域とのレベル差
  int levelDiff[6];
  const cpm_FaceFlag face[6] = {X_MINUS, X_PLUS, Y_MINUS, Y_PLUS, Z_MINUS, Z_PLUS};
  for( int i=0;i<6;i++ )
  {
    levelDiff[face[i]] = GetNeighborLevelDiff(face[i], procGrpNo);
  }

  // 通信バッファを取得
  S_BNDCOMM_BUFFER_LMR *bufInfo = GetBndCommBuffer(procGrpNo);
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }

  //// X方向の通信 ////
  if( dir==X_DIR )
  {
    // リクエスト格納配列
    MPI_Request *reqx = &(req[0]);

    // ポインタ
    T **sendm = (T**)(bufInfo->m_bufSend[X_MINUS]);
    T **recvm = (T**)(bufInfo->m_bufRecv[X_MINUS]);
    T **sendp = (T**)(bufInfo->m_bufSend[X_PLUS]);
    T **recvp = (T**)(bufInfo->m_bufRecv[X_PLUS]);

    // 隣接面
    cpm_FaceFlag faceFlag[2] = {X_MINUS, X_PLUS};

    // 隣接ランクとの通信情報
    int nFace[2];      //隣接ランク数
    int lDiff[2];      //レベル差
    size_t nw_send[2]; //送信サイズ
    size_t nw_recv[2]; //受信サイズ
    int nIDsend[2][4]; //送信先ランクリスト
    int nIDrecv[2][4]; //受信先ランクリスト
    for( int i=0;i<2;i++ )
    {
      // 隣接ランク番号
      for( int j=0;j<4;j++ )
      {
        nIDsend[i][j] = cpm_ParaManager::getRankNull();
        nIDrecv[i][j] = cpm_ParaManager::getRankNull();
      }
      const int *nIDtmp = GetPeriodicRankList(faceFlag[i], nFace[i], procGrpNo );
      if( pm==MINUS2PLUS || pm==BOTH )
      {
        if( i==0 )
        {
          for( int j=0;j<4;j++ )
          {
            nIDsend[i][j] = nIDtmp[j];
          }
        }
        if( i==1 )
        {
          for( int j=0;j<4;j++ )
          {
            nIDrecv[i][j] = nIDtmp[j];
          }
        }
      }
      if( pm==PLUS2MINUS || pm==BOTH )
      {
        if( i==0 )
        {
          for( int j=0;j<4;j++ )
          {
            nIDrecv[i][j] = nIDtmp[j];
          }
        }
        if( i==1 )
        {
          for( int j=0;j<4;j++ )
          {
            nIDsend[i][j] = nIDtmp[j];
          }
        }
      }

      // レベル差とバッファサイズ 
      lDiff[i] = levelDiff[faceFlag[i]];
      if( lDiff[i] == 0 )
      {
        nw_send[i] = size_t(jmax  +2*vc_comm) * size_t(kmax  +2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(jmax  +2*vc_comm) * size_t(kmax  +2*vc_comm) *    vc_comm  * nmax;
      }
      else if( lDiff[i] == 1 )
      {
        nw_send[i] = size_t(jmax/2+2*vc_comm) * size_t(kmax/2+2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(jmax  +4*vc_comm) * size_t(kmax  +4*vc_comm) * (2*vc_comm) * nmax;
      }
      else if( lDiff[i] == -1 )
      {
        nw_send[i] = size_t(jmax  +4*vc_comm) * size_t(kmax  +4*vc_comm) * (2*vc_comm) * nmax;
        nw_recv[i] = size_t(jmax/2+2*vc_comm) * size_t(kmax/2+2*vc_comm) *    vc_comm  * nmax;
      }
      if( nw_send[i] > bufInfo->m_nsend[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
      if( nw_recv[i] > bufInfo->m_nrecv[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for recv error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
    }

    // irecv
    if( (ret = recv_LMR(nIDrecv, nFace, lDiff, nw_recv, recvm, &reqx[0], recvp, &reqx[4], procGrpNo)) != CPM_SUCCESS ) return ret;

    // +側通信面ループ(send)
    if( pm==PLUS2MINUS || pm==BOTH )
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nIDsend[i][j]) )
        {
          // pack
          if( (ret = packPX( array, imax, jmax, kmax, nmax, vc, vc_comm, nIDsend[i][j], j, lDiff[i], sendp[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendp[j], nw_send[i], nIDsend[i][j], &reqx[8+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // -側通信面ループ(send)
    if( pm==MINUS2PLUS || pm==BOTH )
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nIDsend[i][j]) )
        {
          // pack
          if( (ret = packMX( array, imax, jmax, kmax, nmax, vc, vc_comm, nIDsend[i][j], j, lDiff[i], sendm[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendm[j], nw_send[i], nIDsend[i][j], &reqx[12+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // wait(recv)
    if( (ret = Waitall( 8, &reqx[0] )) != CPM_SUCCESS ) return ret;

    // -側unpack
    if( pm==PLUS2MINUS || pm==BOTH )
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nIDrecv[i][j]) )
        {
          // unpack
          if( (ret = unpackMX( array, imax, jmax, kmax, nmax, vc, vc_comm, nIDrecv[i][j], j, lDiff[i], recvm[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // +側unpack
    if( pm==MINUS2PLUS || pm==BOTH)
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nIDrecv[i][j]) )
        {
          // unpack
          if( (ret = unpackPX( array, imax, jmax, kmax, nmax, vc, vc_comm, nIDrecv[i][j], j, lDiff[i], recvp[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // wait(send)
    if( (ret = Waitall( 8, &reqx[8] )) != CPM_SUCCESS ) return ret;

  } // X方向完了


  //// Y方向の通信 ////
  if( dir==Y_DIR )
  {
    // リクエスト格納配列
    MPI_Request *reqy = &(req[16]);

    // ポインタ
    T **sendm = (T**)(bufInfo->m_bufSend[Y_MINUS]);
    T **recvm = (T**)(bufInfo->m_bufRecv[Y_MINUS]);
    T **sendp = (T**)(bufInfo->m_bufSend[Y_PLUS]);
    T **recvp = (T**)(bufInfo->m_bufRecv[Y_PLUS]);

    // 隣接面
    cpm_FaceFlag faceFlag[2] = {Y_MINUS, Y_PLUS};

    // 隣接ランクとの通信情報
    int nFace[2];      //隣接ランク数
    int lDiff[2];      //レベル差
    size_t nw_send[2]; //送信サイズ
    size_t nw_recv[2]; //受信サイズ
    int nIDsend[2][4]; //送信先ランクリスト
    int nIDrecv[2][4]; //受信先ランクリスト
    for( int i=0;i<2;i++ )
    {
      // 隣接ランク番号
      for( int j=0;j<4;j++ )
      {
        nIDsend[i][j] = cpm_ParaManager::getRankNull();
        nIDrecv[i][j] = cpm_ParaManager::getRankNull();
      }
      const int *nIDtmp = GetPeriodicRankList(faceFlag[i], nFace[i], procGrpNo );
      if( pm==MINUS2PLUS || pm==BOTH )
      {
        if( i==0 )
        {
          for( int j=0;j<4;j++ )
          {
            nIDsend[i][j] = nIDtmp[j];
          }
        }
        if( i==1 )
        {
          for( int j=0;j<4;j++ )
          {
            nIDrecv[i][j] = nIDtmp[j];
          }
        }
      }
      if( pm==PLUS2MINUS || pm==BOTH )
      {
        if( i==0 )
        {
          for( int j=0;j<4;j++ )
          {
            nIDrecv[i][j] = nIDtmp[j];
          }
        }
        if( i==1 )
        {
          for( int j=0;j<4;j++ )
          {
            nIDsend[i][j] = nIDtmp[j];
          }
        }
      }

      // レベル差とバッファサイズ 
      lDiff[i] = levelDiff[faceFlag[i]];
      if( lDiff[i] == 0 )
      {
        nw_send[i] = size_t(kmax  +2*vc_comm) * size_t(imax  +2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(kmax  +2*vc_comm) * size_t(imax  +2*vc_comm) *    vc_comm  * nmax;
      }
      else if( lDiff[i] == 1 )
      {
        nw_send[i] = size_t(kmax/2+2*vc_comm) * size_t(imax/2+2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(kmax  +4*vc_comm) * size_t(imax  +4*vc_comm) * (2*vc_comm) * nmax;
      }
      else if( lDiff[i] == -1 )
      {
        nw_send[i] = size_t(kmax  +4*vc_comm) * size_t(imax  +4*vc_comm) * (2*vc_comm) * nmax;
        nw_recv[i] = size_t(kmax/2+2*vc_comm) * size_t(imax/2+2*vc_comm) *    vc_comm  * nmax;
      }
      if( nw_send[i] > bufInfo->m_nsend[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
      if( nw_recv[i] > bufInfo->m_nrecv[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for recv error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
    }

    // irecv
    if( (ret = recv_LMR(nIDrecv, nFace, lDiff, nw_recv, recvm, &reqy[0], recvp, &reqy[4], procGrpNo)) != CPM_SUCCESS ) return ret;

    // +側通信面ループ(send)
    if( pm==PLUS2MINUS || pm==BOTH )
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nIDsend[i][j]) )
        {
          // pack
          if( (ret = packPY( array, imax, jmax, kmax, nmax, vc, vc_comm, nIDsend[i][j], j, lDiff[i], sendp[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendp[j], nw_send[i], nIDsend[i][j], &reqy[8+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // -側通信面ループ(send)
    if( pm==MINUS2PLUS || pm==BOTH )
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nIDsend[i][j]) )
        {
          // pack
          if( (ret = packMY( array, imax, jmax, kmax, nmax, vc, vc_comm, nIDsend[i][j], j, lDiff[i], sendm[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendm[j], nw_send[i], nIDsend[i][j], &reqy[12+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // wait(recv)
    if( (ret = Waitall( 8, &reqy[0] )) != CPM_SUCCESS ) return ret;

    // -側unpack
    if( pm==PLUS2MINUS || pm==BOTH )
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nIDrecv[i][j]) )
        {
          // unpack
          if( (ret = unpackMY( array, imax, jmax, kmax, nmax, vc, vc_comm, nIDrecv[i][j], j, lDiff[i], recvm[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // +側unpack
    if( pm==MINUS2PLUS || pm==BOTH)
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nIDrecv[i][j]) )
        {
          // unpack
          if( (ret = unpackPY( array, imax, jmax, kmax, nmax, vc, vc_comm, nIDrecv[i][j], j, lDiff[i], recvp[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // wait(send)
    if( (ret = Waitall( 8, &reqy[8] )) != CPM_SUCCESS ) return ret;

  } // Y方向完了


  //// Z方向の通信 ////
  if( dir==Z_DIR )
  {
    // リクエスト格納配列
    MPI_Request *reqz = &(req[32]);

    // ポインタ
    T **sendm = (T**)(bufInfo->m_bufSend[Z_MINUS]);
    T **recvm = (T**)(bufInfo->m_bufRecv[Z_MINUS]);
    T **sendp = (T**)(bufInfo->m_bufSend[Z_PLUS]);
    T **recvp = (T**)(bufInfo->m_bufRecv[Z_PLUS]);

    // 隣接面
    cpm_FaceFlag faceFlag[2] = {Z_MINUS, Z_PLUS};

    // 隣接ランクとの通信情報
    int nFace[2];      //隣接ランク数
    int lDiff[2];      //レベル差
    size_t nw_send[2]; //送信サイズ
    size_t nw_recv[2]; //受信サイズ
    int nIDsend[2][4]; //送信先ランクリスト
    int nIDrecv[2][4]; //受信先ランクリスト
    for( int i=0;i<2;i++ )
    {
      // 隣接ランク番号
      for( int j=0;j<4;j++ )
      {
        nIDsend[i][j] = cpm_ParaManager::getRankNull();
        nIDrecv[i][j] = cpm_ParaManager::getRankNull();
      }
      const int *nIDtmp = GetPeriodicRankList(faceFlag[i], nFace[i], procGrpNo );
      if( pm==MINUS2PLUS || pm==BOTH )
      {
        if( i==0 )
        {
          for( int j=0;j<4;j++ )
          {
            nIDsend[i][j] = nIDtmp[j];
          }
        }
        if( i==1 )
        {
          for( int j=0;j<4;j++ )
          {
            nIDrecv[i][j] = nIDtmp[j];
          }
        }
      }
      if( pm==PLUS2MINUS || pm==BOTH )
      {
        if( i==0 )
        {
          for( int j=0;j<4;j++ )
          {
            nIDrecv[i][j] = nIDtmp[j];
          }
        }
        if( i==1 )
        {
          for( int j=0;j<4;j++ )
          {
            nIDsend[i][j] = nIDtmp[j];
          }
        }
      }

      // レベル差とバッファサイズ 
      lDiff[i] = levelDiff[faceFlag[i]];
      if( lDiff[i] == 0 )
      {
        nw_send[i] = size_t(imax  +2*vc_comm) * size_t(jmax  +2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(imax  +2*vc_comm) * size_t(jmax  +2*vc_comm) *    vc_comm  * nmax;
      }
      else if( lDiff[i] == 1 )
      {
        nw_send[i] = size_t(imax/2+2*vc_comm) * size_t(jmax/2+2*vc_comm) *    vc_comm  * nmax;
        nw_recv[i] = size_t(imax  +4*vc_comm) * size_t(jmax  +4*vc_comm) * (2*vc_comm) * nmax;
      }
      else if( lDiff[i] == -1 )
      {
        nw_send[i] = size_t(imax  +4*vc_comm) * size_t(jmax  +4*vc_comm) * (2*vc_comm) * nmax;
        nw_recv[i] = size_t(imax/2+2*vc_comm) * size_t(jmax/2+2*vc_comm) *    vc_comm  * nmax;
      }
      if( nw_send[i] > bufInfo->m_nsend[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
      if( nw_recv[i] > bufInfo->m_nrecv[faceFlag[i]] )
      {
        stmpd_printf("[%d] BndComm buffer length for recv error\n", m_rankNo);
        return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
      }
    }

    // irecv
    if( (ret = recv_LMR(nIDrecv, nFace, lDiff, nw_recv, recvm, &reqz[0], recvp, &reqz[4], procGrpNo)) != CPM_SUCCESS ) return ret;

    // +側通信面ループ(send)
    if( pm==PLUS2MINUS || pm==BOTH )
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nIDsend[i][j]) )
        {
          // pack
          if( (ret = packPZ( array, imax, jmax, kmax, nmax, vc, vc_comm, nIDsend[i][j], j, lDiff[i], sendp[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendp[j], nw_send[i], nIDsend[i][j], &reqz[8+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // -側通信面ループ(send)
    if( pm==MINUS2PLUS || pm==BOTH )
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nIDsend[i][j]) )
        {
          // pack
          if( (ret = packMZ( array, imax, jmax, kmax, nmax, vc, vc_comm, nIDsend[i][j], j, lDiff[i], sendm[j], nw_send[i] )) != CPM_SUCCESS ) return ret;

          // isend
          if( (ret = Isend( sendm[j], nw_send[i], nIDsend[i][j], &reqz[12+j], procGrpNo )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // wait(recv)
    if( (ret = Waitall( 8, &reqz[0] )) != CPM_SUCCESS ) return ret;

    // -側unpack
    if( pm==PLUS2MINUS || pm==BOTH )
    {
      int i=0;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nIDrecv[i][j]) )
        {
          // unpack
          if( (ret = unpackMZ( array, imax, jmax, kmax, nmax, vc, vc_comm, nIDrecv[i][j], j, lDiff[i], recvm[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // +側unpack
    if( pm==MINUS2PLUS || pm==BOTH)
    {
      int i=1;
      for( int j=0;j<nFace[i];j++ )
      {
        if( !IsRankNull(nIDrecv[i][j]) )
        {
          // unpack
          if( (ret = unpackPZ( array, imax, jmax, kmax, nmax, vc, vc_comm, nIDrecv[i][j], j, lDiff[i], recvp[j] )) != CPM_SUCCESS ) return ret;
        }
      }
    }

    // wait(send)
    if( (ret = Waitall( 8, &reqz[8] )) != CPM_SUCCESS ) return ret;

  } // Z方向完了


  // 正常終了
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// １方向(プラス、マイナス)の非同期受信処理
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::recv_LMR( const int nID[2][4], int nFace[2], int levelDiff[2], size_t nw[2]
                            , T* recvm[4], MPI_Request *reqm, T* recvp[4], MPI_Request *reqp
                            , int procGrpNo )
{
  cpm_ErrorCode ret;

  MPI_Request *req[2] = {reqm, reqp};
  T** recvbuf[2] = {recvm, recvp};

  // +/-面ループ
  for( int i=0;i<2;i++ )
  {
    // 通信面ループ
    for( int j=0;j<nFace[i];j++ )
    {
      // ランク番号が有効
      if( !IsRankNull(nID[i][j]) )
      {
        if( (ret = Irecv( recvbuf[i][j], nw[i], nID[i][j], &req[i][j], procGrpNo )) != CPM_SUCCESS )
        {
          return ret;
        }
      }
    }
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)の-X面への送信データのパック(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::packMX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                          , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw )
{
  if( IsRankNull(nID) ) return CPM_SUCCESS;

  // 送信範囲の確定
  int js=0, je=jmax; //j方向範囲
  int ks=0, ke=kmax; //k方向範囲
  int gc=vc_comm;    //通信層数
  if( levelDiff==0 )
  {
    // 同じレベル
    // そのまま
  }
  else if( levelDiff==1 )
  {
    // coarse -> fine
    // 1/4面を送信
    if( faceNo==0 )
    {
      je = jmax/2;
      ke = kmax/2;
    }
    else if( faceNo==1 )
    {
      js = jmax/2;
      ke = kmax/2;
    }
    else if( faceNo==2 )
    {
      je = jmax/2;
      ks = kmax/2;
    }
    else if( faceNo==3 )
    {
      js = jmax/2;
      ks = kmax/2;
    }
  }
  else if( levelDiff==-1 )
  {
    // fine -> coarse
    // 層数を2倍
    gc = vc_comm*2;
  }

  // 送信バッファのサイズ
  int jmaxb = je-js;
  int kmaxb = ke-ks;
  size_t sz = size_t(jmaxb+2*gc) * size_t(kmaxb+2*gc) * size_t(gc) * size_t(nmax);
  if( sz > nw )
  {
    stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // 格納
  for( int n=0;n<nmax;n++){
  for( int k=ks-gc,kk=0-gc;k<ke+gc;k++,kk++ ){
  for( int j=js-gc,jj=0-gc;j<je+gc;j++,jj++ ){
  for( int i=0;i<gc;i++ ){
    sendbuf[_IDXFX(i,jj,kk,n,0,jmaxb,kmaxb,gc)] = array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)];
  }}}}

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)の+X面への送信データのパック(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::packPX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                          , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw )
{
  if( IsRankNull(nID) ) return CPM_SUCCESS;

  // 送信範囲の確定
  int js=0, je=jmax; //j方向範囲
  int ks=0, ke=kmax; //k方向範囲
  int gc=vc_comm;    //通信層数
  if( levelDiff==0 )
  {
    // 同じレベル
    // そのまま
  }
  else if( levelDiff==1 )
  {
    // coarse -> fine
    // 1/4面を送信
    if( faceNo==0 )
    {
      je = jmax/2;
      ke = kmax/2;
    }
    else if( faceNo==1 )
    {
      js = jmax/2;
      ke = kmax/2;
    }
    else if( faceNo==2 )
    {
      je = jmax/2;
      ks = kmax/2;
    }
    else if( faceNo==3 )
    {
      js = jmax/2;
      ks = kmax/2;
    }
  }
  else if( levelDiff==-1 )
  {
    // fine -> coarse
    // 層数を2倍
    gc = vc_comm*2;
  }

  // 送信バッファのサイズ
  int jmaxb = je-js;
  int kmaxb = ke-ks;
  size_t sz = size_t(jmaxb+2*gc) * size_t(kmaxb+2*gc) * size_t(gc) * size_t(nmax);
  if( sz > nw )
  {
    stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // 格納
  for( int n=0;n<nmax;n++){
  for( int k=ks-gc,kk=0-gc;k<ke+gc;k++,kk++ ){
  for( int j=js-gc,jj=0-gc;j<je+gc;j++,jj++ ){
  for( int i=imax-gc;i<imax;i++ ){
    sendbuf[_IDXFX(i,jj,kk,n,imax-gc,jmaxb,kmaxb,gc)] = array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)];
  }}}}

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)の-X面からの受信データの展開(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::unpackMX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                            , const int nID, int faceNo, int levelDiff, T* recvbuf )
{
  if( IsRankNull(nID) ) return CPM_SUCCESS;

  if( levelDiff==0 )
  {
    // 同じレベル
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<0;i++ ){
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = recvbuf[_IDXFX(i,j,k,n,0-vc_comm,jmax,kmax,vc_comm)];
    }}}}
  }
  else if( levelDiff==1 )
  {
    // fine  -> coarse
    // 8cell -> 1cell
    // 1/4面を受信
    int js=0, je=jmax; //j方向範囲
    int ks=0, ke=kmax; //k方向範囲
    int vc_jm=vc_comm, vc_jp=vc_comm; //j方向仮想セル数
    int vc_km=vc_comm, vc_kp=vc_comm; //k方向仮想セル数
    if( faceNo==0 )
    {
      je = jmax/2;
      ke = kmax/2;
      vc_jp = 0;
      vc_kp = 0;
    }
    else if( faceNo==1 )
    {
      js = jmax/2;
      ke = kmax/2;
      vc_jm = 0;
      vc_kp = 0;
    }
    else if( faceNo==2 )
    {
      je = jmax/2;
      ks = kmax/2;
      vc_jp = 0;
      vc_km = 0;
    }
    else if( faceNo==3 )
    {
      js = jmax/2;
      ks = kmax/2;
      vc_jm = 0;
      vc_km = 0;
    }

    // 送信バッファのサイズ
    int gc    = vc_comm*2;

    for( int n=0;n<nmax;n++){
    for( int k=ks-vc_km;k<ke+vc_kp;k++ ){
    for( int j=js-vc_jm;j<je+vc_jp;j++ ){
    for( int i=0-vc_comm;i<0;i++ ){
      int ii=(i   )*2;
      int jj=(j-js)*2;
      int kk=(k-ks)*2;
      T val = recvbuf[_IDXFX(ii  ,jj  ,kk  ,n,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii+1,jj  ,kk  ,n,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii  ,jj+1,kk  ,n,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii+1,jj+1,kk  ,n,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii  ,jj  ,kk+1,n,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii+1,jj  ,kk+1,n,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii  ,jj+1,kk+1,n,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii+1,jj+1,kk+1,n,0-gc,jmax,kmax,gc)];
      val *= 0.125;
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = val;
    }}}}
  }
  else if( levelDiff==-1 )
  {
    // coarse -> fine
    // 1cell  -> 8cell

    // 送信バッファのサイズ
    int jmaxb = jmax/2;
    int kmaxb = kmax/2;

    // 展開する袖数は2倍
    int gc    = vc_comm*2;

    for( int n=0;n<nmax;n++){
    for( int k=0-gc;k<kmax+gc;k++ ){
    for( int j=0-gc;j<jmax+gc;j++ ){
    for( int i=0-gc;i<0;i++ ){
      int ii=(i+2)/2-1;
      int jj=(j+2)/2-1;
      int kk=(k+2)/2-1;
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = recvbuf[_IDXFX(ii,jj,kk,n,0-vc_comm,jmaxb,kmaxb,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)の+X面からの受信データの展開(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::unpackPX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                            , const int nID, int faceNo, int levelDiff, T* recvbuf )
{
  if( IsRankNull(nID) ) return CPM_SUCCESS;

  if( levelDiff==0 )
  {
    // 同じレベル
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=imax;i<imax+vc_comm;i++ ){
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = recvbuf[_IDXFX(i,j,k,n,imax,jmax,kmax,vc_comm)];
    }}}}
  }
  else if( levelDiff==1 )
  {
    // fine  -> coarse
    // 8cell -> 1cell
    // 1/4面を受信
    int js=0, je=jmax; //j方向範囲
    int ks=0, ke=kmax; //k方向範囲
    int vc_jm=vc_comm, vc_jp=vc_comm; //j方向仮想セル数
    int vc_km=vc_comm, vc_kp=vc_comm; //k方向仮想セル数
    if( faceNo==0 )
    {
      je = jmax/2;
      ke = kmax/2;
      vc_jp = 0;
      vc_kp = 0;
    }
    else if( faceNo==1 )
    {
      js = jmax/2;
      ke = kmax/2;
      vc_jm = 0;
      vc_kp = 0;
    }
    else if( faceNo==2 )
    {
      je = jmax/2;
      ks = kmax/2;
      vc_jp = 0;
      vc_km = 0;
    }
    else if( faceNo==3 )
    {
      js = jmax/2;
      ks = kmax/2;
      vc_jm = 0;
      vc_km = 0;
    }

    // 送信バッファのサイズ
    int gc    = vc_comm*2;
    int imax2 = imax*2;

    for( int n=0;n<nmax;n++){
    for( int k=ks-vc_km;k<ke+vc_kp;k++ ){
    for( int j=js-vc_jm;j<je+vc_jp;j++ ){
    for( int i=imax;i<imax+vc_comm;i++ ){
      int ii=(i   )*2;
      int jj=(j-js)*2;
      int kk=(k-ks)*2;
      T val = recvbuf[_IDXFX(ii  ,jj  ,kk  ,n,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii+1,jj  ,kk  ,n,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii  ,jj+1,kk  ,n,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii+1,jj+1,kk  ,n,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii  ,jj  ,kk+1,n,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii+1,jj  ,kk+1,n,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii  ,jj+1,kk+1,n,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(ii+1,jj+1,kk+1,n,imax2,jmax,kmax,gc)];
      val *= 0.125;
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = val;
    }}}}
  }
  else if( levelDiff==-1 )
  {
    // coarse -> fine
    // 1cell  -> 8cell

    // 送信バッファのサイズ
    int jmaxb = jmax/2;
    int kmaxb = kmax/2;
    int imax2 = imax/2;

    // 展開する袖数は2倍
    int gc    = vc_comm*2;

    for( int n=0;n<nmax;n++){
    for( int k=0-gc;k<kmax+gc;k++ ){
    for( int j=0-gc;j<jmax+gc;j++ ){
    for( int i=imax;i<imax+gc;i++ ){
      int ii=i/2;
      int jj=(j+2)/2-1;
      int kk=(k+2)/2-1;
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = recvbuf[_IDXFX(ii,jj,kk,n,imax2,jmaxb,kmaxb,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)の-Y面への送信データのパック(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::packMY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                          , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw )
{
  if( IsRankNull(nID) ) return CPM_SUCCESS;

  // 送信範囲の確定
  int is=0, ie=imax; //i方向範囲
  int ks=0, ke=kmax; //k方向範囲
  int gc=vc_comm;    //通信層数
  if( levelDiff==0 )
  {
    // 同じレベル
    // そのまま
  }
  else if( levelDiff==1 )
  {
    // coarse -> fine
    // 1/4面を送信
    if( faceNo==0 )
    {
      ie = imax/2;
      ke = kmax/2;
    }
    else if( faceNo==1 )
    {
      ie = imax/2;
      ks = kmax/2;
    }
    else if( faceNo==2 )
    {
      is = imax/2;
      ke = kmax/2;
    }
    else if( faceNo==3 )
    {
      is = imax/2;
      ks = kmax/2;
    }
  }
  else if( levelDiff==-1 )
  {
    // fine -> coarse
    // 層数を2倍
    gc = vc_comm*2;
  }

  // 送信バッファのサイズ
  int imaxb = ie-is;
  int kmaxb = ke-ks;
  size_t sz = size_t(imaxb+2*gc) * size_t(kmaxb+2*gc) * size_t(gc) * size_t(nmax);
  if( sz > nw )
  {
    stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // 格納
  for( int n=0;n<nmax;n++){
  for( int k=ks-gc,kk=0-gc;k<ke+gc;k++,kk++ ){
  for( int j=0;j<gc;j++ ){
  for( int i=is-gc,ii=0-gc;i<ie+gc;i++,ii++ ){
    sendbuf[_IDXFY(ii,j,kk,n,imaxb,0,kmaxb,gc)] = array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)];
  }}}}

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)の+Y面への送信データのパック(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::packPY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                          , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw )
{
  if( IsRankNull(nID) ) return CPM_SUCCESS;

  // 送信範囲の確定
  int is=0, ie=imax; //i方向範囲
  int ks=0, ke=kmax; //k方向範囲
  int gc=vc_comm;    //通信層数
  if( levelDiff==0 )
  {
    // 同じレベル
    // そのまま
  }
  else if( levelDiff==1 )
  {
    // coarse -> fine
    // 1/4面を送信
    if( faceNo==0 )
    {
      ie = imax/2;
      ke = kmax/2;
    }
    else if( faceNo==1 )
    {
      ie = imax/2;
      ks = kmax/2;
    }
    else if( faceNo==2 )
    {
      is = imax/2;
      ke = kmax/2;
    }
    else if( faceNo==3 )
    {
      is = imax/2;
      ks = kmax/2;
    }
  }
  else if( levelDiff==-1 )
  {
    // fine -> coarse
    // 層数を2倍
    gc = vc_comm*2;
  }

  // 送信バッファのサイズ
  int imaxb = ie-is;
  int kmaxb = ke-ks;
  size_t sz = size_t(imaxb+2*gc) * size_t(kmaxb+2*gc) * size_t(gc) * size_t(nmax);
  if( sz > nw )
  {
    stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // 格納
  for( int n=0;n<nmax;n++){
  for( int k=ks-gc,kk=0-gc;k<ke+gc;k++,kk++ ){
  for( int j=jmax-gc;j<jmax;j++ ){
  for( int i=is-gc,ii=0-gc;i<ie+gc;i++,ii++ ){
    sendbuf[_IDXFY(ii,j,kk,n,imaxb,jmax-gc,kmaxb,gc)] = array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)];
  }}}}

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)の-Y面からの受信データの展開(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::unpackMY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                            , const int nID, int faceNo, int levelDiff, T* recvbuf )
{
  if( IsRankNull(nID) ) return CPM_SUCCESS;

  if( levelDiff==0 )
  {
    // 同じレベル
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<0;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = recvbuf[_IDXFY(i,j,k,n,imax,0-vc_comm,kmax,vc_comm)];
    }}}}
  }
  else if( levelDiff==1 )
  {
    // fine  -> coarse
    // 8cell -> 1cell
    // 1/4面を受信
    int is=0, ie=imax; //i方向範囲
    int ks=0, ke=kmax; //k方向範囲
    int vc_im=vc_comm, vc_ip=vc_comm; //i方向仮想セル数
    int vc_km=vc_comm, vc_kp=vc_comm; //k方向仮想セル数
    if( faceNo==0 )
    {
      ie = imax/2;
      ke = kmax/2;
      vc_ip = 0;
      vc_kp = 0;
    }
    else if( faceNo==1 )
    {
      ie = imax/2;
      ks = kmax/2;
      vc_ip = 0;
      vc_km = 0;
    }
    else if( faceNo==2 )
    {
      is = imax/2;
      ke = kmax/2;
      vc_im = 0;
      vc_kp = 0;
    }
    else if( faceNo==3 )
    {
      is = imax/2;
      ks = kmax/2;
      vc_im = 0;
      vc_km = 0;
    }

    // 送信バッファのサイズ
    int gc    = vc_comm*2;

    for( int n=0;n<nmax;n++){
    for( int k=ks-vc_km;k<ke+vc_kp;k++ ){
    for( int j=0-vc_comm;j<0;j++ ){
    for( int i=is-vc_im;i<ie+vc_ip;i++ ){
      int ii=(i-is)*2;
      int jj=(j   )*2;
      int kk=(k-ks)*2;
      T val = recvbuf[_IDXFY(ii  ,jj  ,kk  ,n,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(ii+1,jj  ,kk  ,n,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(ii  ,jj+1,kk  ,n,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(ii+1,jj+1,kk  ,n,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(ii  ,jj  ,kk+1,n,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(ii+1,jj  ,kk+1,n,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(ii  ,jj+1,kk+1,n,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(ii+1,jj+1,kk+1,n,imax,0-gc,kmax,gc)];
      val *= 0.125;
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = val;
    }}}}
  }
  else if( levelDiff==-1 )
  {
    // coarse -> fine
    // 1cell  -> 8cell

    // 送信バッファのサイズ
    int imaxb = imax/2;
    int kmaxb = kmax/2;

    // 展開する袖数は2倍
    int gc    = vc_comm*2;

    for( int n=0;n<nmax;n++){
    for( int k=0-gc;k<kmax+gc;k++ ){
    for( int j=0-gc;j<0;j++ ){
    for( int i=0-gc;i<imax+gc;i++ ){
      int ii=(i+2)/2-1;
      int jj=(j+2)/2-1;
      int kk=(k+2)/2-1;
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = recvbuf[_IDXFY(ii,jj,kk,n,imaxb,0-vc_comm,kmaxb,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)の+Y面からの受信データの展開(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::unpackPY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                            , const int nID, int faceNo, int levelDiff, T* recvbuf )
{
  if( IsRankNull(nID) ) return CPM_SUCCESS;

  if( levelDiff==0 )
  {
    // 同じレベル
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=jmax;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = recvbuf[_IDXFY(i,j,k,n,imax,jmax,kmax,vc_comm)];
    }}}}
  }
  else if( levelDiff==1 )
  {
    // fine  -> coarse
    // 8cell -> 1cell
    // 1/4面を受信
    int is=0, ie=imax; //j方向範囲
    int ks=0, ke=kmax; //k方向範囲
    int vc_im=vc_comm, vc_ip=vc_comm; //j方向仮想セル数
    int vc_km=vc_comm, vc_kp=vc_comm; //k方向仮想セル数
    if( faceNo==0 )
    {
      ie = imax/2;
      ke = kmax/2;
      vc_ip = 0;
      vc_kp = 0;
    }
    else if( faceNo==1 )
    {
      ie = imax/2;
      ks = kmax/2;
      vc_ip = 0;
      vc_km = 0;
    }
    else if( faceNo==2 )
    {
      is = imax/2;
      ke = kmax/2;
      vc_im = 0;
      vc_kp = 0;
    }
    else if( faceNo==3 )
    {
      is = imax/2;
      ks = kmax/2;
      vc_im = 0;
      vc_km = 0;
    }

    // 送信バッファのサイズ
    int gc    = vc_comm*2;
    int jmax2 = jmax*2;

    for( int n=0;n<nmax;n++){
    for( int k=ks-vc_km;k<ke+vc_kp;k++ ){
    for( int j=jmax;j<jmax+vc_comm;j++ ){
    for( int i=is-vc_im;i<ie+vc_ip;i++ ){
      int ii=(i-is)*2;
      int jj=(j   )*2;
      int kk=(k-ks)*2;
      T val = recvbuf[_IDXFY(ii  ,jj  ,kk  ,n,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(ii+1,jj  ,kk  ,n,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(ii  ,jj+1,kk  ,n,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(ii+1,jj+1,kk  ,n,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(ii  ,jj  ,kk+1,n,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(ii+1,jj  ,kk+1,n,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(ii  ,jj+1,kk+1,n,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(ii+1,jj+1,kk+1,n,imax,jmax2,kmax,gc)];
      val *= 0.125;
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = val;
    }}}}
  }
  else if( levelDiff==-1 )
  {
    // coarse -> fine
    // 1cell  -> 8cell

    // 送信バッファのサイズ
    int imaxb = imax/2;
    int jmax2 = jmax/2;
    int kmaxb = kmax/2;

    // 展開する袖数は2倍
    int gc    = vc_comm*2;

    for( int n=0;n<nmax;n++){
    for( int k=0-gc;k<kmax+gc;k++ ){
    for( int j=jmax;j<jmax+gc;j++ ){
    for( int i=0-gc;i<imax+gc;i++ ){
      int ii=(i+2)/2-1;
      int jj=j/2;
      int kk=(k+2)/2-1;
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = recvbuf[_IDXFY(ii,jj,kk,n,imaxb,jmax2,kmaxb,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)の-Z面への送信データのパック(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::packMZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                          , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw )
{
  if( IsRankNull(nID) ) return CPM_SUCCESS;

  // 送信範囲の確定
  int is=0, ie=imax; //i方向範囲
  int js=0, je=jmax; //j方向範囲
  int gc=vc_comm;    //通信層数
  if( levelDiff==0 )
  {
    // 同じレベル
    // そのまま
  }
  else if( levelDiff==1 )
  {
    // coarse -> fine
    // 1/4面を送信
    if( faceNo==0 )
    {
      ie = imax/2;
      je = jmax/2;
    }
    else if( faceNo==1 )
    {
      is = imax/2;
      je = jmax/2;
    }
    else if( faceNo==2 )
    {
      ie = imax/2;
      js = jmax/2;
    }
    else if( faceNo==3 )
    {
      is = imax/2;
      js = jmax/2;
    }
  }
  else if( levelDiff==-1 )
  {
    // fine -> coarse
    // 層数を2倍
    gc = vc_comm*2;
  }

  // 送信バッファのサイズ
  int imaxb = ie-is;
  int jmaxb = je-js;
  size_t sz = size_t(imaxb+2*gc) * size_t(jmaxb+2*gc) * size_t(gc) * size_t(nmax);
  if( sz > nw )
  {
    stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // 格納
  for( int n=0;n<nmax;n++){
  for( int k=0;k<gc;k++ ){
  for( int j=js-gc,jj=0-gc;j<je+gc;j++,jj++ ){
  for( int i=is-gc,ii=0-gc;i<ie+gc;i++,ii++ ){
    sendbuf[_IDXFZ(ii,jj,k,n,imaxb,jmaxb,0,gc)] = array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)];
  }}}}

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)の+Z面への送信データのパック(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::packPZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                          , const int nID, int faceNo, int levelDiff, T* sendbuf, size_t nw )
{
  if( IsRankNull(nID) ) return CPM_SUCCESS;

  // 送信範囲の確定
  int is=0, ie=imax; //i方向範囲
  int js=0, je=jmax; //j方向範囲
  int gc=vc_comm;    //通信層数
  if( levelDiff==0 )
  {
    // 同じレベル
    // そのまま
  }
  else if( levelDiff==1 )
  {
    // coarse -> fine
    // 1/4面を送信
    if( faceNo==0 )
    {
      ie = imax/2;
      je = jmax/2;
    }
    else if( faceNo==1 )
    {
      is = imax/2;
      je = jmax/2;
    }
    else if( faceNo==2 )
    {
      ie = imax/2;
      js = jmax/2;
    }
    else if( faceNo==3 )
    {
      is = imax/2;
      js = jmax/2;
    }
  }
  else if( levelDiff==-1 )
  {
    // fine -> coarse
    // 層数を2倍
    gc = vc_comm*2;
  }

  // 送信バッファのサイズ
  int imaxb = ie-is;
  int jmaxb = je-js;
  size_t sz = size_t(imaxb+2*gc) * size_t(jmaxb+2*gc) * size_t(gc) * size_t(nmax);
  if( sz > nw )
  {
    stmpd_printf("[%d] BndComm buffer length for send error\n", m_rankNo);
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // 格納
  for( int n=0;n<nmax;n++){
  for( int k=kmax-gc;k<kmax;k++ ){
  for( int j=js-gc,jj=0-gc;j<je+gc;j++,jj++ ){
  for( int i=is-gc,ii=0-gc;i<ie+gc;i++,ii++ ){
    sendbuf[_IDXFZ(ii,jj,k,n,imaxb,jmaxb,kmax-gc,gc)] = array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)];
  }}}}

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)の-Z面からの受信データの展開(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::unpackMZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                            , const int nID, int faceNo, int levelDiff, T* recvbuf )
{
  if( IsRankNull(nID) ) return CPM_SUCCESS;

  if( levelDiff==0 )
  {
    // 同じレベル
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<0;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = recvbuf[_IDXFZ(i,j,k,n,imax,jmax,0-vc_comm,vc_comm)];
    }}}}
  }
  else if( levelDiff==1 )
  {
    // fine  -> coarse
    // 8cell -> 1cell
    // 1/4面を受信
    int is=0, ie=imax; //i方向範囲
    int js=0, je=jmax; //j方向範囲
    int vc_im=vc_comm, vc_ip=vc_comm; //i方向仮想セル数
    int vc_jm=vc_comm, vc_jp=vc_comm; //j方向仮想セル数
    if( faceNo==0 )
    {
      ie = imax/2;
      je = jmax/2;
      vc_ip = 0;
      vc_jp = 0;
    }
    else if( faceNo==1 )
    {
      is = imax/2;
      je = jmax/2;
      vc_im = 0;
      vc_jp = 0;
    }
    else if( faceNo==2 )
    {
      ie = imax/2;
      js = jmax/2;
      vc_ip = 0;
      vc_jm = 0;
    }
    else if( faceNo==3 )
    {
      is = imax/2;
      js = jmax/2;
      vc_im = 0;
      vc_jm = 0;
    }

    // 送信バッファのサイズ
    int gc    = vc_comm*2;

    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<0;k++ ){
    for( int j=js-vc_jm;j<je+vc_jp;j++ ){
    for( int i=is-vc_im;i<ie+vc_ip;i++ ){
      int ii=(i-is)*2;
      int jj=(j-js)*2;
      int kk=(k   )*2;
      T val = recvbuf[_IDXFZ(ii  ,jj  ,kk  ,n,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(ii+1,jj  ,kk  ,n,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(ii  ,jj+1,kk  ,n,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(ii+1,jj+1,kk  ,n,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(ii  ,jj  ,kk+1,n,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(ii+1,jj  ,kk+1,n,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(ii  ,jj+1,kk+1,n,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(ii+1,jj+1,kk+1,n,imax,jmax,0-gc,gc)];
      val *= 0.125;
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = val;
    }}}}
  }
  else if( levelDiff==-1 )
  {
    // coarse -> fine
    // 1cell  -> 8cell

    // 送信バッファのサイズ
    int imaxb = imax/2;
    int jmaxb = jmax/2;

    // 展開する袖数は2倍
    int gc    = vc_comm*2;

    for( int n=0;n<nmax;n++){
    for( int k=0-gc;k<0;k++ ){
    for( int j=0-gc;j<jmax+gc;j++ ){
    for( int i=0-gc;i<imax+gc;i++ ){
      int ii=(i+2)/2-1;
      int jj=(j+2)/2-1;
      int kk=(k+2)/2-1;
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = recvbuf[_IDXFZ(ii,jj,kk,n,imaxb,jmaxb,0-vc_comm,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)の+Z面からの受信データの展開(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::unpackPZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                            , const int nID, int faceNo, int levelDiff, T* recvbuf )
{
  if( IsRankNull(nID) ) return CPM_SUCCESS;

  if( levelDiff==0 )
  {
    // 同じレベル
    for( int n=0;n<nmax;n++){
    for( int k=kmax;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = recvbuf[_IDXFZ(i,j,k,n,imax,jmax,kmax,vc_comm)];
    }}}}
  }
  else if( levelDiff==1 )
  {
    // fine  -> coarse
    // 8cell -> 1cell
    // 1/4面を受信
    int is=0, ie=imax; //j方向範囲
    int js=0, je=jmax; //j方向範囲
    int vc_im=vc_comm, vc_ip=vc_comm; //j方向仮想セル数
    int vc_jm=vc_comm, vc_jp=vc_comm; //j方向仮想セル数
    if( faceNo==0 )
    {
      ie = imax/2;
      je = jmax/2;
      vc_ip = 0;
      vc_jp = 0;
    }
    else if( faceNo==1 )
    {
      is = imax/2;
      je = jmax/2;
      vc_im = 0;
      vc_jp = 0;
    }
    else if( faceNo==2 )
    {
      ie = imax/2;
      js = jmax/2;
      vc_ip = 0;
      vc_jm = 0;
    }
    else if( faceNo==3 )
    {
      is = imax/2;
      js = jmax/2;
      vc_im = 0;
      vc_jm = 0;
    }

    // 送信バッファのサイズ
    int gc    = vc_comm*2;
    int kmax2 = kmax*2;

    for( int n=0;n<nmax;n++){
    for( int k=kmax;k<kmax+vc_comm;k++ ){
    for( int j=js-vc_jm;j<je+vc_jp;j++ ){
    for( int i=is-vc_im;i<ie+vc_ip;i++ ){
      int ii=(i-is)*2;
      int jj=(j-js)*2;
      int kk=(k   )*2;
      T val = recvbuf[_IDXFZ(ii  ,jj  ,kk  ,n,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(ii+1,jj  ,kk  ,n,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(ii  ,jj+1,kk  ,n,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(ii+1,jj+1,kk  ,n,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(ii  ,jj  ,kk+1,n,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(ii+1,jj  ,kk+1,n,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(ii  ,jj+1,kk+1,n,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(ii+1,jj+1,kk+1,n,imax,jmax,kmax2,gc)];
      val *= 0.125;
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = val;
    }}}}
  }
  else if( levelDiff==-1 )
  {
    // coarse -> fine
    // 1cell  -> 8cell

    // 送信バッファのサイズ
    int imaxb = imax/2;
    int jmaxb = jmax/2;
    int kmax2 = kmax/2;

    // 展開する袖数は2倍
    int gc    = vc_comm*2;

    for( int n=0;n<nmax;n++){
    for( int k=kmax;k<kmax+gc;k++ ){
    for( int j=0-gc;j<jmax+gc;j++ ){
    for( int i=0-gc;i<imax+gc;i++ ){
      int ii=(i+2)/2-1;
      int jj=(j+2)/2-1;
      int kk=k/2;
      array[_IDX_S4D(i,j,k,n,imax,jmax,kmax,vc)] = recvbuf[_IDXFZ(ii,jj,kk,n,imaxb,jmaxb,kmax2,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}



#undef _IDXFX
#undef _IDXFY
#undef _IDXFZ

#endif /* _CPM_PARAMANAGER_BNDCOMM_LMR_H_ */

