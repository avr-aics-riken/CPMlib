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
 * @file   cpm_ParaManager_BndComm.h
 * カーテシアン用パラレルマネージャクラスのインラインヘッダーファイル
 * @date   2015/03/27
 */

#ifndef _CPM_PARAMANAGER_BNDCOMM_H_
#define _CPM_PARAMANAGER_BNDCOMM_H_

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
// 袖通信(Scalar3D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm, int procGrpNo, CPM_PADDING padding )
{
  int sz[3] = {imax, jmax, kmax};
  int pad_size[4] = {0, 0, 0, 0};
  if( padding )
  {
    GetPaddingSize(CPM_ARRAY_S3D, sz, vc, pad_size);
  }
  return BndCommS4D( array, imax, jmax, kmax, 1, vc, vc_comm, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Vector3D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm, int procGrpNo, CPM_PADDING padding )
{
  int sz[3] = {imax, jmax, kmax};
  int pad_size[4] = {0, 0, 0, 0};
  if( padding )
  {
    GetPaddingSize(CPM_ARRAY_V3D, sz, vc, pad_size);
  }
  return BndCommS4D( array, imax, jmax, kmax, 3, vc, vc_comm, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int procGrpNo, CPM_PADDING padding )
{
  int sz[3] = {imax, jmax, kmax};
  int pad_size[4] = {0, 0, 0, 0};
  if( padding )
  {
    GetPaddingSize(CPM_ARRAY_S4D, sz, vc, pad_size, nmax);
  }
  return BndCommS4D( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D版、waitなし)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS3D_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                   , MPI_Request req[48], int procGrpNo, CPM_PADDING padding )
{
  int sz[3] = {imax, jmax, kmax};
  int pad_size[4] = {0, 0, 0, 0};
  if( padding )
  {
    GetPaddingSize(CPM_ARRAY_S3D, sz, vc, pad_size);
  }
  return BndCommS4D_nowait( array, imax, jmax, kmax, 1, vc, vc_comm, req, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Vector3D版、waitなし)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommV3D_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                   , MPI_Request req[48], int procGrpNo, CPM_PADDING padding )
{
  int sz[3] = {imax, jmax, kmax};
  int pad_size[4] = {0, 0, 0, 0};
  if( padding )
  {
    GetPaddingSize(CPM_ARRAY_V3D, sz, vc, pad_size);
  }
  return BndCommS4D_nowait( array, imax, jmax, kmax, 3, vc, vc_comm, req, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4D版、waitなし)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS4D_nowait( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                   , MPI_Request req[48], int procGrpNo, CPM_PADDING padding )
{
  int sz[3] = {imax, jmax, kmax};
  int pad_size[4] = {0, 0, 0, 0};
  if( padding )
  {
    GetPaddingSize(CPM_ARRAY_S4D, sz, vc, pad_size, nmax);
  }
  return BndCommS4D_nowait( array, imax, jmax, kmax, nmax, vc, vc_comm, req, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信のwait、展開(Scalar3D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::wait_BndCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                , MPI_Request req[48], int procGrpNo, CPM_PADDING padding )
{
  int sz[3] = {imax, jmax, kmax};
  int pad_size[4] = {0, 0, 0, 0};
  if( padding )
  {
    GetPaddingSize(CPM_ARRAY_S3D, sz, vc, pad_size);
  }
  return wait_BndCommS4D( array, imax, jmax, kmax, 1, vc, vc_comm, req, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信のwait、展開(Vector3D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::wait_BndCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                , MPI_Request req[48], int procGrpNo, CPM_PADDING padding )
{
  int sz[3] = {imax, jmax, kmax};
  int pad_size[4] = {0, 0, 0, 0};
  if( padding )
  {
    GetPaddingSize(CPM_ARRAY_V3D, sz, vc, pad_size);
  }
  return wait_BndCommS4D( array, imax, jmax, kmax, 3, vc, vc_comm, req, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信のwait、展開(Scalar4D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::wait_BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                , MPI_Request req[48], int procGrpNo, CPM_PADDING padding )
{
  int sz[3] = {imax, jmax, kmax};
  int pad_size[4] = {0, 0, 0, 0};
  if( padding )
  {
    GetPaddingSize(CPM_ARRAY_S4D, sz, vc, pad_size, nmax);
  }
  return wait_BndCommS4D( array, imax, jmax, kmax, nmax, vc, vc_comm, req, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar3D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::PeriodicCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo, CPM_PADDING padding )
{
  int sz[3] = {imax, jmax, kmax};
  int pad_size[4] = {0, 0, 0, 0};
  if( padding )
  {
    GetPaddingSize(CPM_ARRAY_S3D, sz, vc, pad_size);
  }
  return PeriodicCommS4D( array, imax, jmax, kmax, 1, vc, vc_comm, dir, pm, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Vector3D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::PeriodicCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo, CPM_PADDING padding )
{
  int sz[3] = {imax, jmax, kmax};
  int pad_size[4] = {0, 0, 0, 0};
  if( padding )
  {
    GetPaddingSize(CPM_ARRAY_V3D, sz, vc, pad_size);
  }
  return PeriodicCommS4D( array, imax, jmax, kmax, 3, vc, vc_comm, dir, pm, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar4D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::PeriodicCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo, CPM_PADDING padding )
{
  int sz[3] = {imax, jmax, kmax};
  int pad_size[4] = {0, 0, 0, 0};
  if( padding )
  {
    GetPaddingSize(CPM_ARRAY_S4D, sz, vc, pad_size, nmax);
  }
  return PeriodicCommS4D( array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, pad_size, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                               , int pad_size[4], int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 通信バッファを取得
  S_BNDCOMM_BUFFER *bufInfo = GetBndCommBuffer(procGrpNo);
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }

  // 隣接ランクを取得
  const int *nID = GetNeighborRankID(procGrpNo);
  if( !nID )
  {
    return CPM_ERROR_GET_NEIGHBOR_RANK;
  }

  // 通信バッファサイズを計算
  size_t nwX = size_t(jmax+2*vc_comm) * size_t(kmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwY = size_t(kmax+2*vc_comm) * size_t(imax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwZ = size_t(imax+2*vc_comm) * size_t(jmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  if( nwX > bufInfo->m_nwX || nwY > bufInfo->m_nwY || nwZ > bufInfo->m_nwZ )
  {
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // ポインタをセット
  T *sendmx = (T*)(bufInfo->m_bufX[0]);
  T *recvmx = (T*)(bufInfo->m_bufX[1]);
  T *sendpx = (T*)(bufInfo->m_bufX[2]);
  T *recvpx = (T*)(bufInfo->m_bufX[3]);
  T *sendmy = (T*)(bufInfo->m_bufY[0]);
  T *recvmy = (T*)(bufInfo->m_bufY[1]);
  T *sendpy = (T*)(bufInfo->m_bufY[2]);
  T *recvpy = (T*)(bufInfo->m_bufY[3]);
  T *sendmz = (T*)(bufInfo->m_bufZ[0]);
  T *recvmz = (T*)(bufInfo->m_bufZ[1]);
  T *sendpz = (T*)(bufInfo->m_bufZ[2]);
  T *recvpz = (T*)(bufInfo->m_bufZ[3]);

  MPI_Request req[12];
  for( int i=0;i<12;i++ ) req[i] = MPI_REQUEST_NULL;

  //// X face ////
  int nIDmx = nID[X_MINUS];
  int nIDpx = nID[X_PLUS];

  // pack
  if( (ret = packX( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, sendmx, sendpx, nIDmx, nIDpx, procGrpNo )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmx, recvmx, sendpx, recvpx, nwX, &req[0], nIDmx, nIDmx, nIDpx, nIDpx, procGrpNo )) != CPM_SUCCESS ) return ret;

  // wait
  if( (ret = Waitall( 4, &req[0] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackX( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, recvmx, recvpx, nIDmx, nIDpx )) != CPM_SUCCESS ) return ret;

  //// Y face ////
  int nIDmy = nID[Y_MINUS];
  int nIDpy = nID[Y_PLUS];

  // pack
  if( (ret = packY( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, sendmy, sendpy, nIDmy, nIDpy, procGrpNo )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmy, recvmy, sendpy, recvpy, nwY, &req[4], nIDmy, nIDmy, nIDpy, nIDpy, procGrpNo )) != CPM_SUCCESS ) return ret;

  // wait
  if( (ret = Waitall( 4, &req[4] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackY( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, recvmy, recvpy, nIDmy, nIDpy )) != CPM_SUCCESS ) return ret;

  //// Z face ////
  int nIDmz = nID[Z_MINUS];
  int nIDpz = nID[Z_PLUS];

  // pack
  if( (ret = packZ( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, sendmz, sendpz, nIDmz, nIDpz, procGrpNo )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmz, recvmz, sendpz, recvpz, nwZ, &req[8], nIDmz, nIDmz, nIDpz, nIDpz, procGrpNo )) != CPM_SUCCESS ) return ret;

  // wait
  if( (ret = Waitall( 4, &req[8] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackZ( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, recvmz, recvpz, nIDmz, nIDpz )) != CPM_SUCCESS ) return ret;

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4D版、waitなし)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS4D_nowait( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                       , MPI_Request req[12], int pad_size[4], int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  for( int i=0;i<12;i++ )
  {
    req[i] = MPI_REQUEST_NULL;
  }

  // 通信バッファを取得
  S_BNDCOMM_BUFFER *bufInfo = GetBndCommBuffer(procGrpNo);
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }

  // 隣接ランクを取得
  const int *nID = GetNeighborRankID(procGrpNo);
  if( !nID )
  {
    return CPM_ERROR_GET_NEIGHBOR_RANK;
  }

  // 通信バッファサイズを計算
  size_t nwX = size_t(jmax+2*vc_comm) * size_t(kmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwY = size_t(kmax+2*vc_comm) * size_t(imax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwZ = size_t(imax+2*vc_comm) * size_t(jmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  if( nwX > bufInfo->m_nwX || nwY > bufInfo->m_nwY || nwZ > bufInfo->m_nwZ )
  {
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // ポインタをセット
  T *sendmx = (T*)(bufInfo->m_bufX[0]);
  T *recvmx = (T*)(bufInfo->m_bufX[1]);
  T *sendpx = (T*)(bufInfo->m_bufX[2]);
  T *recvpx = (T*)(bufInfo->m_bufX[3]);
  T *sendmy = (T*)(bufInfo->m_bufY[0]);
  T *recvmy = (T*)(bufInfo->m_bufY[1]);
  T *sendpy = (T*)(bufInfo->m_bufY[2]);
  T *recvpy = (T*)(bufInfo->m_bufY[3]);
  T *sendmz = (T*)(bufInfo->m_bufZ[0]);
  T *recvmz = (T*)(bufInfo->m_bufZ[1]);
  T *sendpz = (T*)(bufInfo->m_bufZ[2]);
  T *recvpz = (T*)(bufInfo->m_bufZ[3]);

  //// X face ////
  int nIDmx = nID[X_MINUS];
  int nIDpx = nID[X_PLUS];

  // pack
  if( (ret = packX( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, sendmx, sendpx, nIDmx, nIDpx, procGrpNo )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmx, recvmx, sendpx, recvpx, nwX, &req[0], nIDmx, nIDmx, nIDpx, nIDpx, procGrpNo )) != CPM_SUCCESS ) return ret;

  //// Y face ////
  int nIDmy = nID[Y_MINUS];
  int nIDpy = nID[Y_PLUS];

  // pack
  if( (ret = packY( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, sendmy, sendpy, nIDmy, nIDpy, procGrpNo )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmy, recvmy, sendpy, recvpy, nwY, &req[4], nIDmy, nIDmy, nIDpy, nIDpy, procGrpNo )) != CPM_SUCCESS ) return ret;

  //// Z face ////
  int nIDmz = nID[Z_MINUS];
  int nIDpz = nID[Z_PLUS];

  // pack
  if( (ret = packZ( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, sendmz, sendpz, nIDmz, nIDpz, procGrpNo )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmz, recvmz, sendpz, recvpz, nwZ, &req[8], nIDmz, nIDmz, nIDpz, nIDpz, procGrpNo )) != CPM_SUCCESS ) return ret;

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信のwait、展開(Scalar4D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::wait_BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                    , MPI_Request req[12], int pad_size[4], int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 通信バッファを取得
  S_BNDCOMM_BUFFER *bufInfo = GetBndCommBuffer(procGrpNo);
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }

  // 隣接ランクを取得
  const int *nID = GetNeighborRankID(procGrpNo);
  if( !nID )
  {
    return CPM_ERROR_GET_NEIGHBOR_RANK;
  }

  // 通信バッファサイズを計算
  size_t nwX = size_t(jmax+2*vc_comm) * size_t(kmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwY = size_t(kmax+2*vc_comm) * size_t(imax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwZ = size_t(imax+2*vc_comm) * size_t(jmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  if( nwX > bufInfo->m_nwX || nwY > bufInfo->m_nwY || nwZ > bufInfo->m_nwZ )
  {
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // ポインタをセット
  T *recvmx = (T*)(bufInfo->m_bufX[1]);
  T *recvpx = (T*)(bufInfo->m_bufX[3]);
  T *recvmy = (T*)(bufInfo->m_bufY[1]);
  T *recvpy = (T*)(bufInfo->m_bufY[3]);
  T *recvmz = (T*)(bufInfo->m_bufZ[1]);
  T *recvpz = (T*)(bufInfo->m_bufZ[3]);

  //// X face ////
  int nIDmx = nID[X_MINUS];
  int nIDpx = nID[X_PLUS];

  // wait
  if( (ret = Waitall( 4, &req[0] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackX( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, recvmx, recvpx, nIDmx, nIDpx )) != CPM_SUCCESS ) return ret;

  //// Y face ////
  int nIDmy = nID[Y_MINUS];
  int nIDpy = nID[Y_PLUS];

  // wait
  if( (ret = Waitall( 4, &req[4] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackY( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, recvmy, recvpy, nIDmy, nIDpy )) != CPM_SUCCESS ) return ret;

  //// Z face ////
  int nIDmz = nID[Z_MINUS];
  int nIDpz = nID[Z_PLUS];

  // wait
  if( (ret = Waitall( 4, &req[8] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackZ( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, recvmz, recvpz, nIDmz, nIDpz )) != CPM_SUCCESS ) return ret;

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar4D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::PeriodicCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                    , cpm_DirFlag dir, cpm_PMFlag pm, int pad_size[4], int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 通信バッファを取得
  S_BNDCOMM_BUFFER *bufInfo = GetBndCommBuffer(procGrpNo);
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }

  // 周期境界隣接ランクを取得
  const int *nID = GetPeriodicRankID(procGrpNo);
  if( !nID )
  {
    return CPM_ERROR_GET_PERIODIC_RANK;
  }

  // 通信バッファサイズを計算
  size_t nwX = size_t(jmax+2*vc_comm) * size_t(kmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwY = size_t(kmax+2*vc_comm) * size_t(imax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwZ = size_t(imax+2*vc_comm) * size_t(jmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  if( nwX > bufInfo->m_nwX || nwY > bufInfo->m_nwY || nwZ > bufInfo->m_nwZ )
  {
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // ポインタをセット
  T *sendmx = (T*)(bufInfo->m_bufX[0]);
  T *recvmx = (T*)(bufInfo->m_bufX[1]);
  T *sendpx = (T*)(bufInfo->m_bufX[2]);
  T *recvpx = (T*)(bufInfo->m_bufX[3]);
  T *sendmy = (T*)(bufInfo->m_bufY[0]);
  T *recvmy = (T*)(bufInfo->m_bufY[1]);
  T *sendpy = (T*)(bufInfo->m_bufY[2]);
  T *recvpy = (T*)(bufInfo->m_bufY[3]);
  T *sendmz = (T*)(bufInfo->m_bufZ[0]);
  T *recvmz = (T*)(bufInfo->m_bufZ[1]);
  T *sendpz = (T*)(bufInfo->m_bufZ[2]);
  T *recvpz = (T*)(bufInfo->m_bufZ[3]);

  MPI_Request req[12];
  for( int i=0;i<12;i++ ) req[i] = MPI_REQUEST_NULL;

  //// X face ////
  int nIDsmx = cpm_Base::getRankNull();
  int nIDspx = cpm_Base::getRankNull();
  int nIDrmx = cpm_Base::getRankNull();
  int nIDrpx = cpm_Base::getRankNull();
  if( dir==X_DIR )
  {
    if( pm==MINUS2PLUS || pm==BOTH )
    {
      nIDsmx = nID[X_MINUS];
      nIDrpx = nID[X_PLUS];
    }
    if( pm==PLUS2MINUS || pm==BOTH )
    {
      nIDspx = nID[X_PLUS];
      nIDrmx = nID[X_MINUS];
    }
  }

  // pack
  if( (ret = packX( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, sendmx, sendpx, nIDsmx, nIDspx, procGrpNo )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmx, recvmx, sendpx, recvpx, nwX, &req[0], nIDsmx, nIDrmx, nIDspx, nIDrpx, procGrpNo )) != CPM_SUCCESS ) return ret;

  // wait
  if( (ret = Waitall( 4, &req[0] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackX( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, recvmx, recvpx, nIDrmx, nIDrpx )) != CPM_SUCCESS ) return ret;

  //// Y face ////
  int nIDsmy = cpm_Base::getRankNull();
  int nIDspy = cpm_Base::getRankNull();
  int nIDrmy = cpm_Base::getRankNull();
  int nIDrpy = cpm_Base::getRankNull();
  if( dir==Y_DIR )
  {
    if( pm==MINUS2PLUS || pm==BOTH )
    {
      nIDsmy = nID[Y_MINUS];
      nIDrpy = nID[Y_PLUS];
    }
    if( pm==PLUS2MINUS || pm==BOTH )
    {
      nIDspy = nID[Y_PLUS];
      nIDrmy = nID[Y_MINUS];
    }
  }

  // pack
  if( (ret = packY( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, sendmy, sendpy, nIDsmy, nIDspy, procGrpNo )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmy, recvmy, sendpy, recvpy, nwY, &req[4], nIDsmy, nIDrmy, nIDspy, nIDrpy, procGrpNo )) != CPM_SUCCESS ) return ret;

  // wait
  if( (ret = Waitall( 4, &req[4] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackY( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, recvmy, recvpy, nIDrmy, nIDrpy )) != CPM_SUCCESS ) return ret;

  //// Z face ////
  int nIDsmz = cpm_Base::getRankNull();
  int nIDspz = cpm_Base::getRankNull();
  int nIDrmz = cpm_Base::getRankNull();
  int nIDrpz = cpm_Base::getRankNull();
  if( dir==Z_DIR )
  {
    if( pm==MINUS2PLUS || pm==BOTH )
    {
      nIDsmz = nID[Z_MINUS];
      nIDrpz = nID[Z_PLUS];
    }
    if( pm==PLUS2MINUS || pm==BOTH )
    {
      nIDspz = nID[Z_PLUS];
      nIDrmz = nID[Z_MINUS];
    }
  }

  // pack
  if( (ret = packZ( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, sendmz, sendpz, nIDsmz, nIDspz, procGrpNo )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmz, recvmz, sendpz, recvpz, nwZ, &req[8], nIDsmz, nIDrmz, nIDspz, nIDrpz, procGrpNo )) != CPM_SUCCESS ) return ret;

  // wait
  if( (ret = Waitall( 4, &req[8] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackZ( array, imax, jmax, kmax, nmax, vc, vc_comm, pad_size, recvmz, recvpz, nIDrmz, nIDrpz )) != CPM_SUCCESS ) return ret;

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)のX方向送信バッファのセット
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::packX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int pad_size[4]
                          , T *sendm, T *sendp, int nIDm, int nIDp , int procGrpNo )
{
  //オフセット量のセット
  int is = 0;
  //定義点がVOXELのとき
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FVM ) is = 0;
  //定義点がNODEのとき
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FDM ) is = 1;

  // パディング
  int ip = pad_size[0];
  int jp = pad_size[1];
  int kp = pad_size[2];
  int np = pad_size[3];

  if( !IsRankNull(nIDm) )
  {
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0;i<vc_comm;i++ ){

      sendm[_IDXFX(i,j,k,n,0,jmax,kmax,vc_comm)] = array[_IDX_S4D_PAD(i+is,j,k,n,imax,jmax,kmax,vc,ip,jp,kp)];

    }}}}
  }

  if( !IsRankNull(nIDp) )
  {
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=imax-vc_comm;i<imax;i++ ){

      sendp[_IDXFX(i,j,k,n,imax-vc_comm,jmax,kmax,vc_comm)] = array[_IDX_S4D_PAD(i-is,j,k,n,imax,jmax,kmax,vc,ip,jp,kp)];

    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)のX方向受信バッファを元に戻す
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::unpackX( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int pad_size[4]
                            , T *recvm, T *recvp, int nIDm, int nIDp )
{
  // パディング
  int ip = pad_size[0];
  int jp = pad_size[1];
  int kp = pad_size[2];
  int np = pad_size[3];

  if( !IsRankNull(nIDm) )
  {
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<0;i++ ){
      array[_IDX_S4D_PAD(i,j,k,n,imax,jmax,kmax,vc,ip,jp,kp)] = recvm[_IDXFX(i,j,k,n,0-vc_comm,jmax,kmax,vc_comm)];
    }}}}
  }

  if( !IsRankNull(nIDp) )
  {
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=imax;i<imax+vc_comm;i++ ){
      array[_IDX_S4D_PAD(i,j,k,n,imax,jmax,kmax,vc,ip,jp,kp)] = recvp[_IDXFX(i,j,k,n,imax,jmax,kmax,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)のY方向送信バッファのセット
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::packY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int pad_size[4]
                          , T *sendm, T *sendp, int nIDm, int nIDp, int procGrpNo )
{
  //オフセット量のセット
  int js = 0;
  //定義点がVOXELのとき
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FVM ) js = 0;
  //定義点がNODEのとき
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FDM ) js = 1;

  // パディング
  int ip = pad_size[0];
  int jp = pad_size[1];
  int kp = pad_size[2];
  int np = pad_size[3];

  if( !IsRankNull(nIDm) )
  {
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0;j<vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){

      sendm[_IDXFY(i,j,k,n,imax,0,kmax,vc_comm)] = array[_IDX_S4D_PAD(i,j+js,k,n,imax,jmax,kmax,vc,ip,jp,kp)];

    }}}}
  }

  if( !IsRankNull(nIDp) )
  {
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=jmax-vc_comm;j<jmax;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){

      sendp[_IDXFY(i,j,k,n,imax,jmax-vc_comm,kmax,vc_comm)] = array[_IDX_S4D_PAD(i,j-js,k,n,imax,jmax,kmax,vc,ip,jp,kp)];

    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)のY方向受信バッファを元に戻す
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::unpackY( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int pad_size[4]
                            , T *recvm, T *recvp, int nIDm, int nIDp )
{
  // パディング
  int ip = pad_size[0];
  int jp = pad_size[1];
  int kp = pad_size[2];
  int np = pad_size[3];

  if( !IsRankNull(nIDm) )
  {
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<0;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
      array[_IDX_S4D_PAD(i,j,k,n,imax,jmax,kmax,vc,ip,jp,kp)] = recvm[_IDXFY(i,j,k,n,imax,0-vc_comm,kmax,vc_comm)];
    }}}}
  }

  if( !IsRankNull(nIDp) )
  {
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=jmax;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
      array[_IDX_S4D_PAD(i,j,k,n,imax,jmax,kmax,vc,ip,jp,kp)] = recvp[_IDXFY(i,j,k,n,imax,jmax,kmax,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)のZ方向送信バッファのセット
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::packZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int pad_size[4]
                          , T *sendm, T *sendp, int nIDm, int nIDp, int procGrpNo )
{

  //オフセット量のセット
  int ks = 0;
  //定義点がVOXELのとき
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FVM ) ks = 0;
  //定義点がNODEのとき
  if( GetDefPointType(procGrpNo) == CPM_DEFPOINTTYPE_FDM ) ks = 1;

  // パディング
  int ip = pad_size[0];
  int jp = pad_size[1];
  int kp = pad_size[2];
  int np = pad_size[3];

  if( !IsRankNull(nIDm) )
  {
    for( int n=0;n<nmax;n++){
    for( int k=0;k<vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){

      sendm[_IDXFZ(i,j,k,n,imax,jmax,0,vc_comm)] = array[_IDX_S4D_PAD(i,j,k+ks,n,imax,jmax,kmax,vc,ip,jp,kp)];

    }}}}
  }

  if( !IsRankNull(nIDp) )
  {
    for( int n=0;n<nmax;n++){
    for( int k=kmax-vc_comm;k<kmax;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){

      sendp[_IDXFZ(i,j,k,n,imax,jmax,kmax-vc_comm,vc_comm)] = array[_IDX_S4D_PAD(i,j,k-ks,n,imax,jmax,kmax,vc,ip,jp,kp)];

    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D,4D,Vector3D版)のZ方向受信バッファを元に戻す
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::unpackZ( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int pad_size[4]
                            , T *recvm, T *recvp, int nIDm, int nIDp )
{
  // パディング
  int ip = pad_size[0];
  int jp = pad_size[1];
  int kp = pad_size[2];
  int np = pad_size[3];

  if( !IsRankNull(nIDm) )
  {
    for( int n=0;n<nmax;n++){
    for( int k=0-vc_comm;k<0;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
      array[_IDX_S4D_PAD(i,j,k,n,imax,jmax,kmax,vc,ip,jp,kp)] = recvm[_IDXFZ(i,j,k,n,imax,jmax,0-vc_comm,vc_comm)];
    }}}}
  }

  if( !IsRankNull(nIDp) )
  {
    for( int n=0;n<nmax;n++){
    for( int k=kmax;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
      array[_IDX_S4D_PAD(i,j,k,n,imax,jmax,kmax,vc,ip,jp,kp)] = recvp[_IDXFZ(i,j,k,n,imax,jmax,kmax,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// １方向(プラス、マイナス)の双方向袖通信処理
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::sendrecv( T *sendm, T *recvm, T *sendp, T *recvp, size_t nw, MPI_Request *req
                             , int nIDsm, int nIDrm, int nIDsp, int nIDrp, int procGrpNo )
{
  cpm_ErrorCode ret;

  MPI_Request r0 = MPI_REQUEST_NULL;
  MPI_Request r1 = MPI_REQUEST_NULL;
  MPI_Request r2 = MPI_REQUEST_NULL;
  MPI_Request r3 = MPI_REQUEST_NULL;

  if( !IsRankNull(nIDrm) )
  {
    if( (ret = Irecv( recvm, nw, nIDrm, &r0, procGrpNo )) != CPM_SUCCESS )
    {
      return ret;
    }
  }

  if( !IsRankNull(nIDrp) )
  {
    if( (ret = Irecv( recvp, nw, nIDrp, &r1, procGrpNo )) != CPM_SUCCESS )
    {
      return ret;
    }
  }

  if( !IsRankNull(nIDsp) )
  {
    if( (ret = Isend( sendp, nw, nIDsp, &r2, procGrpNo )) != CPM_SUCCESS )
    {
      return ret;
    }
  }

  if( !IsRankNull(nIDsm) )
  {
    if( (ret = Isend( sendm, nw, nIDsm, &r3, procGrpNo )) != CPM_SUCCESS )
    {
      return ret;
    }
  }

  req[0] = r0;
  req[1] = r1;
  req[2] = r2;
  req[3] = r3;

  return CPM_SUCCESS;
}

#undef _IDXFX
#undef _IDXFY
#undef _IDXFZ

#endif /* _CPM_PARAMANAGER_BNDCOMM_H_ */
