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
 * @file   cpm_ParaManager_BndComm.h
 * パラレルマネージャクラスのインラインヘッダーファイル
 * @date   2012/05/31
 */

#ifndef _CPM_PARAMANAGER_BNDCOMM_H_
#define _CPM_PARAMANAGER_BNDCOMM_H_

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm, int procGrpNo )
{
  return BndCommS4D( array, imax, jmax, kmax, 1, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Vector3D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm, int procGrpNo )
{
  return BndCommS4D( array, imax, jmax, kmax, 3, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm, int procGrpNo )
{
  if( GetDomainType() == CPM_DOMAIN_CARTESIAN )
  {
    cpm_ParaManagerCART *pCART = (cpm_ParaManagerCART*)this;
    return pCART->BndCommS4D(array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  }
  else if( GetDomainType() == CPM_DOMAIN_LMR )
  {
    cpm_ParaManagerLMR *pLMR = (cpm_ParaManagerLMR*)this;
    return pLMR->BndCommS4D(array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  }

  return CPM_ERROR_BNDCOMM;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D版、waitなし)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS3D_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                   , MPI_Request req[48], int procGrpNo )
{
  return BndCommS4D_nowait( array, imax, jmax, kmax, 1, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Vector3D版、waitなし)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommV3D_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                   , MPI_Request req[48], int procGrpNo )
{
  return BndCommS4D_nowait( array, imax, jmax, kmax, 3, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4D版、waitなし)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS4D_nowait( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                   , MPI_Request req[48], int procGrpNo )
{
  if( GetDomainType() == CPM_DOMAIN_CARTESIAN )
  {
    cpm_ParaManagerCART *pCART = (cpm_ParaManagerCART*)this;
    return pCART->BndCommS4D_nowait( array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  }
  else if( GetDomainType() == CPM_DOMAIN_LMR )
  {
    cpm_ParaManagerLMR *pLMR = (cpm_ParaManagerLMR*)this;
    return pLMR->BndCommS4D_nowait( array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  }

  return CPM_ERROR_BNDCOMM;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信のwait、展開(Scalar3D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::wait_BndCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                , MPI_Request req[48], int procGrpNo )
{
  return wait_BndCommS4D( array, imax, jmax, kmax, 1, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信のwait、展開(Vector3D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::wait_BndCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                , MPI_Request req[48], int procGrpNo )
{
  return wait_BndCommS4D( array, imax, jmax, kmax, 3, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信のwait、展開(Scalar4D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::wait_BndCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                , MPI_Request req[48], int procGrpNo )
{
  if( GetDomainType() == CPM_DOMAIN_CARTESIAN )
  {
    cpm_ParaManagerCART *pCART = (cpm_ParaManagerCART*)this;
    return pCART->wait_BndCommS4D( array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  }
  else if( GetDomainType() == CPM_DOMAIN_LMR )
  {
    cpm_ParaManagerLMR *pLMR = (cpm_ParaManagerLMR*)this;
    return pLMR->wait_BndCommS4D( array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  }

  return CPM_ERROR_BNDCOMM;
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar3D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::PeriodicCommS3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  return PeriodicCommS4D( array, imax, jmax, kmax, 1, vc, vc_comm, dir, pm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Vector3D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::PeriodicCommV3D( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  return PeriodicCommS4D( array, imax, jmax, kmax, 3, vc, vc_comm, dir, pm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar4D版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::PeriodicCommS4D( T *array, int imax, int jmax, int kmax, int nmax, int vc, int vc_comm
                                , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  if( GetDomainType() == CPM_DOMAIN_CARTESIAN )
  {
    cpm_ParaManagerCART *pCART = (cpm_ParaManagerCART*)this;
    return pCART->PeriodicCommS4D( array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  }
  else if( GetDomainType() == CPM_DOMAIN_LMR )
  {
    cpm_ParaManagerLMR *pLMR = (cpm_ParaManagerLMR*)this;
    return pLMR->PeriodicCommS4D( array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  }

  return CPM_ERROR_BNDCOMM;
}

#endif /* _CPM_PARAMANAGER_BNDCOMM_H_ */

