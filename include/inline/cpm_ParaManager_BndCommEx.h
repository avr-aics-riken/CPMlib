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
 * @file   cpm_ParaManager_BndCommEx.h
 * パラレルマネージャクラスのインラインヘッダーファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */

#ifndef _CPM_PARAMANAGER_BNDCOMMEX_H_
#define _CPM_PARAMANAGER_BNDCOMMEX_H_

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Vector3DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm, int procGrpNo )
{
  return BndCommS4DEx( array, 3, imax, jmax, kmax, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm, int procGrpNo )
{
  if( GetDomainType() == CPM_DOMAIN_CARTESIAN )
  {
    cpm_ParaManagerCART *pCART = (cpm_ParaManagerCART*)this;
    return pCART->BndCommS4DEx( array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  }
  else if( GetDomainType() == CPM_DOMAIN_LMR )
  {
    cpm_ParaManagerLMR *pLMR = (cpm_ParaManagerLMR*)this;
    return pLMR->BndCommS4DEx( array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  }

  return CPM_ERROR_BNDCOMM;
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Vector3DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommV3DEx_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                    , MPI_Request req[48], int procGrpNo )
{
  return BndCommS4DEx_nowait( array, 3, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Scalar4DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS4DEx_nowait( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                    , MPI_Request req[48], int procGrpNo )
{
  if( GetDomainType() == CPM_DOMAIN_CARTESIAN )
  {
    cpm_ParaManagerCART *pCART = (cpm_ParaManagerCART*)this;
    return pCART->BndCommS4DEx_nowait<T>( array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  }
  else if( GetDomainType() == CPM_DOMAIN_LMR )
  {
    cpm_ParaManagerLMR *pLMR = (cpm_ParaManagerLMR*)this;
    return pLMR->BndCommS4DEx_nowait<T>( array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  }

  return CPM_ERROR_BNDCOMM;
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信のwait、展開(Vector3DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::wait_BndCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                  , MPI_Request req[48], int procGrpNo )
{
  return wait_BndCommS4DEx( array, 3, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信のwait、展開(Scalar4DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::wait_BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                  , MPI_Request req[48], int procGrpNo )
{
  if( GetDomainType() == CPM_DOMAIN_CARTESIAN )
  {
    cpm_ParaManagerCART *pCART = (cpm_ParaManagerCART*)this;
    return pCART->wait_BndCommS4DEx<T>( array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  }
  else if( GetDomainType() == CPM_DOMAIN_LMR )
  {
    cpm_ParaManagerLMR *pLMR = (cpm_ParaManagerLMR*)this;
    return pLMR->wait_BndCommS4DEx<T>( array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  }

  return CPM_ERROR_BNDCOMM;
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Vector3DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::PeriodicCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                  , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  return PeriodicCommS4DEx( array, 3, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar4DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::PeriodicCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                  , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  if( GetDomainType() == CPM_DOMAIN_CARTESIAN )
  {
    cpm_ParaManagerCART *pCART = (cpm_ParaManagerCART*)this;
    return pCART->PeriodicCommS4DEx<T>( array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  }
  else if( GetDomainType() == CPM_DOMAIN_LMR )
  {
    cpm_ParaManagerLMR *pLMR = (cpm_ParaManagerLMR*)this;
    return pLMR->PeriodicCommS4DEx<T>( array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  }

  return CPM_ERROR_PERIODIC;
}

#endif /* _CPM_PARAMANAGER_BNDCOMMEX_H_ */

