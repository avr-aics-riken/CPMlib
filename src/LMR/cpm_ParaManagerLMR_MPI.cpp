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
 * @file   cpm_ParaManager_MPI.cpp
 * パラレルマネージャクラスのMPIインターフェイス関数ソースファイル
 * @date   2012/05/31
 */
#include "stdlib.h"
#include "cpm_ParaManagerLMR.h"

#if !defined(_WIN32) && !defined(WIN32)
#include <unistd.h> // for gethostname() of FX10/K
#endif

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                           , int vc, int vc_comm, int procGrpNo ) 
{
  return BndCommS4D( dtype, array, imax, jmax, kmax, 1, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Vector3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                           , int vc, int vc_comm, int procGrpNo ) 
{
  return BndCommS4D( dtype, array, imax, jmax, kmax, 3, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                           , int vc, int vc_comm, int procGrpNo ) 
{
  if( dtype == MPI_CHAR )
    return BndCommS4D( (char*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_SHORT )
    return BndCommS4D( (short*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_INT )
    return BndCommS4D( (int*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_LONG )
    return BndCommS4D( (long*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return BndCommS4D( (float*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return BndCommS4D( (double*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return BndCommS4D( (long double*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return BndCommS4D( (unsigned char*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return BndCommS4D( (unsigned short*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return BndCommS4D( (unsigned*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_LONG )
    return BndCommS4D( (unsigned long*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return BndCommS4D( (long long int*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return BndCommS4D( (long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return BndCommS4D( (unsigned long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Scalar3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommS3D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                   , int vc, int vc_comm, int procGrpNo ) 
{
  return BndCommS4D_nowait( dtype, array, imax, jmax, kmax, 1, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Vector3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommV3D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                   , int vc, int vc_comm, int procGrpNo ) 
{
  return BndCommS4D_nowait( dtype, array, imax, jmax, kmax, 3, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Scalar4D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommS4D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                                  , int vc, int vc_comm, int procGrpNo ) 
{
  if( dtype == MPI_CHAR )
    return BndCommS4D_nowait( (char*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_SHORT )
    return BndCommS4D_nowait( (short*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_INT )
    return BndCommS4D_nowait( (int*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_LONG )
    return BndCommS4D_nowait( (long*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return BndCommS4D_nowait( (float*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return BndCommS4D_nowait( (double*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return BndCommS4D_nowait( (long double*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return BndCommS4D_nowait( (unsigned char*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return BndCommS4D_nowait( (unsigned short*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return BndCommS4D_nowait( (unsigned*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_LONG )
    return BndCommS4D_nowait( (unsigned long*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return BndCommS4D_nowait( (long long int*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return BndCommS4D_nowait( (long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return BndCommS4D_nowait( (unsigned long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::PeriodicCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  return PeriodicCommS4D( dtype, array, imax, jmax, kmax, 1, vc, vc_comm, dir, pm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Vector3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::PeriodicCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  return PeriodicCommS4D( dtype, array, imax, jmax, kmax, 3, vc, vc_comm, dir, pm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar4D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::PeriodicCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                                , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  if( dtype == MPI_CHAR )
    return PeriodicCommS4D( (char*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_SHORT )
    return PeriodicCommS4D( (short*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_INT )
    return PeriodicCommS4D( (int*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_LONG )
    return PeriodicCommS4D( (long*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return PeriodicCommS4D( (float*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return PeriodicCommS4D( (double*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return PeriodicCommS4D( (long double*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return PeriodicCommS4D( (unsigned char*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return PeriodicCommS4D( (unsigned short*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return PeriodicCommS4D( (unsigned*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_LONG )
    return PeriodicCommS4D( (unsigned long*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return PeriodicCommS4D( (long long int*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return PeriodicCommS4D( (long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return PeriodicCommS4D( (unsigned long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Vector3DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                             , int vc, int vc_comm, int procGrpNo )
{
  return BndCommS4DEx( dtype, array, 3, imax, jmax, kmax, vc, vc_comm, procGrpNo );
}
////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                             , int vc, int vc_comm, int procGrpNo )
{
  if( dtype == MPI_CHAR )
    return BndCommS4DEx( (char*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_SHORT )
    return BndCommS4DEx( (short*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_INT )
    return BndCommS4DEx( (int*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_LONG )
    return BndCommS4DEx( (long*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return BndCommS4DEx( (float*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return BndCommS4DEx( (double*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return BndCommS4DEx( (long double*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return BndCommS4DEx( (unsigned char*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return BndCommS4DEx( (unsigned short*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return BndCommS4DEx( (unsigned*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_LONG )
    return BndCommS4DEx( (unsigned long*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return BndCommS4DEx( (long long int*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return BndCommS4DEx( (long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return BndCommS4DEx( (unsigned long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}
  
////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Vector3DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommV3DEx_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                    , int vc, int vc_comm, int procGrpNo ) 
{
  return BndCommS4DEx_nowait( dtype, array, 3, imax, jmax, kmax, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Scalar4DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommS4DEx_nowait( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                    , int vc, int vc_comm, int procGrpNo ) 
{
  if( dtype == MPI_CHAR )
    return BndCommS4DEx_nowait( (char*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_SHORT )
    return BndCommS4DEx_nowait( (short*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_INT )
    return BndCommS4DEx_nowait( (int*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_LONG )
    return BndCommS4DEx_nowait( (long*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return BndCommS4DEx_nowait( (float*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return BndCommS4DEx_nowait( (double*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return BndCommS4DEx_nowait( (long double*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return BndCommS4DEx_nowait( (unsigned char*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return BndCommS4DEx_nowait( (unsigned short*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return BndCommS4DEx_nowait( (unsigned*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_LONG )
    return BndCommS4DEx_nowait( (unsigned long*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return BndCommS4DEx_nowait( (long long int*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return BndCommS4DEx_nowait( (long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return BndCommS4DEx_nowait( (unsigned long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Vector3DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::PeriodicCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                  , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  return PeriodicCommS4DEx( dtype, array, 3, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar4DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManagerLMR::PeriodicCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                  , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  if( dtype == MPI_CHAR )
    return PeriodicCommS4DEx( (char*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_SHORT )
    return PeriodicCommS4DEx( (short*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_INT )
    return PeriodicCommS4DEx( (int*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_LONG )
    return PeriodicCommS4DEx( (long*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return PeriodicCommS4DEx( (float*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return PeriodicCommS4DEx( (double*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return PeriodicCommS4DEx( (long double*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return PeriodicCommS4DEx( (unsigned char*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return PeriodicCommS4DEx( (unsigned short*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return PeriodicCommS4DEx( (unsigned*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
   else if( dtype == MPI_UNSIGNED_LONG )
    return PeriodicCommS4DEx( (unsigned long*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return PeriodicCommS4DEx( (long long int*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return PeriodicCommS4DEx( (long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return PeriodicCommS4DEx( (unsigned long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}

