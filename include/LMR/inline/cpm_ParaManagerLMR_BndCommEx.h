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
 * @file   cpm_ParaManager_BndCommEx_LMR.h
 * LMR用パラレルマネージャクラスのインラインヘッダーファイル
 * @date   2015/03/27
 */

#ifndef _CPM_PARAMANAGER_BNDCOMMEX_LMR_H_
#define _CPM_PARAMANAGER_BNDCOMMEX_LMR_H_

#define _IDXFX(_N,_I,_J,_K,_NN,_IS,_NJ,_NK,_VC) \
( size_t(_NN) \
* ( size_t(_K+_VC) * size_t(_VC) * size_t(_NJ+2*_VC) \
  + size_t(_J+_VC) * size_t(_VC) \
  + size_t(_I-(_IS)) \
  ) \
+ size_t(_N) \
)

#define _IDXFY(_N,_I,_J,_K,_NN,_NI,_JS,_NK,_VC) \
( size_t(_NN) \
* ( size_t(_K+_VC)   * size_t(_NI+2*_VC) * size_t(_VC) \
  + size_t(_J-(_JS)) * size_t(_NI+2*_VC) \
  + size_t(_I+_VC) \
  ) \
+ size_t(_N) \
)

#define _IDXFZ(_N,_I,_J,_K,_NN,_NI,_NJ,_KS,_VC) \
( size_t(_NN) \
* ( size_t(_K-(_KS)) * size_t(_NI+2*_VC) * size_t(_NJ+2*_VC) \
  + size_t(_J+_VC)   * size_t(_NI+2*_VC) \
  + size_t(_I+_VC) \
  ) \
+ size_t(_N) \
)

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Vector3DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm, int procGrpNo )
{
  return BndCommS4DEx( array, 3, imax, jmax, kmax, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Vector3DEx版、waitなし)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommV3DEx_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                       , int procGrpNo )
{
  return BndCommS4DEx_nowait( array, 3, imax, jmax, kmax, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信のwait、展開(Vector3DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::wait_BndCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                     , int procGrpNo )
{
  return wait_BndCommS4DEx( array, 3, imax, jmax, kmax, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Vector3DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::PeriodicCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                     , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  return PeriodicCommS4DEx( array, 3, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4D版Ex)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm, int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 周期境界フラグ
  bool bPeriodic = false;

  // X方向の通信
  {
    // 面内格子数
    size_t sz_face[2] = {jmax, kmax};

    // 通信マップの取得
    BndCommInfoMap::iterator itM = m_bndCommInfoMapMX.find(procGrpNo);
    if( itM == m_bndCommInfoMapMX.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    BndCommInfoMap::iterator itP = m_bndCommInfoMapPX.find(procGrpNo);
    if( itP == m_bndCommInfoMapPX.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    LeafCommInfoMap &commInfoMapM = itM->second;
    LeafCommInfoMap &commInfoMapP = itP->second;

    // マイナス方向受信処理
    if( (ret = recv_LMR<T>(commInfoMapM, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向受信処理
    if( (ret = recv_LMR<T>(commInfoMapP, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向パックと送信処理
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, X_PLUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // マイナス方向パックと送信処理
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, X_MINUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // ランク内コピー処理
    if( (ret = copy_LMR_Ex(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, commInfoMapP, bPeriodic, X_DIR, BOTH, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // マイナス方向受信待機と展開
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, X_MINUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向受信待機と展開
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, X_PLUS, procGrpNo)) != CPM_SUCCESS )
    {
      return ret;
    }

    // 送信待機
    if( (ret = send_LMR_wait<T>(commInfoMapM)) != CPM_SUCCESS )
    {
      return ret;
    }
    if( (ret = send_LMR_wait<T>(commInfoMapP)) != CPM_SUCCESS )
    {
      return ret;
    }
  } // XDIR


  // Y方向の通信
  {
    // 面内格子数
    size_t sz_face[2] = {imax, kmax};

    // 通信マップの取得
    BndCommInfoMap::iterator itM = m_bndCommInfoMapMY.find(procGrpNo);
    if( itM == m_bndCommInfoMapMY.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    BndCommInfoMap::iterator itP = m_bndCommInfoMapPY.find(procGrpNo);
    if( itP == m_bndCommInfoMapPY.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    LeafCommInfoMap &commInfoMapM = itM->second;
    LeafCommInfoMap &commInfoMapP = itP->second;

    // マイナス方向受信処理
    if( (ret = recv_LMR<T>(commInfoMapM, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向受信処理
    if( (ret = recv_LMR<T>(commInfoMapP, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向パックと送信処理
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, Y_PLUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // マイナス方向パックと送信処理
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, Y_MINUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // ランク内コピー処理
    if( (ret = copy_LMR_Ex(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, commInfoMapP, bPeriodic, Y_DIR, BOTH, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // マイナス方向受信待機と展開
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, Y_MINUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向受信待機と展開
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, Y_PLUS, procGrpNo)) != CPM_SUCCESS )
    {
      return ret;
    }

    // 送信待機
    if( (ret = send_LMR_wait<T>(commInfoMapM)) != CPM_SUCCESS )
    {
      return ret;
    }
    if( (ret = send_LMR_wait<T>(commInfoMapP)) != CPM_SUCCESS )
    {
      return ret;
    }
  } // YDIR

  // Z方向の通信
  {
    // 面内格子数
    size_t sz_face[2] = {imax, jmax};

    // 通信マップの取得
    BndCommInfoMap::iterator itM = m_bndCommInfoMapMZ.find(procGrpNo);
    if( itM == m_bndCommInfoMapMZ.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    BndCommInfoMap::iterator itP = m_bndCommInfoMapPZ.find(procGrpNo);
    if( itP == m_bndCommInfoMapPZ.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    LeafCommInfoMap &commInfoMapM = itM->second;
    LeafCommInfoMap &commInfoMapP = itP->second;

    // マイナス方向受信処理
    if( (ret = recv_LMR<T>(commInfoMapM, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向受信処理
    if( (ret = recv_LMR<T>(commInfoMapP, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向パックと送信処理
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, Z_PLUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // マイナス方向パックと送信処理
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, Z_MINUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // ランク内コピー処理
    if( (ret = copy_LMR_Ex(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, commInfoMapP, bPeriodic, Z_DIR, BOTH, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // マイナス方向受信待機と展開
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, Z_MINUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向受信待機と展開
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, Z_PLUS, procGrpNo)) != CPM_SUCCESS )
    {
      return ret;
    }

    // 送信待機
    if( (ret = send_LMR_wait<T>(commInfoMapM)) != CPM_SUCCESS )
    {
      return ret;
    }
    if( (ret = send_LMR_wait<T>(commInfoMapP)) != CPM_SUCCESS )
    {
      return ret;
    }
  } // ZDIR


  // 正常終了
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx版、waitなし)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::BndCommS4DEx_nowait( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm, int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 周期境界フラグ
  bool bPeriodic = false;

  // X方向の通信
  {
    // 面内格子数
    size_t sz_face[2] = {jmax, kmax};

    // 通信マップの取得
    BndCommInfoMap::iterator itM = m_bndCommInfoMapMX.find(procGrpNo);
    if( itM == m_bndCommInfoMapMX.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    BndCommInfoMap::iterator itP = m_bndCommInfoMapPX.find(procGrpNo);
    if( itP == m_bndCommInfoMapPX.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    LeafCommInfoMap &commInfoMapM = itM->second;
    LeafCommInfoMap &commInfoMapP = itP->second;

    // マイナス方向受信処理
    if( (ret = recv_LMR<T>(commInfoMapM, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向受信処理
    if( (ret = recv_LMR<T>(commInfoMapP, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向パックと送信処理
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, X_PLUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // マイナス方向パックと送信処理
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, X_MINUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // ランク内コピー処理
    if( (ret = copy_LMR_Ex(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, commInfoMapP, bPeriodic, X_DIR, BOTH, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }
  } // XDIR


  // Y方向の通信
  {
    // 面内格子数
    size_t sz_face[2] = {imax, kmax};

    // 通信マップの取得
    BndCommInfoMap::iterator itM = m_bndCommInfoMapMY.find(procGrpNo);
    if( itM == m_bndCommInfoMapMY.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    BndCommInfoMap::iterator itP = m_bndCommInfoMapPY.find(procGrpNo);
    if( itP == m_bndCommInfoMapPY.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    LeafCommInfoMap &commInfoMapM = itM->second;
    LeafCommInfoMap &commInfoMapP = itP->second;

    // マイナス方向受信処理
    if( (ret = recv_LMR<T>(commInfoMapM, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向受信処理
    if( (ret = recv_LMR<T>(commInfoMapP, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向パックと送信処理
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, Y_PLUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // マイナス方向パックと送信処理
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, Y_MINUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // ランク内コピー処理
    if( (ret = copy_LMR_Ex(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, commInfoMapP, bPeriodic, Y_DIR, BOTH, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }
  } // YDIR

  // Z方向の通信
  {
    // 面内格子数
    size_t sz_face[2] = {imax, jmax};

    // 通信マップの取得
    BndCommInfoMap::iterator itM = m_bndCommInfoMapMZ.find(procGrpNo);
    if( itM == m_bndCommInfoMapMZ.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    BndCommInfoMap::iterator itP = m_bndCommInfoMapPZ.find(procGrpNo);
    if( itP == m_bndCommInfoMapPZ.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    LeafCommInfoMap &commInfoMapM = itM->second;
    LeafCommInfoMap &commInfoMapP = itP->second;

    // マイナス方向受信処理
    if( (ret = recv_LMR<T>(commInfoMapM, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向受信処理
    if( (ret = recv_LMR<T>(commInfoMapP, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向パックと送信処理
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, Z_PLUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // マイナス方向パックと送信処理
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, Z_MINUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // ランク内コピー処理
    if( (ret = copy_LMR_Ex(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, commInfoMapP, bPeriodic, Z_DIR, BOTH, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }
  }

  // 正常終了
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信のwait、展開(Scalar4DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::wait_BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm, int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 周期境界フラグ
  bool bPeriodic = false;

  // X方向の通信
  {
    // 面内格子数
    size_t sz_face[2] = {jmax, kmax};

    // 通信マップの取得
    BndCommInfoMap::iterator itM = m_bndCommInfoMapMX.find(procGrpNo);
    if( itM == m_bndCommInfoMapMX.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    BndCommInfoMap::iterator itP = m_bndCommInfoMapPX.find(procGrpNo);
    if( itP == m_bndCommInfoMapPX.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    LeafCommInfoMap &commInfoMapM = itM->second;
    LeafCommInfoMap &commInfoMapP = itP->second;

    // マイナス方向受信待機と展開
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, X_MINUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向受信待機と展開
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, X_PLUS, procGrpNo)) != CPM_SUCCESS )
    {
      return ret;
    }

    // 送信待機
    if( (ret = send_LMR_wait<T>(commInfoMapM)) != CPM_SUCCESS )
    {
      return ret;
    }
    if( (ret = send_LMR_wait<T>(commInfoMapP)) != CPM_SUCCESS )
    {
      return ret;
    }
  } // XDIR


  // Y方向の通信
  {
    // 面内格子数
    size_t sz_face[2] = {imax, kmax};

    // 通信マップの取得
    BndCommInfoMap::iterator itM = m_bndCommInfoMapMY.find(procGrpNo);
    if( itM == m_bndCommInfoMapMY.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    BndCommInfoMap::iterator itP = m_bndCommInfoMapPY.find(procGrpNo);
    if( itP == m_bndCommInfoMapPY.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    LeafCommInfoMap &commInfoMapM = itM->second;
    LeafCommInfoMap &commInfoMapP = itP->second;

    // マイナス方向受信待機と展開
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, Y_MINUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向受信待機と展開
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, Y_PLUS, procGrpNo)) != CPM_SUCCESS )
    {
      return ret;
    }

    // 送信待機
    if( (ret = send_LMR_wait<T>(commInfoMapM)) != CPM_SUCCESS )
    {
      return ret;
    }
    if( (ret = send_LMR_wait<T>(commInfoMapP)) != CPM_SUCCESS )
    {
      return ret;
    }
  } // YDIR

  // Z方向の通信
  {
    // 面内格子数
    size_t sz_face[2] = {imax, jmax};

    // 通信マップの取得
    BndCommInfoMap::iterator itM = m_bndCommInfoMapMZ.find(procGrpNo);
    if( itM == m_bndCommInfoMapMZ.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    BndCommInfoMap::iterator itP = m_bndCommInfoMapPZ.find(procGrpNo);
    if( itP == m_bndCommInfoMapPZ.end() )
    {
      return CPM_ERROR_BNDCOMM_BUFFER;
    }
    LeafCommInfoMap &commInfoMapM = itM->second;
    LeafCommInfoMap &commInfoMapP = itP->second;

    // マイナス方向受信待機と展開
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, Z_MINUS, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }

    // プラス方向受信待機と展開
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, Z_PLUS, procGrpNo)) != CPM_SUCCESS )
    {
      return ret;
    }

    // 送信待機
    if( (ret = send_LMR_wait<T>(commInfoMapM)) != CPM_SUCCESS )
    {
      return ret;
    }
    if( (ret = send_LMR_wait<T>(commInfoMapP)) != CPM_SUCCESS )
    {
      return ret;
    }
  } // ZDIR


  // 正常終了
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar4DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::PeriodicCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                     , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 周期境界フラグ
  bool bPeriodic = true;

  // 各種フラグ
  size_t sz_face[2];
  cpm_FaceFlag facePlus, faceMinus;
  BndCommInfoMap* pBndCommInfoMapM;
  BndCommInfoMap* pBndCommInfoMapP;
  if( dir==X_DIR )
  {
    sz_face[0] = jmax;
    sz_face[1] = kmax;
    facePlus  = X_PLUS;
    faceMinus = X_MINUS;
    pBndCommInfoMapM = &m_bndCommInfoMapMX;
    pBndCommInfoMapP = &m_bndCommInfoMapPX;
  }
  else if( dir==Y_DIR )
  {
    sz_face[0] = imax;
    sz_face[1] = kmax;
    facePlus  = Y_PLUS;
    faceMinus = Y_MINUS;
    pBndCommInfoMapM = &m_bndCommInfoMapMY;
    pBndCommInfoMapP = &m_bndCommInfoMapPY;
  }
  else if( dir==Z_DIR )
  {
    sz_face[0] = imax;
    sz_face[1] = jmax;
    facePlus  = Z_PLUS;
    faceMinus = Z_MINUS;
    pBndCommInfoMapM = &m_bndCommInfoMapMZ;
    pBndCommInfoMapP = &m_bndCommInfoMapPZ;
  }

  // 通信マップの取得
  BndCommInfoMap::iterator itM = pBndCommInfoMapM->find(procGrpNo);
  if( itM == pBndCommInfoMapM->end() )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }
  BndCommInfoMap::iterator itP = pBndCommInfoMapP->find(procGrpNo);
  if( itP == pBndCommInfoMapP->end() )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }
  LeafCommInfoMap &commInfoMapM = itM->second;
  LeafCommInfoMap &commInfoMapP = itP->second;

  // マイナス方向受信処理
  if( pm==PLUS2MINUS || pm==BOTH )
  {
    if( (ret = recv_LMR<T>(commInfoMapM, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }
  }

  // プラス方向受信処理
  if( pm==MINUS2PLUS || pm==BOTH )
  {
    if( (ret = recv_LMR<T>(commInfoMapP, sz_face, nmax, vc_comm, bPeriodic, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }
  }

  // プラス方向パックと送信処理
  if( pm==PLUS2MINUS || pm==BOTH )
  {
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, facePlus, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }
  }

  // マイナス方向パックと送信処理
  if( pm==MINUS2PLUS || pm==BOTH )
  {
    if( (ret = send_LMR_Ex( array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, faceMinus, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }
  }

  // ランク内コピー処理
  if( (ret = copy_LMR_Ex(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, commInfoMapP, bPeriodic, dir, pm, procGrpNo)) != CPM_SUCCESS )
  {
     return ret;
  }

  // マイナス方向受信待機と展開
  if( pm==PLUS2MINUS || pm==BOTH )
  {
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapM, bPeriodic, faceMinus, procGrpNo)) != CPM_SUCCESS )
    {
       return ret;
    }
  }

  // プラス方向受信待機と展開
  if( pm==MINUS2PLUS || pm==BOTH )
  {
    if( (ret = recv_LMR_Ex_wait(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoMapP, bPeriodic, facePlus, procGrpNo)) != CPM_SUCCESS )
    {
      return ret;
    }
  }

  // 送信待機
  if( pm==MINUS2PLUS || pm==BOTH )
  {
    if( (ret = send_LMR_wait<T>(commInfoMapM)) != CPM_SUCCESS )
    {
      return ret;
    }
  }
  if( pm==PLUS2MINUS || pm==BOTH )
  {
    if( (ret = send_LMR_wait<T>(commInfoMapP)) != CPM_SUCCESS )
    {
      return ret;
    }
  }

  // 正常終了
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の１面の送信データのパックと送信
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::send_LMR_Ex( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                               , LeafCommInfoMap &commInfoMap, bool bPeriodic, cpm_FaceFlag face
                               , int procGrpNo)
{
  cpm_ErrorCode ret;

  size_t sz_face[2];
  if( face==X_MINUS || face==X_PLUS )
  {
    sz_face[0] = jmax;
    sz_face[1] = kmax;
  }
  else if( face==Y_MINUS || face==Y_PLUS )
  {
    sz_face[0] = imax;
    sz_face[1] = kmax;
  }
  else if( face==Z_MINUS || face==Z_PLUS )
  {
    sz_face[0] = imax;
    sz_face[1] = jmax;
  }

  for( LeafCommInfoMap::iterator it=commInfoMap.begin();it!=commInfoMap.end();it++ )
  {
    // ランク間通信情報を取得
    int distRank = it->first;
    cpm_LeafCommInfo* pLeafCommInfo = it->second;

    // リクエストをNULLにしておく
    pLeafCommInfo->m_reqSend = MPI_REQUEST_NULL;

    if( m_rankNo == distRank )
    {
      // 同一ランク内なので通信しない
      continue;
    }

    // 送信バッファ
    T* sendbuf = (T*)pLeafCommInfo->GetBndCommSendBufferPtr();

    // 送信サイズの計算とパック
    int commsize = 0;
    T *ptr = sendbuf;
    for( int j=0;j<pLeafCommInfo->m_vecCommInfo.size();j++ )
    {
      cpm_LeafCommInfo::stCommInfo* commInfo = pLeafCommInfo->m_vecCommInfo[j];
      if( commInfo->bPeriodic == bPeriodic )
      {
        // リーフ間の送信サイズ
        size_t csz = commInfo->CalcSendBufferSize(sz_face, vc_comm, nmax);

        // パック
        if( face == X_MINUS )
        {
          if( (ret = packMXEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfo, ptr, csz, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }
        else if( face == X_PLUS )
        {
          if( (ret = packPXEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfo, ptr, csz, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }
        else if( face == Y_MINUS )
        {
          if( (ret = packMYEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfo, ptr, csz, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }
        else if( face == Y_PLUS )
        {
          if( (ret = packPYEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfo, ptr, csz, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }
        else if( face == Z_MINUS )
        {
          if( (ret = packMZEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfo, ptr, csz, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }
        else if( face == Z_PLUS )
        {
          if( (ret = packPZEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfo, ptr, csz, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }

        // 送信サイズを加算
        commsize += (int)csz;

        // バッファのポインタも変更
        ptr += csz;
      }
    }

    // 送信
    if( commsize > 0 )
    {
      if( (ret = Isend( sendbuf, commsize, distRank, &pLeafCommInfo->m_reqSend, procGrpNo )) != CPM_SUCCESS )
      {
        return ret;
      }
//std::cout << "[" << m_rankNo << "] 送信 : dist=" << distRank << " sz=" << commsize << std::endl;
    }
  }

  // 正常終了
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)のランク内コピー処理
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::copy_LMR_Ex( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                               , LeafCommInfoMap &commInfoMapM, LeafCommInfoMap &commInfoMapP
                               , bool bPeriodic, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  cpm_ErrorCode ret;

  size_t sz_face[2];
  if( dir==X_DIR )
  {
    sz_face[0] = jmax;
    sz_face[1] = kmax;
  }
  else if( dir==Y_DIR )
  {
    sz_face[0] = imax;
    sz_face[1] = kmax;
  }
  else if( dir==Z_DIR )
  {
    sz_face[0] = imax;
    sz_face[1] = jmax;
  }

  // ランク内の通信情報を取得
  LeafCommInfoMap::iterator itM = commInfoMapM.find(m_rankNo);
  if( itM == commInfoMapM.end() )
  {
    // ランク内通信は存在しないこともあるため正常終了
    return CPM_SUCCESS;
  }
  LeafCommInfoMap::iterator itP = commInfoMapP.find(m_rankNo);
  if( itP == commInfoMapP.end() )
  {
    // ランク内通信は存在しないこともあるため正常終了
    return CPM_SUCCESS;
  }
  cpm_LeafCommInfo* pLeafCommInfoM = itM->second;
  cpm_LeafCommInfo* pLeafCommInfoP = itP->second;

  // バッファ
  T* copybuf = (T*)pLeafCommInfoM->GetBndCommSendBufferPtr();

  // マイナス側の通信情報を検索
  for( int j=0;j<pLeafCommInfoM->m_vecCommInfo.size();j++ )
  {
    cpm_LeafCommInfo::stCommInfo* commInfoM = pLeafCommInfoM->m_vecCommInfo[j];
    if( commInfoM->bPeriodic == bPeriodic )
    {
      // リーフ間の送信サイズ
      size_t csz = commInfoM->CalcSendBufferSize(sz_face, vc_comm, nmax);

      // 対応するプラス側の通信情報を検索
      cpm_LeafCommInfo::stCommInfo* commInfoP = pLeafCommInfoP->SearchDistCommInfo(commInfoM);
      if( !commInfoP )
      {
        return CPM_ERROR;
      }

      // パックと展開
      if( dir == X_DIR )
      {
        // マイナス側からプラス側へのコピー
        if( pm==MINUS2PLUS || pm==BOTH )
        {
          // マイナス側のパック
          if( (ret = packMXEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoM, copybuf, csz, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
          // プラス側に展開
          if( (ret = unpackPXEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoP, copybuf, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }

        // プラス側からマイナス側へのコピー
        if( pm==PLUS2MINUS || pm==BOTH )
        {
          // プラス側のパック
          if( (ret = packPXEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoP, copybuf, csz, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
          // マイナス側に展開
          if( (ret = unpackMXEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoM, copybuf, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }
      }
      else if( dir == Y_DIR )
      {
        // マイナス側からプラス側へのコピー
        if( pm==MINUS2PLUS || pm==BOTH )
        {
          // マイナス側のパック
          if( (ret = packMYEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoM, copybuf, csz, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
          // プラス側に展開
          if( (ret = unpackPYEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoP, copybuf, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }

        // プラス側からマイナス側へのコピー
        if( pm==PLUS2MINUS || pm==BOTH )
        {
          // プラス側のパック
          if( (ret = packPYEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoP, copybuf, csz, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
          // マイナス側に展開
          if( (ret = unpackMYEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoM, copybuf, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }
      }
      else if( dir == Z_DIR )
      {
        // マイナス側からプラス側へのコピー
        if( pm==MINUS2PLUS || pm==BOTH )
        {
          // マイナス側のパック
          if( (ret = packMZEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoM, copybuf, csz, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
          // プラス側に展開
          if( (ret = unpackPZEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoP, copybuf, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }

        // プラス側からマイナス側へのコピー
        if( pm==PLUS2MINUS || pm==BOTH )
        {
          // プラス側のパック
          if( (ret = packPZEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoP, copybuf, csz, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
          // マイナス側に展開
          if( (ret = unpackMZEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfoM, copybuf, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }
      }
    }
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の１面の受信待機とデータの展開
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::recv_LMR_Ex_wait( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                    , LeafCommInfoMap &commInfoMap, bool bPeriodic, cpm_FaceFlag face
                                    , int procGrpNo)
{
  cpm_ErrorCode ret;

  size_t sz_face[2];
  if( face==X_MINUS || face==X_PLUS )
  {
    sz_face[0] = jmax;
    sz_face[1] = kmax;
  }
  else if( face==Y_MINUS || face==Y_PLUS )
  {
    sz_face[0] = imax;
    sz_face[1] = kmax;
  }
  else if( face==Z_MINUS || face==Z_PLUS )
  {
    sz_face[0] = imax;
    sz_face[1] = jmax;
  }

  for( LeafCommInfoMap::iterator it=commInfoMap.begin();it!=commInfoMap.end();it++ )
  {
    // ランク間通信情報を取得
    int distRank = it->first;
    cpm_LeafCommInfo* pLeafCommInfo = it->second;

    // 同一ランク内のとき、通信しない
    if( m_rankNo == distRank )
    {
      continue;
    }

    // リクエストNULLのとき何もしない
    if( pLeafCommInfo->m_reqRecv == MPI_REQUEST_NULL )
    {
      continue;
    }

    // Wait
    if( (ret = Wait( &pLeafCommInfo->m_reqRecv )) != CPM_SUCCESS )
    {
      return ret;
    }
    pLeafCommInfo->m_reqRecv = MPI_REQUEST_NULL;

    // 受信バッファ
    T* recvbuf = (T*)pLeafCommInfo->GetBndCommRecvBufferPtr();

    // 受信データの展開
    T *ptr = recvbuf;
    for( int j=0;j<pLeafCommInfo->m_vecCommInfo.size();j++ )
    {
      cpm_LeafCommInfo::stCommInfo* commInfo = pLeafCommInfo->m_vecCommInfo[j];
      if( commInfo->bPeriodic == bPeriodic )
      {
        // リーフ間の受信サイズ
        size_t csz = commInfo->CalcRecvBufferSize(sz_face, vc_comm, nmax);

        // 展開
        if( face == X_MINUS )
        {
          if( (ret = unpackMXEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfo, ptr, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }
        else if( face == X_PLUS )
        {
          if( (ret = unpackPXEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfo, ptr, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }
        else if( face == Y_MINUS )
        {
          if( (ret = unpackMYEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfo, ptr, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }
        else if( face == Y_PLUS )
        {
          if( (ret = unpackPYEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfo, ptr, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }
        if( face == Z_MINUS )
        {
          if( (ret = unpackMZEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfo, ptr, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }
        else if( face == Z_PLUS )
        {
          if( (ret = unpackPZEx(array, nmax, imax, jmax, kmax, vc, vc_comm, commInfo, ptr, procGrpNo)) != CPM_SUCCESS )
          {
            return ret;
          }
        }

        // バッファのポインタも変更
        ptr += csz;
      }
    }
  }

  // 正常終了
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の-X面への送信データのパック(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::packMXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                            , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw, int procGrpNo )
{
  // レベル差
  int levelDiff = commInfo->iLevelDiff;

  // face
  int faceNo = commInfo->iFaceIdx;

  // リーフインデクス
  int leafIdx = GetLocalLeafIndex_byID(commInfo->iOwnLeafID, procGrpNo);

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
  for( int k=ks-gc,kk=0-gc;k<ke+gc;k++,kk++ ){
  for( int j=js-gc,jj=0-gc;j<je+gc;j++,jj++ ){
  for( int i=0;i<gc;i++ ){
  for( int n=0;n<nmax;n++){
    sendbuf[_IDXFX(n,i,jj,kk,nmax,0,jmaxb,kmaxb,gc)] = array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)];
  }}}}

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の+X面への送信データのパック(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::packPXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                            , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw, int procGrpNo )
{
  // レベル差
  int levelDiff = commInfo->iLevelDiff;

  // face
  int faceNo = commInfo->iFaceIdx;

  // リーフインデクス
  int leafIdx = GetLocalLeafIndex_byID(commInfo->iOwnLeafID, procGrpNo);

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
  for( int k=ks-gc,kk=0-gc;k<ke+gc;k++,kk++ ){
  for( int j=js-gc,jj=0-gc;j<je+gc;j++,jj++ ){
  for( int i=imax-gc;i<imax;i++ ){
  for( int n=0;n<nmax;n++){
    sendbuf[_IDXFX(n,i,jj,kk,nmax,imax-gc,jmaxb,kmaxb,gc)] = array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)];
  }}}}

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の-Y面への送信データのパック(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::packMYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                            , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw, int procGrpNo )
{
  // レベル差
  int levelDiff = commInfo->iLevelDiff;

  // face
  int faceNo = commInfo->iFaceIdx;

  // リーフインデクス
  int leafIdx = GetLocalLeafIndex_byID(commInfo->iOwnLeafID, procGrpNo);

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
  for( int k=ks-gc,kk=0-gc;k<ke+gc;k++,kk++ ){
  for( int j=0;j<gc;j++ ){
  for( int i=is-gc,ii=0-gc;i<ie+gc;i++,ii++ ){
  for( int n=0;n<nmax;n++){
    sendbuf[_IDXFY(n,ii,j,kk,nmax,imaxb,0,kmaxb,gc)] = array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)];
  }}}}

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の+Y面への送信データのパック(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::packPYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                            , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw, int procGrpNo )
{
  // レベル差
  int levelDiff = commInfo->iLevelDiff;

  // face
  int faceNo = commInfo->iFaceIdx;

  // リーフインデクス
  int leafIdx = GetLocalLeafIndex_byID(commInfo->iOwnLeafID, procGrpNo);

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
  for( int k=ks-gc,kk=0-gc;k<ke+gc;k++,kk++ ){
  for( int j=jmax-gc;j<jmax;j++ ){
  for( int i=is-gc,ii=0-gc;i<ie+gc;i++,ii++ ){
  for( int n=0;n<nmax;n++){
    sendbuf[_IDXFY(n,ii,j,kk,nmax,imaxb,jmax-gc,kmaxb,gc)] = array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)];
  }}}}

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の-Z面への送信データのパック(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::packMZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                            , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw, int procGrpNo )
{
  // レベル差
  int levelDiff = commInfo->iLevelDiff;

  // face
  int faceNo = commInfo->iFaceIdx;

  // リーフインデクス
  int leafIdx = GetLocalLeafIndex_byID(commInfo->iOwnLeafID, procGrpNo);

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
  for( int k=0;k<gc;k++ ){
  for( int j=js-gc,jj=0-gc;j<je+gc;j++,jj++ ){
  for( int i=is-gc,ii=0-gc;i<ie+gc;i++,ii++ ){
  for( int n=0;n<nmax;n++){
    sendbuf[_IDXFZ(n,ii,jj,k,nmax,imaxb,jmaxb,0,gc)] = array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)];
  }}}}

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の+Z面への送信データのパック(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::packPZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                            , cpm_LeafCommInfo::stCommInfo* commInfo, T* sendbuf, size_t nw, int procGrpNo )
{
  // レベル差
  int levelDiff = commInfo->iLevelDiff;

  // face
  int faceNo = commInfo->iFaceIdx;

  // リーフインデクス
  int leafIdx = GetLocalLeafIndex_byID(commInfo->iOwnLeafID, procGrpNo);

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
  for( int k=kmax-gc;k<kmax;k++ ){
  for( int j=js-gc,jj=0-gc;j<je+gc;j++,jj++ ){
  for( int i=is-gc,ii=0-gc;i<ie+gc;i++,ii++ ){
  for( int n=0;n<nmax;n++){
    sendbuf[_IDXFZ(n,ii,jj,k,nmax,imaxb,jmaxb,kmax-gc,gc)] = array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)];
  }}}}

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の-X面からの受信データの展開(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::unpackMXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                              , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo )
{
  // レベル差
  int levelDiff = commInfo->iLevelDiff;

  // face
  int faceNo = commInfo->iFaceIdx;

  // リーフインデクス
  int leafIdx = GetLocalLeafIndex_byID(commInfo->iOwnLeafID, procGrpNo);

  // 展開
  if( levelDiff==0 )
  {
    // 同じレベル
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<0;i++ ){
    for( int n=0;n<nmax;n++){
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = recvbuf[_IDXFX(n,i,j,k,nmax,0-vc_comm,jmax,kmax,vc_comm)];
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

    for( int k=ks-vc_km;k<ke+vc_kp;k++ ){
    for( int j=js-vc_jm;j<je+vc_jp;j++ ){
    for( int i=0-vc_comm;i<0;i++ ){
    for( int n=0;n<nmax;n++){
      int ii=(i   )*2;
      int jj=(j-js)*2;
      int kk=(k-ks)*2;
      T val = recvbuf[_IDXFX(n,ii  ,jj  ,kk  ,nmax,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii+1,jj  ,kk  ,nmax,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii  ,jj+1,kk  ,nmax,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii+1,jj+1,kk  ,nmax,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii  ,jj  ,kk+1,nmax,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii+1,jj  ,kk+1,nmax,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii  ,jj+1,kk+1,nmax,0-gc,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii+1,jj+1,kk+1,nmax,0-gc,jmax,kmax,gc)];
      val *= 0.125;
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = val;
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

    for( int k=0-gc;k<kmax+gc;k++ ){
    for( int j=0-gc;j<jmax+gc;j++ ){
    for( int i=0-gc;i<0;i++ ){
    for( int n=0;n<nmax;n++){
      int ii=(i+2)/2-1;
      int jj=(j+2)/2-1;
      int kk=(k+2)/2-1;
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = recvbuf[_IDXFX(n,ii,jj,kk,nmax,0-vc_comm,jmaxb,kmaxb,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の+X面からの受信データの展開(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::unpackPXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                              , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo )
{
  // レベル差
  int levelDiff = commInfo->iLevelDiff;

  // face
  int faceNo = commInfo->iFaceIdx;

  // リーフインデクス
  int leafIdx = GetLocalLeafIndex_byID(commInfo->iOwnLeafID, procGrpNo);

  // 展開
  if( levelDiff==0 )
  {
    // 同じレベル
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=imax;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = recvbuf[_IDXFX(n,i,j,k,nmax,imax,jmax,kmax,vc_comm)];
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

    for( int k=ks-vc_km;k<ke+vc_kp;k++ ){
    for( int j=js-vc_jm;j<je+vc_jp;j++ ){
    for( int i=imax;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      int ii=(i   )*2;
      int jj=(j-js)*2;
      int kk=(k-ks)*2;
      T val = recvbuf[_IDXFX(n,ii  ,jj  ,kk  ,nmax,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii+1,jj  ,kk  ,nmax,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii  ,jj+1,kk  ,nmax,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii+1,jj+1,kk  ,nmax,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii  ,jj  ,kk+1,nmax,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii+1,jj  ,kk+1,nmax,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii  ,jj+1,kk+1,nmax,imax2,jmax,kmax,gc)]
            + recvbuf[_IDXFX(n,ii+1,jj+1,kk+1,nmax,imax2,jmax,kmax,gc)];
      val *= 0.125;
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = val;
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

    for( int k=0-gc;k<kmax+gc;k++ ){
    for( int j=0-gc;j<jmax+gc;j++ ){
    for( int i=imax;i<imax+gc;i++ ){
    for( int n=0;n<nmax;n++){
      int ii=i/2;
      int jj=(j+2)/2-1;
      int kk=(k+2)/2-1;
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = recvbuf[_IDXFX(n,ii,jj,kk,nmax,imax2,jmaxb,kmaxb,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の-Y面からの受信データの展開(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::unpackMYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                              , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo )
{
  // レベル差
  int levelDiff = commInfo->iLevelDiff;

  // face
  int faceNo = commInfo->iFaceIdx;

  // リーフインデクス
  int leafIdx = GetLocalLeafIndex_byID(commInfo->iOwnLeafID, procGrpNo);

  // 展開
  if( levelDiff==0 )
  {
    // 同じレベル
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<0;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = recvbuf[_IDXFY(n,i,j,k,nmax,imax,0-vc_comm,kmax,vc_comm)];
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

    for( int k=ks-vc_km;k<ke+vc_kp;k++ ){
    for( int j=0-vc_comm;j<0;j++ ){
    for( int i=is-vc_im;i<ie+vc_ip;i++ ){
    for( int n=0;n<nmax;n++){
      int ii=(i-is)*2;
      int jj=(j   )*2;
      int kk=(k-ks)*2;
      T val = recvbuf[_IDXFY(n,ii  ,jj  ,kk  ,nmax,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(n,ii+1,jj  ,kk  ,nmax,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(n,ii  ,jj+1,kk  ,nmax,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(n,ii+1,jj+1,kk  ,nmax,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(n,ii  ,jj  ,kk+1,nmax,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(n,ii+1,jj  ,kk+1,nmax,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(n,ii  ,jj+1,kk+1,nmax,imax,0-gc,kmax,gc)]
            + recvbuf[_IDXFY(n,ii+1,jj+1,kk+1,nmax,imax,0-gc,kmax,gc)];
      val *= 0.125;
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = val;
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

    for( int k=0-gc;k<kmax+gc;k++ ){
    for( int j=0-gc;j<0;j++ ){
    for( int i=0-gc;i<imax+gc;i++ ){
    for( int n=0;n<nmax;n++){
      int ii=(i+2)/2-1;
      int jj=(j+2)/2-1;
      int kk=(k+2)/2-1;
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = recvbuf[_IDXFY(n,ii,jj,kk,nmax,imaxb,0-vc_comm,kmaxb,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の+Y面からの受信データの展開(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::unpackPYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                              , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo )
{
  // レベル差
  int levelDiff = commInfo->iLevelDiff;

  // face
  int faceNo = commInfo->iFaceIdx;

  // リーフインデクス
  int leafIdx = GetLocalLeafIndex_byID(commInfo->iOwnLeafID, procGrpNo);

  // 展開
  if( levelDiff==0 )
  {
    // 同じレベル
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=jmax;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = recvbuf[_IDXFY(n,i,j,k,nmax,imax,jmax,kmax,vc_comm)];
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

    for( int k=ks-vc_km;k<ke+vc_kp;k++ ){
    for( int j=jmax;j<jmax+vc_comm;j++ ){
    for( int i=is-vc_im;i<ie+vc_ip;i++ ){
    for( int n=0;n<nmax;n++){
      int ii=(i-is)*2;
      int jj=(j   )*2;
      int kk=(k-ks)*2;
      T val = recvbuf[_IDXFY(n,ii  ,jj  ,kk  ,nmax,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(n,ii+1,jj  ,kk  ,nmax,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(n,ii  ,jj+1,kk  ,nmax,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(n,ii+1,jj+1,kk  ,nmax,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(n,ii  ,jj  ,kk+1,nmax,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(n,ii+1,jj  ,kk+1,nmax,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(n,ii  ,jj+1,kk+1,nmax,imax,jmax2,kmax,gc)]
            + recvbuf[_IDXFY(n,ii+1,jj+1,kk+1,nmax,imax,jmax2,kmax,gc)];
      val *= 0.125;
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = val;
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

    for( int k=0-gc;k<kmax+gc;k++ ){
    for( int j=jmax;j<jmax+gc;j++ ){
    for( int i=0-gc;i<imax+gc;i++ ){
    for( int n=0;n<nmax;n++){
      int ii=(i+2)/2-1;
      int jj=j/2;
      int kk=(k+2)/2-1;
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = recvbuf[_IDXFY(n,ii,jj,kk,nmax,imaxb,jmax2,kmaxb,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の-Z面からの受信データの展開(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::unpackMZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                              , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo )
{
  // レベル差
  int levelDiff = commInfo->iLevelDiff;

  // face
  int faceNo = commInfo->iFaceIdx;

  // リーフインデクス
  int leafIdx = GetLocalLeafIndex_byID(commInfo->iOwnLeafID, procGrpNo);

  // 展開
  if( levelDiff==0 )
  {
    // 同じレベル
    for( int k=0-vc_comm;k<0;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = recvbuf[_IDXFZ(n,i,j,k,nmax,imax,jmax,0-vc_comm,vc_comm)];
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

    for( int k=0-vc_comm;k<0;k++ ){
    for( int j=js-vc_jm;j<je+vc_jp;j++ ){
    for( int i=is-vc_im;i<ie+vc_ip;i++ ){
    for( int n=0;n<nmax;n++){
      int ii=(i-is)*2;
      int jj=(j-js)*2;
      int kk=(k   )*2;
      T val = recvbuf[_IDXFZ(n,ii  ,jj  ,kk  ,nmax,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(n,ii+1,jj  ,kk  ,nmax,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(n,ii  ,jj+1,kk  ,nmax,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(n,ii+1,jj+1,kk  ,nmax,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(n,ii  ,jj  ,kk+1,nmax,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(n,ii+1,jj  ,kk+1,nmax,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(n,ii  ,jj+1,kk+1,nmax,imax,jmax,0-gc,gc)]
            + recvbuf[_IDXFZ(n,ii+1,jj+1,kk+1,nmax,imax,jmax,0-gc,gc)];
      val *= 0.125;
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = val;
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

    for( int k=0-gc;k<0;k++ ){
    for( int j=0-gc;j<jmax+gc;j++ ){
    for( int i=0-gc;i<imax+gc;i++ ){
    for( int n=0;n<nmax;n++){
      int ii=(i+2)/2-1;
      int jj=(j+2)/2-1;
      int kk=(k+2)/2-1;
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = recvbuf[_IDXFZ(n,ii,jj,kk,nmax,imaxb,jmaxb,0-vc_comm,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)の+Z面からの受信データの展開(通信面毎)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManagerLMR::unpackPZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                              , cpm_LeafCommInfo::stCommInfo* commInfo, T* recvbuf, int procGrpNo )
{
  // レベル差
  int levelDiff = commInfo->iLevelDiff;

  // face
  int faceNo = commInfo->iFaceIdx;

  // リーフインデクス
  int leafIdx = GetLocalLeafIndex_byID(commInfo->iOwnLeafID, procGrpNo);

  // 展開
  if( levelDiff==0 )
  {
    // 同じレベル
    for( int k=kmax;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = recvbuf[_IDXFZ(n,i,j,k,nmax,imax,jmax,kmax,vc_comm)];
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

    for( int k=kmax;k<kmax+vc_comm;k++ ){
    for( int j=js-vc_jm;j<je+vc_jp;j++ ){
    for( int i=is-vc_im;i<ie+vc_ip;i++ ){
    for( int n=0;n<nmax;n++){
      int ii=(i-is)*2;
      int jj=(j-js)*2;
      int kk=(k   )*2;
      T val = recvbuf[_IDXFZ(n,ii  ,jj  ,kk  ,nmax,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(n,ii+1,jj  ,kk  ,nmax,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(n,ii  ,jj+1,kk  ,nmax,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(n,ii+1,jj+1,kk  ,nmax,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(n,ii  ,jj  ,kk+1,nmax,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(n,ii+1,jj  ,kk+1,nmax,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(n,ii  ,jj+1,kk+1,nmax,imax,jmax,kmax2,gc)]
            + recvbuf[_IDXFZ(n,ii+1,jj+1,kk+1,nmax,imax,jmax,kmax2,gc)];
      val *= 0.125;
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = val;
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

    for( int k=kmax;k<kmax+gc;k++ ){
    for( int j=0-gc;j<jmax+gc;j++ ){
    for( int i=0-gc;i<imax+gc;i++ ){
    for( int n=0;n<nmax;n++){
      int ii=(i+2)/2-1;
      int jj=(j+2)/2-1;
      int kk=k/2;
      array[_IDX_S4DEX_LMR(n,i,j,k,leafIdx,nmax,imax,jmax,kmax,vc)] = recvbuf[_IDXFZ(n,ii,jj,kk,nmax,imaxb,jmaxb,kmax2,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

#undef _IDXFX
#undef _IDXFY
#undef _IDXFZ

#endif /* _CPM_PARAMANAGER_BNDCOMMEX_LMR_H_ */
