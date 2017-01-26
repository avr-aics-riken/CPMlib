#ifndef _FFV_LS_H_
#define _FFV_LS_H_
//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   LS.h
 * @brief  LS Class header
 * @author aics
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "omp.h"
#include "FBdefine.h"
#include "kernel_def.h"
#include "FortFunc.h"

using namespace pm_lib;

class LinearSolver {
  
private:
  int size[4];
  int guide;
  int IterationMax;
  REAL_TYPE dh;
  REAL_TYPE coef_ac;
  double eps;
  
  PerfMonitor* PM;   ///< PerfMonitor class
  
  FILE* fp;

  int* vec_tag;
  int* vec_itr;
  int  vec_cnt;

public:
  
  /** デフォルトコンストラクタ */
  LinearSolver() {
    size[0] = size[1] = size[2] = size[3] = 0;
    guide = 0;
    IterationMax = 0;
    dh = 0.0;
    coef_ac = 0.0;
    eps = 0.0;
    
    PM  = NULL;
    fp = NULL;

    vec_tag=NULL;
    vec_itr=NULL;
    vec_cnt=0;
  }
  
  /** コンストラクタ **/
  LinearSolver(const int sz[4],
               const int guide,
               const REAL_TYPE dh,
               const int IterationMax,
               const REAL_TYPE coef_ac,
               const double eps,
               PerfMonitor* PM,
               FILE* fp,
	       int* vec_tag,
	       int* vec_itr,
	       int  vec_cnt
               ) {
    this->guide = guide;
    this->dh = dh;
    this->IterationMax = IterationMax;
    this->coef_ac= coef_ac;
    this->eps = eps;
    
    this->PM     = PM;
    this->fp = fp;
    
    for (int i=0; i<4; i++)
    {
      this->size[i] = sz[i];
    }

    this->vec_tag=vec_tag;
    this->vec_itr=vec_itr;
    this->vec_cnt = vec_cnt;
  }

  /**　デストラクタ */
  ~LinearSolver() {}



protected:


  /**
   * @brief Fdot for 1 array
   * @retval  内積値
   * @param [in]   x   vector1
   */
  double* Fdot1(double* xy, REAL_TYPE* x, int* bcp, int* psz3d, int* psz4dex);


  /**
   * @brief Fdot for 2 arrays
   * @retval  内積値
   * @param [in]   x   vector1
   * @param [in]   y   vector2
   */
  double* Fdot2(double* xy,REAL_TYPE* x, REAL_TYPE* y, int* bcp, int* psz3d, int* psz4dex);

  /**
   * @brief Preconditioner
   * @param [in,out] x  解ベクトル
   * @param [in]     b  RHS vector
   * @param [in]     isPrecond 前処理あり(true)
   */
  void Preconditioner(REAL_TYPE* x,
                      REAL_TYPE* b,
                      REAL_TYPE* exs,
                      int* bcp,
                      int* psz3d,
                      int* psz4dex,
                      const bool isPrecond);

  
  
public:
  
  /**
   * @brief Jacobi relaxation
   * @retval 反復数
   * @param [in,out] x       解ベクトル
   * @param [in]     b       RHS vector
   * @param [in]     wk      work array
   * @param [in]     bcp     BCindex P
   * @param [in]     psz3d   s3dのパディング数
   * @param [in]     psz4dex s4dexのパディング数
   * @param [in]     itrMax  反復最大値
   * @param [in]     converge_check 収束判定あり(true)
   */
  int Jacobi(REAL_TYPE* x,
             REAL_TYPE* b,
             REAL_TYPE* wk,
             REAL_TYPE* exs,
             int* bcp,
             int* psz3d,
             int* psz4dex,
             const int itrMax,
             bool converge_check=true);
  
  /**
   * @brief SOR法
   * @retval 反復数
   * @param [in,out] x      解ベクトル
   * @param [in]     b      RHS vector
   * @param [in]     psz_x  xのパディング数
   * @param [in]     psz_b  bのパディング数
   * @param [in]     itrMax 反復最大値
   * @param [in]     converge_check 収束判定あり(true)
   */
  int PointSOR(REAL_TYPE* x,
               REAL_TYPE* b,
               REAL_TYPE* exs,
               int* bcp,
               int* psz3d,
               int* psz4dex,
               const int itrMax,
               bool converge_check=true);
  
  
  /**
   * @brief 2色オーダリングSORのストライドメモリアクセス版
   * @retval 反復数
   * @param [in,out] x              解ベクトル
   * @param [in]     b              RHS vector
   * @param [in]     itrMax         反復最大値
   * @param [in]     converge_check 収束判定あり(true)
   */
  int SOR2_SMA(REAL_TYPE* x,
               REAL_TYPE* b,
               REAL_TYPE* exs,
               int* bcp,
               int* psz3d,
               int* psz4dex,
               const int itrMax,
               bool converge_check=true);
  
  
  /**
   * @brief 前処理つきBiCGstab
   * @retval 反復数
   * @param [in,out] x       解ベクトル
   * @param [in]     b       RHS vector
   * @param [in]     itrMax  反復最大値
   * @param [in]     isPrecond 前処理あり(true)
   */
  int PBiCGstab(REAL_TYPE* x,
                REAL_TYPE* b,
                REAL_TYPE* exs,
                int* bcp,
                REAL_TYPE** wrk_x,
                REAL_TYPE** wrk_b,
                REAL_TYPE** wrk_q,
                REAL_TYPE** wrk_r,
                REAL_TYPE** wrk_r0,
                REAL_TYPE** wrk_p,
                REAL_TYPE* pcg_p_,
                REAL_TYPE* pcg_s,
                REAL_TYPE* pcg_s_,
                REAL_TYPE* pcg_t_,
                int* psz3d,
                int* psz4dex,
                const int ItrMax,
                const bool isPrecond);
  
  

};

#endif // _FFV_LS_H_
