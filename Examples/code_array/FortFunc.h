/*
 *  FortFunc.h
 *  kernel_test
 *
 *  Created by keno on 11/12/05.
 *  Copyright 2011 __iis__. All rights reserved.
 *
 */

#ifndef _FORTRAN_FUNC_H_
#define _FORTRAN_FUNC_H_

#include "FBdefine.h"

extern "C" {
  // utility.f90
  void fileout_ (int* sz,
                 int* g,
                 int* psz4dex,
                 REAL_TYPE* s,
                 REAL_TYPE* dh,
                 REAL_TYPE* org,
                 char* fname);
  
  void bc_      (int* sz,
                 int* g,
                 int* psz4dex,
                 REAL_TYPE* p,
                 REAL_TYPE* dh,
		 int* vec_tag,
		 int* vec_cnt);
  
  void exact_   (int* sz,
                 int* g,
                 REAL_TYPE* e,
                 REAL_TYPE* dh);
  
  void err_     (int* sz,
                 int* g,
                 REAL_TYPE* dh,
                 double* d,
                 REAL_TYPE* p,
                 REAL_TYPE* e,
		 REAL_TYPE* gosa);

  void set_bcindex_(int* sz,
                    int* g,
                    float* bp);

  void get_vec_sz_(int* sz);

  void change_(int* sz,
	       int* g,
               int* psz4dex,
	       REAL_TYPE* p,
	       int* ch,
	       int* vec_cnt);

  void array_reshape_(int* sz,
		      int* g,
                      int* psz4dex,
		      REAL_TYPE* p0,
		      REAL_TYPE* p1,
		      REAL_TYPE* p,
		      int* vec_cnt,
		      int* table_con,
		      int* num_con,
		      int* table_not_con,
		      int* num_not_con);

  void array_small_(int* sz,
		    int* g,
                    int* psz4dex,
		    REAL_TYPE* p0,
		    REAL_TYPE* p1,
		    int* vec_cnt,
		    int* table_not_con,
		    int* num_not_con);

  void init_array_(int* sz,
		   int* g,
		   REAL_TYPE* p);


  // linear_solver.f90
  void jacobi_        (REAL_TYPE* p,
		       int* sz,
                       int* g,
		       int* psz3d,
		       int* psz4dex,
                       REAL_TYPE* omg,
                       REAL_TYPE* b,
                       int* bp,
                       double* res,
                       REAL_TYPE* wk2,
                       double* flop,
		       int* vec_cnt);
  
  void psor_          (REAL_TYPE* p,
		       int* sz,
                       int* g,
		       int* psz3d,
		       int* psz4dex,
                       REAL_TYPE* omg,
                       REAL_TYPE* b,
                       int* bp,
                       double* res,
                       double* flop,
		       int* vec_cnt);
  
  void psor2sma_core_ (REAL_TYPE* p,
		       int* sz,
                       int* g,
		       int* psz3d,
		       int* psz4dex,
                       int* ip,
                       int* color,
                       REAL_TYPE* omg,
                       REAL_TYPE* b,
                       int* bp,
                       double* res,
                       double* flop,
		       int* vec_cnt);
  
  void src_dirichlet_ (REAL_TYPE* b,
                       int* sz,
                       int* g,
                       int* psz3d,
                       int* psz4dex,
                       int* bp,
                       REAL_TYPE* dh);

  
  // ffv_blas.f90
  void blas_clear_    (REAL_TYPE* x,
                       int* sz,
                       int* g,
		       int* psz4dex,
		       int* vec_cnt);
  
  void blas_copy_     (REAL_TYPE* dst,
                       REAL_TYPE* src,
                       int* sz,
                       int* g,
		       int* psz4dex,
		       int* vec_cnt);
  
  void blas_triad_    (REAL_TYPE* z,
                       REAL_TYPE* x,
                       REAL_TYPE* y,
                       double* a,
                       int* sz,
                       int* g,
		       int* psz4dex,
                       double* flop,
		       int* vec_cnt);
  
  void blas_bicg_1_ (REAL_TYPE* p,
                     REAL_TYPE* r,
                     REAL_TYPE* q,
                     double* beta,
                     double* omg,
                     int* sz,
                     int* g,
		     int* psz4dex,
                     double* flop,
		     int* vec_cnt);
  
  void blas_bicg_2_   (REAL_TYPE* z,
                       REAL_TYPE* x,
                       REAL_TYPE* y,
                       double* a,
                       double* b,
                       int* sz,
                       int* g,
		       int* psz4dex,
                       double* flop,
		       int* vec_cnt);
  
  void blas_dot1_     (double* r,
                       REAL_TYPE* p,
                       int* bp,
                       int* sz,
                       int* g,
		       int* psz3d,
		       int* psz4dex,
                       double* flop,
		       int* vec_cnt);
  
  void blas_dot2_     (double* r,
                       REAL_TYPE* p,
                       REAL_TYPE* q,
                       int* bp,
                       int* sz,
                       int* g,
		       int* psz3d,
		       int* psz4dex,
                       double* flop,
		       int* vec_cnt);
  
  void blas_calc_b_ (double* rhs,
                     REAL_TYPE* b,
                     REAL_TYPE* s_0,
                     int* bp,
                     int* sz,
                     int* g,
                     REAL_TYPE* dh,
                     REAL_TYPE* dt,
                     double* flop);
  
  void blas_calc_rk_  (REAL_TYPE* r,
                       REAL_TYPE* p,
                       REAL_TYPE* b,
                       int* bp,
                       int* sz,
                       int* g,
		       int* psz3d,
		       int* psz4dex,
                       double* flop,
		       int* vec_cnt);
  
  void blas_calc_r2_  (double* res,
                       REAL_TYPE* p,
                       REAL_TYPE* b,
                       int* bp,
                       int* sz,
                       int* g,
                       double* flop);
  
  void blas_calc_ax_  (REAL_TYPE* ap,
                       REAL_TYPE* p,
                       int* bp,
                       int* sz,
                       int* g,
		       int* psz3d,
		       int* psz4dex,
                       double* flop,
		       int* vec_cnt);
}

#endif // _FORTRAN_FUNC_H_
