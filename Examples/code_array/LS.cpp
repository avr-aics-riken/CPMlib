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
 * @file   ffv_LS.C
 * @brief  LS Class
 * @author aics
 */


#include "LS.h"
#define MX_SIZE 256


// #################################################################
double* LinearSolver::Fdot1(double* xy, REAL_TYPE* x, int* bcp, int* psz3d, int* psz4dex)
{
  double flop=0.0;          /// 浮動小数点演算数

  for(int i=0;i<size[3];i++){
    xy[i]=0.0;
  }

  TIMING_start(PM, "Blas_dot1");
  blas_dot1_(xy, x, bcp, size, &guide, psz3d, psz4dex, &flop, &vec_cnt);
  TIMING_stop(PM, "Blas_dot1", flop);

  return xy;
}


// #################################################################
double* LinearSolver::Fdot2(double* xy,REAL_TYPE* x, REAL_TYPE* y, int* bcp, int* psz3d, int* psz4dex)
{
  double flop=0.0;          /// 浮動小数点演算数

  for(int i=0;i<size[3];i++)
    xy[i]=0.0;

  TIMING_start(PM, "Blas_dot2");
  blas_dot2_(xy, x, y, bcp, size, &guide, psz3d, psz4dex, &flop, &vec_cnt);
  TIMING_stop(PM, "Blas_dot2", flop);

  return xy;
}


// #################################################################
int LinearSolver::Jacobi(REAL_TYPE* x,
                         REAL_TYPE* b,
                         REAL_TYPE* wk,
                         REAL_TYPE* exs,
                         int* bcp,
                         int* psz3d,
                         int* psz4dex,
                         const int itrMax,
                         bool converge_check)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  REAL_TYPE omg = coef_ac;
  int lc=0;                       /// ループカウント

  // x     圧力 p^{n+1}
  // b     RHS vector
  // bcp   ビットフラグ

  for (lc=1; lc<=itrMax; lc++)
    {
      double res[MX_SIZE];
      for(int i=0;i<size[3];i++) {res[i]=0.0;}

      // 反復処理
      TIMING_start(PM, "Jacobi_kernel");
      flop_count = 0.0;
      jacobi_(x, size, &guide, psz3d, psz4dex, &omg, b, bcp, res, wk, &flop_count, &vec_cnt);
      TIMING_stop(PM, "Jacobi_kernel", flop_count);

      // 境界条件
      TIMING_start(PM, "BoundaryCondition");
      bc_(size, &guide, psz4dex, x, &dh, vec_tag, &vec_cnt);
      TIMING_stop(PM, "BoundaryCondition");

      // 収束判定
      if ( converge_check )
	{
	  double er[MX_SIZE];
	  for(int i=0;i<size[3];i++) {er[i]=0.0;}
	  
	  //TIMING_start(PM, "Norm_max");
	  //err_(size, &guide, &dh, er, x, exs,gosa);
	  //TIMING_stop(PM, "Norm_max");
	  
	  double res2;
	  double max_res=0.0;
	  int    max_itr=0;
	  fprintf(fp, "%8d ", lc);
	  for(int i=0;i<size[3];i++){
	    res2 = sqrt(res[i]);
	    max_res = max(max_res,res2);
	    max_itr = max(max_itr,vec_itr[i]);
	    fprintf(fp, "%d:%13.6e ",vec_tag[i], res2);
	    //fprintf(fp, "%d:%13.6e %13.6e ",vec_tag[i],er[i], res2);
	  }
	  fprintf(fp, "\n");
	  
          for(int i=1;i<=vec_cnt;i++){//収束したものの探索
            if(sqrt(res[i-1])<eps){//実際はepsで比較する
              int tmp=0;
              for(double r=sqrt(res[vec_cnt-1]);((r<eps)&&(i<vec_cnt));){
		//一番外側のベクトルが収束してないか確認
                vec_cnt--;
		r=sqrt(res[vec_cnt-1]);
              }
              //タグの交換
              tmp=vec_tag[i-1];vec_tag[i-1]=vec_tag[vec_cnt-1];vec_tag[vec_cnt-1]=tmp;
              tmp=res[i-1];res[i-1]=res[vec_cnt-1];res[vec_cnt-1]=tmp;

	      TIMING_start(PM,"Vector_change_p");
              change_(size,&guide,psz4dex,x,&i,&vec_cnt); //配列の交換
	      TIMING_stop(PM,"Vector_change_p");

              TIMING_start(PM,"Vector_change_b");
              change_(size,&guide,psz4dex,b,&i,&vec_cnt); //配列の交換
	      TIMING_stop(PM,"Vector_change_b");

              //TIMING_start(PM,"Vector_change_exs");
              change_(size,&guide,psz4dex,exs,&i,&vec_cnt); //配列の交換
	      //TIMING_stop(PM,"Vector_change_p");
	      vec_cnt--;
            }
	  }
	  //if(max_itr<=lc) break;
	  //if(vec_itr[0]<=lc) break;
	  if (sqrt(res[0]) < eps ) break;
	}

    }

  return lc;
}


// #################################################################
int LinearSolver::PointSOR(REAL_TYPE* x,
                           REAL_TYPE* b,
                           REAL_TYPE* exs,
                           int* bcp,
                           int* psz3d,
                           int* psz4dex,
                           const int itrMax,
                           bool converge_check)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  REAL_TYPE omg = coef_ac;
  int lc=0;                       /// ループカウント
  
  // x     圧力 p^{n+1}
  // b     RHS vector
  // bcp   ビットフラグ
  for (lc=1; lc<=itrMax; lc++)
    {
      double res[MX_SIZE];
      for(int i=0;i<size[3];i++) {res[i]=0.0;}

      // 反復処理
      TIMING_start(PM, "Sor_kernel");
      flop_count = 0.0;
      psor_(x, size, &guide, psz3d, psz4dex, &omg, b, bcp, res, &flop_count, &vec_cnt);
      TIMING_stop(PM, "Sor_kernel", flop_count);

      if ( converge_check )
	{
	  // 境界条件
	  TIMING_start(PM, "BoundaryCondition");
	  bc_(size, &guide, psz4dex, x, &dh, vec_tag, &vec_cnt);
	  TIMING_stop(PM, "BoundaryCondition");
	}
      
      // 収束判定
      if ( converge_check )
	{

	  double er[MX_SIZE];
	  for(int i=0;i<size[3];i++) {er[i]=0.0;}
	  
	  //TIMING_start(PM, "Norm_max");
	  //err_(size, &guide, &dh, er, x, exs,gosa);
	  //TIMING_stop(PM, "Norm_max");
	  double res2;
	  double max_res=0.0;
	  fprintf(fp, "%8d ", lc);
	  for(int i=0;i<size[3];i++){
	    res2 = sqrt(res[i]);
	    max_res = max(max_res,res2);
	    fprintf(fp, "%d:%13.6e ",vec_tag[i], res2);
	    //fprintf(fp, "%d:%13.6e %13.6e ",vec_tag[i],er[i], res2);
	  }
	  fprintf(fp, "\n");
	  
          for(int i=1;i<=vec_cnt;i++){//収束したものの探索
            if(sqrt(res[i-1])<eps){//実際はepsで比較する
              int tmp=0;
              for(double r=sqrt(res[vec_cnt-1]);((r<eps)&&(i<vec_cnt));){
		//一番外側のベクトルが収束してないか確認
                vec_cnt--;
		r=sqrt(res[vec_cnt-1]);
              }
              //タグの交換
              tmp=vec_tag[i-1];vec_tag[i-1]=vec_tag[vec_cnt-1];vec_tag[vec_cnt-1]=tmp;
              tmp=res[i-1];res[i-1]=res[vec_cnt-1];res[vec_cnt-1]=tmp;

	      TIMING_start(PM,"Vector_change_p");
              change_(size,&guide,psz4dex,x,&i,&vec_cnt); //配列の交換
	      TIMING_stop(PM,"Vector_change_p");

              TIMING_start(PM,"Vector_change_b");
              change_(size,&guide,psz4dex,b,&i,&vec_cnt); //配列の交換
	      TIMING_stop(PM,"Vector_change_b");

              //TIMING_start(PM,"Vector_change_exs");
              change_(size,&guide,psz4dex,exs,&i,&vec_cnt); //配列の交換
	      //TIMING_stop(PM,"Vector_change_p");
	      vec_cnt--;
            }
	  }
	  //if(max_itr<=lc) break;
	  //if(vec_itr[0]<=lc) break;
	  if (sqrt(res[0]) < eps ) break;
	}


    }

  return lc;
}


// #################################################################
int LinearSolver::SOR2_SMA(REAL_TYPE* x,
                           REAL_TYPE* b,
                           REAL_TYPE* exs,
                           int* bcp,
                           int* psz3d,
                           int* psz4dex,
                           const int itrMax,
                           bool converge_check)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  int lc=0;                       /// ループカウント
  REAL_TYPE omg = coef_ac;
  
  // x     圧力 p^{n+1}
  // b     RHS vector
  // d_bcp ビットフラグ
  
  for (lc=1; lc<=itrMax; lc++)
    {
      // 2色のマルチカラー(Red&Black)のセットアップ

      // ip = 0 基点(1,1,1)がRからスタート
      //    = 1 基点(1,1,1)がBからスタート
      int ip = 0;

      double res[MX_SIZE];
      for(int i=0;i<size[3];i++) {res[i]=0.0;}

      // 各カラー毎の間に同期, 残差は色間で積算する
      // R - color=0 / B - color=1
      for (int color=0; color<2; color++) {

	TIMING_start(PM, "Sor2sma_kernel");flop_count = 0.0; // 色間で積算しない
	psor2sma_core_(x, size, &guide, psz3d, psz4dex, &ip, &color, &omg, b, bcp, res, &flop_count, &vec_cnt);
	TIMING_stop(PM, "Sor2sma_kernel", flop_count);

	// 境界条件 >> PBiCGstabで前処理のときには境界条件を呼ばない。呼んではならぬ？ bのBCではないから？
	if ( converge_check )
	  {
	    TIMING_start(PM, "BoundaryCondition");
	    bc_(size, &guide, psz4dex, x, &dh, vec_tag, &vec_cnt);
	    TIMING_stop(PM, "BoundaryCondition");
	  }
      }

      // 収束判定
      if ( converge_check )
	{
	  double er[MX_SIZE];
	  for(int i=0;i<size[3];i++) {er[i]=0.0;}
	  
	  //TIMING_start(PM, "Norm_max");
	  //err_(size, &guide, &dh, er, x, exs,gosa);
	  //TIMING_stop(PM, "Norm_max");
	  double res2;
	  double max_res=0.0;
	  int    max_itr=0;
	  fprintf(fp, "%8d ", lc);
	  for(int i=0;i<size[3];i++){
	    res2 = sqrt(res[i]);
	    max_res = max(max_res,res2);
	    max_itr = max(max_itr,vec_itr[i]);
	    fprintf(fp, "%d:%13.6e ",vec_tag[i], res2);
	    //fprintf(fp, "%d:%13.6e %13.6e ",vec_tag[i],er[i], res2);
	  }
	  fprintf(fp, "\n");
	  
          for(int i=1;i<=vec_cnt;i++){//収束したものの探索
            if(sqrt(res[i-1])<eps){//実際はepsで比較する
              int tmp=0;
              for(double r=sqrt(res[vec_cnt-1]);((r<eps)&&(i<vec_cnt));){
		//一番外側のベクトルが収束してないか確認
                vec_cnt--;
		r=sqrt(res[vec_cnt-1]);
              }
              //タグの交換
              tmp=vec_tag[i-1];vec_tag[i-1]=vec_tag[vec_cnt-1];vec_tag[vec_cnt-1]=tmp;
              tmp=res[i-1];res[i-1]=res[vec_cnt-1];res[vec_cnt-1]=tmp;

	      TIMING_start(PM,"Vector_change_p");
              change_(size,&guide,psz4dex,x,&i,&vec_cnt); //配列の交換
	      TIMING_stop(PM,"Vector_change_p");

              TIMING_start(PM,"Vector_change_b");
              change_(size,&guide,psz4dex,b,&i,&vec_cnt); //配列の交換
	      TIMING_stop(PM,"Vector_change_b");

              //TIMING_start(PM,"Vector_change_exs");
              change_(size,&guide,psz4dex,exs,&i,&vec_cnt); //配列の交換
	      //TIMING_stop(PM,"Vector_change_p");
	      vec_cnt--;
            }
	  }
	  //if(max_itr<=lc) break;
	  //if(vec_itr[0]<=lc) break;
	  if (sqrt(res[0]) < eps ) break;
	}

    }

  return lc;

}

// #################################################################
void LinearSolver::Preconditioner(REAL_TYPE* x,
                                  REAL_TYPE* b,
                                  REAL_TYPE* exs,
                                  int* bcp,
                                  int* psz3d,
                                  int* psz4dex,
                                  const bool isPrecond)
{
  // 前処理なし(コピー)
  if ( !isPrecond )
  {
      TIMING_start(PM, "Blas_copy");
      blas_copy_(x, b, size, &guide, psz4dex, &vec_cnt);
      TIMING_stop(PM, "Blas_copy");
    return;
  }
  
  int lc_max = 10;
  
  // 前処理
  //SOR2_SMA(x, b, lc_max, false);
  PointSOR(x, b, exs, bcp, psz3d, psz4dex, lc_max, false);
}


// #################################################################
// PBiCBSTAB 収束判定は残差
// @note 反復回数がマシンによって異なる現象がある．
// Xeon E5では同じ反復回数になるのに対して，Core i7では反復回数が試行毎に異なる．
// 内積のOpenMP並列のためか？
int LinearSolver::PBiCGstab(REAL_TYPE* x,
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
                            const bool isPrecond)
{
  double flop = 0.0;

  TIMING_start(PM, "Blas_copy");
  blas_copy_(wrk_x[0], x, size, &guide, psz4dex, &vec_cnt);
  blas_copy_(wrk_b[0], b, size, &guide, psz4dex, &vec_cnt);
  TIMING_stop(PM, "Blas_copy");

  int tag_con[MX_SIZE];//収束したもののタグ
  int table_con[MX_SIZE];//収束したもののテーブル
  int num_con=0;

  int tag_not_con[MX_SIZE];
  int table_not_con[MX_SIZE];
  int num_not_con=0;

  for(int i=0;i<size[3];i++){
    tag_not_con[i] = vec_tag[i];
  }
  

  int color=0;

  TIMING_start(PM, "Blas_clear");
  blas_clear_(wrk_q[0], size, &guide, psz4dex, &vec_cnt);
  TIMING_stop(PM, "Blas_clear");

  TIMING_start(PM, "Blas_residual");
  flop = 0.0;
  blas_calc_rk_(wrk_r[0], x, b, bcp, size, &guide, psz3d, psz4dex, &flop,&vec_cnt);
  TIMING_stop(PM, "Blas_residual", flop);

  TIMING_start(PM, "Blas_copy");
  blas_copy_(wrk_r0[0], wrk_r[0], size, &guide, psz4dex, &vec_cnt);
  TIMING_stop(PM, "Blas_copy");

  double rho_old[2][MX_SIZE];
  double alpha[2][MX_SIZE]; //1.0とどうちがう？？？
  double omega[2][MX_SIZE];
  double r_omega[2][MX_SIZE];
  int lc=0;             /// ループカウント

  for(int j=0;j<2;j++){
    for(int i=0;i<size[3];i++){
      rho_old[j][i] = 1.0;
      alpha[j][i] = 0.0; //1.0とどうちがう？？？
      omega[j][i]  = 1.0;
      r_omega[j][i] = -omega[j][i];
    }
  }

  for (lc=1; lc<=ItrMax; lc++)
    {
      double rho[MX_SIZE];

      Fdot2(rho,wrk_r[color], wrk_r0[color], bcp, psz3d, psz4dex);
      double max_rho = 0.0;
      for(int i=0;i<size[3];i++){
        max_rho = max(max_rho,fabs(rho[i]));
      }//残差の最大のもの

      if( max_rho < REAL_TYPE_EPSILON )
        {
          lc = 0;
          break;
        }
      
      if( lc == 1 )
        {
          TIMING_start(PM, "Blas_copy");
          blas_copy_(wrk_p[color], wrk_r[color], size, &guide, psz4dex, &vec_cnt);
          TIMING_stop(PM, "Blas_copy");
        }
      else
        {
          double beta[MX_SIZE];
          for(int i=0;i<size[3];i++){
	    beta[i] = rho[i] / rho_old[color][i] * alpha[color][i] / omega[color][i];
	  }

	  TIMING_start(PM, "Blas_bicg1");
	  flop = 0.0;
	  blas_bicg_1_(wrk_p[color], wrk_r[color], wrk_q[color], beta, omega[color], size, &guide, psz4dex, &flop,&vec_cnt);
	  TIMING_stop(PM, "Blas_bicg1", flop);
        }

      TIMING_start(PM, "Blas_clear");
      blas_clear_(pcg_p_, size, &guide, psz4dex, &vec_cnt);
      TIMING_stop(PM, "Blas_clear");

      Preconditioner(pcg_p_, wrk_p[color], exs, bcp, psz3d, psz4dex, isPrecond);

      TIMING_start(PM, "Blas_ax");
      flop = 0.0;
      blas_calc_ax_(wrk_q[color], pcg_p_, bcp, size, &guide, psz3d, psz4dex, &flop,&vec_cnt);
      TIMING_stop(PM, "Blas_ax", flop);

      double r_alpha[MX_SIZE];
      double tmp_fdot2[MX_SIZE];
      Fdot2(tmp_fdot2,wrk_q[color], wrk_r0[color], bcp, psz3d, psz4dex);
      for(int i=0;i<size[3];i++){
        alpha[color][i] = rho[i] / tmp_fdot2[i];
        r_alpha[i] = -alpha[color][i];
      }

      TIMING_start(PM, "Blas_triad");
      flop = 0.0;
      blas_triad_(pcg_s, wrk_q[color], wrk_r[color], r_alpha, size, &guide, psz4dex, &flop,&vec_cnt);
      TIMING_stop(PM, "Blas_triad", flop);

      TIMING_start(PM, "Blas_clear");
      blas_clear_(pcg_s_, size, &guide, psz4dex, &vec_cnt);
      TIMING_stop(PM, "Blas_clear");

      Preconditioner(pcg_s_, pcg_s, exs, bcp, psz3d, psz4dex, isPrecond);

      TIMING_start(PM, "Blas_ax");
      flop = 0.0;
      blas_calc_ax_(pcg_t_, pcg_s_, bcp, size, &guide, psz3d, psz4dex, &flop, &vec_cnt);
      TIMING_stop(PM, "Blas_ax", flop);

      double tmp_fdot1[MX_SIZE];
      Fdot1(tmp_fdot1,pcg_t_, bcp, psz3d, psz4dex);
      Fdot2(tmp_fdot2,pcg_t_, pcg_s, bcp, psz3d, psz4dex);
      for(int i=0;i<size[3];i++){
        omega[color][i] = tmp_fdot2[i] / tmp_fdot1[i];
        r_omega[color][i] = -omega[color][i];
      }

      TIMING_start(PM, "Blas_bicg2");
      flop = 0.0;
      blas_bicg_2_(wrk_x[color], pcg_p_, pcg_s_, alpha[color] , omega[color], size, &guide, psz4dex, &flop, &vec_cnt);
      TIMING_stop(PM, "Blas_bicg2", flop);

      TIMING_start(PM, "Blas_triad");
      flop = 0.0;
      blas_triad_(wrk_r[color], pcg_t_, pcg_s, r_omega[color], size, &guide, psz4dex, &flop, &vec_cnt);
      TIMING_stop(PM, "Blas_triad", flop);
      double res[MX_SIZE];
      Fdot1(res,wrk_r[color], bcp, psz3d, psz4dex);

      for(int i=0;i<size[3];i++){
        rho_old[color][i] = rho[i];
      }

      double er[MX_SIZE];
      for(int i=0;i<size[3];i++){
        er[i] = 0.0;
      }

      //TIMING_start(PM, "Norm_max");
      //err_(size, &guide, &dh, er, x, exs, gosa); 
      //TIMING_stop(PM, "Norm_max");
      
      double res2;
      double max_res=0.0;

      fprintf(fp, "%8d ", lc);
      for(int i=0;i<vec_cnt;i++){
	res2 = sqrt(res[i]);
	max_res = max(max_res,res2);
	fprintf(fp, "%d:%13.6e ",tag_not_con[i], res2);
	//fprintf(fp, "%d:%13.6e %13.6e ",vec_tag[i],er[i], res2);
      }
      fprintf(fp, "\n");
      
      
      num_con=0;
      num_not_con=0;
      for(int i=1;i<=vec_cnt;i++){//収束したものの探索
	if(sqrt(res[i-1])<eps){
	  tag_con[num_con] = tag_not_con[i-1];
	  table_con[num_con] = i;
	  num_con++;
	}
	else{
 	  tag_not_con[num_not_con] = tag_not_con[i-1];
	  table_not_con[num_not_con] = i;
	  num_not_con++;
	}
      }
      //printf("%d,%d,%d\n",num_con,num_not_con,vec_cnt);

      
      if(num_con>0 && num_not_con>0){
	//printf("%d\n",vec_cnt);
	TIMING_start(PM,"Vector_change_exs");
	array_reshape_(size, &guide, psz4dex, wrk_x[color], wrk_x[(color+1)%2], x, &vec_cnt, table_con, &num_con, table_not_con, &num_not_con);
	array_reshape_(size, &guide, psz4dex, wrk_b[color], wrk_b[(color+1)%2], b, &vec_cnt, table_con, &num_con, table_not_con, &num_not_con);
	//printf("%d\n",vec_cnt);
	array_small_(size, &guide, psz4dex, wrk_q[color], wrk_q[(color+1)%2], &vec_cnt, table_not_con, &num_not_con);
	array_small_(size, &guide, psz4dex, wrk_r[color], wrk_r[(color+1)%2], &vec_cnt, table_not_con, &num_not_con);
	array_small_(size, &guide, psz4dex, wrk_r0[color], wrk_r0[(color+1)%2], &vec_cnt, table_not_con, &num_not_con);
	array_small_(size, &guide, psz4dex, wrk_p[color], wrk_p[(color+1)%2], &vec_cnt, table_not_con, &num_not_con);
	TIMING_stop(PM,"Vector_change_exs");

	for(int i=1;i<=num_not_con;i++){
	  rho_old[(color+1)%2][i-1] = rho_old[color][table_not_con[i-1]-1];
	  alpha[(color+1)%2][i-1] = alpha[color][table_not_con[i-1]-1];
	  omega[(color+1)%2][i-1] = omega[color][table_not_con[i-1]-1];
	  r_omega[(color+1)%2][i-1] = r_omega[color][table_not_con[i-1]-1];
	}

	int offset = size[3]-vec_cnt;
	for(int i=0;i<num_con;i++){
	  vec_tag[offset + i] = tag_con[i];
	}

	vec_cnt = num_not_con;
	color = (color + 1)%2;
      }      
	
      if(max_res<=eps) {
	//blas_copy_(x, wrk_x[0], size, &guide,&vec_cnt);

	num_not_con=1;
	TIMING_start(PM,"Vector_change_exs");
	array_reshape_(size, &guide, psz4dex, wrk_x[color], wrk_x[(color+1)%2], x, &vec_cnt, table_con, &num_con, table_not_con, &num_not_con);
	TIMING_stop(PM,"Vector_change_exs");

	int offset = size[3]-vec_cnt;
	for(int i=0;i<num_con;i++){
	  vec_tag[offset + i] = tag_con[i];
	}

	for(int i=0;i<size[3];i++){
	  fprintf(fp, "%d, ",vec_tag[i]);
	}
	fprintf(fp, "\n");
	break;  
      }
    }

  // 境界条件
  TIMING_start(PM, "BoundaryCondition");
  bc_(size, &guide, psz4dex, x, &dh, vec_tag, &vec_cnt);
  TIMING_stop(PM, "BoundaryCondition");

  return lc;
}

