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

#include "cpm_ParaManager.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>

#include "kernel_def.h"
#include "FortFunc.h"
#include "voxinfo.h"
#include "LS.h"

using namespace std;
using namespace pm_lib;

int order_of_PM_key;      ///< PMlib用の登録番号カウンタ < PM_NUM_MAX

// プロトタイプ
void set_timing_label(PerfMonitor* pm);

////////////////////////////////////////////////////////////////////////////////////
// テストプログラムのメイン
////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char **argv )
{
  int ret = 0;
  order_of_PM_key = 0;

  // 並列管理クラスのインスタンスと初期化
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance(argc,argv);
  if( !paraMngr ) return CPM_ERROR_PM_INSTANCE;
  if( paraMngr->GetMyRankID() == 0 )
  {
    cout << "CPMlib Version " << 
    cpm_Base::getVersionInfo() << endl;
    cout << "CPMlib Revision " << 
    cpm_Base::getRevisionInfo() << endl;
  }

  int myrank = paraMngr->GetMyRankID();
  int nrank  = paraMngr->GetNumRank();
  if( nrank != 1 )
  {
    if( myrank == 0 )
    {
      printf("\tthis program is only supportted np=1.\n");
    }
    return 0;
  }

  // プログラム引数の取得
  if( argc != 6 )
  {
    if( myrank == 0 )
    {
      printf("\tUsage: %s array_size linear_solver Iteration vector_num padding_flag\n", argv[0]);
      printf("\t$(ex.) %s 128 sor 1000 2 1\n", argv[0]);
      printf("\n");
      printf("\t  array_size    : array size\n" );
      printf("\t  linear_solver : solver method [jacobi | sor | sor2sma | pbicgstab | bicgstab]\n" );
      printf("\t  Iteration     : number of iteration\n" );
      printf("\t  vector_num    : number of RHS\n" );
      printf("\t  padding_flag  : padding flag [0:off, 1:on]\n" );
    }
    return 9;
  }
  int   dim     = atoi(argv[1]);
  char* method  = argv[2];
  int   ItrMax  = atoi(argv[3]);
  int   nrhs    = atoi(argv[4]);
  bool  padding = (atoi(argv[5])==0 ? false : true);
  if( myrank == 0 )
  {
    printf("\tarray size             : %d^3\n", dim);
    printf("\tparallel vector number : %d\n"  , nrhs);
    printf("\tLinear solver          : %s\n"  , method);
    printf("\tmax iteration          : %d\n"  , ItrMax);
    if( padding )
      printf("\tpadding                : on\n");
    else
      printf("\tpadding                : off\n");
  }

  // 領域情報、パラメータをセット
  int gc = 2; // guide cell
  int sz[3]        = {dim, dim, dim};
  int v_sz[4]      = {dim, dim, dim, nrhs};

  int vec_cnt = v_sz[3];
  int* vec_tag;
  vec_tag = new int[v_sz[3]];
  int* vec_itr;
  vec_itr = new int[v_sz[3]];
  for (int i=0; i<v_sz[3]; i++) {
    *(vec_tag+i) = i+1;
    *(vec_itr+i) = 10000;
  }

  // constant
  double eps        = 1.0e-5;                 // convergence criteria
  REAL_TYPE ac1     = 1.7;                    // acceleration coef. for SOR
  REAL_TYPE ac2     = 0.8;                    // acceleration coef. for jacobi relaxation
  REAL_TYPE dh      = 1.0/(REAL_TYPE)(dim+1); // grid width
  double    drgn[3] = {1.0, 1.0, 1.0};        // region size
  double    dorg[3] = {0.0, 0.0, 0.0};        // original coordinates
  REAL_TYPE rgn[3]  = {(REAL_TYPE)drgn[0], (REAL_TYPE)drgn[1], (REAL_TYPE)drgn[2]};
  REAL_TYPE org[3]  = {(REAL_TYPE)dorg[0], (REAL_TYPE)dorg[1], (REAL_TYPE)dorg[2]};

  // type
  int ls_type = 0;
  REAL_TYPE coef = 1.0;
  char fname[20];
  memset(fname, 0, sizeof(char)*20);
  if( !strcasecmp(method, "jacobi") )
  {
    ls_type = JACOBI;
    strcpy(fname, "jacobi.txt");
    coef = ac2;
  }
  else if( !strcasecmp(method, "sor") )
  {
    ls_type = SOR;
    strcpy(fname, "sor.txt");
    coef = ac1;
  }
  else if( !strcasecmp(method, "sor2sma") )
  {
    ls_type = SOR2SMA;
    strcpy(fname, "sor2sma.txt");
    coef = ac1;
  }
  else if( !strcasecmp(method, "pbicgstab") )
  {
    ls_type = PBICGSTAB;
    strcpy(fname, "pbicgstab.txt");
    coef = ac1;
  }
  else if( !strcasecmp(method, "bicgstab") )
  {
    ls_type = BICGSTAB;
    strcpy(fname, "bicgstab.txt");
    coef = ac1;
  }
  else
  {
    printf("Invalid solver\n");
    exit(0);

  }

  // 領域分割
  if( (ret = paraMngr->VoxelInit( sz, dorg, drgn )) != CPM_SUCCESS )
  {
    cerr << "VoxelInit error : " << ret << endl;
    return ret;
  }

  // 配列領域確保
  REAL_TYPE* p      = NULL;    // pressure
  REAL_TYPE* src0   = NULL;    // source term 0
  REAL_TYPE* wrk    = NULL;    // work array
  REAL_TYPE* exs    = NULL;    // exact solution
  unsigned*  src_bp = NULL;
  int psz3d[3], psz4dex[4];
  if( sizeof(REAL_TYPE) == 4 )
  {
    p    = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
    src0 = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
    exs  = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
    wrk  = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
  }
  else
  {
    p    = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
    src0 = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
    exs  = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
    wrk  = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
  }
  src_bp = (unsigned*)paraMngr->AllocIntS3D(gc, padding, psz3d);
  // 初期化
  {
    size_t sz;
    // s4dex
    sz = v_sz[3]+psz4dex[0];
    sz *= (v_sz[0]+2*gc+psz4dex[1]);
    sz *= (v_sz[1]+2*gc+psz4dex[2]);
    sz *= (v_sz[2]+2*gc+psz4dex[3]);
    paraMngr->InitArray(p   , sz);
    paraMngr->InitArray(src0, sz);
    paraMngr->InitArray(exs , sz);
    paraMngr->InitArray(wrk , sz);
    // s3d
    sz = v_sz[0]+2*gc+psz3d[0];
    sz *= (v_sz[1]+2*gc+psz3d[1]);
    sz *= (v_sz[2]+2*gc+psz3d[2]);
    paraMngr->InitArray((int*)src_bp, sz);
  }
  printf("\t  s4dex padding : %d %d %d %d\n", psz4dex[0], psz4dex[1], psz4dex[2], psz4dex[3]);
  printf("\t  s3d   padding : %d %d %d\n", psz3d[0], psz3d[1], psz3d[2]);

  // for BiCGSTAB
  REAL_TYPE** wrk_x = NULL;
  REAL_TYPE** wrk_b  = NULL;
  REAL_TYPE** wrk_q = NULL;
  REAL_TYPE** wrk_r  = NULL;
  REAL_TYPE** wrk_r0  = NULL; // work for bicgstab
  REAL_TYPE** wrk_p  = NULL;
  REAL_TYPE* pcg_p_ = NULL;
  REAL_TYPE* pcg_s  = NULL;
  REAL_TYPE* pcg_s_ = NULL;
  REAL_TYPE* pcg_t_ = NULL;
  if( ls_type == BICGSTAB || ls_type == PBICGSTAB )
  {
    wrk_x = new REAL_TYPE*[2];
    wrk_b = new REAL_TYPE*[2];
    wrk_q = new REAL_TYPE*[2];
    wrk_r = new REAL_TYPE*[2];
    wrk_r0 = new REAL_TYPE*[2];
    wrk_p = new REAL_TYPE*[2];
    if( sizeof(REAL_TYPE) == 4 )
    {
      wrk_x[0]  = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      wrk_x[1]  = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      wrk_b[0]  = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      wrk_b[1]  = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      wrk_q[0]  = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      wrk_q[1]  = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      wrk_r[0]  = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      wrk_r[1]  = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      wrk_r0[0] = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      wrk_r0[1] = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      wrk_p[0]  = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      wrk_p[1]  = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      pcg_p_ = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      pcg_s  = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      pcg_s_ = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
      pcg_t_ = (REAL_TYPE*)paraMngr->AllocFloatS4DEx(nrhs, gc, padding, psz4dex);
    }
    else
    {
      wrk_x[0]  = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      wrk_x[1]  = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      wrk_b[0]  = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      wrk_b[1]  = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      wrk_q[0]  = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      wrk_q[1]  = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      wrk_r[0]  = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      wrk_r[1]  = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      wrk_r0[0] = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      wrk_r0[1] = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      wrk_p[0]  = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      wrk_p[1]  = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      pcg_p_ = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      pcg_s  = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      pcg_s_ = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
      pcg_t_ = (REAL_TYPE*)paraMngr->AllocDoubleS4DEx(nrhs, gc, padding, psz4dex);
    }
    // 初期化
    {
      size_t sz;
      // s4dex
      sz = v_sz[3]+psz4dex[0];
      sz *= (v_sz[0]+2*gc+psz4dex[1]);
      sz *= (v_sz[1]+2*gc+psz4dex[2]);
      sz *= (v_sz[2]+2*gc+psz4dex[3]);
      paraMngr->InitArray(wrk_x[0], sz);
      paraMngr->InitArray(wrk_x[1], sz);
      paraMngr->InitArray(wrk_b[0], sz);
      paraMngr->InitArray(wrk_b[1], sz);
      paraMngr->InitArray(wrk_q[0], sz);
      paraMngr->InitArray(wrk_q[1], sz);
      paraMngr->InitArray(wrk_r[0], sz);
      paraMngr->InitArray(wrk_r[1], sz);
      paraMngr->InitArray(wrk_r0[0], sz);
      paraMngr->InitArray(wrk_r0[1], sz);
      paraMngr->InitArray(wrk_p[0], sz);
      paraMngr->InitArray(wrk_p[1], sz);
      paraMngr->InitArray(pcg_p_, sz);
      paraMngr->InitArray(pcg_s, sz);
      paraMngr->InitArray(pcg_s_, sz);
      paraMngr->InitArray(pcg_t_, sz);
    }
  }


  // Apply BC
  bc_(v_sz, &gc, psz4dex, p, &dh, vec_tag, &vec_cnt);

  // exact solution
  //exact_(v_sz, &gc, exs, &dh);

  char fname2[20];
  memset(fname2, 0, sizeof(char)*20);
  strcpy(fname2, "exact.sph");
  //fileout_(v_sz, &gc, exs, &dh, org, fname2);

  // setup BC index
  VoxInfo vi;
  vi.setBCIndexP(sz, gc, psz3d, src_bp);

  // history title
  FILE* fp;
  if ( !(fp=fopen(fname, "w")) )
  {
    printf("\tSorry, can't open file.\n");
    assert(0);
  }
  
  fprintf(fp,"Column_Data_00\n");
  fprintf(fp, "Itration          Norm      Residual\n");
    
  // timing
  int num_thread  = omp_get_max_threads();
  PerfMonitor PM;
  PM.initialize( PM_NUM_MAX );
  PM.setRankInfo( myrank );
  PM.setParallelMode("OpenMP", num_thread, nrank);
  set_timing_label(&PM);
  
  // source term
  src_dirichlet_(src0, v_sz, &gc, psz3d, psz4dex, (int *)src_bp, &dh);

  // LS class
  LinearSolver LS(v_sz,
                  gc,
                  dh,
                  ItrMax,
                  coef,
                  eps,
                  &PM,
                  fp,
                  vec_tag,
                  vec_itr,
                  vec_cnt);

  int loop;
  
  // scheme branch
  switch (ls_type) {
    case JACOBI:
      TIMING_start(&PM, "Jacobi");
      loop = LS.Jacobi(p, src0, wrk, exs, (int*)src_bp, psz3d, psz4dex, ItrMax);
      TIMING_stop(&PM, "Jacobi");
      break;
      
    case SOR:
          TIMING_start(&PM, "PointSOR");
          loop = LS.PointSOR(p, src0, exs, (int*)src_bp, psz3d, psz4dex, ItrMax);
          TIMING_stop(&PM, "PointSOR");
      break;

    case SOR2SMA:
          TIMING_start(&PM, "SOR2(SMA)");
          loop = LS.SOR2_SMA(p, src0, exs, (int*)src_bp, psz3d, psz4dex, ItrMax);
          TIMING_stop(&PM, "SOR2(SMA)");
      break;

    case PBICGSTAB:
          TIMING_start(&PM, "BiCGstab_w_precnd");
          loop = LS.PBiCGstab(p, src0, exs, (int*)src_bp, wrk_x, wrk_b, wrk_q, wrk_r, wrk_r0, wrk_p, pcg_p_, pcg_s, pcg_s_, pcg_t_, psz3d, psz4dex, ItrMax, true);
          TIMING_stop(&PM, "BiCGstab_w_precnd");
      break;
      
    case BICGSTAB:
          TIMING_start(&PM, "BiCGstab");
          loop = LS.PBiCGstab(p, src0, exs, (int*)src_bp, wrk_x, wrk_b, wrk_q, wrk_r, wrk_r0, wrk_p, pcg_p_, pcg_s, pcg_s_, pcg_t_, psz3d, psz4dex, ItrMax, false);
          TIMING_stop(&PM, "BiCGstab");
      break;

    default:
      assert(0);
      break;
  }

  // close
  if ( !fp ) fclose(fp);
  
  // file out
  strcpy(fname2, "p.sph");
  //fileout_(v_sz, &gc, psz4dex, p, &dh, org, fname2);

  strcpy(fname2, "e.sph");
  //fileout_(v_sz, &gc, gosa, &dh, org, fname2);

  // profiling
  if ( !(fp=fopen("data-profiling.txt", "w")) )
  {
    printf("\tSorry, can't open 'profiling.txt' file. Write failed.\n");
    assert(0);
  }
  
  // 測定結果の集計(gathreメソッドは全ノードで呼ぶこと)
  PM.gather();
  
  // 結果出力(排他測定のみ)                                                                         
  printf("\n===============================\n");
  PM.print(stdout, "hoge", "foo");
  PM.print(fp, "hoge", "foo");

  if ( !fp ) fclose(fp);
  
  
  return 0;
}

// #################################################################
/**
 * @brief タイミング測定区間にラベルを与えるラッパー
 * @param [in] label     ラベル
 * @param [in] type      測定対象タイプ(COMM or CALC)
 * @param [in] exclusive 排他測定フラグ(ディフォルトtrue)
 */
void set_label(PerfMonitor* pm, const string label, PerfMonitor::Type type, bool exclusive)
{
    // 登録個数のチェック
    order_of_PM_key++;
    
    if ( order_of_PM_key > PM_NUM_MAX )
    {
        printf("\tThe number of labels for Performance monitor goes over limit.\n");
        exit(0);
    }
    
    // 文字数がTM_LABEL_MAX-1を超えるものはカット
    if ( strlen(label.c_str()) > TM_LABEL_MAX-1 )
    {
        printf("\tWarning: Length of timing label must be less than %d\n", TM_LABEL_MAX-1);
    }
    
    // Performance Monitorへの登録
    pm->setProperties(label, type, exclusive);
}

//@brief タイミング測定区間にラベルを与える
void set_timing_label(PerfMonitor* pm)
{
    // 非排他, 計算
    set_label(pm, "Jacobi",           PerfMonitor::CALC, false);
    set_label(pm, "PointSOR",         PerfMonitor::CALC, false);
    set_label(pm, "SOR2(SMA)",        PerfMonitor::CALC, false);
    set_label(pm, "BiCGstab_w_precnd",PerfMonitor::CALC, false);
    set_label(pm, "BiCGstab",         PerfMonitor::CALC, false);
    set_label(pm, "BoundaryCondition",PerfMonitor::CALC, true);
    set_label(pm, "Norm_max",         PerfMonitor::CALC, true);
    
    set_label(pm, "Jacobi_kernel",    PerfMonitor::CALC, true);
    set_label(pm, "Sor_kernel",       PerfMonitor::CALC, true);
    set_label(pm, "Sor2sma_kernel",   PerfMonitor::CALC, true);

    set_label(pm, "Vector_change_p",   PerfMonitor::CALC, true);
    set_label(pm, "Vector_change_b",   PerfMonitor::CALC, true);
    set_label(pm, "Vector_change_exs",   PerfMonitor::CALC, true);
    
    set_label(pm, "Blas_dot1",        PerfMonitor::CALC, true);
    set_label(pm, "Blas_dot2",        PerfMonitor::CALC, true);
    set_label(pm, "Blas_copy",        PerfMonitor::CALC, true);
    set_label(pm, "Blas_clear",       PerfMonitor::CALC, true);
    set_label(pm, "Blas_residual",    PerfMonitor::CALC, true);
    set_label(pm, "Blas_bicg1",       PerfMonitor::CALC, true);
    set_label(pm, "Blas_bicg2",       PerfMonitor::CALC, true);
    set_label(pm, "Blas_ax",          PerfMonitor::CALC, true);
    set_label(pm, "Blas_triad",       PerfMonitor::CALC, true);
        
    // 共通にまとめて利用
    set_label(pm, "Copy_Array",             PerfMonitor::CALC, true);
    set_label(pm, "assign_Const_to_Array",  PerfMonitor::CALC, true);
}

