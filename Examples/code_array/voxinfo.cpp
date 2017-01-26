/*
 *  voxinfo.cpp
 *  kernel_test
 *
 *  Created by keno on 11/12/05.
 *  Copyright 2011 __iis__. All rights reserved.
 *
 */

#include <assert.h>
#include <stdio.h>
#include "voxinfo.h"

void VoxInfo::setBCIndexP(int* m_sz, int gc, int* pad_size, unsigned* bcp)
{
  // 初期化 @note ビットを1に初期化する．初期化範囲はガイドセルを含む全領域．セルフェイスの射影処理で必要．
  int mx = (m_sz[0]+2*gc+pad_size[0]) * (m_sz[1]+2*gc+pad_size[1]) * (m_sz[2]+2*gc+pad_size[2]);
  for (int m=0; m<mx; m++) {
    bcp[m] |= ( 0x3ffff << BC_NDAG_W ); // BC_NDAG_W〜BC_D_Tまで18bitまとめて1に初期化
  }

  // 状態のエンコード
  printf("\tfluid cell = %d\n", encAS(m_sz, gc, pad_size, bcp));
  
  // 外部境界のビットフラグをエンコード
  encPbit_OBC(m_sz, gc, pad_size, X_MINUS, bcp, "Dirichlet", true);
  encPbit_OBC(m_sz, gc, pad_size, X_PLUS,  bcp, "Dirichlet", true);
  encPbit_OBC(m_sz, gc, pad_size, Y_MINUS, bcp, "Dirichlet", true);
  encPbit_OBC(m_sz, gc, pad_size, Y_PLUS,  bcp, "Dirichlet", true);
  encPbit_OBC(m_sz, gc, pad_size, Z_MINUS, bcp, "Dirichlet", true);
  encPbit_OBC(m_sz, gc, pad_size, Z_PLUS,  bcp, "Dirichlet", true);
  
  // 全周Neumannフラグのセルと排他性をチェックし，反復行列の非対角要素/対角要素をエンコードする
  encPbit(m_sz, gc, pad_size, bcp);
}

/**
 @fn void VoxInfo::encPbit_OBC(int face, unsigned* bx, string key, bool dir)
 @brief 外部境界に接するセルにおいて，bx[]に圧力境界条件keyに対応するビットフラグを設定する
 @param face 外部境界面番号
 @param bx BCindex P
 @param key Dirichlet or Neumann
 @param dir 壁面の場合(true)，方向フラグをON
 @retval 固体表面セル数
 @note
 - 流体セルに対してのみ，1-Normal, 0-BC
 - 固体セルに隣接する面のノイマンフラグをゼロにし，方向フラグを立てる
 */
void VoxInfo::encPbit_OBC(int* m_sz, int gc, int* pad_size, int face, unsigned* bx, std::string key, bool dir)
{
  int i,j,k, ix, jx, kx;
  unsigned m;
  unsigned s;
  
  ix = m_sz[0];
  jx = m_sz[1];
  kx = m_sz[2];
  int ip = pad_size[0];
  int jp = pad_size[1];
  int kp = pad_size[2];
  
  switch (face) {
      case X_MINUS:
      i = 1;
      if ("Neumann"==key) {
        for (k=1; k<=kx; k++) {
          for (j=1; j<=jx; j++) {
            m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
            s = bx[m];
            if ( IS_FLUID(s) ) {
              bx[m] = offBit( s, BC_N_W );
            }
          }
        }
      }
      else {
        for (k=1; k<=kx; k++) {
          for (j=1; j<=jx; j++) {
            m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
            s = bx[m];
            if ( IS_FLUID(s) ) {
              bx[m] = offBit( s, BC_D_W );
            }
          }
        }
      }
      break;
      
      case X_PLUS:
      i = ix;
      if ("Neumann"==key) {
        for (k=1; k<=kx; k++) {
          for (j=1; j<=jx; j++) {
            m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
            s = bx[m];
            if ( IS_FLUID(s) ) {
              bx[m] = offBit( s, BC_N_E );
            }
          }
        }
      }
      else {
        for (k=1; k<=kx; k++) {
          for (j=1; j<=jx; j++) {
            m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
            s = bx[m];
            if ( IS_FLUID(s) ) {
              bx[m] = offBit( s, BC_D_E );
            }
          }
        }
      }
      break;
      
      case Y_MINUS:
      j = 1;
      if ("Neumann"==key) {
        for (k=1; k<=kx; k++) {
          for (i=1; i<=ix; i++) {
            m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
            s = bx[m];
            if ( IS_FLUID(s) ) {
              bx[m] = offBit( s, BC_N_S );
            }
          }
        }
      }
      else {
        for (k=1; k<=kx; k++) {
          for (i=1; i<=ix; i++) {
            m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
            s = bx[m];
            if ( IS_FLUID(s) ) {
              bx[m] = offBit( s, BC_D_S );
            }
          }
        }
      }
      break;
      
      case Y_PLUS:
      j = jx;
      if ("Neumann"==key) {
        for (k=1; k<=kx; k++) {
          for (i=1; i<=ix; i++) {
            m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
            s = bx[m];
            if ( IS_FLUID(s) ) {
              bx[m] = offBit( s, BC_N_N );
            }
          }
        }
      }
      else {
        for (k=1; k<=kx; k++) {
          for (i=1; i<=ix; i++) {
            m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
            s = bx[m];
            if ( IS_FLUID(s) ) {
              bx[m] = offBit( s, BC_D_N );
            }
          }
        }
      }
      break;
      
      case Z_MINUS:
      k = 1;
      if ("Neumann"==key) {
        for (j=1; j<=jx; j++) {
          for (i=1; i<=ix; i++) {
            m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
            s = bx[m];
            if ( IS_FLUID(s) ) {
              bx[m] = offBit( s, BC_N_B );
            }
          }
        }
      }
      else {
        for (j=1; j<=jx; j++) {
          for (i=1; i<=ix; i++) {
            m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
            s = bx[m];
            if ( IS_FLUID(s) ) {
              bx[m] = offBit( s, BC_D_B );
            }
          }
        }
      }
      break;
      
      case Z_PLUS:
      k = kx;
      if ("Neumann"==key) {
        for (j=1; j<=jx; j++) {
          for (i=1; i<=ix; i++) {
            m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
            s = bx[m];
            if ( IS_FLUID(s) ) {
              bx[m] = offBit( s, BC_N_T );
            }
          }
        }
      }
      else {
        for (j=1; j<=jx; j++) {
          for (i=1; i<=ix; i++) {
            m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
            s = bx[m];
            if ( IS_FLUID(s) ) {
              bx[m] = offBit( s, BC_D_T );
            }
          }
        }
      }
      break;
  } // end of switch
}

/**
 @fn void VoxInfo::encPbit(int* m_sz, int gc, unsigned* bx)
 @brief ディリクレ条件とノイマン条件の排他性をチェックし，反復行列の非対角要素/対角要素の係数をエンコードする
 @param bx BCindex P
 @note
 - ディリクレ条件とノイマン条件の排他性のチェック
 - 非対角要素と対角要素の係数をエンコードする
 */
void VoxInfo::encPbit(int* m_sz, int gc, int* pad_size, unsigned* bx)
{
  int i, j, k, ix, jx, kx;
  unsigned m, flag, nx;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b, ss;
  unsigned d_e, d_w, d_n, d_s, d_t, d_b;
  unsigned s;
  bool exclusive;
  
  ix = m_sz[0];
  jx = m_sz[1];
  kx = m_sz[2];
  int ip = pad_size[0];
  int jp = pad_size[1];
  int kp = pad_size[2];
  
  // ディリクレ条件とノイマン条件の排他性のチェック
  exclusive = true;
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
        s = bx[m];
        flag = 0;
        
        // ノイマン条件：値がゼロのとき，BCがセットされている
        s_e = BIT_SHIFT(s, BC_N_E);
        s_w = BIT_SHIFT(s, BC_N_W);
        s_n = BIT_SHIFT(s, BC_N_N);
        s_s = BIT_SHIFT(s, BC_N_S);
        s_t = BIT_SHIFT(s, BC_N_T);
        s_b = BIT_SHIFT(s, BC_N_B);
        
        // ディリクレ条件：値がゼロのとき，BCがセットされている
        d_e = BIT_SHIFT(s, BC_D_E);
        d_w = BIT_SHIFT(s, BC_D_W);
        d_n = BIT_SHIFT(s, BC_D_N);
        d_s = BIT_SHIFT(s, BC_D_S);
        d_t = BIT_SHIFT(s, BC_D_T);
        d_b = BIT_SHIFT(s, BC_D_B);
        
        // ノイマンのときに非対角要素の係数をエンコード
        if ( (s_e * d_e) == 0 ) s = offBit( s, BC_NDAG_E );
        if ( (s_w * d_w) == 0 ) s = offBit( s, BC_NDAG_W );
        if ( (s_n * d_n) == 0 ) s = offBit( s, BC_NDAG_N );
        if ( (s_s * d_s) == 0 ) s = offBit( s, BC_NDAG_S );
        if ( (s_t * d_t) == 0 ) s = offBit( s, BC_NDAG_T );
        if ( (s_b * d_b) == 0 ) s = offBit( s, BC_NDAG_B );
        
        bx[m] = s;
        
        if ( (s_e==0) && (d_e==0) ) flag++;
        if ( (s_w==0) && (d_w==0) ) flag++;
        if ( (s_n==0) && (d_n==0) ) flag++;
        if ( (s_s==0) && (d_s==0) ) flag++;
        if ( (s_t==0) && (d_t==0) ) flag++;
        if ( (s_b==0) && (d_b==0) ) flag++;
        
        if ( flag != 0) {
          printf("\tDirichlet and Neumann BC are specified on the same face in cell (%d,%d,%d)\n", i,j,k);
          exclusive = false;
        }
      }
    }
  }
  if ( !exclusive ) assert(0);
  
  // 対角要素の係数のチェックとエンコード
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
        s = bx[m];
        
        s_e = BIT_SHIFT(s, BC_N_E); // 対角要素の係数：0 or 1
        s_w = BIT_SHIFT(s, BC_N_W);
        s_n = BIT_SHIFT(s, BC_N_N);
        s_s = BIT_SHIFT(s, BC_N_S);
        s_t = BIT_SHIFT(s, BC_N_T);
        s_b = BIT_SHIFT(s, BC_N_B);
        
        ss = s_e + s_w + s_n + s_s + s_t + s_b;
        bx[m] = s | (ss<<BC_DIAG);
        
        if ( ss == 0 ) {
          printf("\tError : Coefficient of diagonal element is zero at (%d,%d,%d) : (wesnbt)[%1d %1d %1d %1d %1d %1d]\n", i,j,k,
                 IS_FLUID(bx[ F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i-1, j, k) ]),
                 IS_FLUID(bx[ F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i+1, j, k) ]),
                 IS_FLUID(bx[ F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j-1, k) ]),
                 IS_FLUID(bx[ F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j+1, k) ]),
                 IS_FLUID(bx[ F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k-1) ]),
                 IS_FLUID(bx[ F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k+1) ]) );
        }
        
      }
    }
  }
  
  // ゼロ割防止のためのダミー係数 >> 全領域
  nx = (m_sz[0]+2*gc+ip) * (m_sz[1]+2*gc+kp) * (m_sz[2]+2*gc+kp);
  for (m=0; m<nx; m++) {
    s = bx[m];
    if ( ((s>>BC_DIAG) & 0x7) == 0 ) { // 0x7 = 3 bit
      bx[m] = s | (0x1<<BC_DIAG);
    }
  }
}

/// Active, State bitを設定，内部のみ有効
unsigned VoxInfo::encAS(int* m_sz, int gc, int* pad_size, unsigned* bx)
{
  int i,j,k,ix,jx,kx;
  unsigned m, c=0;
  unsigned s;
  
  ix = m_sz[0];
  jx = m_sz[1];
  kx = m_sz[2];
  int ip = pad_size[0];
  int jp = pad_size[1];
  int kp = pad_size[2];
  
  // all off
  int mx = (ix+2*gc+ip) * (jx+2*gc+jp) * (kx+2*gc+kp);
  for (int m=0; m<mx; m++) {
    s = bx[m];
    s = offBit( s, ACTIVE_BIT ); // Inactive
    s = offBit( s, STATE_BIT  ); // Solid
    bx[m] = s;
  }
  
  // inseide on
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m = F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k);
        s = bx[m];
        s = onBit( s, ACTIVE_BIT ); // Active
        s = onBit( s, STATE_BIT  ); // Fluid
        c++;
        bx[m] = s;
      }
    }
  }
  return c;
}


