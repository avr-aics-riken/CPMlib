/*
 *  FBdefine.h
 *  kernel_test
 *
 *  Created by keno on 11/12/05.
 *  Copyright 2011 __iis__. All rights reserved.
 *
 */

#ifndef _FB_DEFINE_H_
#define _FB_DEFINE_H_

#include "omp.h"
#include <float.h>

//@file FBDefine.h
//@brief FlowBase Definition Header
//@author keno, FSI Team, VCAD, RIKEN

/** 実数型の指定
 * - デフォルトでは、REAL_TYPE=float
 * - コンパイル時オプション-D_REAL_IS_DOUBLE_を付与することで
 *   REAL_TYPE=doubleになる
 */
#ifdef _REAL_IS_DOUBLE_
#define REAL_TYPE double
#define REAL_TYPE_EPSILON DBL_MIN
#else

#define REAL_TYPE float
#define REAL_TYPE_EPSILON FLT_MIN
#endif

// PMlibの登録ラベル個数
#define PM_NUM_MAX 50
#define TM_LABEL_MAX 24

// general
#define NOFACE      6
#define LABEL       64
#define ON          1
#define OFF         0
#define DETAIL      2

// 外部境界条件
#define OBC_MASK      31 // 外部境界と内部境界の識別子

// 面の番号
#define X_MINUS 0
#define X_PLUS  1
#define Y_MINUS 2
#define Y_PLUS  3
#define Z_MINUS 4
#define Z_PLUS  5


// エンコードビット　共通
#define ACTIVE_BIT 31
#define STATE_BIT  30

// エンコードビット　ID
#define TOP_CMP_ID    0  //  コンポーネントの先頭ビット
#define TOP_MATERIAL  6  //  MATERIALの先頭ビット
#define TOP_CELL_ID   12 //  IDの先頭ビット
#define TOP_VF        20 //  Volume Fractionの先頭ビット
#define FORCING_BIT   28 //  外力モデルの識別子

// マスクのビット幅
#define MASK_CMP_ID   0x3f // 6 bit幅
#define MASK_MAT      0x3f // 6 bit幅
#define MASK_CELL_ID  0xff // 8 bit幅
#define MASK_VF       0xff // 8 bit幅
#define MASK_5        0x1f // 5 bit幅

// エンコードビット　P
#define BC_D_T     29
#define BC_D_B     28
#define BC_D_N     27
#define BC_D_S     26
#define BC_D_E     25
#define BC_D_W     24

#define BC_N_T     23
#define BC_N_B     22
#define BC_N_N     21
#define BC_N_S     20
#define BC_N_E     19
#define BC_N_W     18

#define BC_NDAG_T  17
#define BC_NDAG_B  16
#define BC_NDAG_N  15
#define BC_NDAG_S  14
#define BC_NDAG_E  13
#define BC_NDAG_W  12

#define BC_DIAG    9


// state
#define SOLID      0
#define FLUID      1


// 判定マクロ
// BCindex aの状態が流体であればtrueを返す (uint a)
#define IS_FLUID(a) ( ((a >> STATE_BIT) & 0x1) ? true : false )

// aをbだけ右シフトしてデコードする (uint a, b)
#define BIT_SHIFT(a,b) ( (a >> b) & 0x1 )

// コンポーネントエントリを返す (uint a)
#define DECODE_CMP(a) ( (a >> TOP_CMP_ID) & MASK_CMP_ID )

// ID番号を返す (uint a)
#define DECODE_ID(a) ( (a >> TOP_CELL_ID) & MASK_CELL_ID )

// MaterialListへのエントリを返す (uint a)
#define DECODE_MAT(a) ( (a >> TOP_MATERIAL) & MASK_MAT )

// Volume Fraction[0-255]を返す (uint a)
#define DECODE_VF(a) ( (a >> TOP_VF) & MASK_VF )

// BCindex aの第bビットがONかどうかを調べ，ONのとき，trueを返す
#define BIT_IS_SHIFT(a,b) ( ( (a >> b) & 0x1 ) ? true : false )

// BCindex aの第bビットをREALにキャストして返す
#define GET_SHIFT_F(a,b) ( (REAL_TYPE)( (a>>b) & 0x1 ) )

// BCindexにエンコードされたFaceBCのインデクスを返す
#define GET_FACE_BC(a,b) ( (a>>b) & MASK_5 )

// BCindexのセルの6面のいずれかにBCが設定されている場合，true
#define IS_INCLUDE_BC(s) ( (s & 0x3fffffff) != 0 )

// BCindex aの第bビットがONかどうかを調べ，ONのとき，trueを返す
#define TEST_BIT(a,b) ( ( (a >> b) & 0x1 ) ? true : false )

#define F_INDEX_S3D(ix, jx, kx, gc, ip, jp, kp, i, j, k) ( ((ix)+(gc)*2+(ip))*((jx)+(gc)*2+(jp))*((k)+(gc)-1) + ((ix)+(gc)*2+(ip))*((j)+(gc)-1) + (i)+(gc)-1 )

#endif // _FB_DEFINE_H_
