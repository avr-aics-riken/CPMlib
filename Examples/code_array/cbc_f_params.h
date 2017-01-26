!
!  cbc_f_params.h
!1;2c  kernel_test
!
!  Created by keno on 11/12/05.
!  Copyright 2011 __iis__. All rights reserved.
!
!

!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
!
!   *********************************************************
!
!> @file cbc_f_params.h
!! @brief BC Indexのマスク操作のビット定義, FBDefine.hと整合
!! @author keno, FSI Team, VCAD, RIKEN
!<

integer     ::  Active, State, id_fluid, id_solid
integer     ::  bc_diag, bc_mask30, bitw_5
integer     ::  bc_d_T, bc_d_B, bc_d_N, bc_d_S, bc_d_E, bc_d_W
integer     ::  bc_n_T, bc_n_B, bc_n_N, bc_n_S, bc_n_E, bc_n_W
integer     ::  bc_ndag_T, bc_ndag_B, bc_ndag_N, bc_ndag_S, bc_ndag_E, bc_ndag_W
integer     ::  X_plus, X_minus, Y_plus, Y_minus, Z_plus, Z_minus
integer     ::  obc_mask
integer     ::  Vector_num
integer     ::  Grid_size
integer     ::  Gaid_cell

parameter ( bc_mask30 = Z'3fffffff') ! 16進表記，VBCの6面(30bit)をまとめたマスク
parameter ( bitw_5   = 5)  ! FBDefine.h  5bit幅

parameter ( X_minus = 0 ) ! FBDefine.h X_MINUS
parameter ( X_plus  = 1 ) ! FBDefine.h X_PLUS
parameter ( Y_minus = 2 ) ! FBDefine.h Y_MINUS
parameter ( Y_plus  = 3 ) ! FBDefine.h Y_PLUS
parameter ( Z_minus = 4 ) ! FBDefine.h Z_MINUS
parameter ( Z_plus  = 5 ) ! FBDefine.h Z_PLUS

! 状態
parameter ( id_solid  = 0 ) ! FBDefine.h SOLID
parameter ( id_fluid  = 1 ) ! FBDefine.h FLUID

! ビットフラグ共通
parameter ( Active  = 31 ) ! FBDefine.h ACTIVE_BIT
parameter ( State   = 30 ) ! FBDefine.h STATE_BIT

! 圧力ビットフラグ
parameter ( bc_d_T = 29 ) ! FBDefine.h BC_D_T
parameter ( bc_d_B = 28 ) ! FBDefine.h BC_D_B
parameter ( bc_d_N = 27 ) ! FBDefine.h BC_D_N
parameter ( bc_d_S = 26 ) ! FBDefine.h BC_D_S
parameter ( bc_d_E = 25 ) ! FBDefine.h BC_D_E
parameter ( bc_d_W = 24 ) ! FBDefine.h BC_D_W

parameter ( bc_n_T = 23 ) ! FBDefine.h BC_N_T
parameter ( bc_n_B = 22 ) ! FBDefine.h BC_N_B
parameter ( bc_n_N = 21 ) ! FBDefine.h BC_N_N
parameter ( bc_n_S = 20 ) ! FBDefine.h BC_N_S
parameter ( bc_n_E = 19 ) ! FBDefine.h BC_N_E
parameter ( bc_n_W = 18 ) ! FBDefine.h BC_N_W

parameter ( bc_ndag_T = 17 ) ! FBDefine.h BC_NDAG_T
parameter ( bc_ndag_B = 16 ) ! FBDefine.h BC_NDAG_B
parameter ( bc_ndag_N = 15 ) ! FBDefine.h BC_NDAG_N
parameter ( bc_ndag_S = 14 ) ! FBDefine.h BC_NDAG_S
parameter ( bc_ndag_E = 13 ) ! FBDefine.h BC_NDAG_E
parameter ( bc_ndag_W = 12 ) ! FBDefine.h BC_NDAG_W

parameter ( bc_diag = 9 ) ! FBDefine.h BC_DIAG

! 外部境界条件番号
parameter ( obc_mask     = 31) ! FBDefine.h OBC_MASK 内外部境界条件の識別子

parameter ( Gaid_cell   = 2 ) ! FBDefine.h Gaid_cell

parameter ( Grid_size   = 128 ) ! FBDefine.h Grid_size
parameter ( Vector_num   = 1 ) ! FBDefine.h Vector_num

