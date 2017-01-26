!###################################################################################
!
! FFV-C
! Frontflow / violet Cartesian
!
!
! Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
! All rights reserved.
!
! Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!###################################################################################

!> @file   ffv_blas.f90
!! @brief  BLAS routine
!! @author aics, iis
!<

! based on Onishi version

!> ********************************************************************
!! @brief 要素のゼロクリア
!! @param [in,out] x  スカラー
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
  subroutine blas_clear(x, sz, g, psz4dex, vec_cnt)
  implicit none
  integer                                                       ::  i, j, k, l, ix, jx, kx, lx, g, vec_cnt
  integer, dimension(4)                                         ::  sz, psz4dex
  real, dimension(vec_cnt+psz4dex(1), 1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4)) ::  x

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  lx = vec_cnt

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, g, lx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1-g,kx+g+psz4dex(4)
     do j=1-g,jx+g+psz4dex(3)
        do i=1-g,ix+g+psz4dex(2)
           do l=1,lx+psz4dex(1)
              x(l, i, j, k) = 0.0
           end do
        end do
     end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_clear


!> ********************************************************************
!! @brief コピー
!! @param [out]    y  コピー先
!! @param [in]     x  ソース
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!<
  subroutine blas_copy(y, x, sz, g, psz4dex, vec_cnt)
  implicit none
  integer                                                   ::  i, j, k, l, ix, jx, kx, lx, g, vec_cnt
  integer, dimension(4)                                     ::  sz, psz4dex
  real, dimension(vec_cnt+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4)) :: y, x

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  lx = vec_cnt

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, g, lx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1-g,kx+g+psz4dex(4)
     do j=1-g,jx+g+psz4dex(3)
        do i=1-g,ix+g+psz4dex(2)
           do l=1,lx+psz4dex(1)
              y(l, i, j, k) = x(l, i, j, k)
           end do
        end do
     end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_copy


!> ********************************************************************
!! @brief AXPYZ
!! @param [out]    z    ベクトル
!! @param [in]     y    ベクトル
!! @param [in]     x    ベクトル
!! @param [in]     a    係数
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル
!! @param [in]     flop 浮動小数点演算数
!<
subroutine blas_triad(z, x, y, a, sz, g, psz4dex, flop, vec_cnt)
implicit none
integer                                                   ::  i, j, k, l, ix, jx, kx, lx, g, vec_cnt
integer, dimension(4)                                     ::  sz, psz4dex
real, dimension(vec_cnt+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4))    ::  x, y, z
double precision                                          ::  flop
double precision,dimension(sz(4))                         ::  a

ix = sz(1)
jx = sz(2)
kx = sz(3)
lx = vec_cnt

flop = dble(ix) * dble(jx) * dble(kx) * dble(lx) * 2.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, a, lx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
do l=1,lx
   z(l, i, j, k) = a(l) * x(l, i, j, k) + y(l, i, j, k)
end do
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_triad


!> ********************************************************************
!! @brief BiCGstab 2
!! @param [in,out] z    ベクトル
!! @param [in]     y    ベクトル
!! @param [in]     x    ベクトル
!! @param [in]     a    係数
!! @param [in]     b    係数
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル
!! @param [in]     flop 浮動小数点演算数
!<
subroutine blas_bicg_2(z, x, y, a, b, sz, g, psz4dex, flop, vec_cnt)
implicit none
integer                                                       ::  i, j, k, l, ix, jx, kx, lx, g, vec_cnt
integer, dimension(4)                                         ::  sz, psz4dex
real, dimension(vec_cnt+psz4dex(1) ,1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4)) ::  x, y, z
double precision                                              ::  flop
double precision,dimension(sz(4))                             ::  a,b

ix = sz(1)
jx = sz(2)
kx = sz(3)
lx = vec_cnt

flop = dble(ix) * dble(jx) * dble(kx) * dble(lx) * 4.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, a, b, lx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
do l=1,lx
z(l,i, j, k) = a(l) * x(l,i, j, k) + b(l) * y(l,i, j, k) + z(l,i, j, k)
end do
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_bicg_2

!> ********************************************************************
!! @brief 圧力Poissonの定数項bの計算
!! @param [out] rhs  右辺ベクトルbの自乗和
!! @param [out] b    RHS vector b
!! @param [in]  s_u  \sum {u^*}
!! @param [in]  bp   BCindex P
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dh   格子幅
!! @param [in]  dt   時間積分幅
!! @param [out] flop flop count
!<
subroutine blas_calc_b (rhs, b, s_u, bp, sz, g, dh, dt, flop)
implicit none
include 'cbc_f_params.h'
integer                                                     ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                       ::  sz
double precision                                            ::  flop, rhs
real                                                        ::  dh, dt, c1, dv
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  s_u, b
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bp

ix = sz(1)
jx = sz(2)
kx = sz(3)
rhs = 0.0
c1 = dh / dt

flop = flop + dble(ix)*dble(jx)*dble(kx)*4.0d0 + 8.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:rhs) &
!$OMP PRIVATE(dv) &
!$OMP FIRSTPRIVATE(ix, jx, kx, c1)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix
  dv = c1 * s_u(i,j,k) * real(ibits(bp(i,j,k),Active,1))
  b(i,j,k) = dv ! \frac{h^2}{\Delta t} \nabla u^*
  rhs = rhs + dble(dv*dv)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_calc_b


!> ********************************************************************
!! @brief DOT1
!! @param [out] r    内積
!! @param [in]  p    ベクトル
!! @param [in]  bp   BCindex P
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル
!! @param [out] flop flop count
!<
subroutine blas_dot1(r, p, bp, sz, g, psz3d, psz4dex, flop, vec_cnt)
implicit none
include 'cbc_f_params.h'
integer                                                      ::  i, j, k, l, ix, jx, kx, lx, g, vec_cnt 
integer, dimension(4)                                        ::  sz, psz4dex
integer, dimension(3)                                        ::  psz3d
real, dimension(vec_cnt+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4)) ::  p
integer, dimension(1-g:sz(1)+g+psz3d(1), 1-g:sz(2)+g+psz3d(2), 1-g:sz(3)+g+psz3d(3))     ::  bp
double precision                                             ::  flop, q
double precision,dimension(sz(4))                            ::  r

ix = sz(1)
jx = sz(2)
kx = sz(3)
lx = vec_cnt
r  = 0.0

flop = flop + dble(ix)*dble(jx)*dble(kx)*dble(lx)*3.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:r) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx) &
!$OMP PRIVATE(q)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
do l=1,lx
   q = dble(p(l,i, j, k))
   r(l) = r(l) + q*q * dble(ibits(bp(i,j,k),Active,1))
end do
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_dot1


!> ********************************************************************
!! @brief DOT2
!! @param [out] r    内積
!! @param [in]  p    ベクトル
!! @param [in]  q    ベクトル
!! @param [in]  bp   BCindex P
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル
!! @param [out] flop flop count
!<
subroutine blas_dot2(r, p, q, bp, sz, g, psz3d, psz4dex, flop, vec_cnt)
  implicit none
  include 'cbc_f_params.h'
  integer                                                      ::  i, j, k, l, lx, ix, jx, kx, g, vec_cnt 
  integer, dimension(4)                                        ::  sz, psz4dex
  integer, dimension(3)                                        ::  psz3d
  real, dimension(vec_cnt+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4)) ::  p, q
  integer, dimension(1-g:sz(1)+g+psz3d(1), 1-g:sz(2)+g+psz3d(2), 1-g:sz(3)+g+psz3d(3))     ::  bp
  double precision                                             ::  flop
  double precision,dimension(sz(4))                            ::  r

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  lx = vec_cnt
  r  = 0.0
  
  flop = flop + dble(ix)*dble(jx)*dble(kx)*dble(lx)*3.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:r) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1,kx
  do j=1,jx
  do i=1,ix
  do l=1,lx
     r(l) = r(l) + dble(p(l,i, j, k)) * dble(q(l,i, j, k)) * dble(ibits(bp(i,j,k), Active, 1))
  end do
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_dot2


!> ********************************************************************
!! @brief AX
!! @param [out] ap   AX
!! @param [in]  p    解ベクトル
!! @param [in]  bp   BCindexP
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル
!! @param [out] flop flop count
!<
  subroutine blas_calc_ax(ap, p, bp, sz, g, psz3d, psz4dex, flop, vec_cnt)
  implicit none
  include 'cbc_f_params.h'
  integer                                                       ::  i, j, k, l, ix, jx, kx, lx, g, vec_cnt, idx 
  integer, dimension(4)                                         ::  sz, psz4dex
  integer, dimension(3)                                         ::  psz3d
  real                                                          ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
  real                                                          ::  dd, ss, d0, d1, d2, ac
  real, dimension(vec_cnt+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4))  ::  ap, p
  integer, dimension(      1-g:sz(1)+g+psz3d(1), 1-g:sz(2)+g+psz3d(2), 1-g:sz(3)+g+psz3d(3))  ::  bp
  double precision                                              ::  flop

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  lx = vec_cnt

  flop = flop + dble(ix)*dble(jx)*dble(kx)*dble(14.0d0*lx + 4.0d0)

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, d0, d1, d2, ac, idx) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1,kx
  do j=1,jx
  do i=1,ix
     idx = bp(i,j,k)
     ndag_e = real(ibits(idx, bc_ndag_E, 1))  ! e                                         
     ndag_w = real(ibits(idx, bc_ndag_W, 1))  ! w                                         
     ndag_n = real(ibits(idx, bc_ndag_N, 1))  ! n                                         
     ndag_s = real(ibits(idx, bc_ndag_S, 1))  ! s                                         
     ndag_t = real(ibits(idx, bc_ndag_T, 1))  ! t                                         
     ndag_b = real(ibits(idx, bc_ndag_B, 1))  ! b                                         
     
     d0 = real(ibits(idx, bc_diag + 0, 1))
     d1 = real(ibits(idx, bc_diag + 1, 1))
     d2 = real(ibits(idx, bc_diag + 2, 1))
     dd = d2*4.0 + d1*2.0 + d0  ! diagonal                                                
     
     ac = real(ibits(idx, Active, 1))     

     do l=1,lx
        ss =  ndag_e * p(l,i+1,j  ,k  ) &
             + ndag_w * p(l,i-1,j  ,k  ) &
             + ndag_n * p(l,i  ,j+1,k  ) &
             + ndag_s * p(l,i  ,j-1,k  ) &
             + ndag_t * p(l,i  ,j  ,k+1) &
             + ndag_b * p(l,i  ,j  ,k-1)
        ap(l,i, j, k) = (ss - dd * p(l,i, j, k)) * ac
     end do
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_calc_ax


!> ********************************************************************
!! @brief 残差ベクトルの計算
!! @param [out]    r  残差ベクトル
!! @param [in]     p  解ベクトル
!! @param [in]     b  定数項
!! @param [in]     bp BCindexP
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル
!! @param [out] flop flop count
!<
  subroutine blas_calc_rk(r, p, b, bp, sz, g, psz3d, psz4dex, flop, vec_cnt)
  implicit none
  include 'cbc_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx, l, lx, vec_cnt
  integer, dimension(4)                                     ::  sz, psz4dex
  integer, dimension(3)                                     ::  psz3d
  real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
  real                                                      ::  dd, ss, d0, d1, d2, ac
  real, dimension(vec_cnt+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4))    ::  r, p, b
  integer, dimension(1-g:sz(1)+g+psz3d(1), 1-g:sz(2)+g+psz3d(2), 1-g:sz(3)+g+psz3d(3)) ::  bp
  double precision                                          ::  flop

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  lx = vec_cnt

  flop = flop + dble(ix)*dble(jx)*dble(kx)*dble(lx*15.0d0 + 4.0d0)

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, idx, d0, d1, d2, ac) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1,kx
     do j=1,jx
        do i=1,ix
           
           idx = bp(i,j,k)
           ndag_e = real(ibits(idx, bc_ndag_E, 1))  ! e                                         
           ndag_w = real(ibits(idx, bc_ndag_W, 1))  ! w                                         
           ndag_n = real(ibits(idx, bc_ndag_N, 1))  ! n                                         
           ndag_s = real(ibits(idx, bc_ndag_S, 1))  ! s                                         
           ndag_t = real(ibits(idx, bc_ndag_T, 1))  ! t                                         
           ndag_b = real(ibits(idx, bc_ndag_B, 1))  ! b                                         
           
           d0 = real(ibits(idx, bc_diag + 0, 1))
           d1 = real(ibits(idx, bc_diag + 1, 1))
           d2 = real(ibits(idx, bc_diag + 2, 1))
           dd = d2*4.0 + d1*2.0 + d0  ! diagonal                                                
           
           ac = real(ibits(idx, Active, 1))
                      
           do l=1,lx
              ss =  ndag_e * p(l,i+1,j  ,k  ) &
                   + ndag_w * p(l,i-1,j  ,k  ) &
                   + ndag_n * p(l,i  ,j+1,k  ) &
                   + ndag_s * p(l,i  ,j-1,k  ) &
                   + ndag_t * p(l,i  ,j  ,k+1) &
                   + ndag_b * p(l,i  ,j  ,k-1)
              r(l,i, j, k) = (b(l,i, j, k) - (ss - dd * p(l,i, j, k))) * ac
           end do
        end do
     end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine blas_calc_rk


!> ********************************************************************
!! @brief 残差の自乗和のみ
!! @param [out] res  残差の自乗和
!! @param [in]  p    圧力
!! @param [in]  b    RHS vector
!! @param [in]  bp   BCindex P
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [out] flop flop count
!<
subroutine blas_calc_r2 (res, p, b, bp, sz, g, flop)
implicit none
include 'cbc_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop, res
real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
real                                                      ::  dd, ss, dp, d0, d1, d2
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, b
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g,8) ::  bp

ix = sz(1)
jx = sz(2)
kx = sz(3)
res = 0.0

flop = flop + dble(ix)*dble(jx)*dble(kx)*17.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, dp, idx, d0, d1, d2) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix
ndag_e = bp(i,j,k,1)! e, non-diagonal
ndag_w = bp(i,j,k,2)! w
ndag_n = bp(i,j,k,3)! n
ndag_s = bp(i,j,k,4)! s
ndag_t = bp(i,j,k,5)! t
ndag_b = bp(i,j,k,6)! b
dd = bp(i,j,k,7) ! diagonal

ss = ndag_e * p(i+1,j  ,k  ) &
   + ndag_w * p(i-1,j  ,k  ) &
   + ndag_n * p(i  ,j+1,k  ) &
   + ndag_s * p(i  ,j-1,k  ) &
   + ndag_t * p(i  ,j  ,k+1) &
   + ndag_b * p(i  ,j  ,k-1)
dp = ( b(i,j,k) - (ss - dd * p(i,j,k)) ) * real(bp(i,j,k,8))
res = res + dble(dp*dp)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_calc_r2

!> ********************************************************************
!! @brief BiCGstabの部分演算1
!! @param [in,out] p    ベクトル
!! @param [in]     r    ベクトル
!! @param [in]     q    ベクトル
!! @param [in]     beta 係数
!! @param [in]     omg  係数
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル
!! @param [in]     flop 浮動小数点演算数
!<
subroutine blas_bicg_1(p, r, q, beta, omg, sz, g, psz4dex, flop, vec_cnt)
implicit none
integer                                                      ::  i, j, k, l, ix, jx, kx, lx, g, vec_cnt
integer, dimension(4)                                        ::  sz, psz4dex
real, dimension(vec_cnt+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4)) ::  p, r, q
double precision                                             ::  flop
double precision,dimension(sz(4))                            ::  beta, omg

ix = sz(1)
jx = sz(2)
kx = sz(3)
lx = vec_cnt

flop = dble(ix) * dble(jx) * dble(kx) * dble(lx) * 4.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, beta, omg, lx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
do l=1,lx
  p(l,i,j,k) = r(l,i,j,k) + beta(l) * ( p(l,i,j,k) - omg(l) * q(l,i,j,k) )
end do
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine blas_bicg_1
