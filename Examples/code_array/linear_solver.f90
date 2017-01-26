!
!	linear_solver.f90
!	kernel_test
!
!	Created by keno on 11/12/05.
!	Copyright 2011 __iis__. All rights reserved.
!

!> ********************************************************************
!! @brief point SOR法
!! @param [in,out] p    圧力
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     omg  加速係数
!! @param [in]     b    RHS vector
!! @param [in]     bp   BCindex P
!! @param [out]    res  residual
!! @param [out]    flop flop count
!! @note Activeマスクの位置は，固体中のラプラス式を解くように，更新式にはかけず残差のみにする
!<
subroutine psor (p, sz, g, psz3d, psz4dex, omg, b, bp, res, flop, vec_cnt)
implicit none
include 'cbc_f_params.h'
integer                                                      ::  i, j, k, ix, jx, kx, l, lx, g, vec_cnt, idx
integer, dimension(4)                                        ::  sz, psz4dex
integer, dimension(3)                                        ::  psz3d
double precision,dimension(sz(4))                            ::  res
double precision                                             ::  flop, de
real                                                         ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
real                                                         ::  dd, pp, bb, ss, dp, pn, ac, omg, dx, d0, d1, d2
real,    dimension(vec_cnt+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4)) :: p, b
integer, dimension(          1-g:sz(1)+g+psz3d(1), 1-g:sz(2)+g+psz3d(2), 1-g:sz(3)+g+psz3d(3)) :: bp

ix = sz(1)
jx = sz(2)
kx = sz(3)
lx = vec_cnt

res = 0.0
flop = flop + dble(ix)*dble(jx)*dble(kx)*(12.0 + dble(lx)*23.0)

!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b) &
!$OMP PRIVATE(dd, pp, bb, ss, dp, de, pn, ac, dx, d0, d1, d2, idx) &
!$OMP FIRSTPRIVATE(kx,jx,ix,lx,omg)  

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
         
         dx = 1.0 / dd
         
         do l=1,lx
            pp = p(l,i,j,k)
            bb = b(l,i,j,k)
            
            ss = ndag_e * p(l,i+1,j  ,k  ) &
                 + ndag_w * p(l,i-1,j  ,k  ) &
                 + ndag_n * p(l,i  ,j+1,k  ) &
                 + ndag_s * p(l,i  ,j-1,k  ) &
                 + ndag_t * p(l,i  ,j  ,k+1) &
                 + ndag_b * p(l,i  ,j  ,k-1)
            dp = ( (ss - bb)*dx - pp ) * omg
            pn = pp + dp*ac
            p(l,i,j,k) = pn
            
            de  = dble( bb - (ss - pn * dd) )
            res(l) = res(l) + de*de * ac
         end do
         
      end do
   end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine psor

!  **********************************************************************
!> @subroutine jacobi (p, sz, g, omg, b, bp, res, wk2, flop)
!! @brief 緩和Jacobi法
!! @param [in,out] p    圧力
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     omg  加速係数
!! @param [in]     b    RHS vector
!! @param [in]     bp   BCindex P
!! @param [out]    res  residual
!! @param wk2 ワーク用配列
!! @param bp BCindex P
!! @param[out] flop flop count
!<
subroutine jacobi (p, sz, g, psz3d, psz4dex, omg, b, bp, res, wk2, flop, vec_cnt)
implicit none
include 'cbc_f_params.h'
integer                                                       ::  i, j, k, l, ix, jx, kx, lx, g, vec_cnt, idx
integer, dimension(4)                                         ::  sz, psz4dex
integer, dimension(3)                                         ::  psz3d
double precision                                              ::  flop, de
double precision, dimension(sz(4))                            ::  res
real                                                          ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
real                                                          ::  omg, dd, ss, dp, pp, bb, pn, ac, dx, d0, d1, d2
real, dimension(sz(4)+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4)) :: p, b, wk2
integer, dimension(   1-g:sz(1)+g+psz3d(1), 1-g:sz(2)+g+psz3d(2), 1-g:sz(3)+g+psz3d(3)) :: bp


ix = sz(1)
jx = sz(2)
kx = sz(3)
lx = vec_cnt

res = 0.0 ! absolute
flop = flop + dble(ix*jx*kx*(8.0 + lx*23.0))

!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b) &
!$OMP PRIVATE(dd, ss, dp, pp, bb, pn, ac, dx, de, idx, d0, d1, d2) &
!$OMP FIRSTPRIVATE(omg,kx,jx,ix,lx)

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

         dx = 1.0 / dd
         
         do l=1,lx
            pp = p(l,i,j,k)
            bb = b(l,i,j,k)
            
            ss = ndag_e * p(l,i+1,j  ,k  ) &
                 + ndag_w * p(l,i-1,j  ,k  ) &
                 + ndag_n * p(l,i  ,j+1,k  ) &
                 + ndag_s * p(l,i  ,j-1,k  ) &
                 + ndag_t * p(l,i  ,j  ,k+1) &
                 + ndag_b * p(l,i  ,j  ,k-1)
            dp = ( (ss - bb)*dx - pp ) * omg
            pn = pp + dp*ac
            wk2(l,i,j,k) = pn
            
            de  = dble( bb - (ss - pn * dd) )
            res(l) = res(l) + de*de *  ac
         end do
      end do
   end do
end do
!$OMP END DO

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
   do j=1,jx
      do i=1,ix
         do l=1,lx
            p(l,i,j,k)=wk2(l,i,j,k)
         end do
      end do
   end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine jacobi


!> ********************************************************************
!! @brief 2-colored SOR法 stride memory access
!! @param [in,out] p     圧力
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     ip    開始点インデクス
!! @param [in]     color グループ番号
!! @param [in]     omg   加速係数
!! @param [in]     b     RHS vector
!! @param [in]     bp    BCindex P
!! @param [out]    res  residual
!! @param [out]    flop  浮動小数演算数
!! @note resは積算
!<
subroutine psor2sma_core (p, sz, g, psz3d, psz4dex, ip, color, omg, b, bp, res, flop, vec_cnt)
implicit none
include 'cbc_f_params.h'
integer                                                      ::  i, j, k, l, ix, jx, kx, lx, g, vec_cnt, idx
integer, dimension(4)                                        ::  sz, psz4dex
integer, dimension(3)                                        ::  psz3d
double precision                                             ::  flop, de
double precision, dimension(sz(4))                           ::  res
real                                                         ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
real                                                         ::  omg, dd, ss, dp, pp, bb, pn, ac, dx, d0, d1, d2
real, dimension(sz(4)+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4)) :: p, b
integer, dimension(      1-g:sz(1)+g+psz3d(1), 1-g:sz(2)+g+psz3d(2), 1-g:sz(3)+g+psz3d(3)) :: bp
integer                                                      ::  ip, color

ix = sz(1)
jx = sz(2)
kx = sz(3)
lx = vec_cnt

flop = flop + dble(ix)*dble(jx)*dble(kx)*(12.0 + dble(lx)*23.0) * 0.5d0
!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b) &
!$OMP PRIVATE(dd, ss, dp, pp, bb, pn, ac, dx, de, idx, d0, d1, d2) &
!$OMP FIRSTPRIVATE(color, ip, omg,kx,jx,ix,lx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1+mod(k+j+color+ip,2), ix, 2
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

   dx = 1.0 / dd
    
   do l=1,lx
      pp = p(l,i,j,k)
      bb = b(l,i,j,k)
      
      ss = ndag_e * p(l,i+1,j  ,k  ) &
           + ndag_w * p(l,i-1,j  ,k  ) &
           + ndag_n * p(l,i  ,j+1,k  ) &
           + ndag_s * p(l,i  ,j-1,k  ) &
           + ndag_t * p(l,i  ,j  ,k+1) &
           + ndag_b * p(l,i  ,j  ,k-1)
      dp = ( (ss - bb)*dx - pp ) * omg
      pn = pp + dp*ac
      p(l,i,j,k) = pn
      
      de  = dble( bb - (ss - pn * dd) )
      res(l) = res(l) + de*de * ac
   end do

end do 
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine psor2sma_core

!> ********************************************************************
!! @brief Dirichlet source
!! @param [in]     b     RHS vector
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     bp    BCindex P
!! @param [out]    dh    gird width
!! @note resは積算
!<
subroutine src_dirichlet (b, sz, g, psz3d, psz4dex, bp, dh)
implicit none
include 'cbc_f_params.h'
integer                                                      ::  i, j, k, ix, jx, kx, g, l, lx
integer, dimension(4)                                        ::  sz, psz4dex
integer, dimension(3)                                        ::  psz3d
real                                                         ::  dh
real                                                         ::  x, y, pi, z, x2, y2, z2,dh2
real, dimension(sz(4)+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4)) :: b
integer, dimension(1-g:sz(1)+g+psz3d(1), 1-g:sz(2)+g+psz3d(2), 1-g:sz(3)+g+psz3d(3)) :: bp


ix = sz(1)
jx = sz(2)
kx = sz(3)
lx = sz(4)
pi = 2.0*asin(1.0)
dh2 = dh*dh

do l=1,lx
if(mod(l,64)==1) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*sin(pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==2) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.1*sin(pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==3) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*sin(pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==4) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*sin(0.5*pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*sin(0.5*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==5) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.1*sin(0.1*pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*sin(0.1*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==6) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*sin(0.5*pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*sin(0.5*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==7) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*x*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==8) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.1*x*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==9) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*x*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==10) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*x*x*y*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==11) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*2.0*x*x*y*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*2.0*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==12) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*x*x*y*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==13) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.1*(x+y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*(x+y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==14) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.2*(x+y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.2*(x+y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==15) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*abs(sin(4.0*pi*x))
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*abs(sin(4.0*pi*x))
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==16) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*sin(pi*x)*sin(pi*x)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*sin(pi*x)*sin(pi*x)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==17) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==18) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==19) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==20) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*sin(0.5*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==21) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*sin(0.1*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==22) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*sin(0.5*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==23) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==24) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==25) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==26) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==27) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*2.0*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==28) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==29) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.1*(x+y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*(x+y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==30) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.2*(x+y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==31) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*abs(sin(4.0*pi*x))
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==32) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*sin(pi*x)*sin(pi*x)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==33) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*sin(pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*sin(pi*x)*sin(pi*y)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==34) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.1*sin(pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*0.1*sin(pi*x)*sin(pi*y)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.1*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==35) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*sin(pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*0.5*sin(pi*x)*sin(pi*y)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.5*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==36) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*sin(0.5*pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*sin(0.5*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*sin(0.5*pi*x)*sin(pi*y)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*sin(0.5*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==37) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1)  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.1*sin(0.1*pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*sin(0.1*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*0.1*sin(0.1*pi*x)*sin(pi*y)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.1*sin(0.1*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==38) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*sin(0.5*pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*sin(0.5*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*0.5*sin(0.5*pi*x)*sin(pi*y)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.5*sin(0.5*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==39) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*x*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*x*y
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==40) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.1*x*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*0.1*x*y
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.1*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==41) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*x*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*0.5*x*y
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.5*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==42) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*x*x*y*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*x*x*y*y
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==43) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*2.0*x*x*y*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*2.0*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*2.0*x*x*y*y
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*2.0*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==44) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*x*x*y*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*0.5*x*x*y*y
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.5*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==45) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.1*(x+y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*(x+y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*0.1*(x+y)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.1*(x+y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==46) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.2*(x+y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.2*(x+y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*0.2*(x+y)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.2*(x+y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==47) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*abs(sin(4.0*pi*x))
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*abs(sin(4.0*pi*x))
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*0.5*abs(sin(4.0*pi*x))
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.5*abs(sin(4.0*pi*x))
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==48) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*sin(pi*x)*sin(pi*x)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*sin(pi*x)*sin(pi*x)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.5*sin(pi*x)*sin(pi*x)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==49) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*sin(pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==50) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.1*sin(pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.1*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==51) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*sin(pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.5*sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==52) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*sin(0.5*pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*sin(0.5*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*sin(0.5*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==53) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1)  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.1*sin(0.1*pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*sin(0.1*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.1*sin(0.1*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==54) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*sin(0.5*pi*x)*sin(pi*y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*sin(0.5*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.5*sin(0.5*pi*x)*sin(pi*y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==55) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*x*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==56) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.1*x*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.1*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==57) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*x*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.5*x*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==58) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*x*x*y*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==59) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*2.0*x*x*y*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*2.0*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*2.0*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==60) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*x*x*y*y
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.5*x*x*y*y
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==61) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.1*(x+y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.1*(x+y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.1*(x+y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==62) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.2*(x+y)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.2*(x+y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.2*(x+y)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==63) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*abs(sin(4.0*pi*x))
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*abs(sin(4.0*pi*x))
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.5*abs(sin(4.0*pi*x))
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if

if(mod(l,64)==0) then
!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)

b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*0.5*sin(pi*x)*sin(pi*x)
b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*0.5*sin(pi*x)*sin(pi*x)
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL &
!$OMP PRIVATE(x, y) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(k)
b(l,i,1 ,k)  = - (1.0-real(ibits(bp(i,1,k ), bc_d_S, 1)))*0.5*sin(pi*x)*sin(pi*x)
b(l,i,jx,k) = - (1.0-real(ibits(bp(i,jx,k), bc_d_N, 1)))*0.5*sin(pi*x)*sin(pi*x)
end do
end do
!$OMP END DO
!$OMP END PARALLEL
end if


if(l<0) then
!$OMP PARALLEL &                                                                                                                                                
!$OMP PRIVATE(x, y, z) &                                                                                                                                        
!$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)                                                                                                                      
!$OMP DO SCHEDULE(static) COLLAPSE(2)                                                                                                                           
do k=1,kx
do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)
z = dh*real(k)
x2=x*x
y2=y*y
z2=z*z
b(l,i,j,k) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z)
end do
end do
end do
!$OMP END DO                                                                                                                                                    
!$OMP END PARALLEL                                                                                                                                              

do j=1,jx
do i=1,ix
x = dh*real(i)
y = dh*real(j)
z = dh*real(kx)
x2=x*x
y2=y*y
z2=z*z
b(l,i,j,kx) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*sin(x*y)
end do
end do

do k=1,kx
do i=1,ix
x = dh*real(i)
y = dh*real(jx)
z = dh*real(k)
x2=x*x
y2=y*y
z2=z*z
b(l,i,jx,k) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) - (1.0-real(ibits(bp(i ,jx ,k), bc_d_N, 1)))*sin(x*z)
end do
end do

do k=1,kx
do j=1,jx
x = dh*real(ix)
y = dh*real(j)
z = dh*real(k)
x2=x*x
y2=y*y
z2=z*z
b(l,ix,j,k) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) - (1.0-real(ibits(bp(ix,j,k), bc_d_E, 1)))*sin(y*z)
end do
end do

do j=1,jx
x = dh*real(1)
y = dh*real(j)
z = dh*real(kx)
x2=x*x
y2=y*y
z2=z*z
b(l,1,j,kx ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
              -(1.0-real(ibits(bp(1 ,j,kx), bc_d_T, 1)))*sin(x*y)
x = dh*real(ix)
y = dh*real(j)
z = dh*real(1)
x2=x*x
y2=y*y
z2=z*z
b(l,ix,j,1 ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
     -(1.0-real(ibits(bp(ix ,j,1), bc_d_E, 1)))*sin(y*z)
x = dh*real(ix)
y = dh*real(j)
z = dh*real(kx)
x2=x*x
y2=y*y
z2=z*z
b(l,ix,j,kx ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
              -(1.0-real(ibits(bp(ix ,j,kx), bc_d_T, 1)))*sin(x*y)&
              -(1.0-real(ibits(bp(ix ,j,kx), bc_d_E, 1)))*sin(y*z)
end do

do i=1,ix
x = dh*real(i)
y = dh*real(1)
z = dh*real(kx)
x2=x*x
y2=y*y
z2=z*z
b(l,i,1,kx ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
              -(1.0-real(ibits(bp(i ,1,kx), bc_d_T, 1)))*sin(x*y)
x = dh*real(i)
y = dh*real(jx)
z = dh*real(1)
x2=x*x
y2=y*y
z2=z*z
b(l,i,jx,1 ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
              -(1.0-real(ibits(bp(i ,jx,1), bc_d_N, 1)))*sin(x*z)
x = dh*real(i)
y = dh*real(jx)
z = dh*real(kx)
x2=x*x
y2=y*y
z2=z*z
b(l,i,jx,kx ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
     -(1.0-real(ibits(bp(i ,jx,kx), bc_d_T, 1)))*sin(x*y)&
     -(1.0-real(ibits(bp(i ,jx,kx), bc_d_N, 1)))*sin(x*z)
end do

do k=1,kx
x = dh*real(1)
y = dh*real(jx)
z = dh*real(k)
x2=x*x
y2=y*y
z2=z*z
b(l,1,jx,k ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
              -(1.0-real(ibits(bp(1 ,jx,k), bc_d_N, 1)))*sin(x*z)
x = dh*real(ix)
y = dh*real(1)
z = dh*real(k)
x2=x*x
y2=y*y
z2=z*z
b(l,ix,1,k ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
              -(1.0-real(ibits(bp(ix ,1,k), bc_d_E, 1)))*sin(y*z)
x = dh*real(ix)
y = dh*real(jx)
z = dh*real(k)
x2=x*x
y2=y*y
z2=z*z
b(l,ix,jx,k ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
     -(1.0-real(ibits(bp(ix ,jx,k), bc_d_N, 1)))*sin(x*z)&
     -(1.0-real(ibits(bp(ix ,jx,k), bc_d_E, 1)))*sin(y*z)
end do

x = dh*real(1)
y = dh*real(1)
z = dh*real(kx)
x2=x*x
y2=y*y
z2=z*z
b(l,1,1,kx ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
     -(1.0-real(ibits(bp(1 ,1,kx), bc_d_T, 1)))*sin(x*y)

x = dh*real(1)
y = dh*real(jx)
z = dh*real(1)
x2=x*x
y2=y*y
z2=z*z
b(l,1,jx,1 ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
     -(1.0-real(ibits(bp(1 ,jx,1), bc_d_N, 1)))*sin(x*z)

x = dh*real(1)
y = dh*real(jx)
z = dh*real(kx)
x2=x*x
y2=y*y
z2=z*z
b(l,1,jx,kx ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
     -(1.0-real(ibits(bp(1 ,jx,kx), bc_d_N, 1)))*sin(x*z)&
     -(1.0-real(ibits(bp(1 ,jx,kx), bc_d_T, 1)))*sin(x*y)

x = dh*real(ix)
y = dh*real(1)
z = dh*real(1)
x2=x*x
y2=y*y
z2=z*z
b(l,ix,1,1 ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
     -(1.0-real(ibits(bp(ix ,1,1), bc_d_E, 1)))*sin(y*z)

x = dh*real(ix)
y = dh*real(1)
z = dh*real(kx)
x2=x*x
y2=y*y
z2=z*z
b(l,ix,1,kx ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
     -(1.0-real(ibits(bp(ix ,1,kx), bc_d_E, 1)))*sin(y*z)&
     -(1.0-real(ibits(bp(ix ,1,kx), bc_d_T, 1)))*sin(x*y)

x = dh*real(ix)
y = dh*real(jx)
z = dh*real(1)
x2=x*x
y2=y*y
z2=z*z
b(l,ix,jx,1 ) = -dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
     -(1.0-real(ibits(bp(ix ,jx,1), bc_d_E, 1)))*sin(y*z)&
     -(1.0-real(ibits(bp(ix ,jx,1), bc_d_N, 1)))*sin(x*z)

x = dh*real(ix)
y = dh*real(jx)
z = dh*real(kx)
x2=x*x
y2=y*y
z2=z*z
b(l,ix,jx,kx ) = (-dh2*(x2*y2+y2*z2+z2*x2)*sin(x*y*z) &
     -(1.0-real(ibits(bp(ix ,jx,kx), bc_d_E, 1)))*sin(y*z)&
     -(1.0-real(ibits(bp(ix ,jx,kx), bc_d_N, 1)))*sin(x*z)&
     -(1.0-real(ibits(bp(ix ,jx,kx), bc_d_T, 1)))*sin(x*y))

! !$OMP PARALLEL &
! ! !$OMP REDUCTION(+:res) &
! ! !$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t) &
! !$OMP PRIVATE(x, y) &
! !$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)

! !$OMP DO SCHEDULE(static) COLLAPSE(2)

! do j=1,jx
! do i=1,ix
! x = dh*real(i)
! y = dh*real(j)
! end do
! end do
! !$OMP END DO
! !$OMP END PARALLEL
end if

! if(lx>2) then
! l=3
! !$OMP PARALLEL &
! ! !$OMP REDUCTION(+:res) &
! ! !$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t) &
! !$OMP PRIVATE(x, y) &
! !$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)

! !$OMP DO SCHEDULE(static) COLLAPSE(2)
! do j=1,jx
! do i=1,ix
! x = dh*real(i)
! y = dh*real(j)

! b(l,i,j,1 )  = - (1.0-real(ibits(bp(i,j,1 ), bc_d_B, 1)))*sin(pi*x)*sin(pi*y)
! b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*sin(pi*x)*sin(pi*y)
! end do
! end do
! !$OMP END DO
! !$OMP END PARALLEL
! ! !$OMP PARALLEL &
! ! !$OMP REDUCTION(+:res) &
! ! !$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t) &
! ! !$OMP PRIVATE(x, y) &
! ! !$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh, pi)

! ! !$OMP DO SCHEDULE(static) COLLAPSE(2)

! ! do j=1,jx
! ! do i=1,ix
! ! x = dh*real(i)
! ! y = dh*real(j)
! ! b(l,i,j,kx) = - (1.0-real(ibits(bp(i,j,kx), bc_d_T, 1)))*(x+y)
! ! end do
! ! end do
! ! !$OMP END DO
! ! !$OMP END PARALLEL
! !
! end if
end do

return
end subroutine src_dirichlet
