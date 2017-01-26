!
!	utility.f90
!	Mcube
!
!	Created by keno on 11/12/05.
!	Copyright 2011 __iis__. All rights reserved.
!

!   **************************************
subroutine fileout (sz, g, psz4dex, s, dh, org, fname)
  implicit none
  integer                                                       :: nn, ix, jx, kx, i, j, k, g, l, lx
  integer, dimension(4)                                         :: sz, psz4dex
  real, dimension(sz(4)+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(4)):: s
  real                                                          :: rtime, dh
  real, dimension(3)                                            :: org
  character*20                                                  :: fname, tname

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  lx = sz(4)

  nn = 0
  rtime = 0.0

  do l=1,lx
     write(tname,'(i3.3,"_",a10)') l,fname
     open (unit=22, file=tname, form='unformatted')
     write (22) 1, 1
     write (22) ix, jx, kx
     write (22) org(1), org(2), org(3)
     write (22) dh, dh, dh
     write (22) nn, rtime
     write (22) (((s(l,i,j,k),i=1,ix),j=1,jx),k=1,kx)
     close (unit=22)
  end do

  return
end subroutine fileout


!   ****************************
subroutine bc (sz, g, psz4dex, p, dh, vec_tag, vec_cnt)
  implicit none
  integer                                                        :: i, j, k, ix, jx, kx, g, l, lx, vec_cnt
  integer, dimension(4)                                          :: sz, psz4dex
  real, dimension(sz(4)+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3), 1-g:sz(3)+g+psz4dex(3)) :: p
  real                                                           :: pi, x, y, z, dh
  integer, dimension(sz(4))                                      :: vec_tag

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  lx = vec_cnt

  pi = 2.0*asin(1.0)

  do l=1,lx
     if(mod(vec_tag(l),32)==1) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = sin(pi*x)*sin(pi*y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = sin(pi*x)*sin(pi*y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==2) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = 0.1*sin(pi*x)*sin(pi*y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = 0.1*sin(pi*x)*sin(pi*y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==3) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = 0.5*sin(pi*x)*sin(pi*y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = 0.5*sin(pi*x)*sin(pi*y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==4) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = sin(0.5*pi*x)*sin(pi*y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = sin(0.5*pi*x)*sin(pi*y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==5) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = 0.1*sin(0.5*pi*x)*sin(pi*y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = 0.1*sin(0.5*pi*x)*sin(pi*y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==6) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = 0.5*sin(0.5*pi*x)*sin(pi*y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = 0.5*sin(0.5*pi*x)*sin(pi*y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==7) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = x*y
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = x*y
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==8) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = 0.1*x*y
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = 0.1*x*y
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==9) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = 0.5*x*y
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = 0.5*x*y
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==10) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = x*x*y*y
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = x*x*y*y
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==11) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = 2.0*x*x*y*y
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = 2.0*x*x*y*y
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==12) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = 0.5*x*x*y*y
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = 0.5*x*x*y*y
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==13) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = 0.1*(x+y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = 0.1*(x+y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==14) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = 0.2*(x+y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = 0.2*(x+y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==15) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = 0.5*abs(sin(pi*x*4.0))
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = 0.5*abs(sin(pi*x*4.0))
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(mod(vec_tag(l),32)==0) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = 0.5*sin(pi*x)*sin(pi*x)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = 0.5*sin(pi*x)*sin(pi*x)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

       ! XPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,ix+1,j,k) = 0.0 !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0, k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,jx+1,k) = 0.0 !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(vec_tag(l)<0) then
        ! ZMINUS Dirichlet

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, pi, dh) &
        !$OMP PRIVATE(x, y)

        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,0) = 0.0
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! ZPLUS Dirichlet
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx, pi, dh) &
        !$OMP PRIVATE(x, y)

        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do j=1,jx
           do i=1,ix
              x = dh*real(i)
              y = dh*real(j)
              p(l,i,j,kx+1) = sin(x*y)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(jx, kx)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              p(l,0,   j,k) = 0.0 !-p(1, j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! XPLUS

        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx) &
        !$OMP PRIVATE(z, y)        
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do j=1,jx
              z = dh*real(k)
              y = dh*real(j)
              p(l,ix+1,j,k) = sin(y*z) !-p(ix,j,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! YMINUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx)

        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              p(l,i,0,   k) = 0.0 !-p(i,1, k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL


        ! YPLUS
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx) &
        !$OMP PRIVATE(x, z)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1,kx
           do i=1,ix
              x = dh*real(i)
              z = dh*real(k)
              p(l,i,jx+1,k) = sin(x*z) !-p(i,jx,k)
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if
  end do

  return
end subroutine bc

!   *******************************
subroutine exact (sz, g, e, dh)
  implicit none
  integer                                                      ::  i, j, k, ix, jx, kx, g, l, lx, n, m
  integer, dimension(4)                                        ::  sz
  real, dimension(sz(4),1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  e
  real                                                         ::  dh, pi, r2, x, y, z
  double precision                                             :: a,sum,w

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  lx = sz(4)
  r2 = sqrt(2.0)
  pi = 2.0*asin(1.0)

  do l=1,lx
     if(l>0) then
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx) &
        !$OMP PRIVATE(x, y, z)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1, kx
           do j=1, jx
             do i=1, ix
                 x = dh*real(i)
                 y = dh*real(j)
                 z = dh*real(k)
                 e(l,i,j,k) = sin(pi*x)*sin(pi*y) / sinh(r2*pi) * ( sinh(r2*pi*z)-sinh(r2*pi*(z-1.0))  )
              end do
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if

     if(l<0) then
        !$OMP PARALLEL &
        !$OMP FIRSTPRIVATE(ix, jx, kx) &
        !$OMP PRIVATE(x, y, z)
        !$OMP DO SCHEDULE(static) COLLAPSE(2)
        do k=1, kx
           do j=1, jx
              do i=1, ix
                 x = dh*real(i)
                 y = dh*real(j)
                 z = dh*real(k)
                 e(l,i,j,k) = sin(x*y*z)
              end do
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if
  end do
  return
end subroutine exact


!*******************************
subroutine err (sz, g, dh, d, p, e, gosa)
  implicit none
  include "cbc_f_params.h"
  integer                                                      ::  i, j, k, ix, jx, kx, g, l, lx
  integer, dimension(4)                                        ::  sz
  real, dimension(sz(4)+1,1-g:sz(1)+g+1, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p, e, gosa
  real                                                         ::  r2, pi, x, y, z, r, dh
  double precision,dimension(sz(4))                            ::  d

  d=0.0d0
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  lx = sz(4)

  r2 = sqrt(2.0)
  pi = 2.0*asin(1.0)

  !$OMP PARALLEL &              
  !$OMP REDUCTION(max:d) &
  !$OMP PRIVATE(r) &   
  !$OMP FIRSTPRIVATE(ix, jx, kx, lx, dh,pi)       
  !$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1, kx
     do j=1, jx
        do i=1, ix
           do l=1, lx
              r = p(l,i,j,k) -  e(l,i,j,k)
              gosa(l,i,j,k) = r
              d(l) = max(d(l), abs(dble(r)))
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  return
end subroutine err

!   ****************************                                                                       
subroutine set_bcindex (sz, g, bp)
  implicit none
  integer                                                :: i, j, k, ix, jx, kx, g, l, lx
  integer, dimension(3)                                  :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g,8) :: bp

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  lx = 8

  !$OMP PARALLEL &
  !$OMP FIRSTPRIVATE(ix, jx, kx)
  !$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1-g,kx+g
     do j=1-g,jx+g
        do i=1-g,ix+g
           bp(i,j,k,1) = 0!e
           bp(i,j,k,2) = 0!w
           bp(i,j,k,3) = 0!n
           bp(i,j,k,4) = 0!s
           bp(i,j,k,5) = 0!t
           bp(i,j,k,6) = 0!b
           bp(i,j,k,7) = 1!dig
           bp(i,j,k,8) = 0!active
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL


  !$OMP PARALLEL &
  !$OMP FIRSTPRIVATE(ix, jx, kx)
  !$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1,kx
     do j=1,jx
        do i=1,ix
           bp(i,j,k,1) = 1!e
           bp(i,j,k,2) = 1!w
           bp(i,j,k,3) = 1!n
           bp(i,j,k,4) = 1!s
           bp(i,j,k,5) = 1!t
           bp(i,j,k,6) = 1!b
           bp(i,j,k,7) = 6!dig
           bp(i,j,k,8) = 1!active
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  do k=1,kx
     do j=1,jx
        bp(1,j,k,2) = 0!w
        bp(ix,j,k,1) = 0!e
     end do
  end do

  do k=1,kx
     do i=1,ix
        bp(i,1,k,4) = 0!s
        bp(i,jx,k,3) = 0!n
     end do
  end do

  do j=1,jx
     do i=1,ix
        bp(i,j,1,6) = 0!t
        bp(i,j,kx,5) = 0!b
     end do
  end do

  do k=1,kx
     bp(1,1,k,2) = 0!
     bp(1,1,k,4) = 0!

     bp(1,jx,k,2) = 0!
     bp(1,jx,k,3) = 0!

     bp(ix,1,k,1) = 0!
     bp(ix,1,k,4) = 0!

     bp(ix,jx,k,1) = 0!
     bp(ix,jx,k,3) = 0!
  end do

  do j=1,jx
     bp(1,j,1,2) = 0!
     bp(1,j,1,6) = 0!

     bp(1,j,kx,2) = 0!
     bp(1,j,kx,5) = 0!
     bp(ix,j,1,1) = 0!
     bp(ix,j,1,6) = 0!

     bp(ix,j,kx,1) = 0!
     bp(ix,j,kx,5) = 0!
  end do

  do i=1,ix
     bp(i,1,1,4) = 0!
     bp(i,1,1,6) = 0!

     bp(i,1,kx,4) = 0!
     bp(i,1,kx,5) = 0!

     bp(i,jx,1,3) = 0!
     bp(i,jx,1,6) = 0!

     bp(i,jx,kx,3) = 0!
     bp(i,jx,kx,5) = 0!
  end do

  bp(1,1,1,2) = 0!
  bp(1,1,1,4) = 0!
  bp(1,1,1,6) = 0!

  bp(1,1,kx,2) = 0!
  bp(1,1,kx,4) = 0!
  bp(1,1,kx,5) = 0!

  bp(1,jx,1,2) = 0!
  bp(1,jx,1,3) = 0!
  bp(1,jx,1,6) = 0!

  bp(1,jx,kx,2) = 0!
  bp(1,jx,kx,3) = 0!
  bp(1,jx,kx,5) = 0!

  bp(ix,1,1,1) = 0!
  bp(ix,1,1,4) = 0!
  bp(ix,1,1,6) = 0!

  bp(ix,1,kx,1) = 0!
  bp(ix,1,kx,4) = 0!
  bp(ix,1,kx,5) = 0!

  bp(ix,jx,1,1) = 0!
  bp(ix,jx,1,3) = 0!
  bp(ix,jx,1,6) = 0!

  bp(1,jx,kx,2) = 0!
  bp(1,jx,kx,3) = 0!
  bp(1,jx,kx,5) = 0!

  bp(ix,1,1,1) = 0!
  bp(ix,1,1,4) = 0!
  bp(ix,1,1,6) = 0!

  bp(ix,1,kx,1) = 0!
  bp(ix,1,kx,4) = 0!
  bp(ix,1,kx,5) = 0!

  bp(ix,jx,1,1) = 0!
  bp(ix,jx,1,3) = 0!
  bp(ix,jx,1,6) = 0!

  bp(ix,jx,kx,1) = 0!
  bp(ix,jx,kx,3) = 0!
  bp(ix,jx,kx,5) = 0!

  return
end subroutine set_bcindex

!   **************************************
subroutine get_vec_sz (sz)
  implicit none
  include "cbc_f_params.h"
  integer, dimension(4)                                       :: sz

  sz(1) = Grid_size
  sz(2) = Grid_size
  sz(3) = Grid_size
  sz(4) = Vector_num

  return
end subroutine get_vec_sz
!****************************************

subroutine change (sz, g, psz4dex, p, ch, vec_cnt)
  implicit none
  integer                    ::  i, j, k, ix, jx, kx,g, l, lx,ch, vec_cnt
  integer,dimension(4)       ::  sz, psz4dex
  real, dimension(sz(4)+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3),1-g:sz(3)+g+psz4dex(4)) :: p
  real                       ::  tmp

  ix = sz(1)
  jx = sz(2) 
  kx = sz(3)
  lx = sz(4)
  !$OMP PARALLEL &
  !$OMP PRIVATE(tmp) &
  !$OMP FIRSTPRIVATE(ix, jx, kx, ch, vec_cnt)
  !$OMP DO SCHEDULE(static) COLLAPSE(2)
 do k=1-g, kx+g+psz4dex(4)
     do j=1-g, jx+g+psz4dex(3)
        do i=1-g, ix+g+psz4dex(2)
           tmp =  p(ch,i,j,k)
           p(ch,i,j,k) = p(vec_cnt,i,j,k)
           p(vec_cnt,i,j,k) = tmp
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  return
end subroutine change

!-------------------
subroutine array_reshape (sz, g, psz4dex, p0, p1, p, vec_cnt, table_con, num_con, table_not_con, num_not_con)
  implicit none
  integer                    ::  i, j, k, ix, jx, kx, g, l, lx, vec_cnt, offset
  integer,dimension(4)       ::  sz, psz4dex
  integer,dimension(256)     ::  table_con, table_not_con
  integer                    ::  num_con, num_not_con
  real, dimension(vec_cnt    +psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3),1-g:sz(3)+g+psz4dex(4)) ::  p0
  real, dimension(num_not_con+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3),1-g:sz(3)+g+psz4dex(4)) ::  p1
  real, dimension(sz(4)      +psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3),1-g:sz(3)+g+psz4dex(4)) ::  p

  ix = sz(1)
  jx = sz(2) 
  kx = sz(3)
  lx = sz(4)
  offset = sz(4) - vec_cnt 

  !print *,"re",vec_cnt,table_con(1:sz(4)),num_con,table_not_con(1:sz(4)),num_not_con

  !$OMP PARALLEL &
  !$OMP FIRSTPRIVATE(ix, jx, kx, vec_cnt, g, num_con, num_not_con, offset, table_con, table_not_con)
  !$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1-g, kx+g+psz4dex(4)
     do j=1-g, jx+g+psz4dex(3)
        do i=1-g, ix+g+psz4dex(2)
           do l=1,num_not_con
              p1(l,i,j,k)= p0(table_not_con(l),i,j,k)
           end do           
           
           do l=1,num_con
              p(l+offset,i,j,k)=p0(table_con(l),i,j,k)
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  return
end subroutine array_reshape

!------------------
subroutine array_small (sz, g, psz4dex, p0, p1, vec_cnt, table_not_con, num_not_con)
  implicit none
  integer                    ::  i, j, k, ix, jx, kx, g, l, lx, vec_cnt
  integer,dimension(4)       ::  sz, psz4dex
  integer,dimension(256)     ::  table_not_con
  integer                    ::  num_not_con
  real, dimension(vec_cnt    +psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3),1-g:sz(3)+g+psz4dex(4)) ::  p0
  real, dimension(num_not_con+psz4dex(1),1-g:sz(1)+g+psz4dex(2), 1-g:sz(2)+g+psz4dex(3),1-g:sz(3)+g+psz4dex(4)) ::  p1

  ix = sz(1)
  jx = sz(2) 
  kx = sz(3)
  lx = sz(4)
  !print *,"small",vec_cnt,num_not_con
  !print *,"     ",num_not_con,":",table_not_con(1:8)
  
  !$OMP PARALLEL &
  !$OMP FIRSTPRIVATE(ix, jx, kx, g, vec_cnt, table_not_con, num_not_con)
  !$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1-g, kx+g+psz4dex(4)
     do j=1-g, jx+g+psz4dex(3)
        do i=1-g, ix+g+psz4dex(2)
           do l=1,num_not_con
              p1(l,i,j,k)=p0(table_not_con(l),i,j,k)
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  return
end subroutine array_small

!---------------
subroutine init_array (sz, g, p)
  implicit none
  integer                    ::  i, j, k, ix, jx, kx, g, l, lx
  integer,dimension(4)       ::  sz
  real, dimension(sz(4)    +1,1-g:sz(1)+g+1, 1-g:sz(2)+g,1-g:sz(3)+g) ::  p

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  lx = sz(4)
  !print *,"small",vec_cnt,num_not_con
  !print *,"     ",num_not_con,":",table_not_con(1:8)

  !$OMP PARALLEL &
  !$OMP FIRSTPRIVATE(ix, jx, kx, g)
  !$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1-g, kx+g
     do j=1-g, jx+g
        do i=1-g, ix+g+1
           do l=64,sz(4)
              p(l,i,j,k)=1.0
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  return
end subroutine init_array
