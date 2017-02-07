!###################################################################################
!
! CPMlib - Computational space Partitioning Management library
!
! Copyright (c) 2012-2014 Institute of Industrial Science (IIS), The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2014-2016 Advanced Institute for Computational Science (AICS), RIKEN.
! All rights reserved.
!
! Copyright (c) 2016-2017 Research Institute for Information Technology (RIIT), Kyushu University.
! All rights reserved.
!
!###################################################################################

!   *********************************************************
!   ***      Sample code of porting
!   ***      2010 keno@VCAD, RIKEN
!
!   *********************************************************
    module array_def
    real, dimension(:,:,:,:), allocatable, save :: p, q, w
    end module array_def

!   *********************************************************
    subroutine diffusion(iScheme, iStpMax, dDt, dCf, iItrMax, dEps, dOmg, ierr)

    use array_def
    implicit none

    include 'mpif.h'
    include 'cpm_fparam.fi'

    integer :: iScheme, iStpMax, iItrMax, ierr
    real*8  :: dDt, dCf, dEps, dOmg

    integer, parameter :: domainType = 1
    integer            :: nx, ny, nz
    integer            :: n, scheme, itrMax, stp, stpMax
    real               :: dx, dy, dz, ox, oy, oz
    real               :: eps, er, omg, dt, cf
    real               :: pmin, pmax

    integer, parameter :: pg=0
    integer            :: nrank, myrank
    integer            :: sz(3), div(3), pos(3), head(0:2)
    integer            :: nID(4,0:5), nIDNum(0:5)
    integer            :: neiLeafID(4,0:5), neiLeafNum(0:5)
    real*8             :: org8(3), pch8(3), rgn8(3)
    integer            :: i,j
    integer            :: nLeaf, leafID
    character*2        :: chDir(0:5)
    data chDir /"-X","+X","-Y","+Y","-Z","+Z"/

    integer, parameter :: gc = 2

    ! copy param
    scheme = iScheme
    stpMax = iStpMax
    dt     = dDt
    cf     = dCf
    itrMax = iItrMax
    eps    = dEps
    omg    = dOmg

    ! get number of rank
    call cpm_GetNumRank_LMR(nrank, pg, ierr)

    ! get my rank No.
    call cpm_GetMyRankID_LMR(myrank, pg, ierr)

    call flush(5)
    call cpm_Barrier_LMR(pg, ierr)
    if( myrank.eq.0 ) then
      write (*,'(/,A)') "***** MCONV process *****"
      write (*,*) 'number of process=', nrank

      if ( scheme == 1 ) then
        write(*,*) 'Euler Explicit ftcs'
      else
        write(*,*) 'Euler Implicit jacobi'
      endif
    end if

    ! get num leaf
    call cpm_GetLocalNumLeaf_LMR( nLeaf, pg, ierr )

    ! get local voxel size
    call cpm_GetLocalVoxelSize_LMR( sz, pg, ierr )

    ! get local leaf head index, neighbor rank No.
    call flush(5)
    call cpm_Barrier_LMR(pg, ierr)
    do i=1,nrank
      if( i-1.eq.myrank ) then
        do n=1,nLeaf
          call cpm_GetLeafID_LMR( n, leafID, pg, ierr )
          call cpm_GetVoxelHeadIndex_LMR( n, head, pg, ierr )
          do j=0,5
            call cpm_GetNeighborRankList_LMR( n, j, nID(1,j), nIDNum(j), pg, ierr )
            call cpm_GetNeighborLeafList_LMR( n, j, neiLeafID(1,j), neiLeafNum(j), pg, ierr )
          enddo
          call cpm_GetLocalOrigin_LMR( n, org8, pg, ierr )
          call cpm_GetPitch_LMR( n, pch8, pg, ierr )
          call cpm_GetLocalRegion_LMR( n, rgn8, pg, ierr )
          call cpm_GetDivNum_LMR( n, div, pg, ierr )
          call cpm_GetDivPos_LMR( n, pos, pg, ierr )
          nx = sz(1)
          ny = sz(2)
          nz = sz(3)
          dx = pch8(1)
          dy = pch8(2)
          dz = pch8(3)
          ox = org8(1)
          oy = org8(2)
          oz = org8(3)
          write(*,'(A,2(/,A,I4),4(/,A,3(1X,I4)),3(/,A,3(1X,1PE11.4)))') &
            "+--------------------------------------+" &
          , "myrank  = ", myrank &
          , "leaf ID = ", leafID &
          , "  local voxel   =", sz &
          , "  div num       =", div &
          , "  div pos       =", pos &
          , "  head index    =", head &
          , "  local origin  =", ox, oy, oz &
          , "  local pitch   =", dx, dy, dz &
          , "  local region  =", rgn8
          do j=0,5
            if( nIDNum(j).gt.0 ) then
              write(*,'(3A,4(1X,I4))') &
              , "  ", chDir(j), " neighbor rank ID =", nID(1:nIDNum(j),j)
            else
              write(*,'(3A)') &
              , "  ", chDir(j), " neighbor rank ID = none"
            endif
          enddo
          do j=0,5
            if( neiLeafNum(j).gt.0 ) then
              write(*,'(3A,4(1X,I4))') &
              , "  ", chDir(j), " neighbor leaf ID =", neiLeafID(1:neiLeafNum(j),j)
            else
              write(*,'(3A)') &
              , "  ", chDir(j), " neighbor leaf ID = none"
            endif
          enddo
        enddo
      endif
      call flush(5)
      call cpm_Barrier_LMR(pg, ierr)
    enddo

    if( myrank.eq.0 ) then
      write (*,'(/,A)') "***** loop start *****"
    endif

    call allocate_array (nx, ny, nz, nLeaf, gc)

    call pbc (nx, ny, nz, nLeaf, gc, p)

!   time step loop
    do stp=1, stpMax
      if ( scheme == 1 ) then
        call ftcs (nx, ny, nz, nLeaf, gc, p, q, cf, dt, er)

        call pbc (nx, ny, nz, nLeaf, gc, p)

        call minmax(nx, ny, nz, nLeaf, gc, p, pmin, pmax)
        if( myrank == 0 ) write (*,*) stp, er, pmin, pmax

      else
        er = 0.0
        w(0:nx+1, 0:ny+1, 0:nz+1, 1:nLeaf) = p(0:nx+1, 0:ny+1, 0:nz+1, 1:nLeaf)
        ! iteration
        do n=1, itrMax
          call jacobi (nx, ny, nz, nLeaf, gc, p, q, w, cf, dt, omg, er)

          call pbc (nx, ny, nz, nLeaf, gc, p)

          if (er<eps) exit

        end do

        call minmax(nx, ny, nz, nLeaf, gc, p, pmin, pmax)
        if( myrank == 0 ) write (*,*) stp, n, er, pmin, pmax

      end if
    end do

!   Post
    call fileout (nx, ny, nz, nLeaf, gc, p, w)


    stop
    end subroutine diffusion

!   ******************************
    subroutine allocate_array (nx, ny, nz, nLeaf, gc)
    use array_def
    implicit none
    integer :: status
    integer :: nx, ny, nz, nLeaf, gc

    allocate ( p(1-gc:nx+gc, 1-gc:ny+gc, 1-gc:nz+gc, nLeaf) , stat=status )
    allocate ( w(1-gc:nx+gc, 1-gc:ny+gc, 1-gc:nz+gc, nLeaf) , stat=status )
    allocate ( q(1-gc:nx+gc, 1-gc:ny+gc, 1-gc:nz+gc, nLeaf) , stat=status )
    p(1-gc:nx+gc, 1-gc:ny+gc, 1-gc:nz+gc, 1:nLeaf) = 0.d0

    return
    end subroutine allocate_array

!   ****************************************************
    subroutine minmax(nx, ny, nz, nLeaf, gc, p, pmin, pmax)
    implicit none

    include 'cpm_fparam.fi'

    integer                                                 :: nx, ny, nz, nLeaf, gc
    real, dimension(1-gc:nx+gc,1-gc:ny+gc,1-gc:nz+gc,nLeaf) :: p
    real                                                    :: pmin, pmax

    integer                                                 :: i, j, k, n
    real                                                    :: bmin, bmax
    integer, parameter                                      :: pg=0
    integer                                                 :: ierr

    bmin = 1.d+30
    bmax = -1.d+30
    do n=1,nLeaf
    do k=1,nz
    do j=1,ny
    do i=1,nx
      bmin = min(bmin, p(i,j,k,n))
      bmax = max(bmax, p(i,j,k,n))
    enddo
    enddo
    enddo
    enddo

    pmin = bmin
    pmax = bmax
    call cpm_Allreduce_LMR(bmin,pmin,1,CPM_REAL,CPM_MIN,pg,ierr)
    call cpm_Allreduce_LMR(bmax,pmax,1,CPM_REAL,CPM_MAX,pg,ierr)

    end subroutine minmax

!   ****************************************************
    subroutine jacobi (nx, ny, nz, nLeaf, gc, p, q, w, cf, dt, omg, er)
    implicit none

    include 'cpm_fparam.fi'

    integer                                                 :: nx, ny, nz, nLeaf, gc
    real, dimension(1-gc:nx+gc,1-gc:ny+gc,1-gc:nz+gc,nLeaf) :: p, q, w
    real                                                    :: cf, dt, omg, er

    integer                                                 :: i, j, k, n
    real*8                                                  :: pch8(3)
    real                                                    :: dx, dy, dz
    real                                                    :: ddx, ddy, ddz, cpd
    real                                                    :: s0, ss, sx, sy, sz
    integer, parameter                                      :: pg=0
    integer                                                 :: ierr
    real                                                    :: er_buf

    er = 0.0

    do n=1,nLeaf
      call cpm_GetPitch_LMR( n, pch8, pg, ierr )
      dx = pch8(1)
      dy = pch8(2)
      dz = pch8(3)
      ddx = cf*dt/(dx*dx)
      ddy = cf*dt/(dy*dy)
      ddz = cf*dt/(dz*dz)
      cpd=1.0/(1.0+2.0*ddx+2.0*ddy+2.0*ddz)

      do k=1,nz
      do j=1,ny
      do i=1,nx
        sx= p(i+1,j  ,k  ,n) + p(i-1,j  ,k  ,n)
        sy= p(i  ,j+1,k  ,n) + p(i  ,j-1,k  ,n)
        sz= p(i  ,j  ,k+1,n) + p(i  ,j  ,k-1,n)
        s0=ddx*sx+ddy*sy+ddz*sz
        ss=(cpd*(s0+w(i,j,k,n))-p(i,j,k,n))*omg
        q(i,j,k,n)=p(i,j,k,n)+ss
        er = er + ss*ss
      end do
      end do
      end do
    enddo

    p(1:nx, 1:ny, 1:nz, 1:nLeaf) = q(1:nx, 1:ny, 1:nz, 1:nLeaf)

    call cpm_BndCommS3D_LMR(p,nx,ny,nz,gc,1,CPM_REAL,pg,ierr)
    er_buf = er
    call cpm_Allreduce_LMR(er_buf,er,1,CPM_REAL,CPM_SUM,pg,ierr)
    er = sqrt(er)

    return
    end subroutine jacobi

!   ***************************************
    subroutine ftcs (nx, ny, nz, nLeaf, gc, p, q, cf, dt, er)
    implicit none
    include 'cpm_fparam.fi'

    integer                                                 :: nx, ny, nz, nLeaf, gc
    real, dimension(1-gc:nx+gc,1-gc:ny+gc,1-gc:nz+gc,nLeaf) :: p, q
    real                                                    :: cf, dt, er

    integer                                                 :: i, j, k, n
    real*8                                                  :: pch8(3)
    real                                                    :: dx, dy, dz
    real                                                    :: ddx, ddy, ddz
    real                                                    :: ss, sx, sy, sz
    integer, parameter                                      :: pg=0
    integer                                                 :: ierr
    real                                                    :: er_buf

    er = 0.0

    do n=1,nLeaf
      call cpm_GetPitch_LMR( n, pch8, pg, ierr )
      dx = pch8(1)
      dy = pch8(2)
      dz = pch8(3)
      ddx=cf*dt/(dx*dx)
      ddy=cf*dt/(dy*dy)
      ddz=cf*dt/(dz*dz)

      do k=1,nz
      do j=1,ny
      do i=1,nx
        sx= p(i+1,j  ,k  ,n) + p(i-1,j  ,k  ,n) &
          - p(i  ,j  ,k  ,n)*2.0
        sy= p(i  ,j+1,k  ,n) + p(i  ,j-1,k  ,n) &
          - p(i  ,j  ,k  ,n)*2.0
        sz= p(i  ,j  ,k+1,n) + p(i  ,j  ,k-1,n) &
          - p(i  ,j  ,k  ,n)*2.0
        ss=ddx*sx+ddy*sy+ddz*sz
        q(i,j,k,n)=p(i,j,k,n)+ss
        er = er + ss*ss
      end do
      end do
      end do
    end do

    p(1:nx, 1:ny, 1:nz, 1:nLeaf) = q(1:nx, 1:ny, 1:nz, 1:nLeaf)

    call cpm_BndCommS3D_LMR(p,nx,ny,nz,gc,1,CPM_REAL,pg,ierr)
    er_buf = er
    call cpm_Allreduce_LMR(er_buf,er,1,CPM_REAL,CPM_SUM,pg,ierr)
    er = sqrt(er)

    return
    end subroutine ftcs

!   **************************
    subroutine pbc (nx, ny, nz, nLeaf, gc, p)

    implicit none

    include 'cpm_fparam.fi'

    integer                                                 :: nx, ny, nz, nLeaf, gc
    real, dimension(1-gc:nx+gc,1-gc:ny+gc,1-gc:nz+gc,nLeaf) :: p

    real                                                    :: dx, dy, dz
    real*8                                                  :: pch8(3)

    integer                                                 :: head(0:2)
    integer                                                 :: nID(4,0:5), nIDNum(0:5)

    integer                                                 :: i, j, k, n
    real                                                    :: pi, x, y
    integer                                                 :: ierr
    integer, parameter                                      :: pg=0

    pi = 2.0*asin(1.0)

    do n=1,nLeaf
      call cpm_GetPitch_LMR( n, pch8, pg, ierr )
      dx = pch8(1)
      dy = pch8(2)
      dz = pch8(3)
      call cpm_GetVoxelHeadIndex_LMR( n, head, pg, ierr )

      do j=0,5
        call cpm_GetNeighborRankList_LMR( n, j, nID(1,j), nIDNum(j), pg, ierr )
      enddo

      if( nIDNum(Z_MINUS) <= 0 ) then
        do j=1,ny
        do i=1,nx
          x= dx*(i-1+head(X_DIR))+dx*0.5
          y= dy*(j-1+head(Y_DIR))+dy*0.5
          p(i,j,0,n)   =2.0*sin(pi*x)*sin(pi*y)-p(i,j,1,n)
        end do
        end do
      endif
      if( nIDNum(Z_PLUS) <= 0 ) then
        do j=1,ny
        do i=1,nx
          x= dx*(i-1+head(X_DIR))+dx*0.5
          y= dy*(j-1+head(Y_DIR))+dy*0.5
          p(i,j,nz+1,n)=2.0*sin(pi*x)*sin(pi*y)-p(i,j,nz,n)
        end do
        end do
      endif

      if( nIDNum(X_MINUS) <= 0 ) then
        do k=1,nz
        do j=1,ny
          p(0,   j,k,n) = p(1, j,k,n)
        end do
        end do
      endif
      if( nIDNum(X_PLUS) <= 0 ) then
        do k=1,nz
        do j=1,ny
          p(nx+1,j,k,n) = p(nx,j,k,n)
        end do
        end do
      endif

      if( nIDNum(Y_MINUS) <= 0 ) then
        do k=1,nz
        do i=1,nx
          p(i,0,   k,n) = p(i,1, k,n)
        end do
        end do
      endif
      if( nIDNum(Y_PLUS) <= 0 ) then
        do k=1,nz
        do i=1,nx
          p(i,ny+1,k,n) = p(i,ny,k,n)
        end do
        end do
      endif
    end do


    return
    end subroutine pbc

!   ******************************
    subroutine fileout (nx, ny, nz, nLeaf, gc, p, w)
    implicit none
    integer                                                   :: nx, ny, nz, nLeaf, gc
    real  , dimension(1-gc:nx+gc,1-gc:ny+gc,1-gc:nz+gc,nLeaf) :: p
    real*4, dimension(0:nx+1,0:ny+1,0:nz+1)                   :: w

    integer                                                   :: i,j,k,n
    integer                                                   :: step
    real*4                                                    :: org(3), pch(3), time
    character*256                                             :: ename
    integer, parameter                                        :: pg=0
    integer                                                   :: ierr, leafID
    real*8                                                    :: rorg(3), rpch(3)

    step = 0
    time = 0.0

    do n=1,nLeaf
      call cpm_GetLeafID_LMR( n, leafID, pg, ierr )
      write(ename,'("e_id",I4.4,".sph")') leafID

      call cpm_GetLocalOrigin_LMR( n, rorg, pg, ierr )
      call cpm_GetPitch_LMR( n, rpch, pg, ierr )
      pch(1) = rpch(1)
      pch(2) = rpch(2)
      pch(3) = rpch(3)
      org(1) = rorg(1) - pch(1) * 0.5
      org(2) = rorg(2) - pch(2) * 0.5
      org(3) = rorg(3) - pch(3) * 0.5

      do k=0,nz+1
      do j=0,ny+1
      do i=0,nx+1
        w(i,j,k) = real(p(i,j,k,n))
      enddo
      enddo
      enddo

      open (unit=22,file=ename,form='unformatted')
      write (22) 1, 1
      write (22) nx+2,ny+2,nz+2
      write (22) org
      write (22) pch
      write (22) step, time
      write (22) w
      close (unit=22)
    enddo

    return
    end subroutine fileout
