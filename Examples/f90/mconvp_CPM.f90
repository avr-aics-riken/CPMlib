!   *********************************************************
!   ***      Sample code of porting
!   ***      2010 keno@VCAD, RIKEN
!
!   *********************************************************
    module array_def
    real, dimension(:,:,:), allocatable, save :: p, q, w
    end module array_def

!   *********************************************************
    program diffusion

    use array_def
    implicit none

    include 'mpif.h'
    include 'cpm_fparam.fi'

    integer            :: nx, ny, nz
    integer            :: n, scheme, itrMax, stp, stpMax
    real               :: dx, dy, dz
    real               :: eps, er, omg, dt, cf

    integer, parameter :: pg=0
    integer            :: ierr
    integer            :: nrank, myrank, div(3)
    integer            :: iwork(6)
    real               :: rwork(7)
    integer            :: sz(3), head(0:2), nID(0:5)
    real               :: org(3), pch(3), rgn(3)
    real*8             :: org8(3), rgn8(3)


    ! MPI_Init
    call MPI_Init(ierr)

    ! initialize CPM library
    call cpm_Initialize(ierr)

    ! get number of rank
    call cpm_GetNumRank(nrank, pg, ierr)

    ! get my rank No.
    call cpm_GetMyRankID(myrank, pg, ierr)

!   Initialize

    ! user input for rank 0
    if( myrank == 0 ) then
    write (*,*) 'number of process=', nrank

    write (*,*) 'Dimension size'
    write (*,*) 'Input dimension x size ='
    read  (*,*)  nx
    write (*,*) 'Input dimension y size ='
    read  (*,*)  ny
    write (*,*) 'Input dimension z size ='
    read  (*,*)  nz

    write (*,*) 'SELECT time integration method'
    write (*,*) '        1 --- Euler Explicit'
    write (*,*) '        2 --- Euler Implicit'
    read  (*,*) scheme
    
    write (*,*) 'Time step max='
    read  (*,*) stpMax
    
    write (*,*) 'Time increment='
    read  (*,*) dt
    
    write (*,*) 'Coef. of diffusion='
    read  (*,*) cf
      
    if ( scheme == 2 ) then
      write (*,*) 'Iteration max='
      read  (*,*) itrMax
      write (*,*) 'Epsilon ='
      read  (*,*) eps
      write (*,*) 'Relazation coef. ='
      read  (*,*) omg
    end if

    dx = 1.0/real(nx)
    dy = 1.0/real(ny)
    dz = 1.0/real(nz)

    write (*,*) 'nx, ny, nz = ', nx, ny, nz 
    write (*,*) 'dx, dy, dz = ', dx, dy, dz 
    write (*,*) 'dt         = ', dt
    if ( scheme == 1 ) then
      write(*,*) 'Euler Explicit ftcs'
    else
      write(*,*) 'Euler Implicit jacobi'
    end if

    endif

    ! broadcast user input
    iwork(1) = nx
    iwork(2) = ny
    iwork(3) = nz
    iwork(4) = scheme
    iwork(5) = stpMax
    iwork(6) = itrMax
    rwork(1) = dt
    rwork(2) = cf
    rwork(3) = eps
    rwork(4) = omg
    rwork(5) = dx
    rwork(6) = dy
    rwork(7) = dz
    call cpm_Bcast( iwork, 6, CPM_INT , 0, pg, ierr )
    call cpm_Bcast( rwork, 7, CPM_REAL, 0, pg, ierr )
    sz(1)  = iwork(1)
    sz(2)  = iwork(2)
    sz(3)  = iwork(3)
    scheme = iwork(4)
    stpMax = iwork(5)
    itrMax = iwork(6)
    dt     = rwork(1)
    cf     = rwork(2)
    eps    = rwork(3)
    omg    = rwork(4)
    pch(1) = rwork(5)
    pch(2) = rwork(6)
    pch(3) = rwork(7)
    org    = 0.d0

    ! calc region size
    rgn(1) = pch(1) * sz(1)
    rgn(2) = pch(2) * sz(2)
    rgn(3) = pch(3) * sz(3)

    ! Voxel initialize
    org8(1) = org(1)
    org8(2) = org(2)
    org8(3) = org(3)
    rgn8(1) = rgn(1)
    rgn8(2) = rgn(2)
    rgn8(3) = rgn(3)
    call cpm_VoxelInit_nodiv( sz, org8, rgn8, 1, 1, 0, pg, ierr )

    ! get local voxel size, ead index, neighbor rank No.
    call cpm_GetLocalVoxelSize( sz, pg, ierr )
    call cpm_GetVoxelHeadIndex( head, pg, ierr )
    call cpm_GetNeighborRankID( nID, pg, ierr )
    call cpm_GetDivNum( div, pg, ierr )
    nx = sz(1)
    ny = sz(2)
    nz = sz(3)
    dx = pch(1)
    dy = pch(2)
    dz = pch(3)
    write(*,'(A,/,A,I4,/,A,3(1X,I4),/,A,3(1X,I4),/,A,6(1X,I4),/,A,3(1X,I4))') &
      "+--------------------------------------+" &
    , "myrank = ", myrank &
    , "  local voxel   =", sz &
    , "  head index    =", head &
    , "  neighbor rank =", nID &
    , "  div num       =", div
!    stop

    call allocate_array (nx, ny, nz)

    call pbc (nx, ny, nz, p, dx, dy, dz, head, nID)

    
!   time step loop
    do stp=1, stpMax
      if ( scheme == 1 ) then
        call ftcs (nx, ny, nz, p, q, cf, dx, dy, dz, dt, er)

        call pbc (nx, ny, nz, p, dx, dy, dz, head, nID)

        if( myrank == 0 ) &

        write (*,*) stp, er
      else
        er = 0.0
        w(0:nx+1, 0:ny+1, 0:nz+1) = p(0:nx+1, 0:ny+1, 0:nz+1)
        ! iteration
        do n=1, itrMax
          call jacobi (nx, ny, nz, p, q, w, cf, dx, dy, dz, dt, omg, er)

          call pbc (nx, ny, nz, p, dx, dy, dz, head, nID)

          if (er<eps) exit
        end do

        if( myrank == 0 ) &

        write (*,*) stp, n, er
      end if
    end do

!   Post
    call fileout (nx, ny, nz, p, w, dx, dy, dz)


    stop
    end program diffusion

!   ******************************
    subroutine allocate_array (nx, ny, nz)
    use array_def
    implicit none
    integer :: status
    integer :: nx, ny, nz
    
    allocate ( p(0:nx+1, 0:ny+1, 0:nz+1) , stat=status )
    allocate ( w(0:nx+1, 0:ny+1, 0:nz+1) , stat=status )
    allocate ( q(0:nx+1, 0:ny+1, 0:nz+1) , stat=status )
    p(0:nx+1, 0:ny+1, 0:nz+1) = 0.d0

    return
    end subroutine allocate_array

!   ****************************************************
    subroutine jacobi (nx, ny, nz, p, q, w, cf, dx, dy, dz, dt, omg, er)
    implicit none

    include 'cpm_fparam.fi'

    integer                                 :: nx, ny, nz 
    real, dimension (0:nx+1,0:ny+1,0:nz+1)  :: p, q, w
    real                                    :: cf, dx, dy, dz, dt, omg, er

    integer                                 :: i, j, k
    real                                    :: ddx, ddy, ddz, cpd
    real                                    :: s0, ss, sx, sy, sz
    integer, parameter                      :: pg=0
    integer                                 :: ierr
    real                                    :: er_buf
    
    er = 0.0
    ddx = cf*dt/(dx*dx)
    ddy = cf*dt/(dy*dy)
    ddz = cf*dt/(dz*dz)
    cpd=1.0/(1.0+2.0*ddx+2.0*ddy+2.0*ddz)
    
    do k=1,nz
    do j=1,ny
    do i=1,nx
      sx= p(i+1,j  ,k  ) + p(i-1,j  ,k  )
      sy= p(i  ,j+1,k  ) + p(i  ,j-1,k  )
      sz= p(i  ,j  ,k+1) + p(i  ,j  ,k-1)
      s0=ddx*sx+ddy*sy+ddz*sz
      ss=(cpd*(s0+w(i,j,k))-p(i,j,k))*omg
      q(i,j,k)=p(i,j,k)+ss
      er = er + ss*ss
    end do
    end do
    end do


    p(1:nx, 1:ny, 1:nz) = q(1:nx, 1:ny, 1:nz)
    call cpm_BndCommS3D(p,nx,ny,nz,1,1,CPM_REAL,pg,ierr)
    er_buf = er
    call cpm_Allreduce(er_buf,er,1,CPM_REAL,CPM_SUM,pg,ierr)
    er = sqrt(er)
    
    return
    end subroutine jacobi
    
!   ***************************************
    subroutine ftcs (nx, ny, nz, p, q, cf, dx, dy, dz, dt, er)
    implicit none
    include 'cpm_fparam.fi'

    integer                                 :: nx, ny, nz
    real, dimension (0:nx+1,0:ny+1,0:nz+1)  :: p, q
    real                                    :: cf, dx, dy, dz, dt, er

    integer                                 :: i, j, k
    real                                    :: ddx, ddy, ddz
    real                                    :: ss, sx, sy, sz
    integer, parameter                      :: pg=0
    integer                                 :: ierr
    real                                    :: er_buf

    er = 0.0
    ddx=cf*dt/(dx*dx)
    ddy=cf*dt/(dy*dy)
    ddz=cf*dt/(dz*dz)
    
    do k=1,nz
    do j=1,ny
    do i=1,nx
      sx= p(i+1,j  ,k  ) + p(i-1,j  ,k  ) &
        - p(i  ,j  ,k  )*2.0
      sy= p(i  ,j+1,k  ) + p(i  ,j-1,k  ) &
        - p(i  ,j  ,k  )*2.0
      sz= p(i  ,j  ,k+1) + p(i  ,j  ,k-1) &
        - p(i  ,j  ,k  )*2.0
      ss=ddx*sx+ddy*sy+ddz*sz
      q(i,j,k)=p(i,j,k)+ss
      er = er + ss*ss
    end do
    end do
    end do

    p(1:nx, 1:ny, 1:nz) = q(1:nx, 1:ny, 1:nz)


    call cpm_BndCommS3D(p,nx,ny,nz,1,1,CPM_REAL,pg,ierr)
    er_buf = er
    call cpm_Allreduce(er_buf,er,1,CPM_REAL,CPM_SUM,pg,ierr)

    er = sqrt(er)
    
    return
    end subroutine ftcs
    
!   **************************
    subroutine pbc (nx, ny, nz, p, dx, dy, dz, head, nID)

    implicit none

    include 'cpm_fparam.fi'

    integer                               :: nx, ny, nz
    real, dimension(0:nx+1,0:ny+1,0:nz+1) :: p
    real                                  :: dx, dy, dz

    integer                               :: head(0:2), nID(0:5)


    integer                               :: i, j, k
    real                                  :: pi, x, y

    pi = 2.0*asin(1.0)

    if( nID(Z_MINUS) < 0 ) then
      do j=1,ny
      do i=1,nx
        x= dx*(i-1+head(X_DIR))+dx*0.5
        y= dy*(j-1+head(Y_DIR))+dy*0.5
        p(i,j,0)   =2.0*sin(pi*x)*sin(pi*y)-p(i,j,1)
      end do
      end do
    endif
    if( nID(Z_PLUS) < 0 ) then
      do j=1,ny
      do i=1,nx
        x= dx*(i-1+head(X_DIR))+dx*0.5
        y= dy*(j-1+head(Y_DIR))+dy*0.5
        p(i,j,nz+1)=2.0*sin(pi*x)*sin(pi*y)-p(i,j,nz)
      end do
      end do
    endif

    if( nID(X_MINUS) < 0 ) then
      do k=1,nz
      do j=1,ny
        p(0,   j,k) = p(1, j,k)
      end do
      end do
    endif
    if( nID(X_PLUS) < 0 ) then
      do k=1,nz
      do j=1,ny
        p(nx+1,j,k) = p(nx,j,k)
      end do
      end do
    endif

    if( nID(Y_MINUS) < 0 ) then
      do k=1,nz
      do i=1,nx
        p(i,0,   k) = p(i,1, k)
      end do
      end do
    endif
    if( nID(Y_PLUS) < 0 ) then
      do k=1,nz
      do i=1,nx
        p(i,ny+1,k) = p(i,ny,k)
      end do
      end do
    endif

    
    return
    end subroutine pbc

!   ******************************
    subroutine fileout (nx, ny, nz, p, w, dx, dy, dz)
    implicit none
    integer                                 :: nx, ny, nz
    real  , dimension(0:nx+1,0:ny+1,0:nz+1) :: p
    real*4, dimension(0:nx+1,0:ny+1,0:nz+1) :: w
    real                                    :: dx, dy, dz

    integer                                 :: i,j,k
    integer                                 :: step
    real*4                                  :: org(3), pch(3), time
    character*256                           :: ename
    integer, parameter                      :: pg=0
    integer                                 :: ierr, myrank
    real                                    :: rorg(3), rpch(3)


    call cpm_GetMyRankID( myrank, pg, ierr )
    write(ename,'("e_id",I4.4,".sph")') myrank

    step = 0
    time = 0.0

    call cpm_GetLocalOrigin( rorg, pg, ierr )
    call cpm_GetPitch( rpch, pg, ierr )
    pch(1) = rpch(1)
    pch(2) = rpch(2)
    pch(3) = rpch(3)
    org(1) = rorg(1) - pch(1)
    org(2) = rorg(2) - pch(2)
    org(3) = rorg(3) - pch(3)

    do k=0,nz+1
    do j=0,ny+1
    do i=0,nx+1
      w(i,j,k) = real(p(i,j,k))
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

    return
    end subroutine fileout

