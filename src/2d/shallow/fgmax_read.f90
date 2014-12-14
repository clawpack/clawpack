
subroutine fgmax_read(fname,ifg)

    ! Read in data file describing any fixed grids.
    ! The file is assumed to have the form:
    ! 
    ! tstart_max  # start time for monitoring this fgrid.
    ! tend_max    # end time for monitoring this fgrid.
    ! dt_check    # desired maximum time increment between updating max values.
    ! min_level_check # minimum level to check for monitoring values/arrivals
    ! arrival_tol           # tolerance for identifying arrival.
    ! 
    ! point_style   # 0 ==> list of points, 1 ==> 1d transect,  2 ==> 2d grid, 3 ==> quadrilateral
    ! if point_style==0 this is followed by:
    !   npts        # number of grid points
    !   x(1), y(1)  # first grid point
    !   ...
    !   x(npts), y(npts)  # last grid point
    ! if point_style==1:
    !   npts
    !   x1, y1     # first point
    !   x2, y2     # last point
    ! if point_style==2:
    !   nx, ny
    !   x1, y1     # lower left corner of cartesian grid
    !   x2, y2     # upper right corner of cartesian grid
    ! if point_style==3:
    !   n12, n23   # number of points on sides 1-2 (also 3-4) and 2-3 (also 4-1)
    !   x1, y1     # first corner of quadrilateral
    !   x2, y2     # second corner
    !   x3, y3     # third corner
    !   x4, y4     # fourth corner
    

    use fgmax_module
    use amr_module, only: mxnest

    implicit none
    character(150), intent(in) :: fname
    integer, intent(in) :: ifg 
    integer :: k,i,j,point_style,nx,ny,n12,n23
    real(kind=8) :: x1,x2,y1,y2,yj
    real(kind=8) :: x3,x4,y3,y4,x14,y14,x23,y23,xi,eta
    type(fgrid), pointer :: fg
    logical :: foundFile

    open(unit=FG_UNIT,file=trim(fname),status='old')
    inquire(file=trim(fname),exist=foundFile)
    if (.not. foundFile) then
      write(*,*) 'Missing fgmax file...'
      write(*,*) 'Looking for: xx',trim(fname),'xx'
      stop
      endif

    open(unit=FG_UNIT,file=trim(fname),status='old')

    fg => FG_fgrids(ifg)   ! point to next element of array of fgrids
    read(FG_UNIT,*) fg%tstart_max
    read(FG_UNIT,*) fg%tend_max
    read(FG_UNIT,*) fg%dt_check  
    read(FG_UNIT,*) fg%min_level_check
    read(FG_UNIT,*) fg%arrival_tol
    read(FG_UNIT,*) point_style
    if (point_style == 0) then
        read(FG_UNIT,*) fg%npts
        allocate(fg%x(1:fg%npts), fg%y(1:fg%npts))
        do k=1,fg%npts
            read(FG_UNIT,*) fg%x(k), fg%y(k)
            enddo
    else if (point_style == 1) then
        read(FG_UNIT,*) fg%npts
        read(FG_UNIT,*) x1,y1
        read(FG_UNIT,*) x2,y2
        allocate(fg%x(1:fg%npts), fg%y(1:fg%npts))
        do k=1,fg%npts
            fg%x(k) = x1 + (k-1)*(x2-x1)/(fg%npts - 1)
            fg%y(k) = y1 + (k-1)*(y2-y1)/(fg%npts - 1)
            enddo
    else if (point_style == 2) then
        read(FG_UNIT,*) nx, ny
        read(FG_UNIT,*) x1,y1
        read(FG_UNIT,*) x2,y2
        fg%npts = nx*ny
        allocate(fg%x(1:fg%npts), fg%y(1:fg%npts))
        k = 0
        do j=1,ny
            yj = y1 + (j-1)*(y2-y1)/(ny - 1)
            do i=1,nx
                k = k+1
                fg%x(k) = x1 + (i-1)*(x2-x1)/(nx - 1)
                fg%y(k) = yj
                enddo
            enddo
    else if (point_style == 3) then
        read(FG_UNIT,*) n12, n23
        read(FG_UNIT,*) x1,y1
        read(FG_UNIT,*) x2,y2
        read(FG_UNIT,*) x3,y3
        read(FG_UNIT,*) x4,y4
        fg%npts = n12*n23
        allocate(fg%x(1:fg%npts), fg%y(1:fg%npts))
        
        k = 0
        do j=1,n23
            eta = (j-1)/(n23-1.d0)
            x14 = x1 + eta*(x4-x1)
            y14 = y1 + eta*(y4-y1)
            x23 = x2 + eta*(x3-x2)
            y23 = y2 + eta*(y3-y2)
            do i=1,n12
                k = k+1
                xi = (i-1)/(n12-1.d0)
                fg%x(k) = x14 + xi*(x23 - x14)
                fg%y(k) = y14 + xi*(y23 - y14)
                enddo
            enddo
    else
        write(6,*) '*** Unexpected value of point_style in fgmax_read: ',point_style
        write(6,*) '*** Note that format of fgmax file has changed, see documentation'
        stop
        endif


    ! allocate and initialize arrays
    allocate(fg%valuemax(1:FG_NUM_VAL, 1:fg%npts))
    allocate(fg%levelmax(1:fg%npts))
    allocate(fg%aux(1:mxnest, 1:FG_NUM_AUX, 1:fg%npts))
    allocate(fg%auxdone(1:mxnest))
    allocate(fg%tmax(1:FG_NUM_VAL, 1:fg%npts))
    allocate(fg%t_last_updated(1:fg%npts))
    allocate(fg%arrival_time(1:fg%npts))
    fg%valuemax = FG_NOTSET
    fg%levelmax = 0
    fg%aux = FG_NOTSET
    fg%auxdone = .false.
    fg%tmax = FG_NOTSET
    fg%t_last_updated = FG_NOTSET
    fg%arrival_time = FG_NOTSET

    ! Set corners of bounding box.
    fg%x1bb = minval(fg%x)
    fg%x2bb = maxval(fg%x)
    fg%y1bb = minval(fg%y)
    fg%y2bb = maxval(fg%y)

    print *, fg%x1bb,fg%x2bb,fg%y1bb,fg%y2bb

    close(FG_UNIT)


end subroutine fgmax_read
