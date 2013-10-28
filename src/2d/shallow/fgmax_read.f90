
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
    ! npts        # number of grid points
    ! x(1), y(1)  # first grid point
    ! ...
    ! x(npts), y(npts)  # last grid point
    ! <repeat for additional grids from fgridno line>

    use fgmax_module
    use amr_module, only: mxnest

    implicit none
    character(80), intent(in) :: fname
    integer, intent(in) :: ifg 
    integer :: k
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
    read(FG_UNIT,*) fg%npts
    allocate(fg%x(1:fg%npts), fg%y(1:fg%npts))
    do k=1,fg%npts
        read(FG_UNIT,*) fg%x(k), fg%y(k)
        enddo

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
