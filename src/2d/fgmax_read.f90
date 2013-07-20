
subroutine fgmax_read()

    ! Read in data file describing any fixed grids.
    ! The file is assumed to have the form:
    ! 
    ! FG_num_fgrids  # number of fixed grids
    ! 1           # fgridno: number of first fixed grid.  
    ! tstart_max,tend_max # start and end time for monitoring this fgrid.
    ! dt_for_max  # desired maximum time between updating max values.
    ! min_level_for_max # minimum level to use for monitoring max.
    ! tstart_output,tend_output # start and end time for output on this fgrid.
    ! dt_for_max  # desired maximum time between output.
    ! min_level_mor_output # minimum level to output results
    ! 
    ! npts        # number of grid points
    ! x(1), y(1)  # first grid point
    ! ...
    ! x(npts), y(npts)  # last grid point
    ! <repeat for additional grids from fgridno line>

    use fgmax_module
    ! Note: should use mxnest in place of FG_AMR_MAX_LEVELS from above module

    implicit none
    character(30) :: fname = 'setfixedgrids2.data'
    integer :: k,ifg
    type(fgrid), pointer :: fg
    logical :: foundFile

    inquire(file=trim(fname),exist=foundFile)
    if (.not. foundFile) then
      write(*,*) 'No fixedgrid file in new style '
      write(*,*) 'Looking for: ',trim(fname)
      FG_num_fgrids = 0
      return
      endif

    open(unit=FG_UNIT,file=trim(fname),status='old')
    read(FG_UNIT,*) FG_num_fgrids
    print *, 'FG_num_fgrids = ',FG_num_fgrids
    if (FG_num_fgrids > FG_MAXNUM_FGRIDS) then
        write(6,601) FG_num_fgrids
 601    format('*** Too many fixed grids specified: FG_num_fgrids = i3',/, &
               '*** Increase FG_MAXNUM_FGRIDS in fgmax_module.f90')
        stop
        endif

    do ifg = 1,FG_num_fgrids
        fg => FG_fgrids(ifg)   ! point to next element of array of fgrids
        !read(FG_UNIT,*) fg%fgno
        read(FG_UNIT,*) fg%tstart_max, fg%tend_max
        read(FG_UNIT,*) fg%dt_for_max
        read(FG_UNIT,*) fg%min_level_for_max
        read(FG_UNIT,*) fg%tstart_output, fg%tend_output
        read(FG_UNIT,*) fg%dt_for_output
        read(FG_UNIT,*) fg%min_level_for_output
        !read(FG_UNIT,*) fg%num_output
        !if (fg%num_output > 0) then
        !    allocate(fg%t_output(1:fg%num_output))
        !    read(FG_UNIT,*) fg%t_output(1:fg%num_output)
        !    endif
        read(FG_UNIT,*) fg%npts
        allocate(fg%x(1:fg%npts), fg%y(1:fg%npts))
        do k=1,fg%npts
            read(FG_UNIT,*) fg%x(k), fg%y(k)
            enddo

        ! allocate and initialize arrays
        allocate(fg%valuemax(1:FG_NUM_VAL, 1:fg%npts))
        allocate(fg%levelmax(1:fg%npts))
        allocate(fg%aux(1:FG_AMR_MAX_LEVELS, 1:FG_NUM_AUX, 1:fg%npts))
        allocate(fg%tmax(1:FG_NUM_VAL, 1:fg%npts))
        allocate(fg%t_last_updated(1:fg%npts))
        allocate(fg%arrival_time(1:fg%npts))
        fg%valuemax = FG_NOTSET
        fg%levelmax = 0
        fg%aux = FG_NOTSET
        fg%tmax = FG_NOTSET
        fg%t_last_updated = FG_NOTSET
        fg%arrival_time = FG_NOTSET
        !print *, '+++ fg%aux in read: ',fg%aux

        ! Set corners of bounding box.
        fg%x1bb = minval(fg%x)
        fg%x2bb = maxval(fg%x)
        fg%y1bb = minval(fg%y)
        fg%y2bb = maxval(fg%y)

        print *, '++++ bounding box in read:'
        print *, fg%x1bb,fg%x2bb,fg%y1bb,fg%y2bb
        !stop
        enddo

    close(FG_UNIT)


end subroutine fgmax_read
