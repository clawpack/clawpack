module fgmax_module

    save

    ! need to fix this!  Should be mxnest from the computation...
    integer, parameter :: FG_AMR_MAX_LEVELS = 7

    type fgrid

        ! Derived data type for a single fixed grid.
        ! An array of objects of this type is used to keep track of all
        ! fixed grids in a computation.

        ! identifier number for this fixed grid:
        !integer :: fgno ! deprecated

        ! fixed grid points are (x(k),y(k)), for k=1:npts
        integer :: npts
        real(kind=8), allocatable, dimension(:) :: x,y

        ! arrays valuemax(mv,k), aux(level,ma,k), levelmax(k)
        ! where mv = 1:num_values  # number of quantities monitored
        !       ma = 1:num_aux
        !        k = 1:npts
        real(kind=8), allocatable, dimension(:,:,:) :: aux
        real(kind=8), allocatable, dimension(:,:) :: valuemax
        real(kind=8), allocatable, dimension(:,:) :: tmax
        real(kind=8), allocatable, dimension(:) :: t_last_updated
        real(kind=8), allocatable, dimension(:) :: arrival_time
        integer, allocatable, dimension(:) :: levelmax

        ! time range to monitor:
        real(kind=8) :: tstart_max, tend_max
        ! Desired maximum delta t between updating valuemax:
        ! Will only update at start of time step if end of timestep is
        ! to far beyond last update time
        real(kind=8) :: dt_for_max
        ! Mininum level to check when updating valuemax:
        ! Coarser levels will be ignored to save computation time.
        integer :: min_level_for_max


        ! time range for output:
        real(kind=8) :: tstart_output, tend_output
        ! Desired maximum delta t between updating valuemax:
        real(kind=8) :: dt_for_output
        ! Mininum level to create output files
        integer :: min_level_for_output

        !integer :: num_output
        !real(kind=8), allocatable, dimension(:) :: t_output


        ! Coordinates of corners of bounding box.
        ! This will be useful when generalizing to fgrids not aligned with x-y.
        real(kind=8) :: x1bb,x2bb,y1bb,y2bb

        ! keep track of whether all aux arrays have been computed on a given level:
        logical :: auxdone(1:FG_AMR_MAX_LEVELS) = .false.

    end type

    ! declare array fgrids of fixed grids, each of type fgrid.
    ! allow them to be targets of pointers for shorthand in code.
    integer, parameter :: FG_MAXNUM_FGRIDS = 2  ! max number of fixed grids
    type(fgrid), target :: FG_fgrids(FG_MAXNUM_FGRIDS)

    ! special value to flag unset portions of arrays:
    real(kind=8), parameter :: FG_NOTSET = -0.99999d99

    ! unit to use for reading input data and writing fgrid results:
    integer, parameter :: FG_UNIT = 45

    ! number of max vals to monitor
    ! these are specified in fgmax_values
    integer, parameter :: FG_NUM_VAL = 4
    ! number of aux vals to monitor
    integer, parameter :: FG_NUM_AUX = 1

    ! number of fixed grids in use (set by fgmax_read):
    integer :: FG_num_fgrids

    ! turn on debugging output to fort.61:
    logical, parameter :: FG_DEBUG = .false.

    ! depth to use in checking first arrival time:
    real(kind=8), parameter :: FG_HARRIVAL = 0.02d0

    ! level to monitor arrival times 
    ! (should be in type def and read in for each grid!)
    integer, parameter :: FG_arrival_level = 6


end module fgmax_module
