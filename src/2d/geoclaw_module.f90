! ============================================================================
!  File:        geoclaw_mod
! ============================================================================
!    Copyright (C) 2010-04-21 Clawpack Developers http://www.clawpack.org
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

module geoclaw_module

    implicit none
    save
    
    ! ========================================================================
    ! General geoclaw parameters
    ! ========================================================================
    real(kind=8) :: grav,Rearth,R1,R2,pi
    real(kind=8) :: coeffmanning,depthdeep
    real(kind=8) :: frictiondepth
    integer :: icoordsys,maxleveldeep,minlevelwet,maxleveldry
    integer :: icoriolis,ifriction
    integer, parameter :: GEO_PARM_UNIT = 78
    logical            :: varRefTime = .FALSE.
    
    ! Refinement criteria
    real(kind=8), allocatable :: dry_tolerance(:)
    real(kind=8), allocatable :: wave_tolerance(:)
    real(kind=8), allocatable :: speed_tolerance(:)

    ! ========================================================================
    !  Flow grades flagging support
    ! ========================================================================
    double precision, allocatable :: flowgradevalue(:)

    integer, allocatable :: iflowgradevariable(:), iflowgradetype(:)
    integer, allocatable :: iflowgrademinlevel(:)
    integer :: mflowgrades
    
    ! ========================================================================
    !  Multi-layer support
    ! ========================================================================
    integer :: num_layers
    double precision, allocatable :: rho(:)
    double precision, allocatable :: eta_init(:)
    
    ! Multilayer method Parameters
    integer :: eigen_method,inundation_method

    ! Loss of hyperbolicity 
    integer, parameter :: KAPPA_UNIT = 42
    logical :: check_richardson
    double precision :: richardson_tolerance

contains

    ! ========================================================================
    !  set_geo(fname)
    ! ========================================================================
    !  Reads in user parameters from the given file name if provided
    ! ========================================================================
    subroutine set_geo(fname)

        use amr_module, only: mcapa
        implicit none

        ! Input
        character*25, intent(in), optional :: fname

        ! Locals
        integer, parameter :: iunit = 127
        integer :: igrav
        character*25 :: file_name
        logical :: found_file

        ! Common block  !!! replaced by use amr_module above !!!
!       integer :: mcapa
!       common /cmcapa/ mcapa

        open(unit=GEO_PARM_UNIT,file='fort.geo',status="unknown",action="write")

        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'SETGEO:'
        write(GEO_PARM_UNIT,*) '---------'

        ! Earth Radius params needed for latitude longitude grid
        pi = 4.d0*datan(1.d0)

        ! Read user parameters from setgeo.data
        if (present(fname)) then
            file_name = fname
        else
            file_name = 'setgeo.data'
        endif
        inquire(file=file_name,exist=found_file)
        if (.not. found_file) then
            print *, 'You must provide a file ', fname
            stop
        endif

        call opendatafile(iunit, file_name)

        read(iunit,*) igrav

        ! igrav = 1 means next line contains value of g = gravitational const.
        if (igrav /= 1) then
            print *, 'ERROR -- only igrav = 1 is supported in setgeo'
            stop
        endif
        read(iunit,*) grav

        read(iunit,*) icoordsys
        read(iunit,*) icoriolis
        read(iunit,*) Rearth
        read(iunit,*) varRefTime
        close(iunit)

        ! icoordsys = 1 means Cartesian grid in meters
        ! icoordsys = 2 means lat-long grid on sphere
        ! Check that icoordsys is consistent with mcapa:
        if ((icoordsys > 1) .and. (mcapa == 0)) then
            print *, 'ERROR in setgeo:  if icoordsys > 1 then'
            print *, '      mcapa should be nonzero'
            stop
        endif
        if ((icoordsys == 1) .and. (mcapa > 0)) then
            print *, 'ERROR in setgeo:  if icoordsys = 1 then'
            print *, '      mcapa should be zero'
            stop
        endif

        write(GEO_PARM_UNIT,*) '  igravity:',igrav
        write(GEO_PARM_UNIT,*) '   gravity:',grav
        write(GEO_PARM_UNIT,*) '   icoordsys:',icoordsys

        write(GEO_PARM_UNIT,*) '   varRefTime:',varRefTime

    end subroutine set_geo

    ! ========================================================================
    !  set_multilayer(fname)
    ! ========================================================================
    subroutine set_multilayer(fname)
        implicit none
        character(len=25), optional, intent(in) :: fname
        
        ! Locals
        character(len=25) :: file_name
        logical :: found_file
        integer :: i,ios
        integer, parameter :: unit = 124

        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'SET_MULTILAYER:'
        write(GEO_PARM_UNIT,*) '-----------'

        if (present(fname)) then
            file_name = fname
        else
            file_name = 'multilayer.data'
        endif
        inquire(file=file_name,exist=found_file)
        if (.not. found_file) then
            print *,'You must provide a file ', file_name
            stop
        endif

        call opendatafile(unit, file_name)

        read(unit,"(i2)") num_layers
        allocate(rho(num_layers))
        read(unit,*) rho
        allocate(eta_init(num_layers))
        read(unit,*) eta_init
        read(unit,*) check_richardson
        read(unit,"(d16.8)") richardson_tolerance
        read(unit,"(i1)") eigen_method
        read(unit,"(i1)") inundation_method
        close(unit) 
        
        ! Calculate ratios of densities
!         do i=1,num_layers-1
!             r(i) = rho(i) / rho(i+1)
!         enddo
        
        write(GEO_PARM_UNIT,*) '   num_layers:',num_layers
        write(GEO_PARM_UNIT,*) '   rho:',(rho(i),i=1,num_layers)
!         write(GEO_PARM_UNIT,*) '   r (calculated):', (r,i=1,num_layers-1)

        ! Open Kappa output file if num_layers > 1
        ! Open file for writing hyperbolicity warnings if multiple layers
        if (num_layers > 1) then
            open(unit=KAPPA_UNIT, file='fort.kappa', iostat=ios, &
                    status="unknown", action="write")
            if ( ios /= 0 ) stop "Error opening file name fort.kappa"
        endif
        
    end subroutine set_multilayer

    ! ========================================================================
    !  set_shallow(fname)
    ! ========================================================================
    !  Reads in user parameters from the given file name if provided
    ! ========================================================================
    subroutine set_shallow(fname)

        use amr_module, only: mxnest

        implicit none

        ! Input arguments
        character*25, optional, intent(in) :: fname

        ! Locals
        character*25 :: file_name
        logical :: found_file
        integer :: i
        integer, parameter :: iunit = 127

        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'SET_SHALLOW:'
        write(GEO_PARM_UNIT,*) '-----------'

        if (present(fname)) then
            file_name = fname
        else
            file_name = 'setshallow.data'
        endif
        inquire(file=file_name,exist=found_file)
        if (.not. found_file) then
            print *,'You must provide a file ', file_name
            stop
        endif

        call opendatafile(iunit, file_name)

        allocate(dry_tolerance(num_layers))
        read(iunit,*) dry_tolerance
        allocate(wave_tolerance(num_layers))
        read(iunit,*) wave_tolerance
        allocate(speed_tolerance(mxnest))
        read(iunit,*) (speed_tolerance(i),i=1,mxnest)
        read(iunit,*) depthdeep
        read(iunit,*) maxleveldeep
        read(iunit,*) ifriction
        read(iunit,*) coeffmanning
        read(iunit,*) frictiondepth
        close(iunit)

        write(GEO_PARM_UNIT,*) '   dry_tolerance:',dry_tolerance
        write(GEO_PARM_UNIT,*) '   wave_tolerance:',wave_tolerance
        write(GEO_PARM_UNIT,*) '   speed_tolerance:',speed_tolerance
        write(GEO_PARM_UNIT,*) '   maxleveldeep:', maxleveldeep
        write(GEO_PARM_UNIT,*) '   depthdeep:', depthdeep
        write(GEO_PARM_UNIT,*) '   ifriction:', ifriction
        write(GEO_PARM_UNIT,*) '   Manning coefficient:',coeffmanning
        write(GEO_PARM_UNIT,*) '   frictiondepth:',frictiondepth

    end subroutine set_shallow
    ! ========================================================================


    ! ========================================================================
    !  set_flow_grades(fname)
    ! ========================================================================
    subroutine set_flow_grades(fname)

        implicit none

        ! Input arguments
        character*25, optional, intent(in) :: fname

        ! Locals
        integer, parameter :: iunit = 127
        integer :: i
        character*25 :: file_name
        logical :: found_file

        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'SET FLOW GRADES:'
        write(GEO_PARM_UNIT,*) '------------'

        ! Read user parameters from setflowgrades.data
        if (present(fname)) then
            file_name = fname
        else
            file_name = 'setflowgrades.data'
        endif
        inquire(file=file_name,exist=found_file)
        if (.not. found_file) then
            print *, 'You must provide a file ', fname
            print *, 'Or comment out call setflowgrades in setprob'
            stop
        endif

        call opendatafile(iunit, fname)

        read(iunit,*) mflowgrades

        if (mflowgrades == 0) then
            write(GEO_PARM_UNIT,*) '  No flow grades specified'
            return
        endif

        ! Allocate arrays
        allocate(flowgradevalue(mflowgrades),iflowgradevariable(mflowgrades))
        allocate(iflowgradetype(mflowgrades),iflowgrademinlevel(mflowgrades))

        do i=1,mflowgrades
            read(iunit,*) flowgradevalue(i),iflowgradevariable(i), &
                iflowgradetype(i),iflowgrademinlevel(i)
        enddo

        close(iunit)

        write(GEO_PARM_UNIT,*) '   mflowgrades:',  mflowgrades

        do i=1,mflowgrades
            write(GEO_PARM_UNIT,"(d12.3,3i4)") flowgradevalue(i), &
                iflowgradevariable(i),iflowgradetype(i),iflowgrademinlevel(i)

        enddo

    end subroutine set_flow_grades

end module geoclaw_module
