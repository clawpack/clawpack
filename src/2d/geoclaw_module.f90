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
    !  Constants
    ! ========================================================================
    integer, parameter :: GEO_PARM_UNIT = 78
    integer, parameter :: KAPPA_UNIT = 42
    real(kind=8), parameter :: pi = 4.d0*datan(1.d0)
    real(kind=8), parameter :: DEG2RAD = pi / 180.d0
    real(kind=8), parameter :: RAD2DEG = 180.d0 / pi
    
    ! ========================================================================
    !  Physics
    ! ========================================================================
    real(kind=8) :: grav, earth_radius
    integer :: coordinate_system

    ! Rotational velocity of Earth
    real(kind=8), parameter :: omega = 2.0d0 * pi / 86164.2d0
    
    ! Forcing
    logical :: coriolis_forcing ! True then coriolis terms included in src
    real(kind=8) :: theta_0 ! Used if using the beta-plane approximation
    integer :: friction_forcing ! If 0, no forcing, 1 constant N, 2 variable
    real(kind=8) :: manning_coefficient, friction_depth
    
    ! Method parameters    
    real(kind=8), allocatable :: dry_tolerance(:)
    logical :: varRefTime = .FALSE. ! Choose dt refinement automatically
    
    ! ========================================================================
    !  Multi-layer support
    ! ========================================================================
    integer :: num_layers
    real(kind=8), allocatable :: rho(:)
    real(kind=8), allocatable :: eta_init(:)
    
    ! Multilayer method Parameters
    integer :: eigen_method,inundation_method

    ! Loss of hyperbolicity
    logical :: check_richardson
    real(kind=8) :: richardson_tolerance

contains

    ! ========================================================================
    !  set_geo(fname)
    ! ========================================================================
    !  Reads in user parameters from the given file name if provided
    ! ========================================================================
    subroutine set_geo(file_name)

        use amr_module, only: mcapa
        implicit none

        ! Input
        character(len=*), intent(in), optional :: file_name

        ! Locals
        integer, parameter :: unit = 127

        open(unit=GEO_PARM_UNIT,file='fort.geo',status="unknown",action="write")

        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'Physics Parameters:'
        write(GEO_PARM_UNIT,*) '-------------------'

        ! Read user parameters from setgeo.data
        if (present(file_name)) then
            call opendatafile(unit, file_name)
        else
            call opendatafile(unit, 'physics.data')
        endif

        read(unit,*) grav
        read(unit,*) earth_radius
        read(unit,*) coordinate_system
        read(unit,*)
        read(unit,*) num_layers
        allocate(rho(num_layers),eta_init(num_layers),dry_tolerance(num_layers))
        read(unit,*) rho
        read(unit,*) eta_init
        read(unit,*)
        read(unit,*) coriolis_forcing
        if (coordinate_system == 1 .and. coriolis_forcing) then
            read(unit,*) theta_0
            print *,'here'
        else
            theta_0 = 0.d0
        endif
        read(unit,*) friction_forcing
        if (friction_forcing == 0) then
            manning_coefficient = 0.d0
            friction_depth = 1.d10
        else if (friction_forcing == 1) then
            read(unit,*) manning_coefficient
            read(unit,*) friction_depth
        else if ( friction_forcing > 1 ) then
            stop "Friction forcing > 1 unhandled."
        endif
        read(unit,*)
        read(unit,*) dry_tolerance
        read(unit,*) varRefTime
        
        close(unit)

        ! coordinate_system = 1 means Cartesian grid in meters
        ! coordinate_system = 2 means lat-long grid on sphere
        ! Check that coordinate_system is consistent with mcapa:
        if ((coordinate_system > 1) .and. (mcapa == 0)) then
            print *, 'ERROR in setgeo:  if coordinate_system > 1 then'
            print *, '      mcapa should be nonzero'
            stop
        endif
        if ((coordinate_system == 1) .and. (mcapa > 0)) then
            print *, 'ERROR in setgeo:  if coordinate_system = 1 then'
            print *, '      mcapa should be zero'
            stop
        endif

        write(GEO_PARM_UNIT,*) '   gravity:',grav
        write(GEO_PARM_UNIT,*) '   earth_radius:',earth_radius
        write(GEO_PARM_UNIT,*) '   coordinate_system:',coordinate_system
        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '   num_layers:',num_layers
        write(GEO_PARM_UNIT,*) '   rho:',rho
        write(GEO_PARM_UNIT,*) '   eta_init:',eta_init
        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '   coriolis_forcing:',coriolis_forcing
        write(GEO_PARM_UNIT,*) '   theta_0:',theta_0
        write(GEO_PARM_UNIT,*) '   friction_forcing:',friction_forcing
        write(GEO_PARM_UNIT,*) '   manning_coefficient:',manning_coefficient
        write(GEO_PARM_UNIT,*) '   friction_depth:',friction_depth
        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '   dry_tolerance:',dry_tolerance
        write(GEO_PARM_UNIT,*) '   Variable dt Refinement Ratios:',varRefTime

    end subroutine set_geo

    ! ========================================================================
    !  read_multilayer_data(file_name)
    ! ========================================================================
    subroutine read_multilayer_data(file_name)

        implicit none
        
        ! Arguments
        character(len=*), optional, intent(in) :: file_name
        
        ! Locals
        integer, parameter :: unit = 124
        integer :: ios

        ! Only read in this data if we are doing multilayer swe
        if (num_layers > 1) then
            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '--------------------------------------------'
            write(GEO_PARM_UNIT,*) 'Multilayer Parameters:'
            write(GEO_PARM_UNIT,*) '----------------------'

            if (present(file_name)) then
                call opendatafile(unit, file_name)
            else
                call opendatafile(unit, 'multilayer.data')
            endif

            read(unit,*) check_richardson
            read(unit,"(d16.8)") richardson_tolerance
            read(unit,"(i1)") eigen_method
            read(unit,"(i1)") inundation_method
            close(unit) 

            ! Open Kappa output file if num_layers > 1
            ! Open file for writing hyperbolicity warnings if multiple layers
            if (num_layers > 1 .and. check_richardson) then
                open(unit=KAPPA_UNIT, file='fort.kappa', iostat=ios, &
                        status="unknown", action="write")
                if ( ios /= 0 ) stop "Error opening file name fort.kappa"
            endif

            write(GEO_PARM_UNIT,*) '   check_richardson:',check_richardson
            write(GEO_PARM_UNIT,*) '   richardson_tolerance:',richardson_tolerance
            write(GEO_PARM_UNIT,*) '   eigen_method:',eigen_method
            write(GEO_PARM_UNIT,*) '   inundation_method:',inundation_method
        endif
        
    end subroutine read_multilayer_data

    ! ========================================================================
    !  Calculate the coriolis constant f
    !   If coordinate_system == 1 then
    !       A beta-plane approximation is used and y should be in meters
    !   if coordinate_system == 2 then
    !       Grid is in lat-long and y should be in degrees which is converted
    !       to radians
    ! ========================================================================
    real(kind=8) pure function coriolis(y)

        implicit none
        
        ! Input
        real(kind=8), intent(in) :: y
        
        ! Locals
        real(kind=8) :: theta
        
        ! Assume beta plane approximation and y is in meters    
        if (coordinate_system == 1) then
            theta = y / 111d3 * DEG2RAD + theta_0
            coriolis = 2.d0 * omega * (sin(theta_0) + (theta - theta_0)     &
                                                    * cos(theta_0))
        else if (coordinate_system == 2) then        
            theta = y * DEG2RAD
            coriolis = 2.d0 * omega * sin(theta)
        else
            ! Unknown coordinate system, return nothing
            coriolis = 0.d0
        endif
    end function coriolis

end module geoclaw_module
