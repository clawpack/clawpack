! ============================================================================
!  Program:     /Users/mandli/src/local/shallow_water/2d/multi_layer
!  File:        multilayer_module
!  Created:     2010-10-12
!  Author:      Kyle Mandli
! ============================================================================
!      Copyright (C) 2010-10-12 Kyle Mandli <mandli@amath.washington.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

module multilayer_module

    implicit none
    
    ! Physical parameters
    integer :: num_layers
    double precision, allocatable :: rho(:)
    double precision :: r,one_minus_r
    
    ! Algorithm parameters
    integer :: eigen_method,inundation_method
    double precision :: richardson_tolerance
    double precision, allocatable :: wave_tol(:)
    logical :: dry_limit
    
    ! Initial layer depths
!     double precision, allocatable :: eta(:)
!     double precision :: epsilon,sigma,init_location(2),angle
!     integer :: init_type,wave_family
    
    ! Simple bathy states
!     integer :: bathy_type
!     double precision :: bathy_left, bathy_right, bathy_location
!     double precision :: bathy_x0,bathy_x1,bathy_x2,bathy_basin_depth
!     double precision :: bathy_shelf_depth,bathy_shelf_slope,bathy_beach_slope
    
    ! Output files
    integer, parameter :: kappa_file = 42
    
    ! HACK for maux passing
    integer :: ml_maux
    
contains

    subroutine set_multilayer_params(data_file)

        implicit none
        character(len=*), optional, intent(in) :: data_file
        
        integer :: ios
        
        double precision :: A,B,step_height
        
        ! Open file
        if (present(data_file)) then
            open(unit=13,file=data_file,iostat=ios,status='old',action='read',access='sequential')            
            if ( ios /= 0 ) then
                print *,'Error opening "',data_file,'" for reading.'
                stop
            endif
        else
            open(unit=13,file='multilayer.data',iostat=ios,status='old',action='read',access='sequential')            
            if ( ios /= 0 ) then
                print *,'Error opening "multilayer.data" for reading.'
                stop
            endif
        endif
        
        dry_limit = .true.
    
        ! Physics parameters
        read(13,"(i3)") num_layers
        allocate(rho(num_layers))
        allocate(eta_init(num_layers))
        allocate(wave_tol(num_layers))
        read(13,*) rho
        if (num_layers > 1) then
            r = rho(1) / rho(2)
            one_minus_r = 1.d0 - r
        else
            r = -1.d0
            one_minus_r = 0.d0
        endif
        read(13,*)
        
        ! Algorithmic parameters
        read(13,"(i1)") eigen_method
        read(13,"(i1)") inundation_method
        read(13,"(d16.8)") richardson_tolerance
        read(13,*) dry_limit
        read(13,*)
        
        ! Initial conditions
!         read(13,*) eta
!         read(13,*) init_type
!         if (init_type > 0) then
!             read(13,*) epsilon
!             if (init_type <= 2 .or. init_type == 5) then  
!                 read(13,*) init_location
!                 read(13,*) wave_family
!                 if(init_type == 2 .or. init_type == 5) then
!                     read(13,*) angle
!                     read(13,*) sigma
!                 endif
!             else if (init_type == 3) then
!                 read(13,*) init_location
!                 read(13,*) sigma
!             endif
!         endif
!         read(13,*)
        
        ! Bathymetry
!         read(13,*) bathy_type
!         if (bathy_type == 0) then
!             continue
!         else if (bathy_type == 1) then 
!             read(13,"(d16.8)") bathy_location
!             read(13,"(d16.8)") bathy_left
!             read(13,"(d16.8)") bathy_right
!         else if (bathy_type == 2 .or. bathy_type == 3) then
!             read(13,"(d16.8)") bathy_x0
!             read(13,"(d16.8)") bathy_x1
!             read(13,"(d16.8)") bathy_x2
!             read(13,"(d16.8)") bathy_basin_depth
!             read(13,"(d16.8)") bathy_shelf_depth
!             read(13,"(d16.8)") bathy_beach_slope
!             bathy_shelf_slope = (bathy_basin_depth - bathy_shelf_depth) &
!                                         / (bathy_x0 - bathy_x1)
!         else
!             print *,"Invalid bathymetry type ",bathy_type
!             stop
!         endif
        
        close(13)

        ! Open file for writing hyperbolicity warnings
        open(unit=kappa_file, file='fort.kappa', iostat=ios, &
                status="unknown", action="write")
        if ( ios /= 0 ) stop "Error opening file name fort.kappa"
        
    end subroutine set_multilayer_params

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
        close(GEO_PARM_UNIT)
        
    end subroutine read_multilayer_data

end module multilayer_module
