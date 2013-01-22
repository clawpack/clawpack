! ============================================================================
!      Copyright (C) 2010-10-12 Kyle Mandli <kyle.mandli@gmail.com>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

module multilayer_module

    implicit none
    
    ! Physical parameters
    integer :: num_layers
    real(kind=8), allocatable :: rho(:), eta_init(:)
    real(kind=8) :: r,one_minus_r
    
    ! Algorithm parameters
    integer :: eigen_method,inundation_method
    logical :: check_richardson
    real(kind=8) :: richardson_tolerance
    real(kind=8), allocatable :: wave_tol(:)
    
    ! Output files
    integer, parameter :: KAPPA_UNIT = 42
    
contains

    ! ========================================================================
    !  set_multilayer(file_name)
    ! ========================================================================
    subroutine set_multilayer(data_file)

        use geoclaw_module, only: GEO_PARM_UNIT

        implicit none
        character(len=*), optional, intent(in) :: data_file
        
        integer :: ios
        integer, parameter :: IOUNIT = 41
        
        ! Open file
        if (present(data_file)) then
            call opendatafile(IOUNIT, data_file)
        else
            call opendatafile(IOUNIT, 'multilayer.data')
        endif
        
        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'Multilayer Parameters:'
        write(GEO_PARM_UNIT,*) '----------------------'

        ! Read in parameters
        read(IOUNIT,"(i3)") num_layers
        allocate(rho(num_layers))
        allocate(eta_init(num_layers))
        allocate(wave_tol(num_layers))

        read(IOUNIT,*) rho
        if (num_layers > 1) then
            r = rho(1) / rho(2)
            one_minus_r = 1.d0 - r
        else
            r = -1.d0
            one_minus_r = 0.d0
        endif
        read(IOUNIT, *) eta_init
        read(IOUNIT, *)

        ! Algorithmic parameters
        read(IOUNIT, *) check_richardson
        read(IOUNIT, "(d16.8)") richardson_tolerance
        read(IOUNIT, "(i1)") eigen_method
        read(IOUNIT, "(i1)") inundation_method
        
        close(IOUNIT) 

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
        
    end subroutine set_multilayer

end module multilayer_module
