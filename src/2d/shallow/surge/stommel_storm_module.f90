! ==============================================================================
! stommel_storm_module 
!
! Module contains routines for constructing a wind field that enacts the Stommel
! test case.
! ==============================================================================
!  Copyright (C) Kyle T. Mandli <kyle@ices.utexas.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================

module stommel_storm_module

    implicit none
    save

    ! Holland storm type definition
    type stommel_storm_type

        real(kind=8) :: A

    end type stommel_storm_type

    logical, private :: module_setup = .false.
    logical, private, parameter :: DEBUG = .true.

contains

    ! Setup routine for the holland model
    subroutine set_stommel_storm(storm_data_path, storm, log_unit)

        use geoclaw_module, only: coordinate_system

        implicit none

        ! Subroutine I/O
        character(len=*), intent(in), optional :: storm_data_path
        type(stommel_storm_type), intent(in out) :: storm
        integer, intent(in) :: log_unit

        ! Local storage
        integer, parameter :: unit = 701

        if (.not.module_setup) then

            ! Storm type does not work for lat-long coordinate systems
            if (coordinate_system == 2) then
                stop "Stommel storm type does not work on lat-long coordinates."
            endif

            ! Open data file
            if (present(storm_data_path)) then
                call opendatafile(unit,storm_data_path)
            else
                call opendatafile(unit,'storm.data')
            endif

            ! Read in hurricane parameters
            read(unit,*) storm%A

            ! Output to log file
            write(log_unit,*) "Storm Data - Stommel"
            write(log_unit,"('Stormmel Wind Field Amplitude A =',d16.8)") storm%A

            module_setup = .true.
        end if

    end subroutine set_stommel_storm

    ! ==========================================================================
    !  set_stommel_storm_fields()
    ! ==========================================================================
    subroutine set_stommel_storm_fields(maux, mbc, mx, my, &
                                        xlower, ylower, dx, dy, t,aux, &
                                        wind_index, pressure_index, storm)

        use amr_module, only: ylower_domain => ylower, yupper_domain => yupper
        use geoclaw_module, only: pi

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(stommel_storm_type), intent(in out) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Locals
        integer :: j
        real(kind=8) :: y, L

        ! Domain size
        L = yupper_domain - ylower_domain
        
        ! Set wind field
        aux(wind_index + 1,:,:) = 0.d0
        aux(pressure_index,:,:) = 0.d0
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy
            
            ! Set wind
            aux(wind_index,:,j) = -storm%A * cos(pi * y / L)
        enddo
        
    end subroutine set_stommel_storm_fields

end module stommel_storm_module 