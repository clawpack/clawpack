! ==============================================================================
! model_storm_module 
!
! Module contains routines for constructing a wind and pressure field based on a
! provided set of data files.
!
! ==============================================================================
!                   Copyright (C) Clawpack Developers 2017
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================
module data_storm_module

    implicit none
    save

    logical, private :: module_setup = .false.
    logical, private :: DEBUG = .false.

    ! Model storm type definition
    type data_storm_type

        ! Dummy variable
        integer :: mine

    end type data_storm_type

contains

    ! Setup routine for model storms
    subroutine set_storm(storm_data_path, storm, model_type, log_unit)

        implicit none

        ! Subroutine I/O
        character(len=*), optional :: storm_data_path
        type(data_storm_type), intent(in out) :: storm
        integer, intent(in) :: model_type, log_unit

        stop "Data-derived storm are not yet implemented!"

        if (.not. module_setup) then

            module_setup = .true.
        end if

    end subroutine set_storm


    ! ==========================================================================
    !  storm_location(t,storm)
    !    Interpolate location of hurricane in the current time interval
    ! ==========================================================================
    function storm_location(t,storm) result(location)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(data_storm_type), intent(in out) :: storm

        ! Output
        real(kind=8) :: location(2)

        stop "Data-derived storm are not yet implemented!"

    end function storm_location

    ! ==========================================================================
    !  storm_direction
    !   Angle off of due north that the storm is traveling
    ! ==========================================================================
    real(kind=8) function storm_direction(t, storm) result(theta)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(data_storm_type), intent(in) :: storm

        stop "Data-derived storm are not yet implemented!"

    end function storm_direction


    ! ==========================================================================
    !  Use the 1980 Holland model to set the storm fields
    ! ==========================================================================
    subroutine set_HWRF_fields(maux, mbc, mx, my, xlower, ylower,    &
                          dx, dy, t, aux, wind_index,           &
                          pressure_index, storm)

  
        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(data_storm_type), intent(in out) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        stop "HWRF data input is not yet implemented!"

    end subroutine set_HWRF_fields

end module data_storm_module






