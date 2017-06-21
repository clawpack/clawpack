subroutine setprob()

    use adjoint_module, only: read_adjoint_data, set_time_window

    implicit none
    character(len=150) ::  adjointFolder
    real(kind=8) :: t1,t2
    integer :: iunit

    iunit = 7
    call opendatafile(iunit, 'setprob.data')

    read(iunit,*) adjointFolder

    ! time period of interest:
    read(iunit,*) t1
    read(iunit,*) t2

    call set_time_window(t1, t2)                   !# Set time window
    call read_adjoint_data(trim(adjointFolder))    !# Read adjoint solution

end subroutine setprob
