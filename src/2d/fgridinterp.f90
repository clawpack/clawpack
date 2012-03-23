! This routine interpolates q and aux on a computational grid
! to a fgrid not necessarily aligned with the computational grid
! using bilinear interpolation defined on computation grid

subroutine fgridinterp(fgrid,x_low_fg,y_low_fg, &
                       x_hi_fg,y_hi_fg,dx_fg,dy_fg,mx_fg,my_fg,t,num_vars_fg, &
                       q,meqn, &
                       mx_c,my_c,mbc,dx_c,dy_c,num_vars,x_low_c,y_low_c, &
                       maux,aux, &
                       out_arrival_times,out_surface_max,max_check)
     
    use geoclaw_module, only: drytolerance
     
    implicit none
     
    ! Subroutine arguments
    integer, intent(in) :: mx_fg,my_fg,num_vars_fg
    real(kind=8), intent(in) :: x_low_fg,x_hi_fg,y_low_fg,y_hi_fg,dx_fg,dy_fg
    real(kind=8), intent(inout) :: fgrid(num_vars_fg,1:mx_fg,1:my_fg)
     
    integer, intent(in) :: mbc,mx_c,my_c,num_vars,meqn,maux
    real(kind=8), intent(in) :: x_low_c,y_low_c,dx_c,dy_c
    real(kind=8), intent(in) :: q(meqn,1-mbc:mx_c+mbc,1-mbc:my_c+mbc)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx_c+mbc,1-mbc:my_c+mbc)
     
    integer, intent(in) :: out_arrival_times,out_surface_max,max_check
    real(kind=8), intent(in) :: t
     
    ! Locals
    integer :: bathy_index,arrival_index,i,j
    integer :: eta_index,eta_index_min,eta_index_max,eta_index_now
    real(kind=8), parameter :: arrival_tol = 1d-2
    real(kind=8) :: x_hi_c,y_hi_c,tol
     
    ! NaN for filling empty spots
    integer(kind=16) NaN_descriptor
    real(kind=8) :: NaN

    data NaN_descriptor/B'01111111100000100000000000000000'/
    NaN = transfer(NaN,NaN_descriptor)
     
    ! Calculate and extract some important data
    x_hi_c = x_low_c + dx_c * mx_c
    y_hi_c = y_low_c + dy_c * my_c
    
    tol = drytolerance
    
    ! Calculate some of the storage indices
    bathy_index = meqn + 1
    eta_index = meqn + 2
    
    if (max_check > 0) then
        arrival_index = 0
        eta_index_min = 0
        eta_index_max = 0
        eta_index_now = 0
        
        if (out_arrival_times > 0) then
            arrival_index = 1
        endif
        
        if (out_surface_max > 0) then
            eta_index_now = arrival_index + 1
            eta_index_min = arrival_index + 2
            eta_index_max = arrival_index + 3
        endif
    endif
        
    do i=1,mx_fg
        do j=1,my_fg
            
        enddo
    enddo
end subroutine fgridinterp