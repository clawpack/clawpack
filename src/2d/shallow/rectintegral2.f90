
recursive subroutine rectintegral(x1,x2,y1,y2,m,integral)

    real (kind=8), intent(in) :: x1,x2,y1,y2
    integer, intent(in) :: m
    real (kind=8), intent(out) :: integral

    ! Compute the integral of topo over the rectangle (x1,x2) x (y1,y2)
    ! using topo arrays indexed 1 through m (coarse to fine).

    ! The main call to this subroutine has corners of a grid cell for the 
    ! rectangle and m = mtopofiles in order to compute the integral 
    ! over the cell using all topo arrays.

    ! The recursive strategy is to first compute the integral using only topo 
    ! arrays 1 to m-1 and then apply corrections due to adding topo array m.
    ! Corrections are needed if the topo array m intersects the grid cell.
    ! Let the intersection be (x1m,x2m) x (y1m,y2m).
    ! Two corrections are needed, first to subtract out the integral over
    ! the rectangle (x1m,x2m) x (y1m,y2m) computed using
    ! topo arrays 1 to m-1 and then adding in the integral over this 
    ! same region using topo array m.

    ! Assume the function topointegral(x1,x2,y1,y2,m) returns the 
    ! integral over the rectangle based on a single topo array with index m.

    if (m == 1) then
        integral = topointegral(x1,x2,y1,y2,m)
    else
        call rectintegral(x1,y1,x2,y2,m-1,int1)
        call intersection(...)  ! to compute x1m,x2m,y1m,y2m
        call rectintegral(x1m,x2m,y1m,y2m,m-1,int2)
        int3 = topointegral(x1,x2,y1,y2,m)
        integral = int1 - int2 + int3
    endif

end subroutine rectintegral

    

