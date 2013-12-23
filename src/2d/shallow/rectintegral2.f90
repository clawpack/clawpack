
recursive subroutine rectintegral(x1,x2,y1,y2,m,integral)

    ! Compute the integral of topo over the rectangle (x1,x2) x (y1,y2)
    ! using topo arrays indexed mtopoorder(mtopofiles) through mtopoorder(m) 
    ! (coarse to fine).

    ! The main call to this subroutine has corners of a grid cell for the 
    ! rectangle and m = 1 in order to compute the integral over the cell 
    ! using all topo arrays.

    ! The recursive strategy is to first compute the integral using only topo 
    ! arrays mtopoorder(mtopofiles) to mtopoorder(m+1), 
    ! and then apply corrections due to adding topo array mtopoorder(m).
     
    ! Corrections are needed if the new topo array intersects the grid cell.
    ! Let the intersection be (x1m,x2m) x (y1m,y2m).
    ! Two corrections are needed, first to subtract out the integral over
    ! the rectangle (x1m,x2m) x (y1m,y2m) computed using
    ! topo arrays mtopoorder(mtopofiles) to mtopoorder(m+1),
    ! and then adding in the integral over this same region using 
    ! topo array mtopoorder(m).

    ! Note that the function topointegral returns the integral over the 
    ! rectangle based on a single topo array, and that routine calls
    ! bilinearintegral.

    use topo_module, only: xlowtopo,ylowtopo,xhitopo,yhitopo,dxtopo,dytopo, &
        mxtopo,mytopo,i0topo,mtopoorder,mtopofiles,topowork

    implicit none

    ! arguments
    real (kind=8), intent(in) :: x1,x2,y1,y2
    integer, intent(in) :: m
    real (kind=8), intent(out) :: integral

    ! local
    real(kind=8) :: xmlo,xmc,xmhi,ymlo,ymc,ymhi,area,x1m,xmm,x2m, &
        y1m,ymm,y2m, int1,int2,int3
    integer :: mfid, indicator, mp1fid, i0
    real(kind=8), external :: topointegral  


    mfid = mtopoorder(m)
    i0=i0topo(mfid)

    !write(44,*) 'in rectintegral with m = ',m,' mfid = ',mfid
    !write(44,*) 'x1,x2,y1,y2: ',x1,x2,y1,y2
    if (m == mtopofiles) then
         ! innermost step of recursion reaches this point.
         ! only using coarsest topo grid -- compute directly...
         call intersection(indicator,area,xmlo,xmc,xmhi, &
             ymlo,ymc,ymhi, x1,x2,y1,y2, &
             xlowtopo(mfid),xhitopo(mfid),ylowtopo(mfid),yhitopo(mfid))

         !write(44,*) 'x1,x2,y1,y2: ',x1,x2,y1,y2
         !write(44,*) 'indicator, area: ',indicator, area
         if (indicator.eq.1) then
            ! cell overlaps the file
            ! integrate surface over intersection of grid and cell
            !write(44,*) 'calling topointegral with'
            !write(44,*) xmlo,xmc,xmhi,ymlo,ymc, &
            !        ymhi,xlowtopo(mfid),ylowtopo(mfid),dxtopo(mfid), &
            !        dytopo(mfid),mxtopo(mfid),mytopo(mfid)
            integral = topointegral( xmlo,xmc,xmhi,ymlo,ymc, &
                    ymhi,xlowtopo(mfid),ylowtopo(mfid),dxtopo(mfid), &
                    dytopo(mfid),mxtopo(mfid),mytopo(mfid),topowork(i0),1)
            write(44,*) 'innermost     : ',m,integral
         else
            integral = 0.d0
         endif

    else
        ! recursive call to compute area using one fewer topo grids:
        write(44,*) '==== recursive call with m+1 = ',m+1
        call rectintegral(x1,x2,y1,y2,m+1,int1)
        write(44,*) '>>>> return from call with m+1 = ',m+1,'  int1 = ',int1

        ! region of intersection of cell with new topo grid:
        !!mp1fid = mtopoorder(m+1)
        call intersection(indicator,area,x1m,xmm,x2m, &
             y1m,ymm,y2m, x1,x2,y1,y2, &
             xlowtopo(mfid),xhitopo(mfid),ylowtopo(mfid),yhitopo(mfid))
             !!xlowtopo(mp1fid),xhitopo(mp1fid),ylowtopo(mp1fid),yhitopo(mp1fid))

        
        write(44,*) 'area of intersection with grid m = ',m, '  is: ',area
        if (area > 0) then
        
            ! correction to subtract out from previous set of topo grids:
            call rectintegral(x1m,x2m,y1m,y2m,m+1,int2)
            write(44,*) 'computed int2 = ',int2
    
            ! correction to add in for new topo grid:
            !write(44,*) x1m,xmm,x2m, y1m,ymm,y2m, &
            !            xlowtopo(mfid),ylowtopo(mfid),dxtopo(mfid), &
            !            dytopo(mfid),mxtopo(mfid),mytopo(mfid)
            int3 = topointegral(x1m,xmm,x2m, y1m,ymm,y2m, &
                        xlowtopo(mfid),ylowtopo(mfid),dxtopo(mfid), &
                        dytopo(mfid),mxtopo(mfid),mytopo(mfid),topowork(i0),1)
            write(44,*) 'computed int3 = ',int3
    
            ! adjust integral due to corrections for new topo grid:
            write(44,*) 'before update : ',m,int1
            integral = int1 - int2 + int3
            write(44,*) 'after  update : ',m,integral
        else
            integral = int1
        endif
    endif

end subroutine rectintegral

    

