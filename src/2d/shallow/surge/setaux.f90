subroutine setaux(mbc,mx,my,xlow,ylow,dx,dy,maux,aux)
!     ============================================
!
!     # set auxiliary arrays
!
!     aux(1,i,j) = Z(x,y) topography (negative below sea level for topoymetry)
!
!     If coordinate_system=2 then lat-lon coordinates on the sphere and
!        aux(2,i,j) = area ratio (capacity function -- set mcapa = 2)
!        aux(3,i,j) = length ratio for edge
!     
!     aux(4:num_layers + 3,i,j) = Initial layer depths for linearized problem
!     aux(wind_index:wind_index+1, i ,j) = Wind Field (x,y)
!     aux(pressure_index, i , j) = Pressure Field

    use amr_module, only: mcapa, xupper, yupper, xlower, ylower

    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad
    use geoclaw_module, only: sea_level

    use storm_module, only: wind_index, pressure_index, set_storm_fields
    use storm_module, only: ambient_pressure

    use friction_module, only: friction_index, set_friction_field
    
    use topo_module
    
    implicit none
    
    ! Arguments
    integer, intent(in) :: mbc,mx,my,maux
    real(kind=8), intent(in) :: xlow,ylow,dx,dy
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j,m,iint,jint
    real(kind=8) :: x,y,xm,ym,xp,yp,topo_integral
    character(len=*), parameter :: aux_format = "(2i4,4d15.3)"
    
    ! Lat-Long coordinate system in use, check input variables
    if (coordinate_system == 2) then
        if (mcapa /= 2 .or. maux < 3) then
            print *,'ERROR in setaux:  for coordinate_system=2'
            print *,'     need mcapa = 2 and maux >= 3'
            print *,'     have mcapa = ',mcapa,'  maux = ',maux
            stop
        endif
    endif
    
    ! Set default values for aux variables
    aux(1,:,:) = 0.d0 ! Bathymetry
    aux(2,:,:) = 1.d0 ! Grid cell area
    aux(3,:,:) = 1.d0 ! Length ratio for edge
    aux(friction_index,:,:) = 0.d0 ! Friction field

    ! Set these to something non-offensive
    aux(wind_index,:,:) = 0.d0 ! Wind speed x-direction
    aux(wind_index+1,:,:) = 0.d0 ! Wind speed y-direction
    aux(pressure_index,:,:) = ambient_pressure ! Pressure field

    ! Set analytical bathymetry here if requested
    if (test_topography > 0) then
        forall (i=1-mbc:mx+mbc,j=1-mbc:my+mbc)
            aux(1,i,j) = test_topo(xlow + (i - 0.5d0) * dx, &
                                   ylow + (j - 0.5d0) * dy)
        end forall
    endif

    ! Set bathymetry
    do j=1-mbc,my+mbc
        ym = ylow + (j - 1.d0) * dy
        y = ylow + (j - 0.5d0) * dy
        yp = ylow + real(j,kind=8) * dy
        do i=1-mbc,mx+mbc
            xm = xlow + (i - 1.d0) * dx
            x = xlow + (i - 0.5d0) * dx
            xp = xlow + real(i,kind=8) * dx

            ! Set lat-long cell info
            if (coordinate_system == 2) then
                aux(2,i,j) = deg2rad * earth_radius**2 * (sin(yp * deg2rad) - sin(ym * deg2rad)) / dy
                aux(3,i,j) = ym * deg2rad
            endif
            
            ! skip setting aux(1,i,j) in ghost cell if outside physical domain
            ! since topo files may not cover ghost cell, and values
            ! should be extrapolated, which is done in next set of loops.
            if ((y>yupper) .or. (y<ylower) .or. &
                (x>xupper) .or. (x<xlower)) cycle


            ! Use input topography files if available
            if (mtopofiles > 0 .and. test_topography == 0) then
                topo_integral = 0.d0
                call cellgridintegrate(topo_integral,xm,x,xp,ym,y,yp, &
                    xlowtopo,ylowtopo,xhitopo,yhitopo,dxtopo,dytopo, &
                    mxtopo,mytopo,mtopo,i0topo,mtopoorder, &
                    mtopofiles,mtoposize,topowork)

                    aux(1,i,j) = topo_integral / (dx * dy * aux(2,i,j))
            endif
        enddo
    enddo

    ! Copy topo to ghost cells if outside physical domain
    do j=1-mbc,my+mbc
        y = ylow + (j-0.5d0) * dy
        if ((y < ylower) .or. (y>yupper)) then
            do i=1-mbc,mx+mbc
                x = xlow + (i-0.5d0) * dx 
                iint = i + max(0, ceiling((xlower-x)/dx)) &
                         - max(0, ceiling((x-xupper)/dx))
                jint = j + max(0, ceiling((ylower-y)/dy)) &
                         - max(0, ceiling((y-yupper)/dy))
                aux(1,i,j) = aux(1,iint,jint)
            enddo
        endif
    enddo
    do i=1-mbc,mx+mbc
        x = xlow + (i-0.5d0) * dx
        if ((x < xlower) .or. (x > xupper)) then
            do j=1-mbc,my+mbc
                y = ylow + (j-0.5d0) * dy 
                iint = i + max(0, ceiling((xlower-x)/dx)) &
                         - max(0, ceiling((x-xupper)/dx))
                jint = j + max(0, ceiling((ylower-y)/dy)) &
                         - max(0, ceiling((y-yupper)/dy))
                aux(1,i,j) = aux(1,iint,jint)
            enddo
        endif
    enddo

    ! Set friction coefficient based on a set of depth levels
    call set_friction_field(mx,my,mbc,maux,xlow,ylow,dx,dy,aux)

    ! Output bathymetry for debugging
    if (.false.) then
        open(23, file='fort.aux',status='unknown',form='formatted')
        print *,'Writing out aux arrays'
        print *,' '
        do j=1,my
            do i=1,mx
                write(23,aux_format) i,j,(aux(m,i,j),m=1,maux)
            enddo
        enddo
        close(23)
    endif
    
end subroutine setaux