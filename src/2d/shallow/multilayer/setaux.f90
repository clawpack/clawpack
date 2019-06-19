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
!
subroutine setaux(mbc,mx,my,xlow,ylow,dx,dy,maux,aux)

    use amr_module, only: mcapa, xupper, yupper, xlower, ylower, NEEDS_TO_BE_SET
    use amr_module, only: xperdom, yperdom

    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad
    use geoclaw_module, only: sea_level, ambient_pressure

    use storm_module, only: wind_forcing, pressure_forcing
    use storm_module, only: wind_index, pressure_index, set_storm_fields

    use friction_module, only: variable_friction, friction_index
    use friction_module, only: set_friction_field

    use multilayer_module, only: eta_init, num_layers, aux_layer_index

    use topo_module

    use adjoint_module, only : adjoint_flagging,innerprod_index
    
    implicit none
    
    ! Arguments
    integer, intent(in) :: mbc,mx,my,maux
    real(kind=8), intent(in) :: xlow,ylow,dx,dy
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: ii,jj,i,j,m,iint,jint,ilo,jlo
    real(kind=8) :: x,y,xm,ym,xp,yp,topo_integral
    real(kind=8) :: xper,yper,xperm,yperm,xperp,yperp
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

    ! Compute integer indices based off same corner of domain to reduce round
    ! off discrepancies
    ilo = floor((xlow - xlower + .05d0*dx)/dx)
    jlo = floor((ylow - ylower + .05d0*dy)/dy)

    ! Set geometry values
    if (coordinate_system == 1) then
        if (maux == 0 .and. mcapa > 0) then
            print *, "ERROR:  Capacity array requested but number of aux"
            print *, "variables is set to 0."
            stop
        end if
    else if (coordinate_system == 2) then
        do jj = 1 - mbc, my + mbc
            do ii = 1 - mbc, mx + mbc
                ym = ylower + (jlo+jj-1.d0) * dy
                yp = ylower + (jlo+jj) * dy
                aux(2,ii,jj) = deg2rad * earth_radius**2                      &
                                * (sin(yp * deg2rad) - sin(ym * deg2rad)) / dy
                aux(3,ii,jj) = ym * deg2rad
            end do
        end do
    end if

    ! Storm fields if used
    if (wind_forcing) then
        aux(wind_index, :, :) = 0.d0
        aux(wind_index + 1, :, :) = 0.d0
    endif
    if (pressure_forcing) then
        aux(pressure_index, :, :) = ambient_pressure
    endif

    ! Innerproduct field if used
    if (adjoint_flagging) then
        do jj=1-mbc,my+mbc
            do ii=1-mbc,mx+mbc
                if (aux(1,ii,jj) .eq. NEEDS_TO_BE_SET) then
                    aux(innerprod_index,ii,jj) = 0.d0
                endif
            enddo
        enddo
    endif

    ! Initial layer depths for multilayer
    aux(aux_layer_index:num_layers - 1 + aux_layer_index,:,:) = 0.d0

    ! ==========================================================================
    !  Set Bathymetry
    ! Set analytical bathymetry here if requested
    if (test_topography > 0) then
        do jj = 1 - mbc, my + mbc
            ym = ylower + (jlo+jj-1.d0) * dy
            yp = ylower + (jlo+jj) * dy
            y = 0.5d0*(ym+yp)
            do ii = 1 - mbc, mx + mbc
                xm = xlower + (ilo+ii-1.d0) * dx
                xp = xlower + (ilo+ii) * dx
                x = 0.5d0*(xm+xp)
                aux(1,ii,jj) = test_topo(x, y)
            end do
        end do

    ! Set bathymetry based on input files
    else if (mtopofiles > 0) then

        do jj=1-mbc,my+mbc
            ym = ylower + (jlo+jj-1.d0) * dy
            yp = ylower + (jlo+jj) * dy
            y = 0.5d0*(ym+yp)


            do ii=1-mbc,mx+mbc
                xm = xlower + (ilo+ii-1.d0) * dx
                xp = xlower + (ilo+ii) * dx
                x = 0.5d0*(xm+xp)

                ! Parameter NEEDS_TO_BE_SET initialized in amr_module.f90
                ! saves time by otherwise copying instead of reinitializing
                if (aux(1,ii,jj) .ne. NEEDS_TO_BE_SET) then
                    cycle
                endif

                topo_integral = 0.d0
                if ((y>yupper).or.(y<ylower).or.(x>xupper).or.(x<xlower)) then
                    if (.not.(xperdom .or. yperdom)) then
                        ! Skip setting as this cell sticks out of the physical
                        ! domain and we are not setting periodic BCs
                        cycle

                    else
                        ! We evaluate the periodic BC topography by computing
                        ! the appropriate periodic grid cell coordinates and
                        ! again calling cellgridintegrate
                        call wrap_coords(x,y,xperm,xper,xperp,yperm,yper, &
                                         yperp,dx,dy)
                        call cellgridintegrate(topo_integral,  &
                            xperm,xper,xperp,yperm,yper,yperp, &
                            xlowtopo,ylowtopo,xhitopo,yhitopo,dxtopo,dytopo, &
                            mxtopo,mytopo,mtopo,i0topo,mtopoorder, &
                            mtopofiles,mtoposize,topowork)
                    endif
                else
                    ! Cell does not extend outside of physical domain
                    call cellgridintegrate(topo_integral,xm,x,xp,ym,y,yp, &
                            xlowtopo,ylowtopo,xhitopo,yhitopo,dxtopo,dytopo, &
                            mxtopo,mytopo,mtopo,i0topo,mtopoorder, &
                            mtopofiles,mtoposize,topowork)
                endif

                ! Correct for geometry
                if (coordinate_system == 2) then
                    aux(1,ii,jj) = topo_integral / (dx * dy * aux(2,ii,jj))
                else
                    aux(1,ii,jj) = topo_integral / (dx * dy)
                endif
            enddo
        enddo
    else
        print *, "ERROR:  There is no way to set bathymetry!  Either "
        print *, "        provide topography files or request topography "
        print*,  "        defined by a function."
        stop
    end if

    ! Copy topo to ghost cells if outside physical domain and not periodic
    if (.not. yperdom) then
        do jj=1-mbc,my+mbc
            y = ylower + (jlo+jj-.5d0) * dy
            if ((y < ylower) .or. (y>yupper)) then
                do ii=1-mbc,mx+mbc
                    x = xlower + (ilo+ii-.5d0) * dx
                    iint = ii + max(0, ceiling((xlower-x)/dx)) &
                         - max(0, ceiling((x-xupper)/dx))
                    jint = jj + max(0, ceiling((ylower-y)/dy)) &
                         - max(0, ceiling((y-yupper)/dy))
                    aux(1,ii,jj) = aux(1,iint,jint)
                enddo
            endif
        enddo
    endif
    if (.not. xperdom) then
        do ii=1-mbc,mx+mbc
            x =  xlower + (ilo+ii-.5d0) * dx
            if ((x < xlower) .or. (x > xupper)) then
                do jj=1-mbc,my+mbc
                    y = ylower + (jlo+jj-.5d0) * dy
                    iint = ii + max(0, ceiling((xlower-x)/dx)) &
                         - max(0, ceiling((x-xupper)/dx))
                    jint = jj + max(0, ceiling((ylower-y)/dy)) &
                         - max(0, ceiling((y-yupper)/dy))
                    aux(1,ii,jj) = aux(1,iint,jint)
                enddo
            endif
        enddo
    endif

    ! Set friction coefficient based on a set of depth levels
    if (friction_index > 0) then
       call set_friction_field(mx, my, mbc, maux, xlow, ylow, dx, dy, aux)
    endif

    ! Record initial depths if using multiple layers
    do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
            do m=1,num_layers-1
                if (eta_init(m) > aux(1,i,j)) then
                    if (eta_init(m+1) > aux(1,i,j)) then
                        ! There's a layer below this one
                        aux(aux_layer_index + (m - 1), i, j) =      &
                                                 eta_init(m) - eta_init(m+1)
                    else
                        ! This is the last wet layer
                        aux(aux_layer_index + (m - 1), i, j) =      &
                                                    eta_init(m) - aux(1,i,j)
                    endif
                else
                    ! This layer is dry here
                    aux(aux_layer_index + (m - 1), i, j) = 0.d0
                endif
            enddo    
            ! Handle bottom layer seperately
            if (eta_init(num_layers) > aux(1,i,j)) then
                ! Bottom layer is wet here
                aux(aux_layer_index + num_layers - 1,i,j) =             &
                                           eta_init(num_layers) - aux(1,i,j)        
            else
                ! Bottom layer is dry here
                aux(aux_layer_index + num_layers - 1,i,j) = 0.d0
            endif
        enddo
    enddo

contains

    ! Provide wrapper function for providing periodic coordinates
    subroutine wrap_coords(x,y,xperm,xper,xperp,yperm,yper,yperp,dx,dy)

        use amr_module, only: xperdom, yperdom, xupper, yupper, xlower, ylower

        implicit none

        ! Arguments
        real(kind=8), intent(in) :: x,y,dx,dy
        real(kind=8), intent(out) :: xperm,xper,xperp,yperm,yper,yperp
        real(kind=8) :: xdomain, ydomain

        ! need size of domain for wrapping
        xdomain = xupper - xlower
        ydomain = yupper - ylower
        xper = x
        yper = y

        ! test for wrapping of coordinates
        if (x > xupper .and. xperdom) xper = xper - xdomain
        if (x < xlower .and. xperdom) xper = xper + xdomain
        if (y > yupper .and. yperdom) yper = yper - ydomain
        if (y < ylower .and. yperdom) yper = yper + ydomain

        ! adjust rest of variables
        xperm = xper - 0.5d0 * dx
        xperp = xper + 0.5d0 * dx
        yperm = yper - 0.5d0 * dy
        yperp = yper + 0.5d0 * dy

    end subroutine wrap_coords
    
end subroutine setaux
