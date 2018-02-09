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
!

    use amr_module, only: mcapa, xupper, yupper, xlower, ylower

    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad
    use geoclaw_module, only: sea_level, ambient_pressure

    use storm_module, only: wind_forcing, pressure_forcing
    use storm_module, only: wind_index, pressure_index, set_storm_fields

    use friction_module, only: variable_friction, friction_index
    use friction_module, only: set_friction_field

    use multilayer_module, only: eta_init, num_layers, aux_layer_index

    use topo_module
    
    implicit none
    
    ! Arguments
    integer, intent(in) :: mbc,mx,my,maux
    real(kind=8), intent(in) :: xlow,ylow,dx,dy
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: ii,jj,i,j,m,iint,jint,ilo,jlo
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

    ! Check below is new in 5.2.1 -- need to rethink for storm surge
    ! and other applications where other aux arrays might be used?
    if (coordinate_system == 1) then
        if (mcapa > 0) then
            print *,'ERROR in setaux:  for coordinate_system==1'
            print *,'     need mcapa == 0 and maux == 1'
            print *,'     have mcapa = ',mcapa,'  maux = ',maux
            stop
        else if (maux > 1) then
            ! Should not need to set aux(2,:,:) in this case, but 
            ! for some reason it bombs, e.g. in bowl-radial if maux>1.
            aux(2,:,:) = 1.d0 
            !aux(3,:,:) = 1.d0
        endif
    endif

    ! If using a variable friction field initialize the coefficients to 0
    if (variable_friction) then
        aux(friction_index,:,:) = 0.d0
    endif

    ! Storm fields if used
    if (wind_forcing) then
        aux(wind_index, :, :) = 0.d0
        aux(wind_index + 1, :, :) = 0.d0
    endif
    if (pressure_forcing) then
        aux(pressure_index, :, :) = ambient_pressure
    endif

    ! Initial layer depths for multilayer
    aux(aux_layer_index:num_layers - 1 + aux_layer_index,:,:) = 0.d0 

    ! Set analytical bathymetry here if requested
    if (test_topography > 0) then
        forall (i=1-mbc:mx+mbc,j=1-mbc:my+mbc)
            aux(1,i,j) = test_topo(xlow + (i - 0.5d0) * dx, &
                                   ylow + (j - 0.5d0) * dy)
        end forall
    endif

! test:  compute integer indices based off same corner of domain 
!        to reduce round off discrepancies
    ilo = floor((xlow - xlower + .05d0*dx)/dx)
    jlo = floor((ylow - ylower + .05d0*dy)/dy)

    ! Set bathymetry
    do jj=1-mbc,my+mbc

        !ym = ylow + (j - 1.d0) * dy
        !y = ylow + (j - 0.5d0) * dy
        !yp = ylow + real(j,kind=8) * dy

        ym = ylower + (jlo+jj-1.d0) * dy
        yp = ylower + (jlo+jj) * dy
        y = 0.5d0*(ym+yp)

        do ii=1-mbc,mx+mbc
            !xm = xlow + (i - 1.d0) * dx
            !x = xlow + (i - 0.5d0) * dx
            !xp = xlow + real(i,kind=8) * dx

            xm = xlower + (ilo+ii-1.d0) * dx
            xp = xlower + (ilo+ii) * dx
            x = 0.5d0*(xm+xp)

            ! Set lat-long cell info
            if (coordinate_system == 2) then
                aux(2,ii,jj) = deg2rad * earth_radius**2 * (sin(yp * deg2rad) - sin(ym * deg2rad)) / dy
                aux(3,ii,jj) = ym * deg2rad
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


                if (coordinate_system == 2) then
                    aux(1,ii,jj) = topo_integral / (dx * dy * aux(2,ii,jj))
                else
                    aux(1,ii,jj) = topo_integral / (dx * dy)
                endif

            endif
        enddo
    enddo

    ! Copy topo to ghost cells if outside physical domain

   do jj=1-mbc,my+mbc
        y = ylower + (jlo+jj-.5d0) * dy
        if ((y < ylower) .or. (y>yupper)) then
            do ii=1-mbc,mx+mbc
                x = xlower + (ilo+ii-.5d0) * dx
                iint = ii + max(0, ceiling((xlower-x)/dx)) &
                         -  max(0, ceiling((x-xupper)/dx))
                jint = jj + max(0, ceiling((ylower-y)/dy)) &
                         -  max(0, ceiling((y-yupper)/dy))
                aux(1,ii,jj) = aux(1,iint,jint)
            enddo
        endif
    enddo

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

    ! Output for debugging
    if (.false.) then
        open(23, file='fort.aux',status='unknown',form='formatted')
        print *,'Writing out aux arrays'
        print *,' '
        do j=1,my
            do i=1,mx
                write(23,*) i,j,(aux(m,i,j),m=1,maux)
            enddo
        enddo
        close(23)
    endif
    
end subroutine setaux