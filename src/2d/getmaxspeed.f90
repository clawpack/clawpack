real(kind=8) pure function get_max_speed(val,mitot,mjtot,nvar,aux,naux,nghost,hx,hy)

    use geoclaw_module, only: dry_tolerance, coordinate_system
    use geoclaw_module, only: grav, earth_radius, DEG2RAD
      
    implicit none
    
    ! Arguments
    integer, intent(in) :: mitot,mjtot,nvar,naux,nghost
    real(kind=8), intent(in) :: hx,hy
    real(kind=8), intent(in) :: val(nvar,mitot,mjtot), aux(naux,mitot,mjtot)
    
    ! Locals
    integer :: i,j
    real(kind=8) :: ymetric,hyphys,xmetric,hxphys,u,v,sig,sp_over_h


    sp_over_h = 0.d0   ! compute max speed over h, since dx may not equal dy
    if (coordinate_system == 2) then
        do j = nghost+1, mjtot-nghost
            ymetric = earth_radius*deg2rad
            hyphys = ymetric*hy

            do i = nghost+1, mitot-nghost
                xmetric = cos(aux(3,i,j)) * earth_radius * DEG2RAD
                hxphys = xmetric * hx
                if (val(1,i,j) > dry_tolerance) then
                    u  = val(2,i,j) / val(1,i,j)
                    v  = val(3,i,j) / val(1,i,j)
                else
                    u = 0.d0
                    v = 0.d0
                endif
                sig = sqrt(grav*val(1,i,j))
                sp_over_h = max((abs(u)+sig)/hxphys,(abs(v)+sig)/hyphys,sp_over_h)
            end do
        end do
    else  ! speeds in cartesian coords, no metrics needed
        do j = nghost+1, mjtot-nghost
            do i = nghost+1, mitot-nghost
                if (val(1,i,j) > dry_tolerance) then
                    u  = val(2,i,j) / val(1,i,j)
                    v  = val(3,i,j) / val(1,i,j)
                else
                    u = 0.d0
                    v = 0.d0
                endif
                sig = sqrt(grav*val(1,i,j))
                sp_over_h = max((abs(u)+sig)/hx,(abs(v)+sig)/hy,sp_over_h)
            end do
        end do
    endif
      
    get_max_speed = sp_over_h

end function get_max_speed
