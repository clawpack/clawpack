! This routine should be a simplified version of src2
! which applies source terms for a 1-d slice of data along the
! edge of a grid.  This is called only from qad where the conservative
! fix-up is applied and is used to apply source terms over partial
! time steps to the coarse grid cell values used in solving Riemann
! problems at the interface between coarse and fine grids.
subroutine src1d(meqn,mbc,mx1d,q1d,maux,aux1d,t,dt)

    use multilayer_module, only: rho,eta,layers
    use hurricane_module
    use geoclaw_module

    implicit none
    
    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx1d,maux
    double precision, intent(in) :: t,dt
    
    ! Output
    double precision, intent(inout) :: q1d(meqn,mx1d)
    double precision, intent(inout) :: aux1d(maux,mx1d)

    ! Locals
    integer :: i
    double precision :: g,coeff,tol,h(2)
    double precision :: wind_speed,tau,P_atmos_x,P_atmos_y,speed,D

    ! Common block
    double precision dtcom,dxcom,dycom,tcom
    integer icom,jcom
    
    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

    ! Incorporates friction using Manning coefficient
    g=grav
    coeff = coeffmanning
    tol = 1.d-30  !# to prevent divide by zero in gamma

    ! = Friction =============================================================
    if (coeff > 0.d0) then
        do i=1,mx1d
            h(1) = q1d(1,i) / rho(1)
            if (layers > 1) then
                h(2) = q1d(4,i) / rho(2)
            else
                h(2) = 0.d0
            endif
            
            if (h(2) > tol) then
                speed = sqrt(q1d(5,i)**2 + q1d(6,i)**2) / q1d(4,i)
                D = coeff**2 * g * sum(h)**(-7/3) * speed
                
                q1d(5,i) = q1d(5,i) * exp(-D*dt)
                q1d(6,i) = q1d(6,i) * exp(-D*dt)
            else if (h(1) > tol) then
                if (layers > 1) q1d(5:6,i) = 0.d0
                
                speed = sqrt(q1d(2,i)**2 + q1d(3,i)**2) / q1d(1,i)
                
                D = coeff**2 * g * sum(h)**(-7/3) * speed
                
                q1d(2,i) = q1d(2,i) * exp(-D*dt)
                q1d(3,i) = q1d(3,i) * exp(-D*dt)
                
            else
                q1d(2:3,i) = 0.d0
                if (layers > 1)  q1d(5:6,i) = 0.d0
            endif
        enddo          
    endif
    
    ! = Wind Forcing =========================================================
    if (wind_forcing) then
        ! Here we have to take into account geoclaw which we use for the 
        ! single layer case.  It needs to divide by the water's density
        if (layers > 1) then
            do i=1,mx1d
                if (abs(q1d(1,i)) > drytolerance) then
                    wind_speed = sqrt(aux1d(4,i)**2 + aux1d(5,i)**2)
                    if (wind_speed > wind_tolerance) then
                        tau = wind_drag(wind_speed) * rho_air * wind_speed
                        q1d(2,i) = q1d(2,i) + dt * tau * aux1d(4,i)
                        q1d(3,i) = q1d(3,i) + dt * tau * aux1d(5,i)
                    endif
                endif
            enddo
        else
            do i=1,mx1d
                if (abs(q1d(1,i)) > drytolerance) then
                    wind_speed = sqrt(aux1d(4,i)**2 + aux1d(5,i)**2)
                    tau = wind_drag(wind_speed) * rho_air * wind_speed
                    q1d(2,i) = q1d(2,i) + dt * tau * aux1d(4,i) / rho(1)
                    q1d(3,i) = q1d(3,i) + dt * tau * aux1d(5,i) / rho(1)
                endif
            enddo
        endif
    endif
    ! ========================================================================
    
    ! == Pressure Forcing ====================================================
    if (pressure_forcing) then
        stop "Not sure how to proceed, need direction or the right dx or dy."
    endif

end subroutine src1d
