subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)

    use parameters_module

    implicit none
    
    ! Input arguments
    integer, intent(in) :: maxmx,meqn,mbc,mx,maux
    double precision, intent(in) :: xlower,dx
    
    double precision, intent(inout) :: q(1-mbc:maxmx+mbc,meqn)
    double precision, intent(inout) :: aux(1-mbc:maxmx+mbc,maux)

    ! Locals
    integer :: i
    double precision :: x,xmid,eigen_vector(4),gamma,lambda,alpha,h_1,h_2,deta

    if (.not.((0 <= init_type).and.(init_type <= 5))) then
        print "(a,i2)","Invalid initialization type requested, init_type = ",init_type
    endif
    
    do i=1,mx
        x = xlower+(i-0.5)*dx

        ! Set initial perturbation to zero and constant depths
        q(i,1) = aux(i,3) * rho(1)
        q(i,3) = aux(i,4) * rho(2)
        q(i,2) = 0.d0
        q(i,4) = 0.d0
        
        ! Simple Riemann problem specification
        if (init_type == 0) then
            ! Depth is already set above appropriately
            ! Need to set momentum though
            if (x < init_location) then
                q(i,2) = u_left(1) * q(i,1)
                q(i,4) = u_left(2) * q(i,3)
            else
                q(i,2) = u_right(1) * q(i,1)
                q(i,4) = u_right(2) * q(i,3)
            endif

        ! Riemann problem in one wave family
        else if (init_type == 1) then
            ! Calculate wave family for perturbation
            gamma = aux(i,4) / aux(i,3)
            select case(wave_family)
                case(1) ! Shallow water, left-going
                    alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                    lambda = -sqrt(g*aux(i,3)*(1.d0+alpha))
                case(2) ! Internal wave, left-going
                    alpha = 0.5d0 * (gamma - 1.d0 - sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                    lambda = -sqrt(g*aux(i,3)*(1.d0+alpha))
                case(3) ! Internal wave, right-going
                    alpha = 0.5d0 * (gamma - 1.d0 - sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                    lambda = sqrt(g*aux(i,3)*(1.d0+alpha))
                case(4) ! Shallow water, right-going
                    alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                    lambda = sqrt(g*aux(i,3)*(1.d0+alpha))
            end select
            eigen_vector = [1.d0,lambda,alpha,lambda*alpha]
            
            ! Add perturbation
            if ((x < init_location).and.(wave_family >= 3)) then
                q(i,1:2) = q(i,1:2) + rho(1) * epsilon * eigen_vector(1:2)
                q(i,3:4) = q(i,3:4) + rho(2) * epsilon * eigen_vector(3:4)
            else if ((x > init_location).and.(wave_family < 3)) then
                q(i,1:2) = q(i,1:2) + rho(1) * epsilon * eigen_vector(1:2)
                q(i,3:4) = q(i,3:4) + rho(2) * epsilon * eigen_vector(3:4)
            endif
        ! Gaussian hump of water, shallow water style
        else if (init_type == 2) then
            gamma = aux(i,4) / aux(i,3)
            alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
            deta = epsilon * exp(-((x-init_location)/sigma)**2)
            q(i,1) = q(i,1) + rho(1) * deta
            q(i,3) = q(i,3) + rho(2) * alpha * deta
        ! Gaussian on internal layer only
        else if (init_type == 3) then
            deta = epsilon * exp(-((x-init_location)/sigma)**2)
            q(i,1) = q(i,1) - rho(1) * deta
            q(i,3) = q(i,3) + rho(2) * deta
        ! Shelf initial condition from AN paper
        else if (init_type == 4) then
            gamma = aux(i,4) / aux(i,3)
!             alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
            alpha = 0.d0
            xmid = 0.5d0*(-180.e3-80.e3)
            if ((x > -130.e3).and.(x < -80.e3)) then
                deta = epsilon * sin((x-xmid)*PI/(-80.e3-xmid))
                q(i,3) = q(i,3) + rho(2) * alpha * deta
                q(i,1) = q(i,1) + rho(1) * deta * (1.d0 - alpha)
            endif
        endif
    enddo
    
end subroutine qinit