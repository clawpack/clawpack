module qinit_module

    implicit none
    save
    
    ! Type of q initialization
    integer, public :: qinit_type
    
    ! Type of perturbation to add
    integer, private :: wave_family
    real(kind=8), private :: init_location(2),epsilon
    real(kind=8), private :: angle,sigma
    
contains

    subroutine add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
        use geoclaw_module, only: num_layers,eta_init,rho,pi
        use geoclaw_module, only: g => grav
    
        implicit none
    
        ! Subroutine arguments
        integer, intent(in) :: meqn,mbc,mx,my,maux
        real(kind=8), intent(in) :: xlower,ylower,dx,dy
        real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
        
        ! Locals
        integer :: i,j
        real(kind=8) :: x,y,xmid,m,x_c,y_c
        real(kind=8) :: eigen_vector(6),gamma,lambda,alpha,h_1,h_2,deta,r
        
        if (num_layers > 2) then
            print *,"Adding perturbations only supported for num_layers == 2."
            stop
        endif
        
        r = rho(0) / rho(1)
        
        do i=1,mx
            x = xlower + (i - 0.5d0) * dx
            do j=1,my
                y = ylower + (j - 0.5d0) * dy
                
                ! Test perturbations - these only work in the x-direction
                if (qinit_type == 1 .or. qinit_type == 2) then
                    ! Calculate wave family for perturbation
                    gamma = aux(8,i,j) / aux(7,i,j)
                    select case(wave_family)
                        case(1) ! Shallow water, left-going
                            alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                            lambda = -sqrt(g*aux(7,i,j)*(1.d0+alpha))
                        case(2) ! Internal wave, left-going
                            alpha = 0.5d0 * (gamma - 1.d0 - sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                            lambda = -sqrt(g*aux(7,i,j)*(1.d0+alpha))
                        case(3) ! Internal wave, right-going
                            alpha = 0.5d0 * (gamma - 1.d0 - sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                            lambda = sqrt(g*aux(7,i,j)*(1.d0+alpha))
                        case(4) ! Shallow water, right-going
                            alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                            lambda = sqrt(g*aux(7,i,j)*(1.d0+alpha))
                    end select
                    eigen_vector = [1.d0,lambda,0.d0,alpha,lambda*alpha,0.d0]

                    if (qinit_type == 1) then
                        ! Add perturbation
                        if ((x < init_location(1)).and.(wave_family >= 3)) then
                            q(1:3,i,j) = q(1:3,i,j) + rho(1) * epsilon * eigen_vector(1:3)
                            q(4:6,i,j) = q(4:6,i,j) + rho(2) * epsilon * eigen_vector(4:6)
                        else if ((x > init_location(1)).and.(wave_family < 3)) then
                            q(1:2,i,j) = q(1:2,i,j) + rho(1) * epsilon * eigen_vector(1:2)
                            q(4:5,i,j) = q(4:5,i,j) + rho(2) * epsilon * eigen_vector(4:5)
                        endif
                    ! Gaussian wave along a direction on requested wave family
                    else if (qinit_type == 2) then
                        ! Transform back to computational coordinates
                        x_c = x * cos(angle) + y * sin(angle) - init_location(1)
                        deta = epsilon * exp(-(x_c/sigma)**2)
                        q(1,i,j) = q(1,i,j) + rho(1) * deta
                        q(4,i,j) = q(4,i,j) + rho(2) * alpha * deta
                    endif
                ! Symmetric gaussian hump
                else if (qinit_type == 3) then
                    deta = epsilon * exp(-((x-init_location(1))/sigma)**2)  &
                                   * exp(-((y-init_location(2))/sigma)**2)
                    q(1,i,j) = q(1,i,j) + rho(1) * deta
                ! Shelf conditions from AN paper
                else if (qinit_type == 4) then
                    alpha = 0.d0
                    xmid = 0.5d0*(-180.e3 - 80.e3)
                    if ((x > -130.e3).and.(x < -80.e3)) then
                        deta = epsilon * sin((x-xmid)*pi/(-80.e3-xmid))
                        q(4,i,j) = q(4,i,j) + rho(2) * alpha * deta
                        q(1,i,j) = q(1,i,j) + rho(1) * deta * (1.d0 - alpha)
                    endif
                ! Inundation test
                else if (qinit_type == 5) then
                    x_c = (x - init_location(1)) * cos(angle) &
                        + (y - init_location(2)) * sin(angle)
                    deta = epsilon * exp(-(x_c/sigma)**2)
                    q(1,i,j) = q(1,i,j) + rho(1) * deta
                endif
            enddo
        enddo
        
    end subroutine add_perturbation

    subroutine set_qinit(fname)
    
        use geoclaw_module, only: GEO_PARM_UNIT
    
        implicit none
        
        ! Subroutine arguments
        character(len=*), optional, intent(in) :: fname
        
        ! File handling
        integer, parameter :: unit = 7
        character(len=150) :: qinit_fname
        
        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'SETQINIT:'
        write(GEO_PARM_UNIT,*) '-------------'
        
        ! Open the data file
        if (present(fname)) then
            call opendatafile(unit,fname)
        else
            call opendatafile(unit,"setqinit.data")
        endif
        
        read(unit,"(i1)") qinit_type
        if (qinit_type == 0) then
            ! No perturbation specified
            write(GEO_PARM_UNIT,*)  '  qinit_type = 0, no perturbation'
            print *,'  qinit_type = 0, no perturbation'
            return
        endif
        
        if (qinit_type > 0) then
            read(unit,*) epsilon
            if (qinit_type <= 2 .or. qinit_type == 5) then  
                read(unit,*) init_location
                read(unit,*) wave_family
                if(qinit_type == 2 .or. qinit_type == 5) then
                    read(unit,*) angle
                    read(unit,*) sigma
                endif
            else if (qinit_type == 3) then
                read(unit,*) init_location
                read(unit,*) sigma
            endif
        endif
        
        close(unit)
        
    end subroutine set_qinit

end module qinit_module
