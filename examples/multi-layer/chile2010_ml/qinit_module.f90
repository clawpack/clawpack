module qinit_module

    implicit none
    save

    logical, private :: module_setup = .false.
    
    ! Type of q initialization
    integer, public :: qinit_type

    integer, public :: min_level_qinit
    integer, public :: max_level_qinit

    ! Geometry
    real(kind=8), public :: x_low_qinit
    real(kind=8), public :: y_low_qinit
    real(kind=8), public :: t_low_qinit
    real(kind=8), public :: x_hi_qinit
    real(kind=8), public :: y_hi_qinit
    real(kind=8), public :: t_hi_qinit
    real(kind=8), public :: dx_qinit
    real(kind=8), public :: dy_qinit
    
    ! Work array
    real(kind=8), private, allocatable :: qinit(:)

    integer, private :: mx_qinit
    integer, private :: my_qinit

    ! Specifc types of intialization    
    ! Type of perturbation to add
    integer, private :: wave_family
    real(kind=8), private :: init_location(2), epsilon
    real(kind=8), private :: angle, sigma

contains

    subroutine set_qinit(fname)
    
        use geoclaw_module, only: GEO_PARM_UNIT
    
        implicit none
        
        ! Subroutine arguments
        character(len=*), optional, intent(in) :: fname
        
        ! File handling
        integer, parameter :: unit = 7
        character(len=150) :: qinit_fname
        
        if (.not.module_setup) then

            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '--------------------------------------------'
            write(GEO_PARM_UNIT,*) 'SETQINIT:'
            write(GEO_PARM_UNIT,*) '-------------'
            
            ! Open the data file
            if (present(fname)) then
                call opendatafile(unit,fname)
            else
                call opendatafile(unit,"qinit.data")
            endif
            
            read(unit,"(i1)") qinit_type
            if (qinit_type == 0) then
                ! No perturbation specified
                write(GEO_PARM_UNIT,*)  '  qinit_type = 0, no perturbation'
                print *,'  qinit_type = 0, no perturbation'
                return
            else if (qinit_type > 0 .and. qinit_type < 5) then
                read(unit,*) qinit_fname
                read(unit,"(2i2)") min_level_qinit, max_level_qinit

                write(GEO_PARM_UNIT,*) '   min_level, max_level, qinit_fname:'
                write(GEO_PARM_UNIT,*)  min_level_qinit, max_level_qinit, qinit_fname
                
                call read_qinit(qinit_fname)
            else if (qinit_type >= 5) then
                read(unit,*) epsilon
                read(unit,*) init_location
                read(unit,*) wave_family
                read(unit,*) angle
                read(unit,*) sigma

                write(GEO_PARM_UNIT,*) " epsilon = ",  epsilon
                write(GEO_PARM_UNIT,*) " init_location = ",  init_location
                write(GEO_PARM_UNIT,*) " wave_family = ",  wave_family
                write(GEO_PARM_UNIT,*) " angle = ",  angle
                write(GEO_PARM_UNIT,*) " sigma = ",  sigma
            endif

            close(unit)

            module_setup = .true.
        end if

    end subroutine set_qinit



    subroutine add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
        use geoclaw_module, only: sea_level, pi, g => grav
        use multilayer_module, only: aux_layer_index, rho, r
    
        implicit none
    
        ! Subroutine arguments
        integer, intent(in) :: meqn,mbc,mx,my,maux
        real(kind=8), intent(in) :: xlower,ylower,dx,dy
        real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
        
        ! Local
        integer :: i,j
        real(kind=8) :: ximc,xim,x,xc,xip,xipc,yjmc,yjm,y,yc,yjp,yjpc,dq

        real(kind=8) :: xmid,m,x_c,y_c
        real(kind=8) :: eigen_vector(6),gamma,lambda,alpha,h_1,h_2,deta
        
        ! Topography integral function
        real(kind=8) :: topointegral
        
        if (qinit_type > 0 .and. qinit_type < 5) then
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0)*dx
                xim = x - 0.5d0*dx
                xip = x + 0.5d0*dx
                do j=1-mbc,my+mbc
                    y = ylower + (j-0.5d0)*dy
                    yjm = y - 0.5d0*dy
                    yjp = y + 0.5d0*dy

                    ! Check to see if we are in the qinit region at this grid point
                    if ((xip > x_low_qinit).and.(xim < x_hi_qinit).and.  &
                        (yjp > y_low_qinit).and.(yjm < y_hi_qinit)) then

                        xipc=min(xip,x_hi_qinit)
                        ximc=max(xim,x_low_qinit)
                        xc=0.5d0*(xipc+ximc)

                        yjpc=min(yjp,y_hi_qinit)
                        yjmc=max(yjm,y_low_qinit)
                        yc=0.5d0*(yjmc+yjpc)

                        dq = topointegral(ximc,xc,xipc,yjmc,yc,yjpc,x_low_qinit, &
                                          y_low_qinit,dx_qinit,dy_qinit,mx_qinit, &
                                          my_qinit,qinit,1)
                        dq = dq / ((xipc-ximc)*(yjpc-yjmc)*aux(2,i,j))

                        if (qinit_type < 4) then 
                            if (aux(1,i,j) <= sea_level) then
                                q(qinit_type,i,j) = q(qinit_type,i,j) + dq
                            endif
                        else if (qinit_type == 4) then
                            q(1,i,j) = max(dq-aux(1,i,j),0.d0)
                        endif
                    endif
                enddo
            enddo
        else if (qinit_type > 4) then
            do i=1,mx
                x = xlower + (i - 0.5d0) * dx
                do j=1,my
                    y = ylower + (j - 0.5d0) * dy
                    
                    ! Test perturbations - these only work in the x-direction
                    if (qinit_type == 5 .or. qinit_type == 6) then
                        ! Calculate wave family for perturbation
                        gamma = aux(aux_layer_index+1,i,j) / aux(aux_layer_index,i,j)
                        select case(wave_family)
                            case(1) ! Shallow water, left-going
                                alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                                lambda = -sqrt(g*aux(aux_layer_index,i,j)*(1.d0+alpha))
                            case(2) ! Internal wave, left-going
                                alpha = 0.5d0 * (gamma - 1.d0 - sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                                lambda = -sqrt(g*aux(aux_layer_index,i,j)*(1.d0+alpha))
                            case(3) ! Internal wave, right-going
                                alpha = 0.5d0 * (gamma - 1.d0 - sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                                lambda = sqrt(g*aux(aux_layer_index,i,j)*(1.d0+alpha))
                            case(4) ! Shallow water, right-going
                                alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                                lambda = sqrt(g*aux(aux_layer_index,i,j)*(1.d0+alpha))
                        end select
                        eigen_vector = [1.d0,lambda,0.d0,alpha,lambda*alpha,0.d0]

                        if (qinit_type == 5) then
                            ! Add perturbation
                            if ((x < init_location(1)).and.(wave_family >= 3)) then
                                q(1:3,i,j) = q(1:3,i,j) + rho(1) * epsilon * eigen_vector(1:3)
                                q(4:6,i,j) = q(4:6,i,j) + rho(2) * epsilon * eigen_vector(4:6)
                            else if ((x > init_location(1)).and.(wave_family < 3)) then
                                q(1:2,i,j) = q(1:2,i,j) + rho(1) * epsilon * eigen_vector(1:2)
                                q(4:5,i,j) = q(4:5,i,j) + rho(2) * epsilon * eigen_vector(4:5)
                            endif
                        ! Gaussian wave along a direction on requested wave family
                        else if (qinit_type == 6) then
                            ! Transform back to computational coordinates
                            x_c = x * cos(angle) + y * sin(angle) - init_location(1)
                            deta = epsilon * exp(-(x_c/sigma)**2)
                            q(1,i,j) = q(1,i,j) + rho(1) * deta
                            q(4,i,j) = q(4,i,j) + rho(2) * alpha * deta
                        endif
                    ! Symmetric gaussian hump
                    else if (qinit_type == 7) then
                        deta = epsilon * exp(-((x-init_location(1))/sigma)**2)  &
                                       * exp(-((y-init_location(2))/sigma)**2)
                        q(1,i,j) = q(1,i,j) + rho(1) * deta
                    ! Shelf conditions from AN paper
                    else if (qinit_type == 8) then
                        alpha = 0.d0
                        xmid = 0.5d0*(-180.e3 - 80.e3)
                        if ((x > -130.e3).and.(x < -80.e3)) then
                            deta = epsilon * sin((x-xmid)*pi/(-80.e3-xmid))
                            q(4,i,j) = q(4,i,j) + rho(2) * alpha * deta
                            q(1,i,j) = q(1,i,j) + rho(1) * deta * (1.d0 - alpha)
                        endif
                    ! Inundation test
                    else if (qinit_type == 9) then
                        x_c = (x - init_location(1)) * cos(angle) &
                            + (y - init_location(2)) * sin(angle)
                        deta = epsilon * exp(-(x_c/sigma)**2)
                        q(1,i,j) = q(1,i,j) + rho(1) * deta
                    endif
                enddo
            enddo

        endif
        
    end subroutine add_perturbation

        
    ! currently only supports one file type:
    ! x,y,z values, one per line in standard order from NW corner to SE
    ! z is perturbation from standard depth h,hu,hv set in qinit_geo,
    ! if iqinit = 1,2, or 3 respectively.
    ! if iqinit = 4, the z column corresponds to the definition of the 
    ! surface elevation eta. The depth is then set as q(i,j,1)=max(eta-b,0)
    subroutine read_qinit(fname)
    
        use geoclaw_module, only: GEO_PARM_UNIT
        
        implicit none
        
        ! Subroutine arguments
        character(len=150) :: fname
        
        ! Data file opening
        integer, parameter :: unit = 19
        integer :: i,num_points,status
        double precision :: x,y
        
        print *,'  '
        print *,'Reading qinit data from file  ', fname
        print *,'  '

        write(GEO_PARM_UNIT,*) '  '
        write(GEO_PARM_UNIT,*) 'Reading qinit data from'
        write(GEO_PARM_UNIT,*) fname
        write(GEO_PARM_UNIT,*) '  '
        
        open(unit=unit, file=fname, iostat=status, status="unknown", &
             form='formatted',action="read")
        if ( status /= 0 ) then
            print *,"Error opening file", fname
            stop
        endif
        
        ! Initialize counters
        num_points = 0
        mx_qinit = 0
        
        ! Read in first values, determines x_low and y_hi
        read(unit,*) x_low_qinit,y_hi_qinit
        num_points = num_points + 1
        mx_qinit = mx_qinit + 1
        
        ! Sweep through first row figuring out mx
        y = y_hi_qinit
        do while (y_hi_qinit == y)
            read(unit,*) x,y
            num_points = num_points + 1
            mx_qinit = mx_qinit + 1
        enddo
        ! We over count by one in the above loop
        mx_qinit = mx_qinit - 1
        
        ! Continue to count the rest of the lines
        do
            read(unit,*,iostat=status) x,y
            if (status /= 0) exit
            num_points = num_points + 1
        enddo
        if (status > 0) then
            print *,"ERROR:  Error reading qinit file ",fname
            stop
        endif
        
        ! Extract rest of geometry
        x_hi_qinit = x
        y_low_qinit = y
        my_qinit = num_points / mx_qinit
        dx_qinit = (x_hi_qinit - x_low_qinit) / (mx_qinit-1)
        dy_qinit = (y_hi_qinit - y_low_qinit) / (my_qinit-1)
        
        rewind(unit)
        allocate(qinit(num_points))
        
        ! Read and store the data this time
        do i=1,num_points
            read(unit,*) x,y,qinit(i)
        enddo
        close(unit)
        
    end subroutine read_qinit

end module qinit_module