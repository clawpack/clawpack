! ============================================================================
!      Copyright (C) 2010-10-12 Kyle Mandli <kyle.mandli@gmail.com>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

module multilayer_module

    implicit none

    logical, private :: module_setup = .false.

    ! Storage parameters
    integer :: aux_layer_index
    
    ! Physical parameters
    integer :: num_layers
    real(kind=8), allocatable :: eta_init(:)
    real(kind=8) :: r, one_minus_r
    
    ! Algorithm parameters
    integer :: eigen_method,inundation_method
    logical :: check_richardson
    real(kind=8) :: richardson_tolerance
    real(kind=8), allocatable :: wave_tolerance(:), dry_tolerance(:)
    
    ! Output files
    integer, parameter :: KAPPA_UNIT = 42

    ! Eigenspace evaluation functions
    abstract interface
        subroutine eigen_function(h_l, h_r, u_l, u_r, v_l, v_r,                &
                                            n_index, t_index, s, eig_vec)
            implicit none
            real(kind=8), dimension(2), intent(in) :: h_l,h_r,u_l,u_r,v_l,v_r
            integer, intent(in) :: n_index,t_index
            real(kind=8), intent(inout) :: s(6),eig_vec(6,6)
        end subroutine eigen_function
    end interface

    procedure (eigen_function), pointer :: eigen_func
    procedure (eigen_function), pointer :: inundation_eigen_func
    
contains

    ! ========================================================================
    !  set_multilayer(file_name)
    ! ========================================================================
    subroutine set_multilayer(data_file)

        use geoclaw_module, only: GEO_PARM_UNIT, rho
        use geoclaw_module, only: geo_dry_tolerance => dry_tolerance
        use storm_module, only: storm_specification_type

        implicit none
        character(len=*), optional, intent(in) :: data_file
        
        integer :: ios
        integer, parameter :: IOUNIT = 41
            

        if (.not.module_setup) then
            ! Open file
            if (present(data_file)) then
                call opendatafile(IOUNIT, data_file)
            else
                call opendatafile(IOUNIT, 'multilayer.data')
            endif
            
            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '--------------------------------------------'
            write(GEO_PARM_UNIT,*) 'Multilayer Parameters:'
            write(GEO_PARM_UNIT,*) '----------------------'

            ! Read in parameters
            read(IOUNIT,"(i3)") num_layers
            allocate(eta_init(num_layers))
            allocate(wave_tolerance(num_layers))
            allocate(dry_tolerance(num_layers))

            ! read(IOUNIT,*) rho
            ! The densities are in the geoclaw_module, check here to make surge
            ! there are enough values
            if (size(rho, 1) /= num_layers) then
                print *, "The number of values of the water density " //       &
                          "(", size(rho, 1),") does not match the number " //  &
                          "of layers (", num_layers,")."
                stop
            end if
            if (num_layers > 1) then
                r = rho(1) / rho(2)
                one_minus_r = 1.d0 - r
            else
                r = -1.d0
                one_minus_r = 0.d0
            endif
            read(IOUNIT, *) eta_init
            read(IOUNIT, *) wave_tolerance
            read(IOUNIT, *) aux_layer_index
            read(IOUNIT, *)

            ! Algorithmic parameters
            read(IOUNIT, *) check_richardson
            read(IOUNIT, "(d16.8)") richardson_tolerance
            read(IOUNIT, "(i1)") eigen_method
            read(IOUNIT, "(i1)") inundation_method
            
            close(IOUNIT) 

            ! Currently just set dry_tolerance(:) = dry_tolerance
            dry_tolerance = geo_dry_tolerance

            ! Set eigen functions
            select case(eigen_method)
                case(1:2)
                    eigen_func => linearized_eigen
                case(3)
                    eigen_func => vel_diff_eigen
                case(4)
                    eigen_func => lapack_eigen
                case default
                    print *, "Invalid eigenspace method requested: ", eigen_method
                    stop
            end select
            select case(inundation_method)
                case(1:2)
                    inundation_eigen_func => linearized_eigen
                case(3)
                    inundation_eigen_func => vel_diff_eigen
                case(4:5)
                    inundation_eigen_func => lapack_eigen
                case default
                    print *, "Invalid eigenspace method requested: ", eigen_method
                    stop
            end select

            ! Open Kappa output file if num_layers > 1
            ! Open file for writing hyperbolicity warnings if multiple layers
            if (num_layers > 1 .and. check_richardson) then
                open(unit=KAPPA_UNIT, file='fort.kappa', iostat=ios, &
                        status="unknown", action="write")
                if ( ios /= 0 ) stop "Error opening file name fort.kappa"
            endif

            write(GEO_PARM_UNIT,*) '   check_richardson:',check_richardson
            write(GEO_PARM_UNIT,*) '   richardson_tolerance:',richardson_tolerance
            write(GEO_PARM_UNIT,*) '   eigen_method:',eigen_method
            write(GEO_PARM_UNIT,*) '   inundation_method:',inundation_method
            write(GEO_PARM_UNIT,*) '   dry_tolerance:',dry_tolerance

            module_setup = .true.
        end if
        
    end subroutine set_multilayer

    ! ==========================================================================
    !  Eigenspace Evaluation Functions, see interface defined above
    ! ==========================================================================
    subroutine linearized_eigen(h_l,h_r,u_l,u_r,v_l,v_r,n_index,t_index,s,eig_vec)

        use geoclaw_module, only: grav

        implicit none
        
        ! Input
        double precision, dimension(2), intent(in) :: h_l,h_r,u_l,u_r,v_l,v_r
        integer, intent(in) :: n_index,t_index
        
        ! Output
        double precision, intent(inout) :: s(6),eig_vec(6,6)
            
        ! Local
        double precision :: gamma_l,gamma_r,alpha(4),g
        
        g = grav
            
        ! Calculate relevant quantities
        gamma_l = h_l(2) / h_l(1)
        gamma_r = h_r(2) / h_r(1)

        alpha(1) = 0.5d0*(gamma_l-1.d0+sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
        alpha(2) = 0.5d0*(gamma_l-1.d0-sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
        alpha(3) = 0.5d0*(gamma_r-1.d0-sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))
        alpha(4) = 0.5d0*(gamma_r-1.d0+sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))

        s(1) = -sqrt(g*h_l(1)*(1+alpha(1)))
        s(2) = -sqrt(g*h_l(1)*(1+alpha(2)))
        s(3:4) = 0.5d0 * (u_l(:) + u_r(:))
        s(5) = sqrt(g*h_r(1)*(1+alpha(3)))
        s(6) = sqrt(g*h_r(1)*(1+alpha(4)))

        ! Compute eigenspace exactly based on eigenvalues provided
        eig_vec(1,:) = [1.d0,1.d0,0.d0,0.d0,1.d0,1.d0]
        
        eig_vec(n_index,:) = [s(1),s(2),0.d0,0.d0,s(5),s(6)]
        
        eig_vec(t_index,:) = [v_l(1),v_l(1),1.d0,0.d0,v_r(1),v_r(1)]

        eig_vec(4,1:2) = alpha(1:2)
        eig_vec(4,3:4) = 0.d0
        eig_vec(4,5:6) = alpha(3:4)
        
        eig_vec(n_index+3,:) = s * eig_vec(4,:)
        
        eig_vec(t_index+3,1:2) = v_l(2) * alpha(1:2)
        eig_vec(t_index+3,3:4) = [0.d0,1.d0]
        eig_vec(t_index+3,5:6) = v_r(2) * alpha(3:4)

    end subroutine linearized_eigen


    subroutine vel_diff_eigen(h_l,h_r,u_l,u_r,v_l,v_r,n_index,t_index,s,eig_vec)

        use geoclaw_module, only: grav

        implicit none
        
        ! Input
        double precision, dimension(2), intent(in) :: h_l,h_r,u_l,u_r,v_l,v_r
        integer, intent(in) :: n_index,t_index
        
        ! Output
        double precision, intent(inout) :: s(6),eig_vec(6,6)
            
        ! Local
        double precision :: total_depth_l,total_depth_r,mult_depth_l,mult_depth_r
        
        total_depth_l = sum(h_l)
        total_depth_r = sum(h_r)
        mult_depth_l = product(h_l)
        mult_depth_r = product(h_r)
                          
        s(1) = - sqrt(grav*total_depth_l) + 0.5d0 * mult_depth_l/ total_depth_l**(3/2) * one_minus_r
        s(2) = - sqrt(grav * mult_depth_l / total_depth_l * one_minus_r)
        s(3:4) = 0.5d0 * (u_l + u_r)
        s(5) = sqrt(grav * mult_depth_r / total_depth_r * one_minus_r)
        s(6) = sqrt(grav*total_depth_r) - 0.5d0 * mult_depth_r / total_depth_r**(3/2) * one_minus_r

        eig_vec(1,:) = [1.d0,1.d0,0.d0,0.d0,1.d0,1.d0]
        eig_vec(n_index,:) = [s(1),s(2),0.d0,0.d0,s(5),s(6)]
        eig_vec(t_index,:) = [v_l(1),v_l(1),1.d0,0.d0,v_r(1),v_r(1)]
        eig_vec(4,1:2) = ((s(1:2)-u_l(1))**2 - grav*h_l(1)) / (r*grav*h_l(1))
        eig_vec(4,3:4) = 0.d0
        eig_vec(4,5:6) = ((s(5:6)-u_r(1))**2 - grav*h_r(1)) / (r*grav*h_r(1))
        eig_vec(n_index+3,:) = s(:) * eig_vec(4,:)
        eig_vec(t_index+3,1:2) = v_l(2) * eig_vec(4,1:2)
        eig_vec(t_index+3,3:4) = [0.d0,1.d0]
        eig_vec(t_index+3,5:6) = v_r(2) * eig_vec(4,5:6)
        
    end subroutine vel_diff_eigen

#ifdef LAPACK_AVAIL
    subroutine lapack_eigen(h_l,h_r,u_l,u_r,v_l,v_r,n_index,t_index,s,eig_vec)

        use geoclaw_module, only: grav

        implicit none
        
        ! Input
        double precision, dimension(2), intent(in) :: h_l,h_r,u_l,u_r,v_l,v_r
        integer, intent(in) :: n_index,t_index
        
        ! Output
        double precision, intent(inout) :: s(6),eig_vec(6,6)
        
        ! Local
        integer, parameter :: lwork = 6*6
        integer :: i,j,m,info
        double precision :: h_ave(2),u_ave(2),v_ave(2),A(6,6),A_copy(6,6)
        double precision :: imag_evalues(6),empty,work(1,lwork)
        double precision :: g
        
        g = grav

        ! Construct flux matrix
        A = 0.d0
        h_ave(:) = 0.5d0 * (h_l(:) + h_r(:))
        u_ave(:) = 0.5d0 * (u_l(:) + u_r(:))
        v_ave(:) = 0.5d0 * (v_l(:) + v_r(:))
        ! We need to do this since one of the rows cannot be swapped
        if (n_index == 2) then
            A(t_index,1:3) = [-u_ave(1)*v_ave(1),v_ave(1),u_ave(1)]
            A(t_index+3,4:6) = [-u_ave(2)*v_ave(2),v_ave(2),u_ave(2)]
        else
            A(t_index,1:3) = [-u_ave(1)*v_ave(1),u_ave(1),v_ave(1)]
            A(t_index+3,4:6) = [-u_ave(2)*v_ave(2),u_ave(2),v_ave(2)]
        endif
            
        A(1,n_index) = 1.d0
        
        A(n_index,1) = -u_ave(1)**2 + g*h_ave(1)
        A(n_index,n_index) = 2*u_ave(1)
        A(n_index,4) = r*g*h_ave(1)
        
        A(4,n_index+3) = 1.d0
        
        A(n_index+3,1) = g*h_ave(2)
        A(n_index+3,4) = -u_ave(2)**2 + g*h_ave(2)
        A(n_index+3,n_index+3) = 2.d0 * u_ave(2)
        
        A_copy = A
        
        ! Call LAPACK
        call dgeev('N','V',6,A,6,s,imag_evalues,empty,1,eig_vec,6,work,lwork,info)
        if (info < 0) then
            info = -info
            print "(a,i1,a)","The ",info,"th argument had an illegal value."
            stop
        else if (info > 0) then
            print "(a)","The QR algorithm failed to compute all the"
            print "(a)","eigenvalues, and no eigenvectors have been"
            print "(a,i1,a)","computed; elements",i,"+1:4 of WR and WI"
            print "(a)","contain eigenvalues which have converged."
            stop
        endif
        do i=1,6
            if (eig_vec(1,i) /= 0.d0) then
                eig_vec(:,i) = eig_vec(:,i) / eig_vec(1,i)
            endif
            if (abs(imag_evalues(i)) > 0.d0) then
                print "(a,i1,a,d16.8)","Imaginary eigenvalue(",i,") > 0.0",imag_evalues(i)
                stop
            endif
        enddo

    end subroutine lapack_eigen
#else

    subroutine lapack_eigen(h_l,h_r,u_l,u_r,v_l,v_r,n_index,t_index,s,eig_vec)

        use geoclaw_module, only: grav

        implicit none
        
        ! Input
        double precision, dimension(2), intent(in) :: h_l,h_r,u_l,u_r,v_l,v_r
        integer, intent(in) :: n_index,t_index
        
        ! Output
        double precision, intent(inout) :: s(6),eig_vec(6,6)
  
        stop "LAPACK eigensolver was not available at build time."

    end subroutine lapack_eigen

#endif

    ! ==========================================================================
    !  Single layer eigensolver
    !   Note that this routine puts the result in the upper 3x3 matrix, the rest
    !   is left as zeros.
    ! ==========================================================================
    subroutine sl_eigen(h_l,h_r,u_l,u_r,v_l,v_r,n_index,t_index,s,eig_vec)

        use geoclaw_module, only: grav

        implicit none
        
        ! Input
        double precision, dimension(2), intent(in) :: h_l,h_r,u_l,u_r,v_l,v_r
        integer, intent(in) :: n_index,t_index
        
        ! Output
        double precision, intent(inout) :: s(6),eig_vec(6,6)

        s = 0.d0
        s(1) = u_l(1) - sqrt(grav*h_l(1))
        s(2) = 0.5d0 * (u_r(1) + u_l(1))
        s(3) = u_r(1) + sqrt(grav*h_r(1))
        
        eig_vec = 0.d0
        eig_vec(1,1:3) = [1.d0,0.d0,1.d0]
        eig_vec(n_index,1:3) = [s(1),0.d0,s(3)]
        eig_vec(t_index,1:3) = [v_l(1),1.d0,v_r(1)]
        
    end subroutine sl_eigen


end module multilayer_module