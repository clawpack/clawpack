! ============================================================================
!  Reimann Solver for the two-layer shallow water equations
! ============================================================================
subroutine rp1(maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)

    use parameters_module

    implicit none
    
    ! Input arguments
    integer, intent(in) :: maxmx,meqn,mwaves,mbc,mx
    
    double precision, intent(in), dimension(1-mbc:maxmx+mbc, meqn) :: ql,qr
    double precision, intent(in), dimension(1-mbc:maxmx+mbc, *) :: auxl,auxr
    
    ! Output arguments
    double precision, intent(out) :: s(1-mbc:maxmx+mbc, mwaves)
    double precision, intent(out) :: fwave(1-mbc:maxmx+mbc, meqn, mwaves)
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, meqn) :: amdq,apdq
    
    ! Local storage
    integer :: i,j,m,mw
    double precision :: eig_vec(4,4),delta(4),alpha(4),beta(4)
    double precision, dimension(4) :: flux_l,flux_r
    double precision, dimension(2) :: h_l,u_l,hu_l,h_r,u_r,hu_r,h_ave,u_ave
    double precision :: b_l,b_r,gamma_l,gamma_r,tau,w_l,w_r
    double precision :: wind_speed,lambda(4),eta_l(2),eta_r(2),h_hat_l(2),h_hat_r(2)
    logical :: dry_state_l(2), dry_state_r(2), inundation
    
    integer :: layer_index
    double precision :: momentum_transfer(2),flux_transfer_r,flux_transfer_l
    double precision :: inundation_height(2),lambda_l,lambda_r,lambda_hat

    integer :: trans_wave(2-mbc:mx+mbc)
    double precision :: wave_correction(2-mbc:mx+mbc)

    ! Function interfaces
    interface
        function eval_lapack_solve(eig_vec,delta) result(beta)
            implicit none
            double precision, intent(in) :: eig_vec(4,4), delta(4)
            double precision :: beta(4)
        end function eval_lapack_solve
    end interface

    ! Common block
    double precision :: dt,dx,t
    common /comxt/ dt,dx,t

    ! Initialize return variables
    amdq = 0.d0
    apdq = 0.d0

    !          |          |          
    !          |          |          
    !          |          |          
    !----------|----------|----------
    !    i-1         i         i+1
    !  qr(i-1)  ql(i)
    !     Riemann problem at i-1/2 is between qr(i-1) (left state) and ql(i) 
    !     (right state)
    do i=2-mbc,mx+mbc
        inundation = .false.
        dry_state_l = .false.
        dry_state_r = .false.
        
        ! These state variables are the actual depth and momenta
        do j=1,2
            layer_index = 2*(j-1)
            h_l(j) = qr(i-1,layer_index+1) / rho(j)
            h_r(j) = ql(i,layer_index+1) / rho(j)
            hu_l(j) = qr(i-1,layer_index+2) / rho(j)
            hu_r(j) = ql(i,layer_index+2) / rho(j)
            
            h_hat_l(j) = auxr(i-1,j+2)
            h_hat_r(j) = auxl(i,j+2)
            
            ! Check for dry states in this layer
            if (h_l(j) < dry_tolerance) then
                dry_state_l(j) = .true.
                h_l(j) = 0.d0
                u_l(j) = 0.d0
            else
                u_l(j) = qr(i-1,layer_index+2) / qr(i-1,layer_index+1)
            endif
            if (h_r(j) < dry_tolerance) then
                dry_state_r(j) = .true.
                h_r(j) = 0.d0
                u_r(j) = 0.d0
            else
                u_r(j) = ql(i,layer_index+2) / ql(i,layer_index+1)
            endif
        enddo
        
        h_ave = 0.5d0 * (h_l(:) + h_r(:))
        
        b_l = auxr(i-1,1)
        b_r = auxl(i,1)
        w_l = auxr(i-1,2)
        w_r = auxl(i,2)      
            
        ! ====================================================================
        ! Solve Single layer problem seperately
        if (dry_state_r(2).and.dry_state_l(2)) then
            call single_layer_eigen(h_l,h_r,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
            s(i,:) = lambda
            
            delta(1) = rho(1) * (hu_r(1) - hu_l(1))
            flux_r(2) = rho(1) * (h_r(1) * u_r(1)**2 + 0.5d0 * g * h_r(1)**2)
            flux_l(2) = rho(1) * (h_l(1) * u_l(1)**2 + 0.5d0 * g * h_l(1)**2)
            delta(2) = flux_r(2) - flux_l(2) + g * rho(1) * h_ave(1) * (b_r - b_l)
            
            beta(1) = (delta(1)*s(i,4)-delta(2)) / (s(i,4)-s(i,1))
            beta(2) = 0.d0
            beta(3) = 0.d0
            beta(4) = (delta(2) - s(i,1)*beta(1)) / s(i,4)
            
            ! Calculate waves
            forall(mw=1:4)
                fwave(i,:,mw) = eig_vec(:,mw) * beta(mw)
            end forall
            cycle
        endif
        
        ! ====================================================================
        ! Calculate eigen-space values
        ! ====================================================================
        
        ! ====================================================================
        ! Inundation cases
        if (dry_state_r(2).and.(.not.dry_state_l(2)).and.(h_l(2) + b_l > b_r)) then
            print *,"Right inundation problem"
            inundation = .true.
            if (inundation_method == 0) then
                stop "Inundation not allowed."
            else if (inundation_method == 1) then
                ! Linear eigensystem
                inundation_height = [h_r(1),0.d0]
                call linear_eigen(h_l,inundation_height,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                s(i,:) = lambda
                
                ! Corrective wave
                s(i,3) = u_l(2) + 2.d0 * sqrt(g*(1.d0-r)*h_l(2))
                s(i,4) = u_r(1) + sqrt(g*h_r(1))
                alpha(3) = r * g * h_l(2) / ((s(i,3) - u_l(2))**2 - g * h_l(2))
                alpha(4) = 0.d0

                eig_vec(1,3:4) = 1.d0
                eig_vec(2,3:4) = s(i,3:4)
                eig_vec(3,3:4) = alpha(3:4)
                eig_vec(4,3:4) = s(i,3:4)*alpha(3:4)
            else if (inundation_method == 2) then
                inundation_height = [h_r(1),dry_tolerance]
                h_r(2) = dry_tolerance
                call linear_eigen(h_l,inundation_height,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                s(i,:) = lambda
                
!                 s(i,4) = u_r(1) + sqrt(g*h_r(1))
!                 eig_vec(:,4) = [1.d0,s(i,4),0.d0,0.d0]
            else if (inundation_method == 3) then
                inundation_height = [h_r(1),dry_tolerance]
                call velocity_eigen(h_l,inundation_height,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda
                
                ! Correction for fast wave
                s(i,4) = u_r(1) + sqrt(g*h_r(1))
                eig_vec(:,4) = [1.d0,s(i,4),0.d0,0.d0]
            else if (inundation_method == 4) then
                inundation_height = [h_r(1),dry_tolerance]
                call lapack_eigen(h_l,inundation_height,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                s(i,:) = lambda
                
                ! Correction for the fast waves
                s(i,2) = u_r(1) + sqrt(g*h_r(1))
                eig_vec(:,2) = [1.d0,s(i,2),0.d0,0.d0]
            else if (inundation_method == 5) then
                inundation_height = [h_r(1),0.d0]
                call lapack_eigen(h_l,inundation_height,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                s(i,:) = lambda                
            endif   
        else if (dry_state_l(2).and.(.not.dry_state_r(2)).and.(h_r(2) + b_r > b_l)) then
            print *,"Left inundation problem"
            inundation = .true.
            ! Inundation problem eigen
            if (inundation_method == 0) then
                stop "Inundation not allowed."
            else if (inundation_method == 1) then
                ! Linear eigensystem
                inundation_height = [h_l(1),dry_tolerance]
                call linear_eigen(inundation_height,h_r,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                s(i,:) = lambda

                ! Corrections to internal wave
                s(i,2) = u_r(2) - 2.d0 * sqrt(g*(1.d0-r)*h_r(2))
                alpha(2) = r * g * h_r(2) / ((s(i,2) - u_r(2))**2 - g * h_r(2))
                eig_vec(:,2) = [1.d0,s(i,2),alpha(2),alpha(2)*s(i,2)]
                
                ! Correction for the fast waves
                s(i,1) = u_l(1) - sqrt(g*h_l(1))
                eig_vec(:,1) = [1.d0,s(i,1),0.d0,0.d0]
            else if (inundation_method == 2) then
                ! Use linearized eigensystem
                inundation_height = [h_l(1),dry_tolerance]
                call linear_eigen(inundation_height,h_r,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                s(i,:) = lambda
                            
                ! Correction for the fast waves
                s(i,1) = u_l(1) - sqrt(g*h_l(1))
                eig_vec(:,1) = [1.d0,s(i,1),0.d0,0.d0]
            else if (inundation_method == 3) then
                ! Use velocity difference expansion eigensystems
                inundation_height = [h_l(1),dry_tolerance]
                call velocity_eigen(inundation_height,h_r,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                s(i,:) = lambda                
        
                ! Correction for the fast waves
                s(i,1) = u_r(1) - sqrt(g*h_r(1))
                eig_vec(:,1) = [1.d0,s(i,1),0.d0,0.d0]
            else if (inundation_method == 4) then
                ! LAPACK solver with corrective wave and small wet layer
                inundation_height = [h_r(1),dry_tolerance]
                call lapack_eigen(h_l,inundation_height,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                s(i,:) = lambda
                            
                ! Correction for the fast waves
                s(i,1) = u_l(1) - sqrt(g*h_l(1))
                eig_vec(:,1) = [1.d0,s(i,1),0.d0,0.d0]
            else if (inundation_method == 5) then
                ! Use the LAPACK solver with no correction
                inundation_height = [h_l(1),0.d0]
                call lapack_eigen(inundation_height,h_r,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                s(i,:) = lambda
            endif        
          
        ! ====================================================================  
        !  Wall or wet case       
        else
            ! Wall dry state or completely wet case
            if (eigen_method == 1) then
                call linear_eigen(h_hat_l,h_hat_r,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                s(i,:) = lambda
            else if (eigen_method == 2) then
                call linear_eigen(h_l,h_r,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                s(i,:) = lambda
            else if (eigen_method == 3) then 
                call velocity_eigen(h_l,h_r,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                s(i,:) = lambda
            else if (eigen_method == 4) then
                if (dry_state_r(2).and.(.not.dry_state_l(2)).or. &
                        dry_state_l(2).and.(.not.dry_state_r(2))) then
                    call linear_eigen(h_l,h_r,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                    s(i,:) = lambda
                else
                    call lapack_eigen(h_l,h_r,u_l,u_r,b_l,b_r, &
                              &     trans_wave(i),wave_correction(i),lambda,eig_vec)
                    s(i,:) = lambda
                endif
            else
                stop "Invalid eigensystem method requested, method = (1,4)."
            endif
        endif
        
        !  end of eigenspace calculation
        ! ====================================================================
        
        ! ====================================================================
        ! Calculate flux vector to be projected onto e-space
        ! Calculate jumps in fluxes
        if (dry_state_r(2).and.(.not.dry_state_l(2)).and.(.not.inundation)) then
            ! Wall boundary conditions
            h_r(2) = h_l(2)
            hu_r(2) = -hu_l(2)
            u_r(2) = -u_l(2)
            
            ! Top layer eta(2) = b_r - h_l(2) - b_l
            momentum_transfer(1) = g * rho(1) * h_ave(1) * (b_r - h_l(2) - b_l)
            
            momentum_transfer(2) = 0.d0
            flux_transfer_r = 0.d0
            flux_transfer_l = 0.d0
        ! Left state dry, right wet
        else if (dry_state_l(2).and.(.not.dry_state_r(2)).and.(.not.inundation)) then
            ! Wall boundary conditions
            h_l(2) = h_r(2)
            hu_l(2) = -hu_r(2)
            u_l(2) = -u_r(2)
            
            ! Top layer eta(2) = h_r(2) + b_r - b_l
            momentum_transfer(1) = g * rho(1) * h_ave(1) * (b_r + h_r(2) - b_l)
            
            momentum_transfer(2) = 0.d0
            flux_transfer_r = 0.d0
            flux_transfer_l = 0.d0
        ! Fully wet bottom layer or inundation
        else
            momentum_transfer(1) =   g * rho(1) * h_ave(1) * (h_r(2) - h_l(2) + b_r - b_l)
            momentum_transfer(2) = - g * rho(1) * h_ave(1) * (h_r(2) - h_l(2)) + g * rho(2) * h_ave(2) * (b_r - b_l)
            ! Bottom layer momentum transfer flux
            flux_transfer_r = g * rho(1) * product(h_r)
            flux_transfer_l = g * rho(1) * product(h_l)
        endif
        
        ! Flux jumps
        do j=1,2
            layer_index = 2*(j-1)
            flux_r(layer_index+1) = rho(j) * hu_r(j)
            flux_r(layer_index+2) = rho(j) * (h_r(j) * u_r(j)**2 + 0.5d0 * g * h_r(j)**2)
            
            flux_l(layer_index+1) = rho(j) * hu_l(j)
            flux_l(layer_index+2) = rho(j) * (h_l(j) * u_l(j)**2 + 0.5d0 * g * h_l(j)**2)
        enddo
        
        delta = flux_r - flux_l

        ! Bottom layer additional flux
        delta(4) = delta(4) + flux_transfer_r - flux_transfer_l

        ! Momentum source term from layers
        delta(2) = delta(2) + momentum_transfer(1)
        delta(4) = delta(4) + momentum_transfer(2)
        
        ! Wind forcing
        wind_speed = 0.5d0 * (w_l + w_r)
        tau = wind_drag(wind_speed) * rho_air * wind_speed
        delta(2) = delta(2) - tau * wind_speed

        ! ====================================================================
        ! Solve system, solution is stored in delta
        beta = eval_lapack_solve(eig_vec,delta)
        
        ! Calculate waves
        forall(mw=1:4)
            fwave(i,:,mw) = eig_vec(:,mw) * beta(mw)
        end forall
    enddo
    
    ! Calculate amdq and apdq
    ! No entropy fix requested
    if (.not.entropy_fix) then
        do i=2-mbc,mx+mbc
            do mw=1,mwaves
                if (s(i,mw) > 0.d0) then
                    apdq(i,:) = apdq(i,:) + fwave(i,:,mw)
                else                                     
                    amdq(i,:) = amdq(i,:) + fwave(i,:,mw)
                endif
            enddo
        enddo
    else
        do i=2-mbc,mx+mbc
            ! Check to see if an entropy fix is necessary
            if (trans_wave(i) /= 0) then
                print *,"Entropy fix needed: i=",i
                print *,"  s(i,:)=",(s(i,m),m=1,4)
                print *,"  trans_wave=",trans_wave(i)
                print *,"  beta =",wave_correction(i)
                do mw=1,trans_wave(i)-1
                    amdq(i,:) = amdq(i,:) + fwave(i,:,mw)
                enddo
                amdq(i,:) = amdq(i,:) + wave_correction(i) * fwave(i,:,trans_wave(i))
                apdq(i,:) = apdq(i,:) + (1.d0 - wave_correction(i)) * fwave(i,:,trans_wave(i))
                do mw=trans_wave(i)+1,mwaves
                    apdq(i,:) = apdq(i,:) + fwave(i,:,mw)
                enddo
            else
                ! No entropy fix needed                
                do mw=1,mwaves
                    if (s(i,mw) > 0.d0) then
                        apdq(i,:) = apdq(i,:) + fwave(i,:,mw)
                    else                                     
                        amdq(i,:) = amdq(i,:) + fwave(i,:,mw)
                    endif
                enddo
            endif
        enddo
    endif
end subroutine rp1
! ============================================================================



! ============================================================================
!  Eigenvalue routines
! ============================================================================

! ============================================================================
!  Eigenspace reconstruction using linear approximation
subroutine linear_eigen(h_l,h_r,u_l,u_r,b_l,b_r,            &
                         &  transonic_wave,wave_correction,s,eig_vec)

    use parameters_module, only: r,g,dry_tolerance,entropy_fix

    implicit none
    
    ! I/O
    double precision, intent(in) :: h_l(2),h_r(2),u_l(2),u_r(2),b_l,b_r
    
    integer, intent(inout) :: transonic_wave
    double precision, intent(inout) :: wave_correction
    double precision, intent(inout) :: s(4),eig_vec(4,4)
    
    ! Locals
    integer :: mw
    double precision :: alpha(4),speeds(4,2),gamma_l,gamma_r
    double precision :: h_ave(2),u_ave(2)
    double precision :: s_l(4),s_r(4),s_ave(4),work_vec(4,4)
        
    gamma_l = h_l(2) / h_l(1)
    gamma_r = h_r(2) / h_r(1)

    ! Left state alphas
    alpha(1) = 0.5d0*(gamma_l-1.d0+sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
    alpha(2) = 0.5d0*(gamma_l-1.d0-sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
    ! Right state alphas
    alpha(3) = 0.5d0*(gamma_r-1.d0-sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))
    alpha(4) = 0.5d0*(gamma_r-1.d0+sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))

    ! Left state speeds
    speeds(1,1) = u_l(1) - sqrt(g*h_l(1)*(1+alpha(1)))
    speeds(2,1) = u_l(2) - sqrt(g*h_l(1)*(1+alpha(2)))
    speeds(3,1) = u_l(2) + sqrt(g*h_l(1)*(1+alpha(2)))
    speeds(4,1) = u_l(1) + sqrt(g*h_l(1)*(1+alpha(1)))

    ! Right state speeds
    speeds(1,2) = u_r(1) - sqrt(g*h_r(1)*(1+alpha(4)))
    speeds(2,2) = u_r(2) - sqrt(g*h_r(1)*(1+alpha(3)))
    speeds(3,2) = u_r(2) + sqrt(g*h_r(1)*(1+alpha(3)))
    speeds(4,2) = u_r(1) + sqrt(g*h_r(1)*(1+alpha(4)))

    ! Determine wave speeds
    transonic_wave = 0
    wave_correction = 0.d0
    if (entropy_fix) then
        transonic_wave = 0
        wave_correction = 0.d0
        do mw=1,4
            ! Both speeds right going
            if (speeds(mw,1) > 0.d0) then
                s(mw) = speeds(mw,2)
            ! Transonic rarefaction
            ! s_l < 0.0 < s_r
            else if (speeds(mw,2) > 0.d0) then
                ! Assign base speed for rarefaction, should approximate true speed
                s(mw) = 0.5d0 * sum(speeds(mw,:))
                transonic_wave = mw
                wave_correction = abs(speeds(mw,1)) / (abs(speeds(mw,1)) + abs(speeds(mw,2)))
            else
                s(mw) = speeds(mw,1)
            endif
        enddo
    else
        s(:) = [speeds(1,1),speeds(2,1),speeds(3,2),speeds(3,2)]
    endif

    eig_vec(1,:) = 1.d0
    eig_vec(2,:) = s(:)
    eig_vec(3,:) = alpha
    eig_vec(4,:) = s(:)*alpha(:)

end subroutine linear_eigen
! ============================================================================


! ============================================================================
!  Eigenspace reconstruction using the velocity difference expansion
subroutine velocity_eigen(h_l,h_r,u_l,u_r,b_l,b_r,            &
                         &  transonic_wave,wave_correction,s,eig_vec)

    use parameters_module, only: r,g,one_minus_r

    implicit none
    
    ! I/O
    double precision, intent(in) :: h_l(2),h_r(2),u_l(2),u_r(2),b_l,b_r
    
    integer, intent(inout) :: transonic_wave
    double precision, intent(inout) :: wave_correction
    double precision, intent(inout) :: s(4),eig_vec(4,4)
    ! Locals
    double precision :: total_depth_l,total_depth_r,mult_depth_l,mult_depth_r
    double precision :: alpha(4)
    
    total_depth_l = sum(h_l)
    total_depth_r = sum(h_r)
    mult_depth_l = product(h_l)
    mult_depth_r = product(h_r)

    s(1) = (h_l(1)*u_l(1)+h_l(2)*u_l(2)) / total_depth_l &
          - sqrt(g*total_depth_l)
    s(2) = (h_l(2)*u_l(1)+h_l(1)*u_l(2)) / total_depth_l &
          - sqrt(g*one_minus_r*mult_depth_l/total_depth_l &
          * (1-(u_l(1)-u_l(2))**2/(g*one_minus_r*total_depth_l)))
    s(3) = (h_r(2)*u_r(1)+h_r(1)*u_r(2)) / total_depth_r &
          + sqrt(g*one_minus_r*mult_depth_r/total_depth_r &
          * (1-(u_r(1)-u_r(2))**2/(g*one_minus_r*total_depth_r)))
    s(4) = (h_r(1)*u_r(1)+h_r(2)*u_r(2)) / total_depth_r &
            + sqrt(g*total_depth_r)

    ! Determine wave speeds
    transonic_wave = 0
    wave_correction = 0.d0

    alpha(1:2) = ((s(1:2) - u_l(1))**2 - g * h_l(1)) / (g*h_l(1))
    alpha(3:4) = ((s(3:4) - u_r(1))**2 - g * h_r(1)) / (g*h_r(1))

    eig_vec(1,:) = 1.d0
    eig_vec(2,:) = s(:)
    eig_vec(3,:) = alpha
    eig_vec(4,:) = s(:)*alpha(:)

end subroutine velocity_eigen
! ============================================================================


! ============================================================================
!  Eigenspace reconstruction using LAPACK
subroutine lapack_eigen(h_l,h_r,u_l,u_r,b_l,b_r,            &
                         &  transonic_wave,wave_correction,s,eig_vec)

    use parameters_module, only: entropy_fix

    implicit none
    
    ! I/O
    double precision, intent(in) :: h_l(2),h_r(2),u_l(2),u_r(2),b_l,b_r
    
    integer, intent(inout) :: transonic_wave
    double precision, intent(inout) :: wave_correction
    double precision, intent(inout) :: s(4),eig_vec(4,4)
    
    ! Local
    integer :: j
    double precision :: h_ave(2),u_ave(2)
    double precision :: s_l(4),s_r(4),vec_work(4,4)
    
    ! Solve eigenvalue problem
    h_ave(:) = 0.5d0 * (h_l(:) + h_r(:))
    u_ave(:) = 0.5d0 * (u_l(:) + u_r(:))
    call eval_lapack_eigen(h_ave,u_ave,s,eig_vec)
    
    transonic_wave = 0
    wave_correction = 0.d0
    
    if (entropy_fix) then
        ! Check to see if we may be at a transonic rarefaction
        call eval_lapack_eigen(h_l,u_l,s_l,vec_work)        
        call eval_lapack_eigen(h_r,u_r,s_r,vec_work)
        
        ! Check each wave for a transonic problem
        do j=1,4
            if (s_l(j) < 0.d0 .and. 0.d0 < s_r(j)) then
                print *,"Transonic wave detected in wave family ",j,"."
                transonic_wave = j
                wave_correction = (s_l(j) - s(j)) / (s_l(j) - s_r(j))
            endif
        enddo
    endif
    
end subroutine lapack_eigen
! ============================================================================


! ============================================================================
! Single layer eigenspace reconstruction
subroutine single_layer_eigen(h_l,h_r,u_l,u_r,b_l,b_r,            &
                         &  transonic_wave,wave_correction,s,eig_vec)

    use parameters_module, only: g

    implicit none
    
    ! I/O
    double precision, intent(in) :: h_l(2),h_r(2),u_l(2),u_r(2),b_l,b_r
    
    integer, intent(inout) :: transonic_wave
    double precision, intent(inout) :: wave_correction
    double precision, intent(inout) :: s(4),eig_vec(4,4)
    
    transonic_wave = 0
    wave_correction = 0.d0
    
    s(1) = u_l(1) - sqrt(g*h_l(1))
    s(2) = 0.d0
    s(3) = 0.d0
    s(4) = u_r(1) + sqrt(g*h_r(1))
    
    eig_vec(1,:) = 1.d0
    eig_vec(2,:) = s(:)
    eig_vec(3,:) = 0.d0
    eig_vec(4,:) = 0.d0

end subroutine single_layer_eigen
! ============================================================================




! ============================================================================
!  Helper routines
! ============================================================================

! ============================================================================
!  Evaluate eigenvalues using LAPACK's DGEEV function
subroutine eval_lapack_eigen(h,u,lambda,vec)

    use parameters_module, only: r,g

    implicit none
    double precision, intent(in) :: h(2),u(2)
    double precision, intent(inout) :: lambda(4),vec(4,4)
    
    integer, parameter :: lwork = 16
    integer :: info
    double precision :: A(4,4),real_lambda(4),imaginary_lambda(4)
    double precision :: empty,work(1,lwork)

    ! Quasi-linear matrix
    A(1,:) = [0.d0,1.d0,0.d0,0.d0]
    A(2,:) = [-u(1)**2 + g*h(1),2.d0*u(1),g*h(1),0.d0]
    A(3,:) = [0.d0,0.d0,0.d0,1.d0]
    A(4,:) = [g*r*h(2),0.d0,-u(2)**2 + g*h(2),2.d0*u(2)]
    
    ! Call LAPACK eigen solver
    call dgeev('N','V',4,A,4,real_lambda,imaginary_lambda,empty,1,vec,4,work,lwork,info)
    if (info < 0) then
        info = -info
        print "(a,i1,a)","The ",info,"th argument had an illegal value."
        stop
    else if (info > 0) then
        print "(a)","The QR algorithm failed to compute all the"
        print "(a)","eigenvalues, and no eigenvectors have been"
        print "(a,i1,a)","computed; elements",info,"+1:4 of WR and WI"
        print "(a)","contain eigenvalues which have converged."
        stop
    endif
    ! Only need to check if there is a positive imaginary eigenvalue as they 
    ! should come in pairs
    if (any(imaginary_lambda > 0.d0)) then
        print "(a,4d16.8)","Imaginary eigenvalues computed: ",imaginary_lambda
        stop
    endif

end subroutine eval_lapack_eigen
! ============================================================================

! ============================================================================
!  Solve R beta = delta using LAPACK's DGESV function
function eval_lapack_solve(eig_vec,delta) result(beta)

    implicit none

    ! Matrix and RHS
    double precision, intent(in) :: eig_vec(4,4), delta(4)
    
    ! Output
    double precision :: beta(4)
    
    ! Local storage
    double precision :: R(4,4)
    integer :: ipiv(4),info,i,m

    ! Move input data into work arrays
    R = eig_vec
    beta = delta

    ! Call LAPACK dense linear solver
    call dgesv(4,1,R,4,ipiv,beta,4,info)
    if (.not.(info == 0)) then 
        print *, "Error solving R beta = delta, ",info
        print *, "  R=",(eig_vec(1,m),m=1,4)
        do i=2,4
            print *, "    ",(eig_vec(i,m),m=1,4)
        enddo
        print *, "  delta=",(delta(m),m=1,4)
        stop
    endif

end function eval_lapack_solve
! ============================================================================


