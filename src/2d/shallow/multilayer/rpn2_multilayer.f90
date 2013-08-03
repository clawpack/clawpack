subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! ============================================================================
!  Solves normal Riemann problem for the multilayer shallow water equations in
!  2D with topography and wind forcing:
!    (h_1)_t + (h_1 u_1)_x + (h_1 v_1)_y = 0
!    (h_1 u_1)_t + (h_1 u_1^2 + 1/2 g h_1^2)_x + (h_1 u_1 v_1)_y = -gh_1(r(h_2)_x + B_x)
!    (h_1 v_1)_t + (h_1 u_1 v_1)_x + (h_1 v_1^2 + 1/2 g h_1^2)_y = -gh_1(r(h_2)_y + B_y)
!    (h_2)_t + (h_2 u_2)_x + (h_2 v_2)_y = 0
!    (h_2 u_2)_t + (h_2 u_2^2 + 1/2 g h_2^2)_x + (h_2 u_2 v_2)_y = -gh_2(h_1 + B)_x + Tau |W| W_x
!    (h_2 v_2)_t + (h_2 u_2 v_2)_x + (h_2 v_2^2 + 1/2 g h_2^2)_y = -gh_2(h_1 + B)_y + Tau |W| W_y
!
!  On input, ql contains the state vector at the left edge of each cell and qr
!  contains the state vector at the right edge of each cell
!
!           |            |          
!    qr(i-1)|ql(i)  qr(i)|ql(i+1)   
! ----------|------------|----------
!    i-1          i         i+1
!
!  The i-1/2 Riemann problem has left state qr(i-1) and right state ql(i)
!
!  If ixy == 1 then the sweep direction is x, ixy == 2 implies the y direction
!
!  Kyle T. Mandli (10-11-2010)
! ============================================================================

    use geoclaw_module
    use amr_module, only: mcapa
    use hurricane_module
    use multilayer_module

    implicit none

    ! Input arguments
    integer, intent(in) :: ixy,maxm,meqn,mwaves,mbc,mx
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql,qr
    double precision, dimension(ml_maux,1-mbc:maxm+mbc), intent(in) :: auxl,auxr

    ! Output arguments
    double precision, dimension(meqn, mwaves, 1-mbc:maxm+mbc), intent(out) :: fwave
    double precision, dimension(mwaves, 1-mbc:maxm+mbc), intent(out) :: s
    double precision, dimension(meqn, 1-mbc:maxm+mbc), intent(out) :: apdq, amdq

    ! Counters
    integer :: i,j,m,mw,k,maxiter,info
    integer :: n_index,t_index,layer_index
    
    ! Physics
    double precision :: g,dxdc
    
    ! State variables
    double precision, dimension(2) :: h_l,h_r,hu_l,hu_r,hv_l,hv_r
    double precision, dimension(2) :: u_l,u_r,v_l,v_r,advected_speed
    double precision :: b_l,b_r,w_normal,w_transverse,kappa_l,kappa_r
    double precision :: eta_l(2),eta_r(2),h_ave(2),momentum_transfer(2)
    double precision :: h_hat_l(2),h_hat_r(2),gamma_l,gamma_r
    double precision :: flux_transfer_l,flux_transfer_r,lambda(6)
    double precision :: temp_depth(2),temp_u(2),temp_v(2)

    ! Solver variables
    double precision, dimension(6) :: delta,flux_r,flux_l,pivot
    double precision, dimension(6,6) :: eig_vec,A
    double precision :: beta(6),alpha(4)
    logical :: dry_state_l(2), dry_state_r(2)
    
    ! Single layer locals
    integer, parameter :: max_iterations = 1
    double precision :: wall(3),fw(3,3),sw(3),phi_r(2),phi_l(2)
    double precision :: s_l,s_r,s_roe(2),s_E(2),u_hat,c_hat,sm(2)
    double precision :: h_star,h_star_test,h_star_HLL,s_l_test,s_r_test
    logical :: rare(2)

    ! Common block variables
    integer :: icom,jcom
    double precision :: dtcom,dxcom,dycom,tcom
    
    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

    external dgesv
    
    ! Initialize output variables
    amdq = 0.d0
    apdq = 0.d0
    
    ! Set normal direction
    if (ixy == 1) then
        n_index = 2
        t_index = 3
    else
        n_index = 3
        t_index = 2
    endif

    g = grav
    
    ! ========================================================================
    ! Loop through Riemann problems
    ! ========================================================================
    do i=2-mbc,mx+mbc
        dry_state_l = .false.
        dry_state_r = .false.
        
        ! Parse states and set appropriate zeros
        ! Note that the "u-direction" is the direction of sweeping which 
        ! could actually be the x or y-directions depending on ixy
        
        do j=1,2
            layer_index = 3*(j-1)
            h_l(j) = qr(layer_index+1,i-1) / rho(j)
            hu_l(j) = qr(layer_index+n_index,i-1) / rho(j)
            hv_l(j) = qr(layer_index+t_index,i-1) / rho(j)
            
            h_r(j) = ql(layer_index+1,i) / rho(j)
            hu_r(j) = ql(layer_index+n_index,i) / rho(j)
            hv_r(j) = ql(layer_index+t_index,i) / rho(j)
            
            h_hat_l(j) = auxr(j+6,i-1)
            h_hat_r(j) = auxl(j+6,i)
            
            h_ave(:) = 0.5d0 * (h_l(:) + h_r(:))
            
            ! Check for dry states
            if (h_l(j) < drytolerance) then
                dry_state_l(j) = .true.
                hu_l(j) = 0.d0
                hv_l(j) = 0.d0
                u_l(j) = 0.d0
                v_l(j) = 0.d0
            else
                u_l(j) = hu_l(j) / h_l(j)
                v_l(j) = hv_l(j) / h_l(j)
            endif
            if (h_r(j) < drytolerance) then
                dry_state_r(j) = .true.
                hu_r(j) = 0.d0
                hv_r(j) = 0.d0
                u_r(j) = 0.d0
                v_r(j) = 0.d0
            else
                u_r(j) = hu_r(j) / h_r(j)
                v_r(j) = hv_r(j) / h_r(j)
            endif
        enddo
        
        b_l = auxr(1,i-1)
        b_r = auxl(1,i)
        
        ! Calculate wind stress
!         w_normal = 0.5d0 * (auxr(i-1,n_index+2) + auxl(i,n_index+2))
!         w_transverse = 0.5d0 * (auxr(i-1,t_index+2) + auxl(i,t_index+2))
!         wind_speed = sqrt(w_normal**2 + w_transverse**2)
!         tau = wind_drag(wind_speed) * rho_air * wind_speed
!         if (ixy == 1) then
!             wind_stress = tau * w_normal
!         else if (ixy == 2) then
!             wind_stress = tau * w_normal
!         endif

        ! ====================================================================
        !  Top layer only
        ! ====================================================================            
        if (dry_state_l(2).and.dry_state_r(2)) then
            wall = 1.d0
            
            ! Completely dry cell
            if (dry_state_l(1).and.dry_state_r(1)) then
                s(:,i) = 0.d0
                fwave(:,:,i) = 0.d0
                cycle
            endif
            
            ! Calculate momentum fluxes
            phi_l(1) = 0.5d0 * g * h_l(1)**2 + h_l(1) * u_l(1)**2
            phi_r(1) = 0.5d0 * g * h_r(1)**2 + h_r(1) * u_r(1)**2
             
            ! Check for dry state to right
            if (h_r(1) <= drytolerance) then
                call riemanntype(h_l(1),h_l(1),u_l(1),-u_l(1),h_star, &
                                 sm(1),sm(2),rare(1),rare(2),1,drytolerance,g)
                h_star_test = max(h_l(1),h_star)
                ! Right state should become ghost values that mirror left for wall problem
                if (h_star_test + b_l < b_r) then 
                    wall(2:3)=0.d0
                    h_r(1) = h_l(1)
                    hu_r(1) = -hu_l(1)
                    b_r = b_l
                    phi_r(1) = phi_l(1)
                    u_r(1) = -u_l(1)
                    v_r(1) = v_l(1)
                elseif (h_l(1) + b_l < b_r) then
                    b_r = h_l(1) + b_l
                endif
            ! Check for drystate to left, i.e right surface is lower than left topo
            else if (h_l(1) <= drytolerance) then 
                call riemanntype(h_r(1),h_r(1),-u_r(1),u_r(1),h_star, &
                                 sm(1),sm(2),rare(1),rare(2),1,drytolerance,g)
                h_star_test = max(h_r(1),h_star)
                ! Left state should become ghost values that mirror right
                if (h_star_test + b_r < b_l) then  
                   wall(1:2) = 0.d0
                   h_l(1) = h_r(1)
                   hu_l(1) = -hu_r(1)
                   b_l = b_r
                   phi_l(1) = phi_r(1)
                   u_l(1) = -u_r(1)
                   v_l(1) = v_r(1)
                elseif (h_r(1) + b_r < b_l) then
                   b_l = h_r(1) + b_r
                endif
             endif

             ! Determine wave speeds
             s_l = u_l(1) - sqrt(g*h_l(1)) ! 1 wave speed of left state
             s_r = u_r(1) + sqrt(g*h_r(1)) ! 2 wave speed of right state
             
             u_hat = (sqrt(g*h_l(1))*u_l(1) + sqrt(g*h_r(1))*u_r(1)) &
                        / (sqrt(g*h_r(1))+sqrt(g*h_l(1))) ! Roe average
             c_hat = sqrt(g*0.5d0*(h_r(1)+h_l(1))) ! Roe average
             s_roe(1) = u_hat - c_hat ! Roe wave speed 1 wave
             s_roe(2) = u_hat + c_hat ! Roe wave speed 2 wave

             s_E(1) = min(s_l,s_roe(1)) ! Eindfeldt speed 1 wave
             s_E(2) = max(s_r,s_roe(2)) ! Eindfeldt speed 2 wave

             ! Solve Riemann problem
             call riemann_aug_JCP(max_iterations,3,3,h_l(1),h_r(1),hu_l(1), &
                                    hu_r(1),hv_l(1),hv_r(1),b_l,b_r,u_l(1), &
                                    u_r(1),v_l(1),v_r(1),phi_l(1),phi_r(1), &
                                    s_E(1),s_E(2),drytolerance,g,sw,fw)
            
            ! Eliminate ghost fluxes for wall
            do mw=1,3
                sw(mw)=sw(mw)*wall(mw)
                do m=1,3
                   fw(m,mw)=fw(m,mw)*wall(mw)
                enddo
            enddo

            ! Update speeds and waves
            ! Note that we represent all the waves in the first three arrays
            ! so it does not directly correspond to the two-layer case's wave
            ! structure
            s(:,i) = 0.d0
            fwave(:,:,i) = 0.d0
            
            s(1:3,i) = sw(:)
            fwave(1,1:3,i) = fw(1,:) * rho(1)
            fwave(n_index,1:3,i) = fw(2,:) * rho(1)
            fwave(t_index,1:3,i) = fw(3,:) * rho(1)
            
            ! Go on to next cell, lat-long and fluctuation calculations are 
            ! outside of this loop
            cycle
        endif

        
        ! ====================================================================
        !  Calculate Eigenstructure
        !   The parameter eigen_method if either a completely wet problem or a
        !   wall dry state problem exists.  Otherwise the inundation_method
        !   controls the method for the eigenspace calculation.  rare(1) is
        !   set to true if inundation occurs into the right state and rare(2)
        !   if inundation occurs in the left state.
        rare = .false.
        
        ! Dry state only to right
        if (dry_state_r(2).and.(.not.dry_state_l(2))) then
            ! Inundation occurs
            if (h_l(2) + b_l > b_r) then
                rare(1) = .true.
                if (inundation_method == 1) then
                    ! Linearized static eigenspace with zero depth
                    temp_depth = [h_r(1),0.d0]
                    call linearized_eigen(h_l,temp_depth,u_l,u_r,v_l,v_r, &
                        n_index,t_index,lambda,eig_vec)
                else if (inundation_method == 2) then
                    ! Linearized with eigenspace with small depth
                    temp_depth = [h_r(1),drytolerance]
                    call linearized_eigen(h_l,temp_depth,u_l,u_r,v_l,v_r, &
                        n_index,t_index,lambda,eig_vec)
                else if (inundation_method == 3) then
                    ! Velocity difference with eigenspace with small depth
                    temp_depth = [h_r(1),drytolerance]
                    call vel_diff_eigen(h_l,temp_depth,u_l,u_r,v_l,v_r, &
                        n_index,t_index,lambda,eig_vec)
                else if (inundation_method == 4) then
                    ! LAPACK with zero depth
                    temp_depth = [h_r(1),0.d0]
                    call lapack_eigen(h_l,temp_depth,u_l,u_r,v_l,v_r, &
                        n_index,t_index,lambda,eig_vec)
                else if (inundation_method == 5) then
                    ! LAPACK with small depth
                    temp_depth = [h_r(1),drytolerance]
                    call lapack_eigen(h_l,temp_depth,u_l,u_r,v_l,v_r, &
                        n_index,t_index,lambda,eig_vec)
                endif
                s(:,i) = lambda
                
                ! Internal wave correction
                if (inundation_method == 1) then
                    s(5,i) = u_l(2) + 2.d0 * sqrt(g*(1.d0-r)*h_l(2))
                    alpha(3) = r * g * h_l(2) / ((s(5,i) - u_l(2))**2 - g * h_l(2))
                    
                    eig_vec(1,5) = 1.d0
                    eig_vec(n_index,5) = s(i,5)
                    eig_vec(t_index,5) = v_r(1)
                    eig_vec(4,5) = alpha(3)
                    eig_vec(n_index,5) = s(i,5) * alpha(3)
                    eig_vec(t_index,6) = v_r(2) * alpha(3)
                endif
                ! Fast wave correction
                if (inundation_method /= 5) then
                    s(6,i) = u_r(1) + sqrt(g*h_r(1))
                    eig_vec(1,6) = 1.d0
                    eig_vec(n_index,6) = s(6,i)
                    eig_vec(t_index,6) = v_r(1)
                    eig_vec(4:6,6) = 0.d0
                endif
            ! Wall boundary RP
            else                
                ! Wall state
                temp_depth = [h_r(1),0.d0]
                temp_u = [u_r(1),-u_l(2)]
                temp_v = [v_r(1),v_l(2)]
                if (eigen_method == 1) then
                    call linearized_eigen(h_hat_l,h_hat_r,u_l,temp_u,v_l,temp_v,n_index,t_index,lambda,eig_vec)
                else if (eigen_method == 2 .or. eigen_method == 4) then
                    call linearized_eigen(h_l,temp_depth,u_l,temp_u,v_l,temp_v,n_index,t_index,lambda,eig_vec)
                else if (eigen_method == 3) then
                    call vel_diff_eigen(h_l,temp_depth,u_l,temp_u,v_l,temp_v,n_index,t_index,lambda,eig_vec)
                endif
                s(:,i) = lambda
            endif
        else if (dry_state_l(2).and.(.not.dry_state_r(2))) then
            ! Inundation
            if (h_r(2) + b_r > b_l) then
                rare(2) = .true.
                if (inundation_method == 1) then
                    temp_depth = [h_l(1),0.d0]
                    call linearized_eigen(temp_depth,h_r,u_l,u_r,v_l,v_r,n_index,t_index,lambda,eig_vec)
                else if (inundation_method == 2) then
                    temp_depth = [h_l(1),drytolerance]
                    call linearized_eigen(temp_depth,h_r,u_l,u_r,v_l,v_r,n_index,t_index,lambda,eig_vec)
                else if (inundation_method == 3) then
                    temp_depth = [h_l(1),drytolerance]
                    call vel_diff_eigen(temp_depth,h_r,u_l,u_r,v_l,v_r,n_index,t_index,lambda,eig_vec)
                else if (inundation_method == 4) then
                    temp_depth = [h_l(1),drytolerance]
                    call lapack_eigen(temp_depth,h_r,u_l,u_r,v_l,v_r,n_index,t_index,lambda,eig_vec)
                else if (inundation_method == 5) then
                    temp_depth = [h_l(1),0.d0]
                    call lapack_eigen(temp_depth,h_r,u_l,u_r,v_l,v_r,n_index,t_index,lambda,eig_vec)
                endif
                s(:,i) = lambda    
                            
                ! Internal wave correction
                if (inundation_method == 1) then
                    s(2,i) = u_r(2) - 2.d0 * sqrt(g*(1.d0-r)*h_r(2))
                    alpha(2) = r * g * h_r(2) / ((s(2,i) - u_r(2))**2 - g*h_r(2))
                    eig_vec(1,2) = 1.d0
                    eig_vec(n_index,2) = s(2,i) 
                    eig_vec(t_index,2) = v_l(1)
                    eig_vec(4,2) = alpha(2)
                    eig_vec(n_index,2) = alpha(2)*s(2,i)
                    eig_vec(t_index,2) = alpha(2)*v_l(2)
                endif
                ! Fast wave correction
                if (inundation_method /= 5) then
                    s(1,i) = u_l(1) - sqrt(g*h_l(1))
                    eig_vec(1,1) = 1.d0
                    eig_vec(n_index,1) = s(1,i)
                    eig_vec(t_index,1) = v_l(1)
                    eig_vec(4:6,1) = 0.d0
                endif
            ! Wall boundary
            else
                ! Wall state
                temp_depth = [h_l(1),0.d0]
                temp_u = [u_l(1),-u_r(2)]
                temp_v = [v_l(1),v_r(2)]
                if (eigen_method == 1) then
                    call linearized_eigen(h_hat_l,h_hat_r,temp_u,u_r,temp_v,v_r,n_index,t_index,lambda,eig_vec)
                else if (eigen_method == 2 .or. eigen_method == 4) then
                    call linearized_eigen(h_l,temp_depth,temp_u,u_r,temp_v,v_r,n_index,t_index,lambda,eig_vec)
                else if (eigen_method == 3) then
                    call vel_diff_eigen(h_l,temp_depth,temp_u,u_r,temp_v,v_r,n_index,t_index,lambda,eig_vec)
                endif
                s(:,i) = lambda
            endif
        ! Completely wet state
        else            
            if (eigen_method == 1) then
                call linearized_eigen(h_hat_l,h_hat_r,u_l,u_r,v_l,v_r, &
                    n_index,t_index,lambda,eig_vec)
            else if (eigen_method == 2) then
                call linearized_eigen(h_l,h_r,u_l,u_r,v_l,v_r,n_index,t_index, &
                    lambda,eig_vec)
            else if (eigen_method == 3) then
                call vel_diff_eigen(h_l,h_r,u_l,u_r,v_l,v_r,n_index,t_index, &
                    lambda,eig_vec)
            else if (eigen_method == 4) then
                call lapack_eigen(h_l,h_r,u_l,u_r,v_l,v_r,n_index,t_index, &
                    lambda,eig_vec)
            endif
            s(:,i) = lambda
        endif
        
        ! ====================================================================
        ! Compute jump in fluxes
        ! Dry state, bottom layer to right
        if(dry_state_r(2).and.(.not.dry_state_l(2)).and.(.not.rare(1))) then
            h_r(2) = h_l(2)
            hu_r(2) = -hu_l(2)
            u_r(2) = -u_l(2)
            hv_r(2) = hv_l(2)
            v_r(2) = v_l(2)
        
            flux_transfer_r = 0.d0
            flux_transfer_l = 0.d0
            momentum_transfer(1) = g * rho(1) * h_ave(1) * (b_r - h_l(2) - b_l)
            momentum_transfer(2) = 0.d0
        ! ====================================================================
        ! Dry state, bottom layer to left
        else if(dry_state_l(2).and.(.not.dry_state_r(2)).and.(.not.rare(2))) then    
            h_l(2) = h_r(2)
            hu_l(2) = -hu_r(2)
            u_l(2) = -u_r(2)
            hv_l(2) = hv_r(2)
            v_l(2) = v_r(2)
        
            flux_transfer_r = 0.d0
            flux_transfer_l = 0.d0
            momentum_transfer(1) = g * rho(1) * h_ave(1) * (b_r + h_r(2) - b_l)
            momentum_transfer(2) = 0.d0
        ! ====================================================================
        ! Full two layer case
        else
            momentum_transfer(1) =  g * rho(1) * h_ave(1) * (h_r(2) - h_l(2) + b_r - b_l)
            momentum_transfer(2) = -g * rho(1) * h_ave(1) * (h_r(2) - h_l(2)) + g * rho(2) * h_ave(2) * (b_r - b_l)
            flux_transfer_r = g * rho(1) * h_r(1) * h_r(2)
            flux_transfer_l = g * rho(1) * h_l(1) * h_l(2)
        endif
        
        do j=1,2
            layer_index = 3*(j-1)
            flux_r(layer_index+1) = rho(j) * hu_r(j)
            flux_r(layer_index+n_index) = rho(j) * (h_r(j) * u_r(j)**2 + 0.5d0 * g * h_r(j)**2)
            flux_r(layer_index+t_index) = rho(j) * h_r(j) * u_r(j) * v_r(j)
            
            flux_l(layer_index+1) = rho(j) * hu_l(j)
            flux_l(layer_index+n_index) = rho(j) * (h_l(j) * u_l(j)**2 + 0.5d0 * g * h_l(j)**2)
            flux_l(layer_index+t_index) = rho(j) * h_l(j) * u_l(j) * v_l(j)
        enddo
        ! Add extra flux terms
        flux_r(3 + n_index) = flux_r(3 + n_index) + flux_transfer_r
        flux_l(3 + n_index) = flux_l(3 + n_index) + flux_transfer_l
        
        delta = flux_r - flux_l
            
        ! Momentum transfer and bathy terms
        delta(n_index) = delta(n_index) + momentum_transfer(1)
        delta(n_index+3) = delta(n_index+3) + momentum_transfer(2)
        
        ! ====================================================================
        ! Project jump in fluxes - Use LAPACK's dgesv routine
        !    N - (int) - Number of linear equations (6)
        !    NRHS - (int) - Number of right hand sides (1)
        !    A - (dp(6,6)) - Coefficient matrix, in this case eig_vec
        !    LDA - (int) - Leading dimension of A (6)
        !    IPIV - (int(N)) - Pivot indices
        !    B - (dp(LDB,NRHS)) - RHS of equations (delta)
        !    LDB - (int) - Leading dimension of B (6)
        !    INFO - (int) - Status of result
        !  Note that the solution (betas) are in delta after the call
        A = eig_vec ! We need to do this as the return matrix is modified and
                    ! we have to use eig_vec again to compute fwaves
        call dgesv(6,1,A,6,pivot,delta,6,info)
        if (.not.(info == 0)) then
            if (dry_state_l(2)) then
                print *,"left dry"
            endif
            if (dry_state_r(2)) then
                print *,"right dry"
            endif
            print *,"        left            |             right"
            print *,"====================================================="
            print *,h_l(1),h_r(1)
            print *,hu_l(1),hu_r(1)
            print *,hv_l(1),hv_r(1)
            print *,h_l(2),h_r(2)
            print *,hu_l(2),hu_r(2)
            print *,hv_l(2),hv_r(2)
            print *,b_l,b_r
            print *,""
            print "(a,i2)","In normal solver: ixy=",ixy
            print "(a,i3)","  Error solving R beta = delta,",info
            print "(a,i3,a,i3)","  Location: ",icom," ",jcom
            print "(a,6d16.8)","  Eigenspeeds: ",s(i,:)
            print "(a)","  Eigenvectors:"
            do j=1,6
                print "(a,6d16.8)","  ",(eig_vec(j,mw),mw=1,6)
            enddo
            stop
        endif
        beta = delta

        ! ====================================================================
        ! Compute fwaves
        forall(mw=1:mwaves)
            fwave(:,mw,i) = eig_vec(:,mw) * beta(mw)
        end forall
            
    enddo
    ! == End of Riemann Solver Loop per grid cell ============================
    
    ! ========================================================================
    ! Capacity for mapping from latitude longitude to physical space
    if (mcapa > 0) then
        do i=2-mbc,mx+mbc
            if (ixy == 1) then
                dxdc=(earth_radius*deg2rad)
            else
                dxdc=auxl(i,3)
            endif

            do mw=1,mwaves
    	        s(mw,i)=dxdc*s(mw,i)
    	        fwave(:,mw,i)=dxdc*fwave(:,mw,i)
            enddo
        enddo
    endif

    ! ========================================================================
    !  Compute fluctuations 
    do i=2-mbc,mx+mbc
        do mw=1,mwaves
            if (s(mw,i) > 0.d0) then
                apdq(:,i) = apdq(:,i) + fwave(:,mw,i)
            else
                amdq(:,i) = amdq(:,i) + fwave(:,mw,i)
            endif
            h_r(1) = ql(1,i) / rho(1)
            h_l(1) = qr(1,i-1) / rho(1)
            h_r(2) = ql(4,i) / rho(2)
            h_l(2) = qr(4,i-1) / rho(2)
            dry_state_r(2) = h_r(2) < drytolerance
            dry_state_l(2) = h_l(2) < drytolerance
            rare(1) = h_l(2) + b_l > b_r
            rare(2) = h_r(2) + b_r > b_l
            if (dry_state_r(2).and.(.not.dry_state_l(2)).and.(.not.rare(1)) &
             .or.(dry_state_l(2).and.(.not.dry_state_r(2)).and.(.not.rare(2)))) then
                do m=4,6
                    if (apdq(i,m) /= 0.d0) then
                        print *,"========================"
                        print *,"Wave ",mw," equation ",m
                        print *,"s = ",s(mw,i)
                        print *,"f = ",fwave(m,mw,i)
                        print *,"amdq = ",(amdq(m,i))
                        print *,"apdq = ",(apdq(m,i))
                        stop "Flux non-zero going into a wall, aborting calculation."
                    endif
                enddo
                apdq(4:6,i) = 0.d0
            endif
        enddo
    enddo

end subroutine rpn2