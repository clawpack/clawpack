!======================================================================
       subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx, &
                       ql,qr,auxl,auxr,wave,s,amdq,apdq)
!======================================================================
!
! Solves normal Riemann problems for the adjoint of the linearized
! 2D shallow water equations with topography.
! Uses q-waves, not f-waves, since the forward solve is in conservation
! form using f-waves.
!
! The full nonlinear forward problem in the x-direction is:
!     h_t + (hu)_x = 0
!     (hu)_t + (h*u^2 + 0.5*g*h^2)_x = -g*h*B_x(x)
! where h is the depth, hu the momentum, and B(x) the topography/bathymetry.
!
! The linearized forward problem is:
!     eta_t + mom_x = 0
!     mom_t + [g*h0(x) * eta]_x = 0
! where (eta, mom) are the surface displacement eta and linearized momentum,
! and h0(x) is the background depth (ocean at rest), so h0(x) = max(-B(x), 0).
!
! Note the forward problem is in conservation form q_t + (A(x)q)_x = 0, 
! so the adjoint equation is non-conservative and has the form
!     q_hat_t + A(x)^T q_hat_x = 0
! with coefficient matrix is
!     A(x)^T = [0  g*h0(x)]
!              [1     0   ]
!
! We need to solve the adjoint backwards in time, but we convert this into
! a problem going forward in time in order to apply GeoClaw, by solving
!     q_hat_t - A(x)^T q_hat_x = 0.
!
! Hence the Riemann solver uses the eigenvalues/vectors of -A(x)^T, which are:
!     lambda_1 = -c0,  r_1 = [+c0, 1]^T,
!     lambda_2 = +c0,  r_2 = [-c0, 1]^T.
!
! If h0 < drytol on one side of the interface, then "solid wall" BCs
! are used to find the single wave propagating into the ocean with mom=0 behind.
!
! On input, ql contains the state vector at the left edge of each cell
!     qr contains the state vector at the right edge of each cell
!
! This data is along a slice in the x-direction if ixy=1
!     or the y-direction if ixy=2.

!  Note that the i'th Riemann problem has left state qr(i-1,:)
!     and right state ql(i,:)
!  From the basic clawpack routines, this routine is called with
!     ql = qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           !
!      # This Riemann solver is for the shallow water adjoint               !
!             equations. It is modified from rpn2_geoclaw.f.                !
!                                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use geoclaw_module, only: g => grav, drytol => dry_tolerance
      use geoclaw_module, only: earth_radius, deg2rad
      use amr_module, only: mcapa, use_fwaves

      implicit none

      !input
      integer, intent(in) :: maxm,meqn,maux,mwaves,mbc,mx,ixy

      real(kind=8), intent(in) ::  ql(meqn, 1-mbc:maxm+mbc)
      real(kind=8), intent(in) ::  qr(meqn, 1-mbc:maxm+mbc)
      real(kind=8), intent(in) ::  auxl(maux,1-mbc:maxm+mbc)
      real(kind=8), intent(in) ::  auxr(maux,1-mbc:maxm+mbc)
      real(kind=8), intent(out) ::  wave(meqn, mwaves, 1-mbc:maxm+mbc)
      real(kind=8), intent(out) ::  s(mwaves, 1-mbc:maxm+mbc)
      real(kind=8), intent(out) ::  apdq(meqn,1-mbc:maxm+mbc)
      real(kind=8), intent(out) ::  amdq(meqn,1-mbc:maxm+mbc)

      ! local 
      integer :: m,i,mu,nv
      real(kind=8) :: delta(2)
      real(kind=8) :: h0R,h0L,huR,huL,etaR,etaL
      real(kind=8) :: hR,hL,bR,bL,cL,cR,alf1,alf2,dxdc
      logical dryL,dryR

      if (mwaves .ne. 2) then
         write(6,*) '*** Adjoint solver requires num_waves = 2'
         stop
         endif

      if (use_fwaves) then
         write(6,*) '*** Adjoint solver requires use_fwaves=False'
         stop
         endif

      ! initialize for dry areas:
      wave(:,:,:) = 0.d0
      s(:,:) = 1.d0
      amdq(:,:) = 0.d0
      apdq(:,:) = 0.d0

      !loop through Riemann problems at each grid cell
      do i=2-mbc,mx+mbc

         ! background depth based on topography:
         bL = auxr(1,i-1)
         bR = auxl(1,i)
         h0L = max(-bL, 0.d0)
         h0R = max(-bR, 0.d0)

         !skip problem if in a completely dry area
         if (h0L <= drytol .and. h0R <= drytol) then
             cycle ! go to next i
         endif
         
         !set normal direction
         if (ixy.eq.1) then
            mu=2
            nv=3
         else
            mu=3
            nv=2
         endif

        ! linearized eta:
        hL = qr(1,i-1)
        hR = ql(1,i)
        etaL = hL + bL
        etaR = hR + bR

        ! linearized momenta:
        huL = qr(mu,i-1)
        huR = ql(mu,i)

        ! wave speeds
        cL = sqrt(g*h0L) ! 1 wave speed is -cL
        cR = sqrt(g*h0R) ! 2 wave speed is +cR


        ! Check for wet/dry boundary
        if (hR <= drytol) then
            ! set "ghost cell" values for zero momentum BC:
            etaR = etaL
            huR = -huL
            dryR = .true.
            cR = cL
        else
            dryR = .false.
        endif

        if (hL <= drytol) then
            ! set "ghost cell" values for zero momentum BC:
            etaL = etaR
            huL = -huR
            dryL = .true.
            cL = cR
        else
            dryL = .false.
        endif

        
        ! Roe average speeds 
        ! More consistent with forward solver, but doesn't make much difference
        !cL = sqrt(g*0.5d0*(h0L+h0R)) ! 1 wave speed is -cL
        !cR = sqrt(g*0.5d0*(h0L+h0R)) ! 2 wave speed is +cR
        
        if ((cL<drytol) .and. (cR<drytol)) then
            ! should not happen
            write(6,*) '*** g,h0L,h0r: ',g,h0L,h0r
            endif

        !--------------------end initializing----------
        ! solve Riemann problem.

        ! q-wave splitting
        delta(1) = etaR - etaL
        delta(2) = huR - huL

        alf1 = (delta(1) + cR*delta(2))/(cL + cR)
        alf2 = (-delta(1) + cL*delta(2))/(cL + cR)

        ! eliminate one wave if dry to either side:
        if (dryL) alf1 = 0.d0
        if (dryR) alf2 = 0.d0

        ! Compute the waves.
        wave(1,1,i) = cL*alf1
        wave(mu,1,i) = alf1
        wave(nv,1,i) = 0.d0
        s(1,i) = -cL

        wave(1,2,i) = -cR*alf2
        wave(mu,2,i) = alf2
        wave(nv,2,i) = 0.d0
        s(2,i) = cR
        
!====== Adjust for mapping from latitude longitude to physical space====
        if (mcapa.gt.0) then
            if (ixy.eq.1) then
               dxdc=(earth_radius*deg2rad)
            else
               dxdc=earth_radius*cos(auxl(3,i))*deg2rad
            endif
            s(1,i) = dxdc * s(1,i)
            s(2,i) = dxdc * s(2,i)
        endif

        ! Fluctuations: 
        ! 1-wave is always left-going, 2-wave is right-going
        amdq(1:meqn,i) = s(1,i) * wave(1:meqn,1,i)
        apdq(1:meqn,i) = s(2,i) * wave(1:meqn,2,i)

      enddo

      return
      end subroutine rpn2
