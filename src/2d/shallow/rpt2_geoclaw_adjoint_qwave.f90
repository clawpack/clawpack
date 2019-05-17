! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr, &
                aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           !
!      # This Riemann solver is for the shallow water adjoint               !
!             equations. It is similar to rpt2_vc_acoustics.f90
!                                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     # Split asdq into down-going flux bmasdq and up-going flux bpasdq.
!     # imp=1  means  asdq=amdq,    imp=2 means asdq=apdq

    use geoclaw_module, only: g => grav, drytol => dry_tolerance
    use geoclaw_module, only: coordinate_system,earth_radius,deg2rad

    implicit none

    integer, intent(in) :: ixy,imp,maxm,meqn,mwaves,maux,mbc,mx
    real(kind=8), intent(in) ::    ql(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) ::    qr(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) ::    asdq(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) ::   aux1(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) ::   aux2(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) ::   aux3(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(out) :: bmasdq(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(out) :: bpasdq(meqn, 1-mbc:maxm+mbc)

!   local:
    integer :: i,mu,mv,i1
    real(kind=8) :: h0D,h00,h0U,bD,b0,bU,c0,cU,cD,a1,a2,dxdcU,dxdcD,s1,s2

!   # initialize for dry state cells:
    bmasdq = 0.d0
    bpasdq = 0.d0
    
    if (ixy == 1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    do i = 2-mbc, mx+mbc

         if (maxval(abs(asdq(:,i))) == 0.d0) then
            cycle  ! go to next i
            endif
    
!        # imp is used to flag whether wave is going to left or right,
!        # since material properties are different on the two sides
    
        if (imp == 1) then
!            # asdq = amdq, moving to left
            i1 = i-1
        else
!            # asdq = apdq, moving to right
            i1 = i
        endif
        
    
!       # The flux difference asdq is split into downward moving part
!       # traveling at speed -cD relative to the medium below and
!       # an upward moving part traveling
!       # at speed +cU relative to the medium above.
   
!       # Note that the sum of these parts does not give all of asdq
!       # since there is also reflection at the interfaces which decreases
!       # the flux.
    
!       # gravity wave speed in each row of cells:

        ! background depth based on topography:
        bD = aux1(1,i1)
        b0 = aux2(1,i1)
        bU = aux3(1,i1)
        h0U = max(-bU, 0.d0)
        h00 = max(-b0, 0.d0)
        h0D = max(-bD, 0.d0)
        cD = sqrt(g*h0D)
        c0 = sqrt(g*h00)
        cU = sqrt(g*h0U)

        if (h00 < drytol) then
            cycle  ! no waves in this cell
            endif

!       # transmitted part of down-going wave:
        if (h0D < drytol) then
            a1 = 0.d0  ! no transmitted wave
        else
            a1 = (asdq(1,i) + asdq(mv,i)*c0) / (c0 + cD)
        endif

!       # transmitted part of up-going wave:
        if (h0U < drytol) then
            a2 = 0.d0  ! no transmitted wave
        else
            a2 = (-asdq(1,i) + asdq(mv,i)*c0) / (c0 + cU)
        endif
        
        !====== Adjust for mapping from latitude longitude to physical space====
        if (coordinate_system.eq.2) then
           if (ixy.eq.2) then
               dxdcU=(earth_radius*deg2rad)
               dxdcD = dxdcU
           else
               dxdcU = earth_radius*cos(aux3(3,i1))*deg2rad
               dxdcD = earth_radius*cos(aux1(3,i1))*deg2rad
           endif
        else
            dxdcU = 1.d0
            dxdcD = 1.d0
        endif

    
!       # The down-going flux difference bmasdq is the product  -cD * wave1
!       # down-going eigenvector of -A(x)^T is r1 = [cD, 1]^T
    
        s1 = -cD * dxdcD  ! adjusted down-going wave speed
        bmasdq(1,i) = s1 * a1*cD
        bmasdq(mu,i) = 0.d0
        bmasdq(mv,i) = s1 * a1
    
!       # The up-going flux difference bpasdq is the product  cU * wave2
!       # up-going eigenvector of -A(x)^T is r2 = [-cU, 1]^T
    
        s2 = cU * dxdcU  ! adjusted up-going wave speed
        bpasdq(1,i) = s2 * a2*(-cU)
        bpasdq(mu,i) = 0.d0
        bpasdq(mv,i) = s2 * a2
    
        enddo

    return
    end subroutine rpt2
