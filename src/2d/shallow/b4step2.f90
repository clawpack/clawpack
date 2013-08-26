! ============================================
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux)
! ============================================
! 
! # called before each call to step
! # use to set time-dependent aux arrays or perform other tasks.
! 
! This particular routine sets negative values of q(1,i,j) to zero,
! as well as the corresponding q(m,i,j) for m=1,meqn.
! This is for problems where q(1,i,j) is a depth.
! This should occur only because of rounding error.
! 
! Also calls movetopo if topography might be moving.

    use geoclaw_module, only: dry_tolerance
    use geoclaw_module, only: g => grav
    use topo_module
    use dtopo_module
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower, ylower, dx, dy, t, dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Local storage
    integer :: index,i,j,k
    real(kind=8) :: h,u,v

    ! Check for NaNs in the solution
    call check4nans(meqn,mbc,mx,my,q,t,1)

    ! check for h < 0 and reset to zero
    ! check for h < drytolerance
    ! set hu = hv = 0 in all these cells
    forall(i=1-mbc:mx+mbc, j=1-mbc:my+mbc,q(1,i,j) < dry_tolerance)
        q(1,i,j) = max(q(1,i,j),0.d0)
        q(2:3,i,j) = 0.d0
    end forall

    ! Move the topography if needed
    ! write(26,*) 'B4STEP2: t, num_dtopo: ', t,num_dtopo

    do i=1,num_dtopo
        call movetopo(mbc,mx,my,                                  &
                      xlower,ylower,dx,dy,t,dt,maux,aux,                      &
                      dtopowork(i0dtopo(i):i0dtopo(i)+mdtopo(i)-1),           &
                      xlowdtopo(i),ylowdtopo(i),xhidtopo(i),yhidtopo(i),      &
                      t0dtopo(i),tfdtopo(i),dxdtopo(i),dydtopo(i),dtdtopo(i), &
                      mxdtopo(i),mydtopo(i),mtdtopo(i),mdtopo(i),             &
                      minleveldtopo(i),maxleveldtopo(i),topoaltered(i))
    enddo

end subroutine b4step2
