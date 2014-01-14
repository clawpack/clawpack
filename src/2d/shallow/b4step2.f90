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
    use topo_module, only: num_dtopo,topotime
    use topo_module, only: tfdtopo,t0dtopo,topo_finalized
    use topo_module, only: xlowdtopo,xhidtopo,ylowdtopo,yhidtopo

    use amr_module, only: xlowdomain => xlower
    use amr_module, only: ylowdomain => ylower
    use amr_module, only: xhidomain => xupper
    use amr_module, only: yhidomain => yupper
    use amr_module, only: xperdom,yperdom,spheredom

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn
    integer, intent(inout) :: mbc,mx,my,maux
    real(kind=8), intent(inout) :: xlower, ylower, dx, dy, t, dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Local storage
    integer :: index,i,j,k,dummy
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

    ! update topography if needed
    if ((minval(xlowdtopo)<= xlower + real(mx+mbc,kind=8)*dx).and. &
            (minval(ylowdtopo)<= ylower + real(my+mbc,kind=8)*dy).and. &
            (maxval(xhidtopo) >= xlower).and. &
            (maxval(yhidtopo) >= ylower) ) then
      call topo_update(t)
      call setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
!      call bc2amr(q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc), &
!               aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc), &
!               mx+mbc,my+mbc,meqn,maux,dx,dy,1,t, &
!               xlower-real(mbc,kind=8)*dx,xlower + real(mx+mbc,kind=8)*dx, &
!               ylower-real(mbc,kind=8)*dy,ylower + real(mbc,kind=8)*dy, &
!               xlowdomain,ylowdomain,xhidomain,yhidomain,xperdom,yperdom,spheredom)
   endif
   ! unfortunately if setaux is not called after topo_finalized in the above loop, it's possible
   ! that even though all topofiles are finalized in time,
   ! some levels of grids may not have their aux values updated to current time
   ! the only fix I can come to is to always call setaux when overlapping dtopo
   !
   ! Also, should we require the intersection clause even for the above topo_update?
   ! that seems safe...at the moment
   ! the rationale for requiring an intersection to call topo_update is if
   ! there are very fine grids taking many timesteps away from moving topography
   ! such as in the case of highly resolved near-shore grids, topo might only need to
   ! get reset at much coarser time intervals of coarser grids overlapping dtopo
   ! Dave George, Jan 10 2014


end subroutine b4step2
