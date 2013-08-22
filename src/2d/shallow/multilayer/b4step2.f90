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

    use geoclaw_module, only:num_layers,dry_tolerance,rho
    use geoclaw_module, only:KAPPA_UNIT,check_richardson,richardson_tolerance
    use geoclaw_module, only:g => grav
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
    real(kind=8) :: h(num_layers),u(num_layers),v(num_layers)
    real(kind=8) :: kappa,one_minus_r
    logical :: dry_state(num_layers)

    ! Format strings
    character(len=*), parameter :: hyp_warning = '("Hyperbolicity may have failed at index (",i4,",",i4,")")'
    character(len=*), parameter :: hyp_info = '("  layer = ",i2,"  kappa = ",d16.8)'

    ! Check for NaNs in the solution
    call check4nans(meqn,mbc,mx,my,q,t,1)

    ! check for h < 0 and reset to zero
    ! check for h < drytolerance
    ! set hu = hv = 0 in all these cells
    forall(i=1-mbc:mx+mbc, j=1-mbc:my+mbc, k=1:num_layers, &
           q(3*(k-1)+1,i,j) / rho(k) < dry_tolerance(k))
        q(3*(k-1)+2:3*(k-1)+3,i,j) = 0.d0
    end forall
!     do i=1-mbc,mx+mbc
!         do j=1-mbc,my+mbc
!             do k=1,num_layers
!                 if (q(3*(k-1)+1,i,j) / rho(k) < dry_tolerance(k)) then
!                     print *,3*(k-1)
!                     q(3*(k-1)+1:3*(k-1)+2,i,j) = 0.d0
!                 endif
!             enddo
!         enddo
!     enddo

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

    ! Check Richardson number -- Only implemented for 2 layers
    if (num_layers == 2 .and. check_richardson) then
        do i=1,mx
            do j=1,my
                dry_state = .false.
                do k=1,num_layers
                    index = 3*(k-1)
                    h(k) = q(index+1,i,j) / rho(k)
                    if (h(k) > dry_tolerance(k)) then
                        u(k) = q(index+2,i,j) / q(index+1,i,j)
                        v(k) = q(index+3,i,j) / q(index+1,i,j)
                    else
                        dry_state(k) = .true.
                        u(k) = 0.d0
                        v(k) = 0.d0
                    endif
                enddo
                ! Calculate for each layer pairing
                do k=1,num_layers-1
                    one_minus_r = 1.d0 - rho(k) / rho(k+1)
                    if (sum(h(k:k+1)) > dry_tolerance(k)) then
                        kappa = (u(k) - u(k+1))**2 / (g*one_minus_r*sum(h(k:k+1)))
                        if (kappa > richardson_tolerance) then
                            write(KAPPA_UNIT,hyp_warning) i,j
                            write(KAPPA_UNIT,hyp_info) k,kappa
                            print hyp_warning, i,j
                            print hyp_info, k,kappa
                        endif
                        kappa = (v(k) - v(k+1))**2 / (g*one_minus_r*sum(h(k:k+1)))
                        if (kappa > richardson_tolerance) then
                            write(KAPPA_UNIT,hyp_warning) i,j
                            write(KAPPA_UNIT,hyp_info) k,kappa
                            print hyp_warning, i,j
                            print hyp_info, k,kappa
                        endif
                    endif
                enddo
            enddo
        enddo
    endif

end subroutine b4step2
