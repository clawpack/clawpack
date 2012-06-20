c     ============================================
      subroutine b4step2(maxmx,maxmy,mbc,mx,my,meqn,q,
     &            xlower,ylower,dx,dy,t,dt,maux,aux)
c     ============================================
c
c     # called before each call to step
c     # use to set time-dependent aux arrays or perform other tasks.
c
c     This particular routine sets negative values of q(i,j,1) to zero,
c     as well as the corresponding q(i,j,m) for m=1,meqn.
c     This is for problems where q(i,j,1) is a depth.
c     This should occur only because of rounding error.

c     Also calls movetopo if topography might be moving.

      use geoclaw_module
      use topo_module
      use dtopo_module

      implicit double precision (a-h,o-z)

      dimension q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
      dimension aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
      
      integer :: layer,layer_index
      double precision :: h(2),u(2),v(2),g,kappa
      logical :: dry_state(2)

c=====================Parameters===========================================


c     # check for NANs in solution:
      call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,1)

c     # check for h < 0 and reset to zero
c     # check for h < drytolerance
c     # set hu = hv = 0 in all these cells

      do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
            do k=1,layers
                if (q(1,i,j) / rho(k) < drytolerance) then
                    q(2:meqn,i,j) = 0.d0
                endif
            enddo
        enddo
      enddo

      write(26,*) 'B4STEP2: t, num_dtopo: ', t,num_dtopo
      do i=1,num_dtopo
          call movetopo(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,t,dt,maux,aux,
     &      dtopowork(i0dtopo(i):i0dtopo(i)+mdtopo(i)-1),
     &      xlowdtopo(i),ylowdtopo(i),xhidtopo(i),yhidtopo(i),
     &      t0dtopo(i),tfdtopo(i),dxdtopo(i),dydtopo(i),dtdtopo(i),
     &      mxdtopo(i),mydtopo(i),mtdtopo(i),mdtopo(i),
     &      minleveldtopo(i),maxleveldtopo(i),topoaltered(i))
      enddo

      ! Check Richardson number -- Only implemented for 2 layers
      if (layers == 2 .and. check_richardson) then
      do i=1,mx
          do j=1,my
              dry_state = .false.
              do k=1,layers
                  index = 3*(k-1)
                  h(k) = q(index+1,i,j) / rho(k)
                  if (h(k) > drytolerance) then
                      u(k) = q(index+2,i,j) / q(index+1,i,j)
                      v(k) = q(index+3,i,j) / q(index+1,i,j)
                  else
                      dry_state(k) = .true.
                      u(k) = 0.d0
                      v(k) = 0.d0
                  endif
              enddo
              
              ! Calculate for each layer pairing
              do k=1,layers-1
                  if (sum(h(k:k+1)) > drytolerance) then
                      one_minus_r = rho(k) / rho(k+1)
                      kappa = (u(k) - u(k+1))**2
     &                          / (g*one_minus_r*sum(h(k:k+1)))
                      if (kappa > richardson_tolerance) then
                          write(kappa_file,100) i,j
                          write(kappa_file,101) k,kappa
                          print 100, i,j
                          print 101, k,kappa
                      endif
                      kappa = (v(k) - v(k+1))**2
     &                          / (g*one_minus_r*sum(h(k:k+1)))
                      if (kappa > richardson_tolerance) then
                          write(kappa_file,100) i,j
                          write(kappa_file,101) k,kappa
                          print 100, i,j
                          print 101, k,kappa
                      endif
                  endif
              enddo
          enddo
      enddo
      endif
      
100   format ("Hyperbolicity may have failed at index (",i4,",",i4,")")
101   format("  layer = ",i2,"  kappa = ",d16.8)

      return
      end
