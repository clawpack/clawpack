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

      use multilayer_module
      use hurricane_module
      use geoclaw_module
      use topo_module
      use dtopo_module

      implicit double precision (a-h,o-z)

      dimension q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
      dimension vel(2,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
      double precision :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
      
      integer :: layer,layer_index
      double precision :: h(2),u(2),v(2),g,kappa
      logical :: dry_state(2)
      
      double precision, pointer, dimension(:,:,:) :: wind
      double precision, pointer, dimension(:,:) :: pressure
      
      g = grav
      
c=====================Parameters===========================================


c     # check for NANs in solution:
      call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,1)

c     # check for h < 0 and reset to zero
c     # check for h < drytolerance
c     # set hu = hv = 0 in all these cells

      do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
            do m=1,layers
              index = 3*(m-1)
              if (q(index+1,i,j) / rho(m) < drytolerance) then
                 q(index+1,i,j) = max(q(index+1,i,j),0.d0)
                 q(index+2,i,j) = 0.d0
                 q(index+3,i,j) = 0.d0
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
    
      ! Set wind and pressure aux variables for this grid
      write(26,*) "B4STEP2:  Setting aux array for wind and pressure"
      call hurricane_wind(maxmx,maxmy,maux,mbc,mx,my,xlower,ylower,dx,
     &                      dy,t,aux)
      call hurricane_pressure(maxmx,maxmy,maux,mbc,mx,my,xlower,ylower,
     &                         dx,dy,t,aux)
     
      ! Check Richardson number
      if (layers > 1) then
      do i=1,mx
          do j=1,my
              dry_state = .false.
              do layer=1,2
                  m = 3*(layer-1)
                  h(layer) = q(m+1,i,j)
                  if (h(layer) > drytolerance) then
                      u(layer) = q(m+2,i,j) / q(m+1,i,j)
                      v(layer) = q(m+3,i,j)/ q(m+1,i,j)
                  else
                      dry_state(layer) = .true.
                      u(layer) = 0.d0
                      v(layer) = 0.d0
                  endif
              enddo
              if (sum(h) > drytolerance) then
                  kappa=(u(1) - u(2))**2 / (g*one_minus_r*sum(h))
                  if ((kappa > richardson_tolerance)
     &                  .and.(.not.dry_state(2))) then
                      write(kappa_file,100) i,j,kappa
                      print 100,i,j,kappa
                  endif
                  aux(10,i,j)=(v(1) - v(2))**2 / (g*one_minus_r*sum(h))
                  if ((kappa > richardson_tolerance)
     &                  .and.(.not.dry_state(2))) then
                      write(kappa_file,100),i,j,kappa
                      print 100,i,j,kappa
                  endif
               endif
          enddo
      enddo
      endif
      
100   format ("Hyperbolicity may have failed (",i4,",",i4,") = ",d16.8)
       
       ! Calculate vorticity
C        aux(9,:,:) = 0.d0
C        vel = 0.d0
C        do i=1,mx
C            do j=1,my
C                if (abs(q(1,i,j)) > drytolerance) then
C                    vel(1,i,j) = q(2,i,j) / q(1,i,j)
C                    vel(2,i,j) = q(3,i,j) / q(1,i,j)
C                endif
C            enddo
C        enddo
C        do i=1,mx
C            do j=1,my
C                aux(9,i,j) = (vel(2,i+1,j) - vel(2,i-1,j)) / (2.d0*dx) -         
C       &                     (vel(1,i,j+1) - vel(1,i,j-1)) / (2.d0*dy)
C            enddo
C        enddo
      
      return
      end
    