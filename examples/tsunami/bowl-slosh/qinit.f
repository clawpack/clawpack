

c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c      # Set initial sea level flat unless iqinit = 1, in which case
c      # an initial perturbation of the q(i,j,1) is specified and has
c      # been strored in qinitwork.


       use geoclaw_module
       use qinit_module

       implicit double precision (a-h,o-z)
       dimension q(meqn, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
       dimension aux(maux, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)


c      grav = 9.81d0   !# should come from module
       a = 1.d0
       sigma = 0.5d0
       h0 = 0.1d0
       omega = dsqrt(2.d0*grav*h0) / a
       write(6,*) 'omega = ',omega

       do i=1-mbc,mx+mbc
          x = xlower + (i-0.5d0)*dx
          do j=1-mbc,my+mbc
             y = ylower + (j-0.5d0)*dy
             eta = sigma*h0/a**2 * (2.*x - sigma) 
             q(1,i,j) = dmax1(0.d0,eta-aux(1,i,j))
             q(2,i,j) = 0.d0
             q(3,i,j) = sigma*omega * q(1,i,j)
             enddo
          enddo


       return
       end
