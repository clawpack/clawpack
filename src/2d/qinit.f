

c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c      # Set initial sea level flat unless iqinit = 1, in which case
c      # an initial perturbation of the q(i,j,1) is specified and has
c      # been strored in qinitwork.


       use geoclaw_module

       implicit double precision (a-h,o-z)
       dimension q(meqn,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
       dimension aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)

       include 'qinit.i'


       do i=1-mbc,mx+mbc
          x = xlower + (i-0.5d0)*dx
          do j=1-mbc,my+mbc
             y = ylower + (j-0.5d0)*dy
             q(1,i,j)=dmax1(0.d0,sealevel-aux(1,i,j))
             q(2,i,j)=0.d0
             q(3,i,j)=0.d0
             enddo
          enddo

       if (iqinit .gt. 0) then

         do i=1-mbc,mx+mbc
           x = xlower + (i-0.5d0)*dx
           xim = x - 0.5d0*dx
           xip = x + 0.5d0*dx
           do j=1-mbc,my+mbc
             y = ylower + (j-0.5d0)*dy
             yjm = y - 0.5d0*dy
             yjp = y + 0.5d0*dy


             if (xip.gt.xlowqinit.and.xim.lt.xhiqinit
     &          .and.yjp.gt.ylowqinit.and.yjm.lt.yhiqinit) then

                   xipc=min(xip,xhiqinit)
                   ximc=max(xim,xlowqinit)
                   xc=0.5d0*(xipc+ximc)

                   yjpc=min(yjp,yhiqinit)
                   yjmc=max(yjm,ylowqinit)
                   yc=0.5d0*(yjmc+yjpc)

                   dq = topointegral(ximc,xc,xipc,yjmc,yc,yjpc,
     &                    xlowqinit,ylowqinit,dxqinit,dyqinit,
     &                    mxqinit,myqinit,qinitwork,1)
                   dq=dq/((xipc-ximc)*(yjpc-yjmc)*aux(2,i,j))

                   if (iqinit.lt.4) then 
                      if (aux(1,i,j).le.0.d0) then
                          q(iqinit,i,j) = q(iqinit,i,j) + dq
                      endif
                   elseif (iqinit.eq.4) then
                      q(1,i,j) = max(dq-aux(1,i,j),0.d0)
                   endif
             endif
            enddo
           enddo
         endif

       return
       end
