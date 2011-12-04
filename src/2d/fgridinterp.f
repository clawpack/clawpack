c======================================================================
      subroutine fgridinterp(fgrid,xlowfg,ylowfg,
     &       xhifg,yhifg,dxfg,dyfg,mxfg,myfg,t,mvarsfg,q,meqn,
     &       mxc,myc,mbc,dxc,dyc,nvar,xlowc,ylowc,maux,aux,
     &       ioutarrivaltimes,ioutsurfacemax,maxcheck) 
c======================================================================

      use geoclaw_module

      implicit double precision (a-h,o-z)

      dimension q(meqn,1-mbc:mxc+mbc,1-mbc:myc+mbc)
      dimension aux(maux,1-mbc:mxc+mbc,1-mbc:myc+mbc)
      dimension fgrid(mvarsfg,1:mxfg,1:myfg)

c=====================FGRIDINTERP=======================================
c         # This routine interpolates q and aux on a computational grid
c         # to a fgrid not necessarily aligned with the computational grid
c         # using bilinear interpolation defined on computation grid
c=======================================================================

      xhic=xlowc + dxc*mxc  
      yhic=ylowc + dyc*myc    

      tol=drytolerance
      arrivaltol=1.d-2

      indb=meqn+1
      indeta=meqn+2

      if (maxcheck.gt.0) then 
          indarrive=0
          indetamin=0
          indetamax=0
          indetanow=0

        if (ioutarrivaltimes.gt.0) then
           indarrive = 1
        endif

        if (ioutsurfacemax.gt.0) then
          indetanow = indarrive+1
          indetamin = indarrive+2
          indetamax = indarrive+3
        endif
      endif


      do ifg=1,mxfg
         xfg=xlowfg + (ifg-1)*dxfg
      do jfg=1,myfg
         yfg=ylowfg + (jfg-1)*dyfg

         if (.not.((xfg.lt.xlowc.or.xfg.gt.xhic).or.
     &         (yfg.lt.ylowc.or.yfg.gt.yhic))) then
             

c        # find where xfg,yfg is in the computational grid

         ic1 = int((xfg-(xlowc+0.5d0*dxc))/(dxc))+1
         jc1 = int((yfg-(ylowc+0.5d0*dyc))/(dyc))+1
         if (ic1.eq.mxc) ic1=mxc-1
         if (jc1.eq.myc) jc1=myc-1 
         ic2= ic1+1
         jc2= jc1+1

         xc1= xlowc + dxc*(ic1-0.5d0)
         yc1= ylowc + dyc*(jc1-0.5d0)
         xc2= xlowc + dxc*(ic2-0.5d0)
         yc2= ylowc + dyc*(jc2-0.5d0)
 
c        # interpolate bilinear used to interpolate to xfg,yfg

c        # define constant parts of bilinear
         xterm=(xfg-xc1)/dxc
         yterm=(yfg-yc1)/dyc
         xyterm= xterm*yterm

         if (maxcheck.eq.0) then 
         do m=1,meqn
            z11=q(m,ic1,jc1)
            z21=q(m,ic2,jc1)
            z12=q(m,ic1,jc2)
            z22=q(m,ic2,jc2)
            a=z21-z11
            b=z12-z11
            d=z11
            c=z22-(a+b+d)
            fgrid(m,ifg,jfg) = a*xterm + b*yterm + c*xyterm + d

         enddo
         z11=aux(1,ic1,jc1)
         z21=aux(1,ic2,jc1)
         z12=aux(1,ic1,jc2)
         z22=aux(1,ic2,jc2) 
         a=z21-z11
         b=z12-z11
         d=z11
         c=z22-(a+b+d)
         fgrid(indb,ifg,jfg) = a*xterm + b*yterm + c*xyterm + d
         endif
c        # This next output variable is the surface using bilinear interpolation,
c        # using a surface that only uses the wet eta points near the shoreline

         z11=aux(1,ic1,jc1)+q(1,ic1,jc1)
         z21=aux(1,ic2,jc1)+q(1,ic2,jc1)
         z12=aux(1,ic1,jc2)+q(1,ic1,jc2)
         z22=aux(1,ic2,jc2)+q(1,ic2,jc2)

         h11=q(1,ic1,jc1)
         h21=q(1,ic2,jc1)
         h12=q(1,ic1,jc2)
         h22=q(1,ic2,jc2)
         depthindicator= min(h11,h12,h21,h22)
         totaldepth= h11+h22+h21+h12

         if (depthindicator.lt.tol.and.totaldepth.gt.4.d0*tol) then !near shoreline
             if (h11.lt.tol) then
                z11w=  (h12*z12 + h21*z21 + h22*z22)/totaldepth
                z11=z11w
             endif
             if (h12.lt.tol) then
                z12w=  (h11*z11 + h21*z21 + h22*z22)/totaldepth
                z12=z12w
             endif
             if (h21.lt.tol) then
                z21w=  (h11*z11 + h12*z12 + h22*z22)/totaldepth
                z21=z21w
             endif
             if (h22.lt.tol) then
                z22w=  (h12*z12 + h21*z21 + h11*z11)/totaldepth
                z22=z22w
             endif            
         endif
         if (totaldepth.le.4.d0*tol) then
            z22=d_nan()
         endif

         a=z21-z11
         b=z12-z11
         d=z11
         c=z22-(a+b+d)


c        #If eta max/min are saved on this grid initialized if necessary
         if (ioutsurfacemax.gt.0.and.maxcheck.eq.2) then 
            if (.not.(fgrid(indetamin,ifg,jfg).eq.
     &        fgrid(indetamin,ifg,jfg))) fgrid(indetamin,ifg,jfg)=0.d0
            if (.not.(fgrid(indetamax,ifg,jfg).eq.
     &        fgrid(indetamax,ifg,jfg))) fgrid(indetamax,ifg,jfg)=0.d0
             
         endif

c        # check wich task to perform
         if (maxcheck.eq.0) then 
           fgrid(indeta,ifg,jfg) = a*xterm + b*yterm + c*xyterm + d
           fgrid(mvarsfg,ifg,jfg) = t
         elseif (maxcheck.eq.1.and.ioutsurfacemax.gt.0) then
           fgrid(indetanow,ifg,jfg) = a*xterm + b*yterm + c*xyterm + d   
         elseif (maxcheck.eq.2.and.ioutsurfacemax.gt.0) then
            fgrid(indetamin,ifg,jfg) = 
     &          min(fgrid(indetamin,ifg,jfg),fgrid(indetanow,ifg,jfg))
            fgrid(indetamax,ifg,jfg) = 
     &          max(fgrid(indetamax,ifg,jfg),fgrid(indetanow,ifg,jfg))            
         endif

c        #If arrival times are saved on this grid
         if (maxcheck.eq.1.and.ioutarrivaltimes.gt.0) then
            check=fgrid(indarrive,ifg,jfg)
            if (.not.(check.eq.check)) then !# check=NaN: Waves haven't arrived previously 
            if (dabs(fgrid(indeta,ifg,jfg)).gt.arrivaltol) then
                fgrid(indarrive,ifg,jfg)= t
            endif
            endif
         endif



         endif
      enddo
      enddo
      
      return

      end
