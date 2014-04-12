c     ============================================
      subroutine movetopo(mbc,mx,my,
     &     xlow,ylow,dx,dy,t,dt,maux,aux,dtopo,
     &     xlowdtopo,ylowdtopo,xhidtopo,yhidtopo,t0dtopo,tfdtopo,
     &     dxdtopo,dydtopo,dtdtopo,mxdtopo,mydtopo,mtdtopo,mdtopo,
     &     minleveldtopo,maxleveldtopo,topoaltered)
c     ============================================
c
c     This routine changes the topography dynamically.
c     It also changes the actual topo array after the motion has ended.
c     Topography is stored in aux(i,j,1).
c
c
      use geoclaw_module
      use amr_module

      use topo_module, only: topowork
      use topo_module, only: mtopofiles,mtoposize
      use topo_module, only: xlowtopo, ylowtopo, tlowtopo
      use topo_module, only: xhitopo, yhitopo, thitopo
      use topo_module, only: dxtopo, dytopo
      use topo_module, only: topotime
      use topo_module, only: mxtopo, mytopo

      use topo_module, only: i0topo, mtopo, mtopoorder
      use topo_module, only: topoID


      implicit double precision (a-h,o-z)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension dtopo(mxdtopo,mydtopo,mtdtopo)

      logical topoaltered

      dimension auxorig(maux,-1:mx+mbc,-1:my+mbc)

        t0=t  !# start of coming timestep
        tf=t+dt !# end of coming timestep
        t2=t0+0.5d0*dt
        tt=t2 !# this is the time in the fault file that will be used for this timestep

c       write(26,*) 'MOVETOPO: tt, dt, t0dtopo: ',tt,dt,t0dtopo

        if (topoaltered.and.t0.le.tfdtopo) then
           write(*,*) 'MOVETOPO: WARNING!'
           write(*,*) 'FAULT MOTION HAS COMPLETED PREMATURELY'
        endif
        if (topoaltered) return
        if (tf.lt.t0dtopo) return
        if (t0.gt.tfdtopo) go to 28

c     # change topography

        write(*,*) 'MOVETOPO: setting dtopo at time = ',t
c       write(26,*) 'MOVETOPO: setting dtopo at time = ',t

        if (tt.ge.tfdtopo) then
              kkm = mtdtopo
              kkp = mtdtopo
              tau = 1.d0
              go to 14
        endif

        if (tt.le.t0dtopo) then
              kkm = 1
              kkp = 1
              tau= 1.d0
              go to 14
        endif

        dkk = (tt-t0dtopo)/dtdtopo
        kkm = dint(dkk)+1
        kkp = dint(dkk)+2

        kkm = max(kkm,1)
        kkm = min(kkm,mtdtopo)

        kkp = max(kkp,1)
        kkp = min(kkp,mtdtopo)

        tdtopom = t0dtopo+dtdtopo*(kkm-1)
        tdtopop = tdtopom+dtdtopo
        tau = 1.d0-max(0.d0,((tt-tdtopom)/dtdtopo))
        tau = max(tau,0.d0)

 14     continue


c     # recreate original topography:
      auxorig(1,:,:) = NEEDS_TO_BE_SET
      call setaux(mbc,mx,my,xlow,ylow,dx,dy,
     &                  maux,auxorig)

c=======loop through the computational grid row by row====================
        do j=1-mbc,my+mbc
            ycell = ylow + (j-0.5d0)*dy
            yjm =   ylow + (j-1.d0)*dy
            yjp =   ylow + dble(j)*dy

c           if (yjm.ge.ylowdtopo.and.yjp.le.yhidtopo) then
           if (yjp.gt.ylowdtopo.and.yjm.lt.yhidtopo) then
              yjpc=min(yjp,yhidtopo)
              yjmc=max(yjm,ylowdtopo)
              dyc=yjpc-yjmc
              ycellc=0.5d0*(yjpc+yjmc)
c=========sweep through the columns of the computational row================
           do i=1-mbc,mx+mbc
              xcell =  xlow + (i- 0.5d0)*dx
              xim =    xlow + (i- 1.0d0)*dx
              xip =   xlow + dble(i)*dx


c              if (xim.ge.xlowdtopo.and.xip.le.xhidtopo) then
              if (xip.gt.xlowdtopo.and.xim.lt.xhidtopo) then
                 xipc=min(xip,xhidtopo)
                 ximc=max(xim,xlowdtopo)
                 dxc=xipc-ximc
                 xcellc=0.5d0*(xipc+ximc)
c==========alter aux(i,j,1) if it is in the dtopo region=====================

                 aux(1,i,j) = auxorig(1,i,j)

c==================dynamically alter the topography to the interpolated value=
                 dauxijm=0.d0
                 dauxijm = topointegral(ximc,xcellc,xipc,yjmc,
     &              ycellc,yjpc,xlowdtopo,ylowdtopo,dxdtopo,dydtopo,
     &              mxdtopo,mydtopo,dtopo(1,1,kkm),1)
                 dauxijm = dauxijm/(dxc*dyc*aux(2,i,j))

                 dauxijp=0.d0
                 dauxijp = topointegral(ximc,xcellc,xipc,yjmc,
     &              ycellc,yjpc,xlowdtopo,ylowdtopo,dxdtopo,dydtopo,
     &              mxdtopo,mydtopo,dtopo(1,1,kkp),1)
                 dauxijp = dauxijp/(dxc*dyc*aux(2,i,j))

                 aux(1,i,j)=aux(1,i,j)+tau*dauxijm+(1.d0-tau)*dauxijp

              endif
           enddo
           endif
         enddo

c============= end of altering topography dynamically========================


c        # if the aux cell is outside the physical domain,
c        # boundary conditions must be reset

         if (xlowdtopo.lt.xlower.and.
     &      xlow+(0.5d0-mbc)*dx.lt.xlower) then
            do j=1-mbc,my+mbc
               ycell = ylow +(j-0.5d0)*dy
               yjm = ylow  + (j-1.0d0)*dy
               yjp = ylow + dble(j)*dy

               xcell =  xlow + (1- 0.5d0)*dx
               xim =    xlow
               xip =    xlow + dx
               if ((xim.ge.xlowdtopo.and.xip.le.xhidtopo).and.
     &                     (yjm.ge.ylowdtopo.and.yjp.le.yhidtopo)) then

                   do i=1-mbc,0
                      aux(1,i,j) = aux(1,1,j)
                   enddo
               endif
            enddo
          endif

          if (xhidtopo.gt.xupper.and.
     &      xlow+(mx+mbc-0.5d0)*dx.gt.xupper) then
            do j=1-mbc,my+mbc
               ycell = ylow+(j-0.5d0)*dy
               yjm = ylow+(j-1.0d0)*dy
               yjp = ylow + dble(j)*dy

               xcell =  xlow + (mx-0.5d0)*dx
               xim =  xlow + (mx-1.d0)*dx
               xip =  xlow + dble(mx)*dx
               if ((xim.ge.xlowdtopo.and.xip.le.xhidtopo).and.
     &                     (yjm.ge.ylowdtopo.and.yjp.le.yhidtopo)) then

                   do i=mx+1,mx+mbc
                      aux(1,i,j) = aux(1,mx,j)
                   enddo
               endif
            enddo
          endif

          if (ylowdtopo.lt.ylower.and.
     &             ylow+(.5d0-mbc)*dy.lt.ylower) then
                  do i=1-mbc,mx+mbc
               ycell = ylow+(1-0.5d0)*dy
               yjm = ylow
               yjp = ylow + dy

               xcell =  xlow + (i- 0.5d0)*dx
               xim =  xlow + (i- 1.d0)*dx
               xip = xlow + dble(i)*dx

               if ((xim.ge.xlowdtopo.and.xip.le.xhidtopo).and.
     &                     (yjm.ge.ylowdtopo.and.yjp.le.yhidtopo)) then

                   do j=1-mbc,0
                      aux(1,i,j) = aux(1,i,1)
                   enddo
               endif
            enddo
          endif

          if (yhidtopo.gt.yupper.and.
     &             ylow+(my+mbc-0.5d0)*dy.gt.yupper) then
            do i=1-mbc,mx+mbc
               ycell = ylow+(my-0.5d0)*dy
               yjm = ylow+(my-1.d0)*dy
               yjp = ylow+ dble(my)*dy

               xcell =  xlow + (i- 0.5d0)*dx
               xim =   xlow + (i- 1.d0)*dx
               xip =  xlow + dble(i)*dx

               if ((xim.ge.xlowdtopo.and.xip.le.xhidtopo).and.
     &                      (yjm.ge.ylowdtopo.and.yjp.le.yhidtopo)) then

                   do j=my+1,my+mbc
                      aux(1,i,j) = aux(1,i,my)
                   enddo
               endif
            enddo
          endif
c===========================B.C.'s reset========================================



c============= statically change topowork array for next call to setaux====
 28     continue
        if (tt-dt.gt.tfdtopo) then
c       if (.false.) then
          topoaltered=.true.  !# so this procedure only done once
          write(*,*) 'MOVETOPO: Altering topography arrays at t=',t
          write(*,*) 'MOVETOPO: Topography Motion Complete'
          do m=1,mtopofiles
            if (xlowtopo(m).lt.xhidtopo.and.xhitopo(m).gt.xlowdtopo
     &                .and.ylowtopo(m).lt.yhidtopo.and.
     &                yhitopo(m).gt.ylowdtopo) then

            do ib=1,mxtopo(m)
              xcell=xlowtopo(m) + (ib-.5d0)*dxtopo(m)
              xim=xlowtopo(m) + (ib-1.d0)*dxtopo(m)
              xip=xlowtopo(m) + ib*dxtopo(m)

             if (xim.lt.xhidtopo.and.xip.gt.xlowdtopo) then
                ximc=max(xim,xlowdtopo)
                xipc=min(xip,xhidtopo)
                dxc=xipc-ximc
                xcellc=0.5d0*(ximc+xipc)

               do jb=1,mytopo(m)
                     ycell=ylowtopo(m) + (jb-.5d0)*dytopo(m)
                     yjm=ylowtopo(m) + (jb-1.d0)*dytopo(m)
                     yjp = ylowtopo(m) + jb*dytopo(m)

                if (yjm.lt.yhidtopo.and.yjp.gt.ylowdtopo) then !orig changed to match 4-x
C                 if (yjm.lt.yhidtopo.and.yjp.gt.ylowdtopo) then
                    yjmc=max(yjm,ylowdtopo)
                    yjpc=min(yjp,yhidtopo)
                    dyc = yjpc-yjmc
                    ycellc= 0.5d0*(yjmc+yjpc)

                    ztopoij = 0.d0
                    ztopoij = topointegral(ximc,xcellc,xipc,yjmc,ycellc,
     &                   yjpc,xlowdtopo,ylowdtopo,dxdtopo,dydtopo,
     &                   mxdtopo,mydtopo,dtopo(1,1,mtdtopo),1)

c                  # topointegral actually integrates,
c                  #so need to divide by area
c                  #physical area of cell depends on coordinate system
                   ztopoij=ztopoij/(dxc*dyc)
                   if (coordinate_system.eq.2) then
                     capac_area = deg2rad*earth_radius**2*
     &               (sin(yjp*deg2rad)-sin(yjm*deg2rad))/dyc
                     ztopoij=ztopoij/capac_area
                   endif

                   jbr = mytopo(m)-jb+1
                   ij = i0topo(m) + (jbr-1)*mxtopo(m) +ib -1
                   topowork(ij) = topowork(ij)+ztopoij
                endif
               enddo
            endif
           enddo
          endif
          enddo
        endif
c=================================================================================


      return
      end
