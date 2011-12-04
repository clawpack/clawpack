c======================================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  fwave,s,amdq,apdq)
c======================================================================
c
c Solves normal Riemann problems for the 2D SHALLOW WATER equations
c     with topography:
c     #        h_t + (hu)_x + (hv)_y = 0                           #
c     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
c     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #

c On input, ql contains the state vector at the left edge of each cell
c     qr contains the state vector at the right edge of each cell
c
c This data is along a slice in the x-direction if ixy=1
c     or the y-direction if ixy=2.

c  Note that the i'th Riemann problem has left state qr(i-1,:)
c     and right state ql(i,:)
c  From the basic clawpack routines, this routine is called with
c     ql = qr
c
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           !
!      # This Riemann solver is for the shallow water equations.            !
!                                                                           !
!       It allows the user to easily select a Riemann solver in             !
!       riemannsolvers_geo.f. this routine initializes all the variables    !
!       for the shallow water equations, accounting for wet dry boundary    !
!       dry cells, wave speeds etc.                                         !
!                                                                           !
!           David George, Vancouver WA, Feb. 2009                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*     ! this version rpn2_ez_fast.f is faster by using a simpler fwave
*        solver in the deep ocean.

      use geoclaw_module

      implicit none

      !input
      integer maxm,meqn,mwaves,mbc,mx,ixy

      double precision  fwave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision  s(1-mbc:maxm+mbc, mwaves)
      double precision  ql(1-mbc:maxm+mbc, meqn)
      double precision  qr(1-mbc:maxm+mbc, meqn)
      double precision  apdq(1-mbc:maxm+mbc, meqn)
      double precision  amdq(1-mbc:maxm+mbc, meqn)
      double precision  auxl(1-mbc:maxm+mbc, *)
      double precision  auxr(1-mbc:maxm+mbc, *)

      double precision drytol,g

      !local only
      integer m,i,mw,maxiter,mu,nv
      double precision wall(3)
      double precision fw(3,3)
      double precision sw(3)

      double precision hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL
      double precision bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      double precision s1m,s2m,s1,s2
      double precision hstar,hstartest,hstarHLL,sLtest,sRtest
      double precision tw,dxdc

      logical rare1,rare2
      integer mcapa

      common /cmcapa/  mcapa

      g=grav
      drytol=drytolerance

      !loop through Riemann problems at each grid cell
      do i=2-mbc,mx+mbc

!-----------------------Initializing-----------------------------------
         !inform of a bad riemann problem from the start
         if((qr(i-1,1).lt.0.d0).or.(ql(i,1) .lt. 0.d0)) then
            write(*,*) 'Negative input: hl,hr,i=',qr(i-1,1),ql(i,1),i
         endif

         !Initialize Riemann problem for grid interface
         do mw=1,mwaves
              s(i,mw)=0.d0
              do m=1,meqn
                 fwave(i,m,mw)=0.d0
              enddo
         enddo

c        !set normal direction
         if (ixy.eq.1) then
            mu=2
            nv=3
         else
            mu=3
            nv=2
         endif

         !zero (small) negative values if they exist
         if (qr(i-1,1).lt.0.d0) then
            do m=1,meqn
               qr(i-1,m)=0.d0
            enddo
         endif

         if (ql(i,1).lt.0.d0) then
            do m=1,meqn
               ql(i,1)=0.d0
            enddo
         endif

         !skip problem if in a completely dry area
         if (qr(i-1,1).le.drytol.and.ql(i,1).le.drytol) then
            go to 30
         endif

         !Riemann problem variables
         hL = qr(i-1,1)
         hR = ql(i,1)
         huL = qr(i-1,mu)
         huR = ql(i,mu)
         bL = auxr(i-1,1)
         bR = auxl(i,1)

         hvL=qr(i-1,nv)
         hvR=ql(i,nv)

         !check for wet/dry boundary
         if (hR.gt.drytol) then
            uR=huR/hR
            vR=hvR/hR
            phiR = 0.5d0*g*hR**2 + huR**2/hR
         else
            hR = 0.d0
            huR = 0.d0
            hvR = 0.d0
            uR = 0.d0
            vR = 0.d0
            phiR = 0.d0
         endif

         if (hL.gt.drytol) then
            uL=huL/hL
            vL=hvL/hL
            phiL = 0.5d0*g*hL**2 + huL**2/hL
         else
            hL=0.d0
            huL=0.d0
            hvL=0.d0
            uL=0.d0
            vL=0.d0
            phiL = 0.d0
         endif

!        check if a simple solver will do the trick
         if ((abs(hL-hR).lt.0.001d0*min(hL,hR)).and.
     &            (max(uL**2,uR**2)/max(g*hL,g*hR).lt.0.01d0)) then
            s1 = 0.5d0*(uR+uL) - sqrt(g*0.5d0*(hL+hR))
            s2 = 0.5d0*(uR+uL) + sqrt(g*0.5d0*(hL+hR))
            call riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,
     &         bL,bR,uL,uR,vL,vR,phiL,phiR,s1,s2,drytol,g,sw,fw)
            go to 25
         endif

         wall(1) = 1.d0
         wall(2) = 1.d0
         wall(3) = 1.d0
         if (hR.le.drytol) then
            call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,
     &                                  rare1,rare2,1,drytol,g)
            hstartest=max(hL,hstar)
            if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
c                bR=hstartest+bL
               wall(2)=0.d0
               wall(3)=0.d0
               hR=hL
               huR=-huL
               bR=bL
               phiR=phiL
               uR=-uL
               vR=vL
            elseif (hL+bL.lt.bR) then
               bR=hL+bL
            endif
         elseif (hL.le.drytol) then ! right surface is lower than left topo
            call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,
     &                                  rare1,rare2,1,drytol,g)
            hstartest=max(hR,hstar)
            if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
c               bL=hstartest+bR
               wall(1)=0.d0
               wall(2)=0.d0
               hL=hR
               huL=-huR
               bL=bR
               phiL=phiR
               uL=-uR
               vL=vR
            elseif (hR+bR.lt.bL) then
               bL=hR+bR
            endif
         endif

         !determine wave speeds
         sL=uL-sqrt(g*hL) ! 1 wave speed of left state
         sR=uR+sqrt(g*hR) ! 2 wave speed of right state

         uhat=(sqrt(g*hL)*uL + sqrt(g*hR)*uR)/(sqrt(g*hR)+sqrt(g*hL)) ! Roe average
         chat=sqrt(g*0.5d0*(hR+hL)) ! Roe average
         sRoe1=uhat-chat ! Roe wave speed 1 wave
         sRoe2=uhat+chat ! Roe wave speed 2 wave

         sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
         sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

         !--------------------end initializing...finally----------
         !solve Riemann problem.

         maxiter = 1

         call riemann_aug_JCP(maxiter,3,3,hL,hR,huL,
     &        huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,
     &                                    drytol,g,sw,fw)

c         call riemann_ssqfwave(maxiter,meqn,mwaves,hL,hR,huL,huR,
c     &     hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

c          call riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,
c     &      bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

c        !eliminate ghost fluxes for wall
         do mw=1,3
            sw(mw)=sw(mw)*wall(mw)
            do m=1,meqn
               fw(m,mw)=fw(m,mw)*wall(mw)
            enddo
         enddo

 25      continue

         do mw=1,mwaves
            s(i,mw)=sw(mw)
            fwave(i,1,mw)=fw(1,mw)
            fwave(i,mu,mw)=fw(2,mw)
            fwave(i,nv,mw)=fw(3,mw)
         enddo

 30      continue
      enddo


c==========Capacity for mapping from latitude longitude to physical space====

        if (mcapa.gt.0) then
         do i=2-mbc,mx+mbc
          if (ixy.eq.1) then
            dxdc=(Rearth*pi/180.d0)
          else
	    dxdc=auxl(i,3)
          endif

          do mw=1,mwaves
c             if (s(i,mw) .gt. 316.d0) then
c               # shouldn't happen unless h > 10 km!
c                write(6,*) 'speed > 316: i,mw,s(i,mw): ',i,mw,s(i,mw)
c                endif
	     s(i,mw)=dxdc*s(i,mw)
             do m=1,meqn
               fwave(i,m,mw)=dxdc*fwave(i,m,mw)
             enddo
          enddo
         enddo
        endif

c===============================================================================


c============= compute fluctuations=============================================
         do i=1-mbc,mx+mbc
            do m=1,meqn
               amdq(i,m)=0.0d0
               apdq(i,m)=0.0d0
               do  mw=1,mwaves
                  if (s(i,mw).lt.0.d0) then
                     amdq(i,m)=amdq(i,m) + fwave(i,m,mw)
                  elseif (s(i,mw).gt.0.d0) then
                     apdq(i,m)=apdq(i,m) + fwave(i,m,mw)
                  else
		            amdq(i,m) = amdq(i,m) + .5d0*fwave(i,m,mw)
		            apdq(i,m) = apdq(i,m) + .5d0*fwave(i,m,mw)
                  endif
               enddo
            enddo
         enddo


      return
      end subroutine









