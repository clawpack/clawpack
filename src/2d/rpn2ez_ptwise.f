
c
c
c     =====================================================================
      subroutine rpn2_ptwise(ixy, meqn, mwaves, ldw, maux, 
     &           qLeft, qRight, auxLeft, auxRight, 
     &           wave, s, amdq, apdq, rpfwave)
c     =====================================================================
c
c Solves normal Riemann problems for the 2D SHALLOW WATER equations
c     with topography:
c     #        h_t + (hu)_x + (hv)_y = 0                           #
c     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
c     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #
c     
c     # This is a new-style Riemann solver that works pointwise, with
c     # Left state qLeft and right state qRight.  (Note change in notation!!)
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

      use geoclaw_module

      implicit none
      integer meqn,mwaves,ixy,maux,ldw
      double precision qLeft(meqn), qRight(meqn) 
      double precision auxLeft(maux), auxRight(maux) 
      double precision wave(ldw, mwaves) 
      double precision s(mwaves)
      double precision amdq(meqn), apdq(meqn)
      logical rpfwave

      double precision drytol,g

      !local only
      integer m,i,mw,maxiter,mu,nv
      double precision wall(3)
      double precision fw(3,3)
      double precision sw(3)

      double precision hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL
      double precision bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      double precision s1m,s2m
      double precision hstar,hstartest,hstarHLL,sLtest,sRtest
      double precision tw,dxdc

      logical rare1,rare2
      integer mcapa

      common /cmcapa/  mcapa

      g=grav
      drytol=drytolerance
c
c     # This solver returns fwaves:
      rpfwave = .true.


!-----------------------Initializing-----------------------------------
         !inform of a bad riemann problem from the start
         if((qLeft(1).lt.0.d0).or.(qRight(1) .lt. 0.d0)) then
            write(*,*) 'Negative input: hl,hr=',qLeft(1),qRight(1)
         endif

         !Initialize Riemann problem for grid interface
         do mw=1,mwaves
              s(mw)=0.d0
              do m=1,meqn
                 wave(m,mw)=0.d0
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
         if (qLeft(1).lt.0.d0) then
            do m=1,meqn
               qLeft(m)=0.d0
            enddo
         endif

         if (qRight(1).lt.0.d0) then
            do m=1,meqn
               qRight(1)=0.d0
            enddo
         endif

         !skip problem if in a completely dry area
         if (qLeft(1).le.drytol.and.qRight(1).le.drytol) then
            go to 30
         endif

         !Riemann problem variables
         hL = qLeft(1)
         hR = qRight(1)
         huL = qLeft(mu)
         huR = qRight(mu)
         bL = auxLeft(1)
         bR = auxRight(1)

         hvL=qLeft(nv)
         hvR=qRight(nv)

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

c         call riemann_ssqwave(maxiter,meqn,mwaves,hL,hR,huL,huR,
c     &     hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

c          call riemann_wave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,
c     &      bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

c        !eliminate ghost fluxes for wall
         do mw=1,3
            sw(mw)=sw(mw)*wall(mw)
            do m=1,meqn
               fw(m,mw)=fw(m,mw)*wall(mw)
            enddo
         enddo

         do mw=1,mwaves
            s(mw)=sw(mw)
            wave(1,mw)=fw(1,mw)
            wave(mu,mw)=fw(2,mw)
            wave(nv,mw)=fw(3,mw)
         enddo

 30      continue


c==========Capacity for mapping from latitude longitude to physical space====

        if (mcapa.gt.0) then
          if (ixy.eq.1) then
             dxdc=(Rearth*pi/180.d0)
          else
	     dxdc=auxRight(3)
          endif

          do mw=1,mwaves
c             if (s(mw) .gt. 316.d0) then
c               # shouldn't happen unless h > 10 km!
c                write(6,*) 'speed > 316: i,mw,s(mw): ',i,mw,s(mw)
c                endif
	           s(mw)=dxdc*s(mw)
             do m=1,meqn
               wave(m,mw)=dxdc*wave(m,mw)
             enddo
          enddo
        endif

c===============================================================================


c============= compute fluctuations=============================================
            do m=1,meqn
               amdq(m)=0.0d0
               apdq(m)=0.0d0
               do  mw=1,mwaves
                  if (s(mw).lt.0.d0) then
                     amdq(m)=amdq(m) + wave(m,mw)
                  elseif (s(mw).gt.0.d0) then
                     apdq(m)=apdq(m) + wave(m,mw)
                  else
   	             amdq(m) = amdq(m) + .5d0*wave(m,mw)
 	             apdq(m) = apdq(m) + .5d0*wave(m,mw)
                  endif
               enddo
            enddo


      return
      end subroutine

