c
c ---------------------------------------------------------------
c
        recursive subroutine filrecur(level,nvar,valbig,aux,naux,
     1                      time,mitot,mjtot,
     2                      nrowst,ncolst,ilo,ihi,jlo,jhi)

c :::::::::::::::::::::::::::: FILPATCH :::::::::::::::::::::::::;
c
c  fill the portion of valbig from rows  nrowst
c                             and  cols  ncolst
c  the patch can also be described by the corners (xlp,ybp) by (xrp,ytp).
c  vals are needed at time time, and level level,
c
c  first fill with  values obtainable from the level level
c  grids. if any left unfilled, then enlarge remaining rectangle of
c  unfilled values by 1 (for later linear interp), and recusively
c  obtain the remaining values from  coarser levels.
c
c :::::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::;

      use geoclaw_module

      implicit double precision (a-h,o-z)

      include  "call.i"

      logical   set, sticksout
      dimension valbig(nvar,mitot,mjtot), aux(naux,mitot,mjtot)

c     use stack-based scratch arrays instead of alloc, since dont really
c     need to save beyond these routines, and to allow dynamic memory resizing
c
c     use 1d scratch arrays that are potentially the same size as
c     current grid, since may not coarsen.
c     need to make it 1d instead of 2 and do own indexing, since
c     when pass it in to subroutines they treat it as having different
c     dimensions than the max size need to allocate here
c
      dimension valcrse((ihi-ilo+2)*(jhi-jlo+2)*nvar)  ! NB this is a 1D array
      dimension auxcrse((ihi-ilo+2)*(jhi-jlo+2)*naux)  ! the +2 is to expand on coarse grid to enclose fine
c
      dimension flaguse(ihi-ilo+1,jhi-jlo+1)

      logical reloop
      logical fineflag((ihi-ilo+2)*(jhi-jlo+2)*nvar)
      double precision finemass((ihi-ilo+2)*(jhi-jlo+2))
      double precision etacrse((ihi-ilo+2)*(jhi-jlo+2))
      double precision velmax((ihi-ilo+2)*(jhi-jlo+2))
      double precision velmin((ihi-ilo+2)*(jhi-jlo+2))
      double precision slopex((ihi-ilo+2)*(jhi-jlo+2))
      double precision slopey((ihi-ilo+2)*(jhi-jlo+2))
      integer icount((ihi-ilo+2)*(jhi-jlo+2))

c OLD INDEXING
c$$$      ivalc(i,j,ivar) = i + nrowc*(j - 1)
c$$$     &                    + nrowc*ncolc*(ivar-1)
c$$$      icrse(i,j) = i + nrowc*(j-1)
c$$$c
c$$$c     # index into first component of aux = topo:
c$$$      iauxc(i,j) =  i + nrowc*(j-1)

c NEW INDEXING - ORDER SWITCHED
      ivalc(ivar,i,j) = ivar + nvar*(i-1) + nvar*nrowc*(j-1)
      icrse(i,j) = i + nrowc*(j-1)
c
c     # index into first component of aux = topo:
      iauxc(i,j) = 1 + naux*(i-1) + naux*nrowc*(j-1)
c
      sticksout(iplo,iphi,jplo,jphi)  =
     &            (iplo .lt. 0 .or. jplo .lt. 0 .or.
     &             iphi .ge. iregsz(levc) .or. jphi .ge. jregsz(levc))


!--      write(*,*)" entering filrecur with level ",level
!--      write(*,*)"     and patch indices ilo,ihi,jlo,jhi ",
!--     &             ilo,ihi,jlo,jhi

c         write(*,*)" in filrecur for level ",level,mitot,mjtot
c
c        We begin by filling values for grids at level level. If all values can be
c        filled in this way, we return;

        nrowp   = ihi - ilo + 1
        ncolp   = jhi - jlo + 1
c        locuse  = igetsp(nrowp*ncolp)
        hxf     = hxposs(level)
        hyf     = hyposs(level)
        xlp     = xlower + ilo*hxf
        xrp     = xlower + (ihi+1)*hxf
        ybp     = ylower + jlo*hyf
        ytp     = ylower + (jhi+1)*hyf

        call intfil
     &  (valbig,mitot,mjtot,time,flaguse,nrowst,ncolst,
     &   ilo,ihi,jlo,jhi,level,nvar,naux)
c     &  (valbig,mitot,mjtot,time,locuse,nrowst,ncolst,

c
c Trimbd returns set = true if all of the entries are filled (=1.).
c set = false, otherwise. If set = true, then no other levels are
c are required to interpolate, and we return.
c
c Note that the used array is filled entirely in intfil, i.e. the
c marking done there also takes into account the points filled by
c the boundary conditions. bc2amr will be called later, after all 4
c boundary pieces filled.

c        call trimbd(alloc(locuse),nrowp,ncolp,set,il,ir,jb,jt)
        call trimbd(flaguse,nrowp,ncolp,set,il,ir,jb,jt)

        if (set) go to 90 ! all done except for bcs
c
c otherwise make recursive calls to coarser levels to fill remaining unset points
c
        if (level .eq. 1) then
           write(outunit,*)" error in filrecur - level 1 not set"
           write(outunit,900) nrowst,ncolst
           write(*,*)" error in filrecur - level 1 not set"
           write(*,*)" should not need more recursion "
           write(*,*)" to set patch boundaries"
           write(*,900) nrowst,ncolst
900        format("start at row: ",i4," col ",i4)
           stop
        endif

c set = false. we will have to interpolate some values from coarser
c levels. We begin by initializing the level level arrays, so that we can use
c purely recursive formulation for interpolating.


        levc = level - 1
        hxc  = hxposs(levc)
        hyc  = hyposs(levc)

        isl  = il + ilo - 1
        isr  = ir + ilo - 1
        jsb  = jb + jlo - 1
        jst  = jt + jlo - 1
c
c       coarsen
        lratiox = intratx(levc)
        lratioy = intraty(levc)
        iplo   = (isl-lratiox  +nghost*lratiox)/lratiox - nghost
        jplo   = (jsb-lratioy  +nghost*lratioy)/lratioy - nghost
        iphi   = (isr+lratiox  )/lratiox
        jphi   = (jst+lratioy  )/lratioy

        xlc  =  xlower + iplo*hxc
        ybc  =  ylower + jplo*hyc
        xrc  =  xlower + (iphi+1)*hxc
        ytc  =  ylower + (jphi+1)*hyc

        nrowc   =  iphi - iplo + 1
        ncolc   =  jphi - jplo + 1
        ntot    = nrowc*ncolc*(nvar+naux)
c        write(*,876) nrowc,ncolc, ihi-ilo+2,jhi-jlo+2
 876    format(" needed coarse grid size ",2i5," allocated ",2i5)
        if (nrowc .gt. ihi-ilo+2 .or. ncolc .gt. jhi-jlo+2) then
            write(*,*)" did not make big enough work space in filrecur "
            write(*,*)" need coarse space with nrowc,ncolc ",nrowc,ncolc
            write(*,*)" coarser level is ",levc," nghost ",nghost
            write(*,*)" with ratios lratiox lratioy ",lratiox,lratioy
            write(*,*)" made space for ilo,ihi,jlo,jhi ",ilo,ihi,jlo,jhi
            write(*,*)" isl,isr,jsb,jst ",isl,isr,jsb,jst
            write(*,*)"iplo,iphi,jplo,jphi ",iplo,iphi,jplo,jphi
            write(*,*)" orig il,ir,jb,jt ",il,ir,jb,jt
            write(*,*)"filpatch called with mitot,mjtot ",mitot,mjtot
            call outtre(1,.false.,nvar,naux)
            stop
        endif
c        loccrse = igetsp(ntot)
c        locauxc = loccrse + nrowc*ncolc*nvar
        if (naux.gt.0) then
              maxmx = nrowc - 2*nghost
              mx = maxmx
              maxmy = ncolc - 2*nghost
              my = maxmy
              xl = xlc + nghost*hxc
              yb = ybc + nghost*hyc
              call setaux(maxmx,maxmy,nghost,mx,my,xl,yb,hxc,hyc,
     &                    naux,auxcrse)
c     &                    naux,alloc(locauxc))
        endif

        if ((xperdom .or. (yperdom .or. spheredom)) .and.
     &       sticksout(iplo,iphi,jplo,jphi)) then
            call prefilrecur(levc,nvar,valcrse,auxcrse,
     1                    naux,time,nrowc,ncolc,1,1,
     2                    iplo,iphi,jplo,jphi)
        else
c          call filpatch2(levc,nvar,alloc(loccrse),alloc(locauxc),naux,
          call filrecur(levc,nvar,valcrse,auxcrse,naux,
     1                   time,nrowc,ncolc,1,1,
     2                   iplo,iphi,jplo,jphi)
        endif

c       interpolate back up

20      continue


*       !loop through coarse cells determining intepolation slopes
*       !these will be saved for fine grid loop
*       !prevents finding the same slope possibly lratiox*lratioy times
*       !all fine gid depths will be found before any momentum
         reloop = .false.
         toldry= drytolerance
         do ic  = 2, nrowc-1
         do jc  = 2, ncolc-1
            icount(icrse(ic,jc)) = 0
            finemass(icrse(ic,jc)) = 0.d0
            do ivar=1,nvar
               fineflag(ivalc(ivar,ic,jc)) = .false.
            enddo

*           !find interpolation slope for eta = q(1,:)+ aux(1,:)
            do i=-1,1
               etacrse(icrse(ic+i,jc)) = valcrse(ivalc(1,ic+i,jc))
     &            +  auxcrse(iauxc(ic+i,jc))
               if (valcrse(ivalc(1,ic+i,jc)).lt.toldry) then
                  etacrse(icrse(ic+i,jc)) = sealevel
                  endif
               enddo
            s1 = etacrse(icrse(ic,jc))- etacrse(icrse(ic-1,jc))
            s2 = etacrse(icrse(ic+1,jc))- etacrse(icrse(ic,jc))
            if (s1*s2.le.0) then
               slopex(icrse(ic,jc))= 0.d0
            else
               slopex(icrse(ic,jc))=dmin1(dabs(s1),dabs(s2))*dsign(1.d0,
     &               etacrse(icrse(ic+1,jc))- etacrse(icrse(ic-1,jc)))
               endif
            do j=-1,1
               etacrse(icrse(ic,jc+j)) = valcrse(ivalc(1,ic,jc+j))
     &            +  auxcrse(iauxc(ic,jc+j))
               if (valcrse(ivalc(1,ic,jc+j)).lt.toldry) then
                  etacrse(icrse(ic,jc+j)) = sealevel
                  endif
               enddo
            s1 = etacrse(icrse(ic,jc))  - etacrse(icrse(ic,jc-1))
            s2 = etacrse(icrse(ic,jc+1))- etacrse(icrse(ic,jc))
            if (s1*s2.le.0) then
               slopey(icrse(ic,jc))= 0.d0
            else
               slopey(icrse(ic,jc))=dmin1(dabs(s1),dabs(s2))*dsign(1.d0,
     &               etacrse(icrse(ic,jc+1))- etacrse(icrse(ic,jc-1)))
               endif

            end do
            end do

*        !loop through patch: note this includes multiple coarse cells
         do iff = 1,nrowp
            ic = 2 + (iff - (isl - ilo) - 1)/lratiox
            eta1 = (-0.5d0+dble(mod(iff-1,lratiox)))/dble(lratiox)
         do jf  = 1,ncolp
            jc = 2 + (jf - (jsb - jlo) - 1)/lratioy
            eta2 = (-0.5d0+dble(mod(jf -1,lratioy)))/dble(lratioy)
            cc = icrse(ic,jc)
c           flag = alloc(iadflag(iff,jf))
            flag = flaguse(iff,jf)
            if (flag .eq. 0.0) then
*              !interp. from coarse cells to fine grid to find eta
               icount(icrse(ic,jc)) = icount(icrse(ic,jc)) + 1
               etafine =  etacrse(icrse(ic,jc))
     &            + eta1*slopex(icrse(ic,jc))+eta2*slopey(icrse(ic,jc))
               hfine = max(etafine - aux(1,iff+nrowst-1,jf+ncolst-1),
     &            0.d0)
               valbig(1,iff+nrowst-1,jf+ncolst-1) = hfine

               finemass(icrse(ic,jc)) = finemass(icrse(ic,jc)) +
     &                        valbig(1,iff+nrowst-1,jf+ncolst-1)
               if (valbig(1,iff+nrowst-1,jf+ncolst-1).lt.toldry) then
                  fineflag(ivalc(1,ic,jc)) = .true.
                  reloop = .true.
                  endif
               endif
            enddo
            enddo

c        ! determine momentum
         do ivar = 2,nvar
*           !find interpolation slope for momentum = q(:,ivar)
            do ic  = 2, nrowc-1
            do jc  = 2, ncolc-1

               s1 = valcrse(ivalc(ivar,ic,jc))
     &               - valcrse(ivalc(ivar,ic-1,jc))
               s2 = valcrse(ivalc(ivar,ic+1,jc))
     &               - valcrse(ivalc(ivar,ic,jc))
               if (s1*s2.le.0) then
                  slopex(icrse(ic,jc))= 0.d0
               else
                  slopex(icrse(ic,jc))=dmin1(dabs(s1),dabs(s2))
     &               *dsign(1.d0, valcrse(ivalc(ivar,ic+1,jc))
     &                  - valcrse(ivalc(ivar,ic-1,jc)))
                  endif
               s1 = valcrse(ivalc(ivar,ic,jc))
     &               - valcrse(ivalc(ivar,ic,jc-1))
               s2 = valcrse(ivalc(ivar,ic,jc+1))
     &               - valcrse(ivalc(ivar,ic,jc))
               if (s1*s2.le.0) then
                  slopey(icrse(ic,jc))= 0.d0
               else
                  slopey(icrse(ic,jc))=dmin1(dabs(s1),dabs(s2))
     &               *dsign(1.d0, valcrse(ivalc(ivar,ic,jc+1))
     &                  - valcrse(ivalc(ivar,ic,jc-1)))
                  endif

               if (valcrse(ivalc(1,ic,jc)).gt.toldry) then
                  velmax(icrse(ic,jc)) = valcrse(ivalc(ivar,ic,jc))
     &                                       /valcrse(ivalc(1,ic,jc))
                  velmin(icrse(ic,jc)) =  valcrse(ivalc(ivar,ic,jc))
     &                                       /valcrse(ivalc(1,ic,jc))
               else
                  velmax(icrse(ic,jc)) = 0.d0
                  velmin(icrse(ic,jc)) = 0.d0
                  endif

*              !look for bounds on velocity to avoid generating new extrema
*              !necessary since interpolating momentum linearly
*              !yet depth is not interpolated linearly
               do ii = -1,1,2
                  if (valcrse(ivalc(1,ic+ii,jc)).gt.toldry) then
                     velmax(icrse(ic,jc)) = max(velmax(icrse(ic,jc))
     &                  ,valcrse(ivalc(ivar,ic+ii,jc))
     &                  /valcrse(ivalc(1,ic+ii,jc)))
                     velmin(icrse(ic,jc)) = min(velmin(icrse(ic,jc))
     &                  ,valcrse(ivalc(ivar,ic+ii,jc))
     &                  /valcrse(ivalc(1,ic+ii,jc)))
                     endif
                  if (valcrse(ivalc(1,ic,jc+ii)).gt.toldry) then
                     velmax(icrse(ic,jc)) = max(velmax(icrse(ic,jc))
     &                  ,valcrse(ivalc(ivar,ic,jc+ii))
     &                  /valcrse(ivalc(1,ic,jc+ii)))
                     velmin(icrse(ic,jc)) = min(velmin(icrse(ic,jc))
     &                  ,valcrse(ivalc(ivar,ic,jc+ii))
     &                  /valcrse(ivalc(1,ic,jc+ii)))
                     endif
                  enddo

               end do
               end do

*           !determine momentum in fine cells
            do iff = 1,nrowp
               ic = 2 + (iff - (isl - ilo) - 1)/lratiox
               eta1 = (-0.5d0+dble(mod(iff-1,lratiox)))/dble(lratiox)
            do jf  = 1,ncolp
               jc = 2 + (jf  - (jsb - jlo) - 1)/lratioy
               eta2 = (-0.5d0+dble(mod(jf -1,lratioy)))/dble(lratioy)

               flag = flaguse(iff,jf)
               if (flag .eq. 0.0) then
                  if (.not.(fineflag(ivalc(1,ic,jc)))) then
*                    !this is a normal wet cell. intepolate normally
                     hvf = valcrse(ivalc(ivar,ic,jc))
     &                   + eta1*slopex(icrse(ic,jc))
     &                   + eta2*slopey(icrse(ic,jc))
                     vf = hvf/valbig(1,iff+nrowst-1,jf+ncolst-1)
                     if (vf.lt.velmin(icrse(ic,jc)).or.
     &                        vf.gt.velmax(icrse(ic,jc))) then
                        fineflag(ivalc(ivar,ic,jc))=.true.
                        reloop = .true.
                     else
                        valbig(ivar,iff+nrowst-1,jf+ncolst-1) = hvf
                        endif
                     endif
                  endif
               enddo
               enddo


*           !loop again and reset momentum to conserve momentum
*           !in the case of gained momentum, or if velocity bounds violated
            if (reloop) then
               do iff = 1,nrowp
                  ic = 2 + (iff - (isl - ilo) - 1)/lratiox
                  eta1 =(-0.5d0+dble(mod(iff-1,lratiox)))/dble(lratiox)
               do jf  = 1,ncolp
                  jc = 2 + (jf  - (jsb - jlo) - 1)/lratioy
                  eta2 = (-0.5d0+dble(mod(jf -1,lratioy)))/dble(lratioy)
                  flag = flaguse(iff,jf)
                  if (flag.eq.0.0) then
                     if (fineflag(ivalc(1,ic,jc))
     &               .or.fineflag(ivalc(ivar,ic,jc))) then
                        if (finemass(icrse(ic,jc)).gt.toldry) then
                           hcrse = valcrse(ivalc(1,ic,jc))
                           hcnt = dble(icount(icrse(ic,jc)))
                           hfineave = finemass(icrse(ic,jc))/hcnt
                           dividemass = max(hcrse,hfineave)
                           hfine = valbig(1,iff+nrowst-1,jf+ncolst-1)
                           Vnew = valcrse(ivalc(ivar,ic,jc))/dividemass
                           valbig(ivar,iff+nrowst-1,jf+ncolst-1) =
     &                           Vnew*valbig(1,iff+nrowst-1,jf+ncolst-1)
                        else
                           valbig(ivar,iff+nrowst-1,jf+ncolst-1)=0.d0
                           endif
                        endif
                     endif
                  enddo
                  enddo
               endif

            enddo




c        call reclam(loccrse,ntot)

 90       continue
c
c  set bcs, whether or not recursive calls needed. set any part of patch that stuck out
c
        call bc2amr(valbig,aux,mitot,mjtot,nvar,naux,
     1              hxf,hyf,level,time,
     2              xlp,xrp,ybp,ytp,
     3              xlower,ylower,xupper,yupper,
     4              xperdom,yperdom,spheredom)


c        call reclam(locuse,nrowp*ncolp)

        return
        end
