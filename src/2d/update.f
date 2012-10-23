
c
c -----------------------------------------------------------
c
      subroutine update (level, nvar, naux)
      
      use geoclaw_module, only: rho, dry_tolerance
c
c     # modified for shallow water on topography to use surface level eta
c     # rather than depth h = q(i,j,1)
c     # eta = q(i,j,1) + aux(i,j,1)
c
      use amr_module
      implicit double precision (a-h,o-z)


      integer listgrids(numgrids(level))

      iadd(ivar,i,j)  = loc    + ivar-1 + nvar*((j-1)*mitot+i-1)
      iaddf(ivar,i,j) = locf   + ivar-1 + nvar*((j-1)*mi+i-1)
      iaddfaux(i,j)   = locfaux + mcapa-1 + naux*((j-1)*mi + (i-1))
      iaddcaux(i,j)   = loccaux + mcapa-1 + naux*((j-1)*mitot+(i-1))
c  topo in aux component 1
      iaddftopo(i,j)   = locfaux +  naux*((j-1)*mi + (i-1))
      iaddctopo(i,j)   = loccaux +  naux*((j-1)*mitot+(i-1))

c
c :::::::::::::::::::::::::: UPDATE :::::::::::::::::::::::::::::::::
c update - update all grids at level 'level'.
c          this routine assumes cell centered variables.
c          the update is done from 1 level finer meshes under it.
c input parameter:
c    level  - ptr to the only level to be updated. levels coarser than
c             this will be at a diffeent time.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      lget = level
      if (uprint) write(outunit,100) lget
100   format(19h    updating level ,i5)
c     need to set up data structure for parallel distrib of grids
      call prepgrids(listgrids,numgrids(level),level)
c
c  grid loop for each level
c
      dt     = possk(lget)

c      mptr = lstart(lget)
c 20   if (mptr .eq. 0) go to 85

!$OMP PARALLEL DO PRIVATE(ng,mptr,loc,loccaux,nx,ny,mitot,mjtot,
!$OMP&                    ilo,jlo,ihi,jhi,locuse,mkid,iclo,ichi,
!$OMP&                    jclo,jchi,mi,mj,locf,locfaux,iplo,iphi,
!$OMP&                    jplo,jphi,iff,jff,totrat,i,j,ivar,capac,
!$OMP&                    capa,bc,etasum,hsum,husum,hvsum,drytol,
!$OMP&                    newt,ico, jco,hf,bf,huf,hvf,
!$OMP&                    etaf,etaav,hav,nwet,hc,huc,hvc),
!$OMP&            SHARED(lget,numgrids,listgrids,level,intratx,intraty,
!$OMP&                   nghost,uprint,nvar,naux,mcapa,node,listsp,
!$OMP&                   alloc,lstart,dry_tolerance,rho),
!$OMP&            DEFAULT(none)
      do ng = 1, numgrids(lget)
c        mptr    = mget(ng, level)
         mptr    = listgrids(ng)
         loc     = node(store1,mptr)
         loccaux = node(storeaux,mptr)
         nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         mitot   = nx + 2*nghost
         mjtot   = ny + 2*nghost
         ilo     = node(ndilo,mptr)
         jlo     = node(ndjlo,mptr)
         ihi     = node(ndihi,mptr)
         jhi     = node(ndjhi,mptr)
c
         if (node(cfluxptr,mptr) .eq. 0) go to 25
c         locuse = igetsp(mitot*mjtot)
         call upbnd(alloc(node(cfluxptr,mptr)),alloc(loc),nvar,naux,
     1              mitot,mjtot,listsp(lget),mptr) ! took out next to last arg
c     1              mitot,mjtot,listsp(lget),alloc(locuse),mptr)
c         call reclam(locuse,mitot*mjtot)
c
c  loop through all intersecting fine grids as source updaters.
c
 25      mkid = lstart(lget+1)
 30        if (mkid .eq. 0) go to 80
           iclo   = node(ndilo,mkid)/intratx(lget)
           jclo   = node(ndjlo,mkid)/intraty(lget)
           ichi   = node(ndihi,mkid)/intratx(lget)
           jchi   = node(ndjhi,mkid)/intraty(lget)

           mi      = node(ndihi,mkid)-node(ndilo,mkid) + 1 + 2*nghost
           mj      = node(ndjhi,mkid)-node(ndjlo,mkid) + 1 + 2*nghost
           locf    = node(store1,mkid)
           locfaux = node(storeaux,mkid)
c
c  calculate starting and ending indices for coarse grid update, if overlap
c
         iplo = max(ilo,iclo)
         jplo = max(jlo,jclo)
         iphi = min(ihi,ichi)
         jphi = min(jhi,jchi)

         if (iplo .gt. iphi .or. jplo .gt. jphi) go to 75
c
c  calculate starting index for fine grid source pts.
c
         iff    = iplo*intratx(lget) - node(ndilo,mkid) + nghost + 1
         jff    = jplo*intraty(lget) - node(ndjlo,mkid) + nghost + 1
         totrat = intratx(lget) * intraty(lget)

         do 71 i = iplo-ilo+nghost+1, iphi-ilo+nghost+1
         do 70 j = jplo-jlo+nghost+1, jphi-jlo+nghost+1
           if (uprint) then
              write(outunit,101) i,j,mptr,iff,jff,mkid
 101          format(' updating pt. ',2i4,' of grid ',i3,' using ',2i4,
     1               ' of grid ',i4)
              write(outunit,102)(alloc(iadd(ivar,i,j)),ivar=1,nvar)
 102          format(' old vals: ',4e12.4)
           endif
c
c
c  update using intrat fine points in each direction
c
c
c
c     # modified by dlg 11/01/05 below:
c
c     Strategy:
c     average eta,h,hu,hv in wet cells only.
c     set hcoarse = ave(eta_wet) - b_coarse, or zero if that is negative

c     set hucoarse = (min(hcoarse,h_average))/h_average)* hu_average
c     note: this conserves mass and momentum as long as all fine cells are wet.
c     Momentum is reduced by the same factor as mass in the coarse cell
c     and is never increased given an increase in mass

      if (mcapa .eq. 0) then
         capac=1.0d0
      else
         capac=alloc(iaddcaux(i,j))
         endif

      bc = alloc(iaddctopo(i,j))

      etasum = 0.d0
      hsum = 0.d0
      husum = 0.d0
      hvsum = 0.d0

      nwet=0

      do jco  = 1, intraty(lget)
         do ico  = 1, intratx(lget)
            if (mcapa .eq. 0) then
               capa=1.0d0
            else
               capa=alloc(iaddfaux(iff+ico-1,jff+jco-1))
               endif

            hf = alloc(iaddf(1,iff+ico-1,jff+jco-1))*capa / rho(1)
            bf = alloc(iaddftopo(iff+ico-1,jff+jco-1))*capa
            huf= alloc(iaddf(2,iff+ico-1,jff+jco-1))*capa / rho(1)
            hvf= alloc(iaddf(3,iff+ico-1,jff+jco-1))*capa / rho(1)

            if (hf > dry_tolerance(1)) then
               etaf = hf+bf
               nwet=nwet+1
            else
               etaf = 0.d0
               huf=0.d0
               hvf=0.d0
               endif

               hsum   = hsum + hf
               husum  = husum + huf
               hvsum  = hvsum + hvf
               etasum = etasum + etaf
            enddo
         enddo

      if (nwet.gt.0) then
         etaav=etasum/dble(nwet)
         hav= hsum/dble(nwet)
*         hc=max(etaav-bc*capac,0.d0) !tsunamiclaw method
         hc=min(hav,(max(etaav-bc*capac,0.d0)))
         huc=(min(hav,hc)/hsum)*husum
         hvc=(min(hav,hc)/hsum)*hvsum
      else
         hc=0.d0
         huc=0.d0
         hvc=0.d0
         endif

c     # set h on coarse grid based on surface, not conservative near shoreline

      alloc(iadd(1,i,j)) = hc / capac * rho(1)
      alloc(iadd(2,i,j)) = huc / capac * rho(1)
      alloc(iadd(3,i,j)) = hvc / capac * rho(1)
c
      if (uprint) write(outunit,103)(alloc(iadd(ivar,i,j)),
     .     ivar=1,nvar)
 103  format(' new vals: ',4e12.4)
c
      jff = jff + intraty(lget)
 70   continue
      iff = iff + intratx(lget)
      jff    = jplo*intraty(lget) - node(ndjlo,mkid) + nghost + 1
 71   continue
c
 75   mkid = node(levelptr,mkid)
      go to 30
c
 80   continue
      end do
!$OMP END PARALLEL DO
c
c 80   mptr = node(levelptr, mptr)
c      go to 20
c
c 85   continue
c
 99   return
      end
