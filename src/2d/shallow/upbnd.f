c
c ------------------------------------------------------------
c
       subroutine upbnd(listbc,val,nvar,maux,mitot,mjtot,
     1                  maxsp,mptr)
c     1                  maxsp,iused,mptr)
 
      use geoclaw_module, only: dry_tolerance
      use amr_module
      implicit double precision (a-h,o-z)

 
       dimension val(nvar,mitot,mjtot),listbc(5,maxsp),
     1           iused(mitot,mjtot)

c  OLD INDEXING
c      iaddaux(i,j) = locaux + i-1 +  mitot*(j-1) 
c    1                 + mitot*mjtot*(mcapa-1)
c NEW INDEXING ORDER SWITCHED
       iaddaux(i,j) = locaux + mcapa-1 +  maux*(i-1) 
     1                 + maux*mitot*(j-1)
 
c
c :::::::::::::::::::::::::::: UPBND :::::::::::::::::::::::::::::
c We now correct the coarse grid with the flux differences stored
c with each of the fine grids. We use an array   iused
c to indicate whether the flux has been updated or not for that zone.
c iused(i,j) = sum from (l=1,4) i(l)*2**(l-1), where i(l) = 1 if the
c flux for the  l-th side of the (i,j)-th cell has already been
c updated, and i(l) = 0 if not.
 
c if there is a capacity fn. it needs to be included in update formula
c indicated by mcapa not zero (is index of capacity fn.)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
 
      do 10 j=1,mjtot
      do 10 i=1,mitot
         iused(i,j) = 0.
 10   continue
 
      locaux = node(storeaux,mptr)
      levc   = node(nestlevel,mptr)
      area   = hxposs(levc)*hyposs(levc)


      if (uprint) write(outunit,*)" upbnding grid ",mptr

      do 40 ispot = 1,maxsp
         icrse = listbc(1,ispot)
         if (icrse.eq.0) go to 99
 
         jcrse = listbc(2,ispot)
         iside = listbc(3,ispot)
c        continue to use iside1/norm for debugging, but should soon delete         
c        this if/then/else block needed due to new categories corresponding
c        to mapped bcs. should still only have one update per side of coarse cell though
         if (iside .lt. 5) then   
           iside1 = iside
         elseif (iside .eq. 5) then
           iside1 = 2
         else  ! iside is 6
           iside1 = 4
         endif
         norm = 2**(iside1-1)
         iflag =iused(icrse,jcrse)/norm
         if (mod(iflag,2).eq.1) then
            write(6,*)" ***  double flux update CAN happen in upbnd ***"
            go to 40
         endif
         mkid = listbc(4,ispot)
         kidlst = node(ffluxptr,mkid)
         lkid = listbc(5,ispot)
c         if (mod(iside,4).gt.1) then
c         modified to include other side options
          if (iside .eq. 2 .or. iside .eq. 3 .or. iside .eq. 6) then
c           (iside .eq. 2 .or. iside .eq. 3)
            sgnm = -1.
         else
c           (iside .eq. 4 .or. iside .eq. 1)
            sgnm = 1.
         endif

c        ## debugging output
         if (uprint) then
           write(outunit,101) icrse,jcrse,
     .         (val(ivar,icrse,jcrse),ivar=1,nvar)
 101       format(" old ",1x,2i4,4e15.7)
         endif

         if (mcapa .gt. 0) then
c            # capacity array:  need to divide by capa in each cell.
c            # modify sgnm which is reset for each grid cell.
c            # Note capa is stored in aux(icrse,jcrse,mcapa)
             sgnm = sgnm / alloc(iaddaux(icrse,jcrse))
         endif

c        # If coarse cell is dry then don't include updates from fine
c        # grid fluxes that might give unphysical wetting because coarse
c        # cell may be much higher than some fine cells:
         if (val(1,icrse,jcrse) > dry_tolerance) then
           do 20 ivar = 1,nvar
            val(ivar,icrse,jcrse) = val(ivar,icrse,jcrse) +
     1      sgnm*alloc(kidlst+nvar*(lkid-1)+ivar-1)/area
 20        continue
         else
c          write(6,*) '+++ throw out ',alloc(kidlst+nvar*(lkid-1))
           continue
         endif
         
c        # Reset small h to zeros once again:
         if (val(1,icrse,jcrse) < dry_tolerance) then
               val(:,icrse,jcrse) = 0.d0
         endif
         iused(icrse,jcrse) = iused(icrse,jcrse) + norm

c        ## debugging output
         if (uprint) then
           write(outunit,102) mkid,
     .         (val(ivar,icrse,jcrse),ivar=1,nvar)
 102       format(" new ","(grid",i3,")",4e15.7)
         endif

 40   continue
c
 99   return
      end
