c
c ------------------------------------------------------------------
c
      subroutine filval(val,mitot,mjtot,hx,hy,lev,time,
     1                  valc,auxc,mic,mjc,
     2                  xleft,xright,ybot,ytop,nvar,
     3                  mptr,ilo,ihi,jlo,jhi,aux,naux,locflip,
     4                  sp_over_h)

      use geoclaw_module

      implicit double precision (a-h,o-z)

      include "call.i"

      dimension   val(nvar,mitot,mjtot), valc(nvar,mic,mjc)
      dimension   aux(naux,mitot,mjtot), auxc(naux,mic,mjc)

      double precision coarseval(3)
      logical fineflag(3)


c
c :::::::::::::::::::::::::::::: FILVAL ::::::::::::::::::::::::::
c
c create and fill coarser (lev-1) patch with one extra coarse cell all
c around, plus the ghost cells . will interpolate from this patch to grid mptr
c without needing special boundary code.
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     # indext into eta array for surface values:
c       iaddeta(i,j) = loceta + i-1 + mic*(j-1)

      levc    = lev - 1
      lratiox = intratx(levc)
      lratioy = intraty(levc)
      hxcrse  = hx*lratiox
      hycrse  = hy*lratioy
      xl      = xleft  - hxcrse
      xr      = xright + hxcrse
      yb      = ybot   - hycrse
      yt      = ytop   + hycrse
c
c     set integer indices for coarser patch enlarged by 1 cell
c     (can stick out of domain). proper nesting will insure this one
c     call is sufficient.
      iclo   = ilo/lratiox - 1
      jclo   = jlo/lratioy - 1
      ichi   = (ihi+1)/lratiox - 1 + 1
      jchi   = (jhi+1)/lratioy - 1 + 1
      ng     = 0

c    :::  mcapa  is the capacity function index

      if (naux .eq. 0) then
c     if (mcapa .eq. 0) then
        if (xperdom .or. yperdom .or. spheredom) then
          call preintcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,levc,
     &                    locflip)
        else
          call intcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,levc,1,1)
        endif
      else  ! intersect grids and copy all (soln and aux)
        if (xperdom .or. yperdom .or. spheredom) then
          call preicall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,
     &                  jchi,levc,locflip)
        else
          call icall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi,
     &               levc,1,1)
        endif
      endif
      call bc2amr(valc,auxc,mic,mjc,nvar,naux,
     1            hxcrse,hycrse,levc,time,
     2            xl,xr,yb,yt,
     3            xlower,ylower,xupper,yupper,
     4            xperdom,yperdom,spheredom)

c-----------------------------
c     # for shallow water over topography,
c     # in coarse cells convert from h,
c     # to eta,  before interpolating:
c-----------------------------
      toldry = drytolerance
c     #prepare slopes - use min-mod limiters
      do j=2, mjc-1
      do i=2, mic-1
         fineflag(1) = .false.
*        !interpolate eta to find depth---------------------------------------
         do ii=-1,1
            coarseval(2+ii) = valc(1,i+ii,j)+ auxc(1,i+ii,j)
            if (valc(1,i+ii,j).le.toldry) then
               coarseval(2+ii)=sealevel
               endif
            enddo
         s1p=coarseval(3)-coarseval(2)
         s1m=coarseval(2)-coarseval(1)
         slopex=dmin1(dabs(s1p),dabs(s1m))*dsign(1.d0,
     &      coarseval(3)-coarseval(1))
         if (s1m*s1p.le.0.d0) slopex=0.d0
         do jj=-1,1
            coarseval(2+jj) = valc(1,i,j+jj)+ auxc(1,i,j+jj)
            if (valc(1,i,j+jj).le.toldry) then
               coarseval(2+jj)=sealevel
               endif
            enddo
         s1p=coarseval(3)-coarseval(2)
         s1m=coarseval(2)-coarseval(1)
         slopey=dmin1(dabs(s1p),dabs(s1m))*dsign(1.d0,
     &      coarseval(3)-coarseval(1))
         if (s1m*s1p.le.0.d0) slopey=0.d0
*       !interp. from coarse cells to fine grid to find depth
         finemass = 0.d0
         do ico = 1,lratiox
            do jco = 1,lratioy
               yoff = (float(jco) - .5)/lratioy - .5
               xoff = (float(ico) - .5)/lratiox - .5
               jfine = (j-2)*lratioy + nghost + jco
               ifine = (i-2)*lratiox + nghost + ico
               val(1,ifine,jfine) = coarseval(2)
     &            + xoff*slopex + yoff*slopey
               val(1,ifine,jfine)=max(0.d0,
     &            val(1,ifine,jfine)-aux(1,ifine,jfine))
               finemass = finemass + val(1,ifine,jfine)
               if (val(1,ifine,jfine).le.toldry) then
                  fineflag(1) = .true.
                  val(2,ifine,jfine)=0.d0
                  val(3,ifine,jfine)=0.d0
                  endif
               enddo
            enddo
* !------determine momentum----------------------------------
*        !finemass is the total mass in all new fine grid cells
*        !all fine mass has been determined for this coarse grid cell
*        !if all fine cells are dry, momentum has already been set
         if (finemass.ge.toldry) then
            do ivar = 2,nvar
               fineflag(ivar)=.false.
               s1p=valc(ivar,i+1,j)-valc(ivar,i,j)
               s1m=valc(ivar,i,j)-valc(ivar,i-1,j)
               slopex=dmin1(dabs(s1p),dabs(s1m))*dsign(1.d0,
     &            valc(ivar,i+1,j)-valc(ivar,i-1,j))
               if (s1m*s1p.le.0.d0) slopex=0.d0
               s1p=valc(ivar,i,j+1)-valc(ivar,i,j)
               s1m=valc(ivar,i,j)-valc(ivar,i,j-1)
               slopey=dmin1(dabs(s1p),dabs(s1m))*dsign(1.d0,
     &            valc(ivar,i,j+1)-valc(ivar,i,j-1))
               if (s1m*s1p.le.0.d0) slopey=0.d0
               if (valc(1,i,j).gt.toldry) then
                  velmax = valc(ivar,i,j)/valc(1,i,j)
                  velmin = valc(ivar,i,j)/valc(1,i,j)
               else
                  velmax = 0.d0
                  velmin = 0.d0
                  endif
               do ii = -1,1,2
                  if (valc(1,i+ii,j).gt.toldry) then
                     vel = valc(ivar,i+ii,j)/valc(1,i+ii,j)
                     velmax = max(vel,velmax)
                     velmin = min(vel,velmin)
                     endif
                  if (valc(1,i,j+ii).gt.toldry) then
                     vel = valc(ivar,i,j)/valc(1,i,j+ii)
                     velmax = max(vel,velmax)
                     velmin = min(vel,velmin)
                     endif
                  enddo

*              !try to set momentum
               do ico = 1,lratiox
                  if (fineflag(1).or.fineflag(ivar)) exit
                  do jco = 1,lratioy
                     jfine = (j-2)*lratioy + nghost + jco
                     ifine = (i-2)*lratiox + nghost + ico
                     yoff = (float(jco) - .5)/lratioy - .5
                     xoff = (float(ico) - .5)/lratiox - .5
                     hvf = valc(ivar,i,j)+ xoff*slopex + yoff*slopey
                     vf = hvf/val(1,ifine,jfine)
                     if (vf.gt.velmax.or.vf.lt.velmin) then
                        fineflag(ivar)=.true.
                        exit
                     else
                        val(ivar,ifine,jfine) = hvf
                        endif
                     enddo
                  enddo

*              !momentum is set to preserve old momentum or not violate
*              !generating new extrema in velocities
               if (fineflag(1).or.fineflag(ivar)) then !more mass now, conserve momentum
                  area = dble(lratiox*lratioy)
                  dividemass = max(finemass,valc(1,i,j))
                  Vnew = area*valc(ivar,i,j)/dividemass

                  do ico = 1,lratiox
                     do jco = 1,lratioy
                        jfine = (j-2)*lratioy + nghost + jco
                        ifine = (i-2)*lratiox + nghost + ico
                        val(ivar,ifine,jfine) = Vnew*val(1,ifine,jfine)
                        enddo
                     enddo
                  endif

               enddo
            endif

         enddo !end of coarse loop
         enddo !end of coarse loop

c
c      if (mcapa .ne. 0) then
c        call fixcapaq(val,aux,mitot,mjtot,valc,auxc,mic,mjc,
c     &                nvar,naux,levc)
c      endif
c
c  overwrite interpolated values with fine grid values, if available.
c
      call intcopy(val,mitot,mjtot,nvar,ilo-nghost,ihi+nghost,
     &             jlo-nghost,jhi+nghost,lev,1,1)

c
c    scan for max wave speed on newly created grid. this will be used to set appropriate
c    time step and appropriate refinement in time. For this app not nec to refine by same
c    amount as refinement in space since refinement at shores where h is shallow has lower
c    speeds.
c
      if (varRefTime) then   ! keep consistent with setgrd_geo and qinit_geo
         sp_over_h = get_max_speed(val,mitot,mjtot,nvar,aux,naux,nghost,
     &                           hx,hy)
      endif

 99   return
      end
