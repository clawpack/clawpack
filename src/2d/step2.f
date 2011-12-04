c
c
c     ==========================================================
      subroutine step2(maxm,maxmx,maxmy,meqn,maux,mbc,mx,my,
     &                 qold,aux,dx,dy,dt,cflgrid,
     &                 fm,fp,gm,gp,
     &                 faddm,faddp,gaddm,gaddp,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,rpn2,rpt2)
c     ==========================================================
c
c     # clawpack routine ...  modified for AMRCLAW
c
c     # Take one time step, updating q.
c     # On entry, qold gives
c     #    initial data for this step
c     #    and is unchanged in this version.
c
c     # fm, fp are fluxes to left and right of single cell edge
c     # See the flux2 documentation for more information.
c

c     # modified again for GeoClaw

c------------------step2_geo.f-----------------------
c     This version of step2 relimits the fluxes in order to
c     maintain positivity.

c     to do so set relimit=.true.

c     The only modification is the 101 loop

c------------last modified 12/30/04--------------------------
c
      use geoclaw_module
    
      implicit double precision (a-h,o-z)


      dimension qold(meqn,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      dimension   fm(meqn,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      dimension   fp(meqn,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      dimension   gm(meqn,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      dimension   gp(meqn,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      dimension  q1d(meqn,1-mbc:maxm+mbc)
      dimension faddm(meqn,1-mbc:maxm+mbc)
      dimension faddp(meqn,1-mbc:maxm+mbc)
      dimension gaddm(meqn,1-mbc:maxm+mbc, 2)
      dimension gaddp(meqn,1-mbc:maxm+mbc, 2)
      dimension aux(maux,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      dimension aux1(maux,1-mbc:maxm+mbc)
      dimension aux2(maux,1-mbc:maxm+mbc)
      dimension aux3(maux,1-mbc:maxm+mbc)
      dimension dtdx1d(1-mbc:maxm+mbc)
      dimension dtdy1d(1-mbc:maxm+mbc)
      dimension work(mwork)

      logical relimit
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
      include "call.i"
c

c     # to relimit fluxes to maintain positivity
      relimit=.false.
c     # for relimiting to keep values very near zero positive.

c
c     # store mesh parameters that may be needed in Riemann solver but not
c     # passed in...
      dxcom = dx
      dycom = dy
      dtcom = dt
c
c
c     # partition work array into pieces needed for local storage in
c     # flux2 routine.  Find starting index of each piece:
c
      i0wave = 1
      i0s = i0wave + (maxm+2*mbc)*meqn*mwaves
      i0amdq = i0s + (maxm+2*mbc)*mwaves
      i0apdq = i0amdq + (maxm+2*mbc)*meqn
      i0cqxx = i0apdq + (maxm+2*mbc)*meqn
      i0bmadq = i0cqxx + (maxm+2*mbc)*meqn
      i0bpadq = i0bmadq + (maxm+2*mbc)*meqn
      iused = i0bpadq + (maxm+2*mbc)*meqn - 1
c
      if (iused.gt.mwork) then
c        # This shouldn't happen due to checks in claw2
         write(outunit,*) 'not enough work space in step2'
         write(*      ,*) 'not enough work space in step2'
         stop
         endif
c
c
      cflgrid = 0.d0
      dtdx = dt/dx
      dtdy = dt/dy
c
      do 10 j=1-mbc,my+mbc
         do 10 i=1-mbc,mx+mbc
            do 10 m=1,meqn
               fm(m,i,j) = 0.d0
               fp(m,i,j) = 0.d0
               gm(m,i,j) = 0.d0
               gp(m,i,j) = 0.d0
   10          continue
c
      if (mcapa.eq.0) then
c        # no capa array:
         do 5 i=1-mbc,maxm+mbc
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
    5       continue
         endif

c
c
c     # perform x-sweeps
c     ==================
c
      do 50 j = 0,my+1
c
c        # copy data along a slice into 1d arrays:
         do 20 i = 1-mbc, mx+mbc
           do 20 m=1,meqn
               q1d(m,i) = qold(m,i,j)
   20          continue
c
         if (mcapa.gt.0)  then
           do 21 i = 1-mbc, mx+mbc
               dtdx1d(i) = dtdx / aux(mcapa,i,j)
   21          continue
           endif
c
         if (maux .gt. 0)  then
             do 22 i = 1-mbc, mx+mbc
               do 22 ma=1,maux
                 aux1(ma,i) = aux(ma,i,j-1)
                 aux2(ma,i) = aux(ma,i,j  )
                 aux3(ma,i) = aux(ma,i,j+1)
   22            continue
           endif
c
c
c        # Store the value of j along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         jcom = j
c
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(1,maxm,meqn,maux,mbc,mx,
     &              q1d,dtdx1d,aux1,aux2,aux3,
     &              faddm,faddp,gaddm,gaddp,cfl1d,
     &              work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &              work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2)
         cflgrid = dmax1(cflgrid,cfl1d)
c
c        # update fluxes for use in AMR:
c
          do 25 i=1,mx+1
           do 25 m=1,meqn
               fm(m,i,j) = fm(m,i,j) + faddm(m,i)
               fp(m,i,j) = fp(m,i,j) + faddp(m,i)
               gm(m,i,j) = gm(m,i,j) + gaddm(m,i,1)
               gp(m,i,j) = gp(m,i,j) + gaddp(m,i,1)
               gm(m,i,j+1) = gm(m,i,j+1) + gaddm(m,i,2)
               gp(m,i,j+1) = gp(m,i,j+1) + gaddp(m,i,2)
   25          continue
   50    continue
c
c
c
c     # perform y sweeps
c     ==================
c
c
      do 100 i = 0, mx+1
c
c        # copy data along a slice into 1d arrays:
         do 70 j = 1-mbc, my+mbc
           do 70 m=1,meqn
               q1d(m,j) = qold(m,i,j)
   70          continue
c
         if (mcapa.gt.0)  then
           do 71 j = 1-mbc, my+mbc
               dtdy1d(j) = dtdy / aux(mcapa,i,j)
   71          continue
           endif
c
         if (maux .gt. 0)  then
             do 72 j = 1-mbc, my+mbc
               do 72 ma=1,maux
                 aux1(ma,j) = aux(ma,i-1,j)
                 aux2(ma,j) = aux(ma,i,  j)
                 aux3(ma,j) = aux(ma,i+1,j)
   72            continue
           endif
c
c
c        # Store the value of i along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         icom = i
c
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(2,maxm,meqn,maux,mbc,my,
     &              q1d,dtdy1d,aux1,aux2,aux3,
     &              faddm,faddp,gaddm,gaddp,cfl1d,
     &              work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &              work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2)
c
         cflgrid = dmax1(cflgrid,cfl1d)
c
c        #
c        # update fluxes for use in AMR:
c
         do 75 m=1,meqn
            do 75 j=1,my+1
               gm(m,i,j)   = gm(m,i,j)  + faddm(m,j)
               gp(m,i,j)   = gp(m,i,j)  + faddp(m,j)
               fm(m,i,j)   = fm(m,i,j)  + gaddm(m,j,1)
               fp(m,i,j)   = fp(m,i,j)  + gaddp(m,j,1)
               fm(m,i+1,j) = fm(m,i+1,j)+ gaddm(m,j,2)
               fp(m,i+1,j) = fp(m,i+1,j)+ gaddp(m,j,2)
   75          continue
  100    continue

c     # relimit correction fluxes if they drive a cell negative
         if (relimit) then
         dtdxij=dtdx
         dtdyij=dtdy
         do 101 i=1,mx
            do 101 j=1,my
                  if (mcapa.gt.0) then
                     dtdxij=dtdx/aux(mcapa,i,j)
                     dtdyij=dtdy/aux(mcapa,i,j)
                  endif

                  p   = dmax1(0.d0,dtdxij*fm(1,i+1,j))
     &                + dmax1(0.d0,dtdyij*gm(1,i,j+1))
     &                - dmin1(0.d0,dtdxij*fp(1,i,j))
     &                - dmin1(0.d0,dtdyij*gp(1,i,j))
                  phi = dmin1(1.d0,
     &                 dabs(qold(1,i,j)/(p+drytolerance)))

                  if (phi.lt.1.d0) then
                     do m=1,meqn
                        if (fp(1,i,j).lt.0.d0) then
                           cm=fp(m,i,j)-fm(m,i,j)
                           fm(m,i,j)=phi*fm(m,i,j)
                           fp(m,i,j)=fm(m,i,j)+cm
                        endif
                        if (gp(1,i,j).lt.0.d0) then
                           cm=gp(m,i,j)-gm(m,i,j)
                           gm(m,i,j)=phi*gm(m,i,j)
                           gp(m,i,j)=gm(m,i,j)+cm
                        endif
                        if (fm(1,i+1,j).gt.0.d0) then
                           cm=fp(m,i+1,j)-fm(m,i+1,j)
                           fp(m,i+1,j)=phi*fp(m,i+1,j)
                           fm(m,i+1,j)=fp(m,i+1,j)-cm
                        endif
                        if (gm(1,i,j+1).gt.0.d0) then
                           cm=gp(m,i,j+1)-gm(m,i,j+1)
                           gp(m,i,j+1)=phi*gp(m,i,j+1)
                           gm(m,i,j+1)=gp(m,i,j+1)-cm
                        endif
                     enddo
                  endif
 101     continue
         endif

c
      return
      end
