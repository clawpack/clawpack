c
c --------------------------------------------------------------
c Modified from errf1.f in amrclaw to use aux arrays
c and correctly account for wet/dry cells.
c --------------------------------------------------------------
c
      subroutine errf1(rctfine,nvar,rctcrse,mptr,mi2tot,mj2tot,
     2                 mitot,mjtot,rctflg,mibuff,mjbuff,auxfine,
     2                 naux,auxcrse,nx,ny,mask_selecta)

      use adjoint_module, only: innerprod_index,
     .        totnum_adjoints, adjoints, trange_start, trange_final,
     .        levtol, eptr, errors
      use innerprod_module, only: calculate_innerproduct
      use amr_module, only: rnode,node,hxposs,hyposs,possk,tol,nghost
      use amr_module, only: iorder,cornylo,cornxlo,edebug,eprint,goodpt
      use amr_module, only: outunit,timemult,badpt,nestlevel,mcapa
      use amr_module, only: t0,tfinal,numcells
      use refinement_module, only: wave_tolerance
      use geoclaw_module, only:dry_tolerance, sea_level
      implicit none

      integer, intent(in) :: nvar, mptr, mi2tot, mj2tot,nx,ny
      integer, intent(in) :: mitot, mjtot,mibuff,mjbuff, naux
      real(kind=8), intent(in) :: rctfine(nvar,mitot,mjtot)
      real(kind=8), intent(inout) :: rctcrse(nvar,mi2tot,mj2tot)
      real(kind=8), intent(inout) :: auxfine(naux,mitot,mjtot)
      real(kind=8), intent(in) :: auxcrse(naux,mi2tot,mj2tot)
      real(kind=8) :: aux_temp(1:nx,1:ny)
      real(kind=8), intent(inout) :: rctflg(mibuff,mjbuff)
      real(kind=8) :: est(nvar,mitot,mjtot)
      logical :: mask_selecta(totnum_adjoints)

      logical :: allowflag
      external allowflag

c     Local variables
      integer :: i, j, ifine, jfine, jj, ii, levm, m, k
      integer :: ico, jco, nwet
      real(kind=8) :: errmax, err2, order
      real(kind=8) :: hx, hy, dt, rflag, time
      real(kind=8) :: xofi, ybot, xleft, yofj, capa, capacrse
      real(kind=8) :: hcrse, hucrse, hvcrse, etacrse
      real(kind=8) :: bcrse(mitot,mjtot)
      real(kind=8) :: hf, huf, hvf, etaf, bf, etac1
      real(kind=8) :: hc, huc, hvc, etac, hav, etaav, bav
      real(kind=8) :: etasum, hsum, husum, hvsum, bsum
      real(kind=8) :: etaerr, herr, huerr, hverr, innerprod
      real(kind=8) :: tol_exact
c
c
c ::::::::::::::::::::::::::::: ERRF1 ::::::::::::::::::::::::::::::::
c
c  Richardson error estimator:  Used when flag_richardson is .true.
c  Compare error estimates in rctfine, rctcrse, 
c  A point is flagged if the error estimate is greater than tol
c  later we check if its in a region where its allowed to be flagged
c  or alternatively required.
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

c
c     Local variables
c
      time  = rnode(timemult, mptr)
      xleft = rnode(cornxlo,mptr)
      levm  = node(nestlevel, mptr)
      hx    = hxposs(levm)
      ybot  = rnode(cornylo,mptr)
      hy    = hyposs(levm)
      dt    = possk(levm)
 
      errmax = 0.0d0
      err2   = 0.0d0
      order  = dble(2**(iorder+1) - 2)

      auxfine(innerprod_index,:,:) = 0.0d0

c     Calculating correct tol for this level
c     --------------------
c     Total error allowed in this time step
      tol_exact = tol*dt/tfinal
c     Error allowed at this level
      tol_exact = tol_exact/(2**(levm - 1))
c     Error allowed per cell at this level
      tol_exact = tol_exact/(numcells(levm)*hx*hy)

      if (t0+possk(levm) .eq. time) levtol(levm) = tol_exact

c
c     Print values for debugging
c
      if (.not. (edebug)) go to 20
         write(outunit,107) mptr
 107     format(//,' coarsened grid values for grid ',i4)
         do 10 jj = nghost+1, mj2tot-nghost
            j = mj2tot + 1 - jj
            write(outunit,101) (rctcrse(1,i,j),
     .                          i = nghost+1, mi2tot-nghost)
10       continue
         write(outunit,108) mptr
 108     format(//, ' fine grid values for grid ',i4)
         do 15 jj = nghost+1, mjtot-nghost
            j = mjtot + 1 - jj
            write(outunit,101) (rctfine(1,i,j),i=nghost+1,mitot-nghost)
15       continue
101      format(' ',10e15.7)

c
c Main loop: flag points and set extra aux values
c
 20   continue
      jfine = nghost+1
      do j = nghost+1,mj2tot-nghost
          ifine = nghost+1
          do i = nghost+1,mi2tot-nghost

              bcrse(ifine,jfine) = auxcrse(1,i,j)
              bcrse(ifine+1,jfine) = auxcrse(1,i,j)
              bcrse(ifine+1,jfine+1) = auxcrse(1,i,j)
              bcrse(ifine,jfine+1) = auxcrse(1,i,j)

              ifine = ifine + 2
          enddo
          jfine = jfine + 2
      enddo

      jfine = nghost+1
      do 35  j = nghost+1, mj2tot-nghost
          yofj  = ybot + (dble(jfine) - .5d0)*hy
          ifine = nghost+1
c
          do 30  i  = nghost+1, mi2tot-nghost
              rflag = goodpt
              xofi  = xleft + (dble(ifine) - .5d0)*hx

              herr  = 0.d0
              huerr = 0.d0
              hverr = 0.d0
              etaerr= 0.d0
c
c             Calculate error for each value of q and eta
c             if course grid is in wet state
c
              if (mcapa .eq. 0) then
                  capacrse=1.0d0
              else
                  capacrse=auxcrse(2,i,j)
              endif

              hcrse = rctcrse(1,i,j)
              hucrse= rctcrse(2,i,j)
              hvcrse= rctcrse(3,i,j)

              etacrse = hcrse+bcrse(ifine,jfine)

              if (hcrse > dry_tolerance) then

c               # Calculating which values to use for error
                etasum= 0.d0
                hsum  = 0.d0
                husum = 0.d0
                hvsum = 0.d0
                bsum = 0.d0

                nwet=0

                do jco = 1,2
                  do ico = 1,2
                    if (mcapa .eq. 0) then
                      capa=1.0d0
                    else
                      capa=auxfine(2,ifine+ico-1,jfine+jco-1)
                    endif

                    hf =rctfine(1,ifine+ico-1,jfine+jco-1)*capa
                    huf=rctfine(2,ifine+ico-1,jfine+jco-1)*capa
                    hvf=rctfine(3,ifine+ico-1,jfine+jco-1)*capa

                    bf =auxfine(1,ifine+ico-1,jfine+jco-1)*capa

                    if (hf > dry_tolerance) then
                        etaf = hf+bf
                        nwet=nwet+1
                    else
                        etaf = 0.d0
                        huf=0.d0
                        hvf=0.d0
                        bf = 0.d0
                    endif

                    hsum   = hsum + hf
                    husum  = husum + huf
                    hvsum  = hvsum + hvf
                    etasum = etasum + etaf
                    bsum = bsum + bf
                  enddo
                enddo

c               # Finding coarsened values
                if (nwet .gt. 0) then
                    hav  = hsum/dble(nwet)
                    etaav = etasum/dble(nwet)
                    bav = bsum/dble(nwet)

c                    hc=max(etaav-bcrse(ifine,jfine)*capacrse,0.d0) !tsunamiclaw method
                    hc = min(hav,
     .                   (max(etaav-bcrse(ifine,jfine)*capacrse,0.d0)))
                    if (hsum .ne. 0) then
                        huc = (min(hav,hc)/hsum)*husum
                        hvc = (min(hav,hc)/hsum)*hvsum
                    else
                        huc = 0.d0
                        hvc = 0.d0
                    endif

                else
                    hc  = 0.d0
                    huc = 0.d0
                    hvc = 0.d0
                endif

                hc = hc / capacrse
                huc = huc / capacrse
                hvc = hvc / capacrse
                etac = hc + bcrse(ifine,jfine)

c               # Finding errors
                herr  = dabs((hc-hcrse)/ order)
                huerr = dabs((huc-hucrse)/ order)
                hverr = dabs((hvc-hvcrse)/ order)
                etaerr= dabs((etac-etacrse)/ order)

                est(1,ifine,jfine) = dabs((etac-etacrse)/ order)
                est(2,ifine,jfine) = dabs((huc-hucrse)/ order)
                est(3,ifine,jfine) = dabs((hvc-hvcrse)/ order)

c               retaining directionality of the wave
                est(1,ifine,jfine) = sign(est(1,ifine,jfine),etacrse)
                est(2,ifine,jfine) = sign(est(2,ifine,jfine),hucrse)
                est(3,ifine,jfine) = sign(est(3,ifine,jfine),hvcrse)

                do k = 1,nvar
                    est(k,ifine + 1,jfine) = est(k,ifine,jfine)
                    est(k,ifine + 1,jfine + 1) = est(k,ifine,jfine)
                    est(k,ifine,jfine + 1) = est(k,ifine,jfine)
                enddo
              endif

              ifine = ifine + 2
 30       continue

          jfine = jfine + 2
 35    continue

      do 12 k = 1,totnum_adjoints
c         ! Consider only snapshots that are within the desired time range
          if (mask_selecta(k)) then
c             set innerproduct for fine grid
              aux_temp(:,:) =
     .            calculate_innerproduct(time,est,k,nx,ny,
     .            xleft,ybot,hx,hy,nvar,nghost,bcrse)

              do 22  i  = 1, nx
                  do 23  j  = 1, ny
                     auxfine(innerprod_index,i,j) =
     .                max(auxfine(innerprod_index,i,j),
     .                    aux_temp(i,j))

                     if (auxfine(innerprod_index,i,j) .ge. tol) then
c                     ## never set rctflg to good, since flag2refine may
c                     ## have previously set it to bad
c                     ## can only add bad pts in this routine
                         rctflg(i,j)    = badpt
                     endif

 23           continue
 22           continue
          endif
 12   continue

c
c     Print out values for debugging
c
      if (eprint) then
         write(outunit,118)
 118     format(' on fine grid (no ghost cells) flagged points are')
         if (edebug) then
          do 56 jj = nghost+1, mjtot-nghost
           j = mjtot + 1 - jj
 106       format(1h ,80i1)
           write(outunit,106)
     &      (nint(rctflg(i,j)),i=nghost+1,mitot-nghost)
 56       continue
        endif
      endif

      return
      end
