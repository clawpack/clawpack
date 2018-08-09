c
c --------------------------------------------------------------
c Modified from errf1.f in amrclaw to use aux arrays
c and correctly account for wet/dry cells.
c --------------------------------------------------------------
c
      subroutine errf1(rctfine,nvar,rctcrse,mptr,mi2tot,mj2tot,
     2                 mitot,mjtot,rctflg,mibuff,mjbuff,auxfine,
     2                 naux,auxcrse)

      use amr_module, only: rnode,node,hxposs,hyposs,possk,tol,nghost
      use amr_module, only: iorder,cornylo,cornxlo,edebug,eprint
      use amr_module, only: outunit,timemult,nestlevel,mcapa
      use amr_module, only: DONTFLAG,DOFLAG,UNSET
      use refinement_module, only: wave_tolerance
      use geoclaw_module, only:dry_tolerance, sea_level
      implicit none

      integer, intent(in) :: nvar, mptr, mi2tot, mj2tot
      integer, intent(in) :: mitot, mjtot,mibuff,mjbuff, naux
      real(kind=8), intent(in) :: rctfine(nvar,mitot,mjtot)
      real(kind=8), intent(inout) :: rctcrse(nvar,mi2tot,mj2tot)
      real(kind=8), intent(inout) :: auxfine(naux,mitot,mjtot)
      real(kind=8), intent(in) :: auxcrse(naux,mi2tot,mj2tot)
      real(kind=8), intent(inout) :: rctflg(mibuff,mjbuff)

      logical :: allowflag
      external allowflag

c     Local variables
      integer :: i, j, ifine, jfine, jj, ii, levm, m
      integer :: ico, jco, nwet
      real(kind=8) :: errmax, err2, order
      real(kind=8) :: hx, hy, dt, rflag, time
      real(kind=8) :: xofi, ybot, xleft, yofj, capa, capacrse
      real(kind=8) :: hcrse, hucrse, hvcrse, etacrse, bcrse
      real(kind=8) :: hf, huf, hvf, etaf, bf, etac1
      real(kind=8) :: hc, huc, hvc, etac, hav, etaav, bav
      real(kind=8) :: etasum, hsum, husum, hvsum, bsum
      real(kind=8) :: etaerr, herr, huerr, hverr
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
      do 35  j = nghost+1, mj2tot-nghost
          yofj  = ybot + (dble(jfine) - .5d0)*hy
          ifine = nghost+1
c
          do 30  i  = nghost+1, mi2tot-nghost
              rflag = DONTFLAG

! Only check errors if flag hasn't been set yet.
! If flag == DONTFLAG then refinement is forbidden by a region,
! if flag == DOFLAG checking is not needed

c Note: here rctcrse is being used as a temporary flag
c the fine grid amrflags array is stored in rctflg, and will be
c updated based on rctcrse at the end of this routine
           if(rctflg(ifine,jfine) == UNSET 
     .        .or. rctflg(ifine+1,jfine) == UNSET 
     .        .or. rctflg(ifine,jfine+1) == UNSET 
     .        .or. rctflg(ifine+1,jfine+1) == UNSET) then

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
              bcrse = auxcrse(1,i,j)
              etacrse = hcrse+bcrse

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

c                    hc=max(etaav-bcrse*capacrse,0.d0) !tsunamiclaw method
                    hc = min(hav,(max(etaav-bcrse*capacrse,0.d0)))
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
                etac = hc + bcrse

c               # Finding errors
                herr  = dabs((hc-hcrse)/ order)
                huerr = dabs((huc-hucrse)/ order)
                hverr = dabs((hvc-hvcrse)/ order)
                etaerr= dabs((etac-etacrse)/ order)
              endif

c             Flag point if error is large enough
              if (etaerr .ge. tol) then
                  rflag  = DOFLAG
              endif

            endif

            rctcrse(1,i,j) = rflag

            ifine = ifine + 2
 30       continue

          jfine = jfine + 2
 35    continue

c
c  print out intermediate flagged rctcrse (for debugging)
c
      if (eprint) then
         err2 = dsqrt(err2/dble((mi2tot-2*nghost)*(mj2tot-2*nghost)))
         write(outunit,103) mptr, levm, errmax, err2
 103     format(' grid ',i4,' level ',i4,
     .          ' max. error = ',e15.7,' err2 = ',e15.7)
         if (edebug) then
           write(outunit,*) ' flagged points on coarsened grid ',
     .                      '(no ghost cells) for grid ',mptr
           do 45 jj = nghost+1, mj2tot-nghost
              j = mj2tot + 1 - jj
              write(outunit,106) (nint(rctcrse(1,i,j)),
     .                            i=nghost+1,mi2tot-nghost)
106           format(1h ,80i1)
45         continue
         endif
      endif

c
c     If coarse point is flagged, set point in flag array
c
      jfine   = nghost+1
      do 70 j = nghost+1, mj2tot-nghost
          ifine   = nghost+1
          do 60 i = nghost+1, mi2tot-nghost
             if (rctcrse(1,i,j) .eq. DOFLAG) then
c                ## never set rctflg to DONTFLAG, since flag2refine or
c                ## flag2refine may have previously set it to DOFLAG
c                ## can only add DOFLAG pts in this routine
                 rctflg(ifine,jfine)    = DOFLAG
                 rctflg(ifine+1,jfine)  = DOFLAG
                 rctflg(ifine,jfine+1)  = DOFLAG
                 rctflg(ifine+1,jfine+1)= DOFLAG
             endif

             ifine   = ifine + 2
 60        continue

           jfine   = jfine + 2
 70    continue

c
c     Print out values for debugging
c
      if (eprint) then
         write(outunit,118)
 118     format(' on fine grid (no ghost cells) flagged points are')
         if (edebug) then
          do 56 jj = nghost+1, mjtot-nghost
           j = mjtot + 1 - jj
           write(outunit,106)
     &      (nint(rctflg(i,j)),i=nghost+1,mitot-nghost)
 56       continue
        endif
      endif

      return
      end
