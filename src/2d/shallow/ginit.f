c
c -------------------------------------------------------------
c
      subroutine ginit(msave, first, nvar, naux, start_time)
c
      use amr_module
      implicit double precision (a-h,o-z)
      logical first


c ::::::::::::::::::::::::::::: GINIT ::::::::::::::::::::::::
c
c  initializes soln on all grids at 'level'  by calling qinit
c  if first = true, (first call to init), then allocate the
c  soln storage area too, else was already allocated.
c
c :::::::::::::::::::::::::::::::::::::::;::::::::::::::::::::

      if (msave .eq. 0) go to 99

      level = node(nestlevel,msave)
      hx    = hxposs(level)
      hy    = hyposs(level)
      mptr  = msave
 
 10       nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          mitot   = nx + 2*nghost
          mjtot   = ny + 2*nghost
          corn1   = rnode(cornxlo,mptr)
          corn2   = rnode(cornylo,mptr)
          if(.not. (first)) go to 20
              loc                 = igetsp(mitot*mjtot*nvar)
              node(store1,mptr)   = loc
              if (naux .gt. 0) then
                locaux              = igetsp(mitot*mjtot*naux)
                do k = 1, mitot*mjtot*naux, naux  ! only set first component of aux to signal
                   alloc(locaux+k-1) = NEEDS_TO_BE_SET ! new system checks this val before setting
                end do
                
                call setaux(nghost,nx,ny,corn1,corn2,hx,hy,
     &                    naux,alloc(locaux))
              else 
                locaux = 1
              endif
              node(storeaux,mptr) = locaux
              if (level .lt. mxnest) then
                loc2              = igetsp(mitot*mjtot*nvar)
                node(store2,mptr) = loc2
              endif
              rnode(timemult, mptr) = start_time
              go to 30
 20       continue
c
c  if 2nd time through, put initial values in store2 so finer grids
c  can be advanced with interpolation of their boundary values.
c  new time soln should still be in location store1.
c
          loc     = node(store2,mptr)
          locaux  = node(storeaux,mptr)
c
   30     continue
          call qinit(nvar,nghost,nx,ny,corn1,corn2,hx,hy,
     &               alloc(loc),naux,alloc(locaux))

c
          mptr  = node(levelptr, mptr)
      if (mptr .ne. 0) go to 10
c
c
 99   continue
      return
      end
