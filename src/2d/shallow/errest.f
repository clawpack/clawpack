c
c -------------------------------------------------------------
c Modified from errest.f in amrclaw to pass aux arrays into errf1.f
c -------------------------------------------------------------
c
      subroutine errest (nvar,naux,lcheck,mptr,nx,ny)
c
      use amr_module, only: node,store1,store2,storeaux,alloc
      use amr_module, only: nghost,tempptr,storeflags,ibuff
      use adjoint_module, only : adjoint_flagging
      use adjointsup_module, only: errf1a
      implicit none

      integer, intent(in) :: nvar, naux, lcheck, mptr, nx, ny

c
c   ### changed to stack based storage 2/23/13
c   ### and broken into smaller routines to minimize 
c   ### stack space

      real(kind=8) :: valbgc(nvar,nx/2+2*nghost,ny/2+2*nghost)
      real(kind=8) :: auxbgc(naux,nx/2+2*nghost,ny/2+2*nghost)

c   # Other local variables
      integer :: locamrflags,locaux,locbig,locnew,locold
      integer :: mbuff,mi2tot,mibuff,midub,mitot
      integer :: mj2tot,mjbuff,mjdub,mjtot

c :::::::::::::::::::::::::: ERREST :::::::::::::::::::::::::::::::::::
c for this grid at level lcheck:
c  estimate the error by taking a large (2h,2k) step based on the
c  values in the old storage loc., then take one regular (and for
c  now wasted) step based on the new info.   compare using an
c  error relation for a pth order  accurate integration formula.
c  flag error plane as either bad (needs refinement), or good.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
       mitot  = nx + 2*nghost
       mjtot  = ny + 2*nghost
       locnew = node(store1,mptr)
       locold = node(store2,mptr)
       locaux = node(storeaux,mptr)
       mi2tot = nx/2  + 2*nghost
       mj2tot = ny/2  + 2*nghost
c
c     prepare double the stencil size worth of boundary values,
c            then coarsen them for the giant step integration.
       midub = nx+4*nghost
       mjdub = ny+4*nghost
c
       call prepbigstep(nvar,naux,lcheck,mptr,nx,ny,midub,
     .                    mjdub,valbgc,auxbgc,mi2tot,mj2tot)

c
c  the one giant step based on old values is done. now take
c  one regular step based on new values. 
c  boundary values already in locbig, (see subr. flagger)
c
      locbig = node(tempptr,mptr)
      locaux = node(storeaux,mptr)
      call prepregstep(nvar,naux,lcheck,mptr,nx,ny,mitot,mjtot,
     .                 alloc(locbig),alloc(locaux))
c
c     ## locamrflags allocated in flagger. may previously have been used 
c     ## by flag2refine so make sure not to overwrite
      locamrflags = node(storeflags, mptr)    
      mbuff = max(nghost,ibuff+1)  
      mibuff = nx + 2*mbuff 
      mjbuff = ny + 2*mbuff
      if (adjoint_flagging) then
          call errf1a(alloc(locbig),nvar,valbgc,mptr,mi2tot,mj2tot,
     1           mitot,mjtot,alloc(locamrflags),mibuff,mjbuff,
     1           alloc(locaux),naux,auxbgc)
      else
          call errf1(alloc(locbig),nvar,valbgc,mptr,mi2tot,mj2tot,
     1           mitot,mjtot,alloc(locamrflags),mibuff,mjbuff,
     1           alloc(locaux),naux,auxbgc)
      endif

c
      return
      end
