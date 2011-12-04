c
c -----------------------------------------------------
c
      subroutine valout (lst, lend, time, nvar, naux)
c
      implicit double precision (a-h,o-z)
      character*10  matname1, matname2

      include  "call.i"

c OLD INDEXING
c     iadd(i,j,ivar) = loc + i - 1 + mitot*((ivar-1)*mjtot+j-1)
c     iaddaux(i,j,iaux) = locaux + i - 1 + mitot*((iaux-1)*mjtot+j-1)
c NEW INDEXING ORDER SWITCHED
      iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
      iaddaux(iaux,i,j) = locaux + iaux-1 + naux*(i-1) +
     .                                      naux*mitot*(j-1)



c ::::::::::::::::::::::::::: VALOUT ::::::::::::::::::::::::::::::::::;
c graphics output of soln values for contour or surface plots.
c modified for GeoClaw to output the surface level along with q.
c    surface = q(i,j,1) + aux(i,j,1)
c :::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::::::::::;

c
c     ### MATLAB graphics output
c
      if (matlabout) then
c        ###  make the file names and open output files
         matname1 = 'fort.qxxxx'
         matname2 = 'fort.txxxx'
         matunit1 = 50
         matunit2 = 60
         nstp     = matlabu
         do 55 ipos = 10, 7, -1
            idigit = mod(nstp,10)
            matname1(ipos:ipos) = char(ichar('0') + idigit)
            matname2(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
 55      continue
         open(unit=matunit1,file=matname1,status='unknown',
     .       form='formatted')

         level = lst
         ngrids = 0
 65      if (level .gt. lfine) go to 90
            mptr = lstart(level)
 70         if (mptr .eq. 0) go to 80
              ngrids  = ngrids + 1
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              write(matunit1,1001) mptr, level, nx, ny
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx',/,
     &       i5,'                 my')


c  old        xcorn = rnode(cornxlo,mptr) - .5d0*hxposs(level)
c  old        ycorn = rnode(cornylo,mptr) - .5d0*hyposs(level)
              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              write(matunit1,1002)
     &              xlow,ylow,hxposs(level),hyposs(level)
 1002 format(e18.8,'    xlow', /,
     &       e18.8,'    ylow', /,
     &       e18.8,'    dx', /,
     &       e18.8,'    dy',/)

      do j = nghost+1, mjtot-nghost
         do i = nghost+1, mitot-nghost
            do ivar=1,nvar
               if (dabs(alloc(iadd(ivar,i,j))) .lt. 1d-90) then
                  alloc(iadd(ivar,i,j)) = 0.d0
               endif
            enddo
            surface = alloc(iadd(1,i,j)) + alloc(iaddaux(1,i,j))
            write(matunit1,109) (alloc(iadd(ivar,i,j)), ivar=1,nvar),
     &         surface
         enddo
         write(matunit1,*) ' '
      enddo
  109       format(4e26.16)


            mptr = node(levelptr, mptr)
            go to 70
 80      level = level + 1
         go to 65

 90     continue

        open(unit=matunit2,file=matname2,status='unknown',
     .       form='formatted')

c     # nvar+1 variable printed since surface also printed

      write(matunit2,1000) time,nvar+1,ngrids,3,2
 1000 format(e18.8,'    time', /,
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 naux'/,
     &       i5,'                 ndim'/,/)
c

      write(6,601) matlabu,time
  601 format('GeoClaw: Frame ',i4,
     &       ' output files done at time t = ', d12.6,/)

      matlabu = matlabu + 1

      close(unit=matunit1)
      close(unit=matunit2)
      endif

      return
      end
