c
c -----------------------------------------------------
c
      subroutine valout (lst, lend, time, nvar, naux)
c
      use amr_module

      use storm_module, only: wind_index, pressure_index
      use storm_module, only: output_storm_location

      use friction_module, only: friction_index

      implicit double precision (a-h,o-z)
      character*10  matname1, matname2, matname3

      real(kind=8) :: h, hu, hv, eta
      integer, parameter :: NUM_FIELDS = 6
      real(kind=8) :: storm_field(NUM_FIELDS)

c     # Output the results for a general system of conservation laws
c     # in 2 dimensions
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw2.m or Python tools
c
c     # set outaux = .true. to also output the aux arrays to fort.a<iframe>

      logical outaux

      iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
      iaddaux(iaux,i,j) = locaux + iaux-1 + naux*(i-1) +
     .                                      naux*mitot*(j-1)
c
      outaux = .false.

c     ASCII Output
      if (output_format == 1) then
         ! Output the storm location at this time
         call output_storm_location(time)

c        ###  make the file names and open output files
         matname1 = 'fort.qxxxx'
         matname2 = 'fort.txxxx'
         matname3 = 'fort.axxxx'
         matunit1 = 50
         matunit2 = 60
         matunit3 = 70
         nstp     = matlabu
         do 55 ipos = 10, 7, -1
            idigit = mod(nstp,10)
            matname1(ipos:ipos) = char(ichar('0') + idigit)
            matname2(ipos:ipos) = char(ichar('0') + idigit)
            matname3(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
 55      continue
         open(unit=matunit1,file=matname1,status='unknown',
     .       form='formatted')

         level = lst
         ngrids = 0
c65      if (level .gt. lfine) go to 90
 65      if (level .gt. lend) go to 90
            mptr = lstart(level)
 70         if (mptr .eq. 0) go to 80
              ngrids  = ngrids + 1
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              if (ny.gt.1) then
                  write(matunit1,1001) mptr, level, nx, ny
                else
c                 # output in 1d format if ny=1:
                  write(matunit1,1003) mptr, level, nx
                endif
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx',/,
     &       i5,'                 my')
 1003 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx')


              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              if (ny.gt.1) then
                  write(matunit1,1002)
     &              xlow,ylow,hxposs(level),hyposs(level)
                else
                  write(matunit1,1004)
     &              xlow,hxposs(level)
                endif
 1002 format(e18.8,'    xlow', /,
     &       e18.8,'    ylow', /,
     &       e18.8,'    dx', /,
     &       e18.8,'    dy',/)
 1004 format(e18.8,'    xlow', /,
     &       e18.8,'    dx', /)


        do j = nghost+1, mjtot-nghost
          do i = nghost+1, mitot-nghost
            do ivar=1,nvar
              if (abs(alloc(iadd(ivar,i,j))) < 1d-90) then
                alloc(iadd(ivar,i,j)) = 0.d0
              endif
            enddo

            ! Extract all but bottom layer depth and momenta
            h = alloc(iadd(1,i,j))
            hu = alloc(iadd(2,i,j))
            hv = alloc(iadd(3,i,j))
            eta = h + alloc(iaddaux(1,i,j))
            if (abs(eta) < 1d-90) then
              eta = 0.d0
            end if

            ! Storm fields and location
            storm_field(1) = alloc(iaddaux(friction_index,i,j))
            ! Actually output the wind field stress term
            storm_field(2) = alloc(iaddaux(wind_index,i,j))
            storm_field(3) = alloc(iaddaux(wind_index+1,i,j))
            storm_field(4) = alloc(iaddaux(pressure_index,i,j))
            storm_field(5) = alloc(iaddaux(pressure_index+1,i,j))
            storm_field(6) = alloc(iaddaux(pressure_index+2,i,j))
            forall (k=1:NUM_FIELDS,abs(storm_field(k)) < 1d-90)
              storm_field(k) = 0.d0
            end forall

            write(matunit1,109) h,hu,hv,eta,
     &                         (storm_field(k),k=1,NUM_FIELDS)
            
          enddo
          write(matunit1,*) ' '
        enddo
  109       format(4e26.16)


            mptr = node(levelptr, mptr)
            go to 70
 80      level = level + 1
         go to 65

 90     continue

        if (outaux) then
c        # output aux array to fort.aXXXX
         open(unit=matunit3,file=matname3,status='unknown',
     .       form='formatted')
         level = lst
         ngrids = 0
 165     if (level .gt. lfine) go to 190
            mptr = lstart(level)
 170        if (mptr .eq. 0) go to 180
              ngrids  = ngrids + 1
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              if (ny.gt.1) then
                  write(matunit3,1001) mptr, level, nx, ny
                else
c                 # output in 1d format if ny=1:
                  write(matunit3,1003) mptr, level, nx
                endif
              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              if (ny.gt.1) then
                  write(matunit3,1002)
     &              xlow,ylow,hxposs(level),hyposs(level)
                else
                  write(matunit3,1004)
     &              xlow,hxposs(level)
                endif

         do j = nghost+1, mjtot-nghost
            do i = nghost+1, mitot-nghost
               do ivar=1,naux
                  if (abs(alloc(iaddaux(ivar,i,j))) .lt. 1d-90) then
                     alloc(iaddaux(ivar,i,j)) = 0.d0
                  endif
               enddo
               write(matunit3,109) (alloc(iaddaux(ivar,i,j)), 
     &                              ivar=1,naux)
            enddo
            write(matunit3,*) ' '
         enddo

            mptr = node(levelptr, mptr)
            go to 170
 180     level = level + 1
         go to 165

 190    continue
        close(unit=matunit3)
        endif !# end outputting aux array

      open(unit=matunit2,file=matname2,status='unknown',
     .       form='formatted')
      if (ny.gt.1) then 
          ndim = 2
        else
c         # special case where 2d AMR is used for a 1d problem
c         # and we want to use 1d plotting routines
          ndim = 1
        endif

      write(matunit2,1000) time,4 + NUM_FIELDS,ngrids,naux,ndim
 1000 format(e18.8,'    time', /,
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 naux'/,
     &       i5,'                 ndim'/,/)
c

      write(6,601) matlabu,time
  601 format('AMRCLAW: Frame ',i4,
     &       ' output files done at time t = ', d13.6,/)

      matlabu = matlabu + 1

      close(unit=matunit1)
      close(unit=matunit2)
      endif

      return
      end
