c :::::::::::::::::::::::::::: RELOAD ::::::::::::::::::::::::::::::::
c read back in the fort.b*, fort.q*, and fort.t* files
c
c Note: This assumes that the binary output format was used
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c ---------------------------------------------------------
c
      subroutine reload(adjfile, k)
c
      use adjoint_module

      implicit double precision (a-h,o-z)

      integer, intent(in) :: k
      integer :: mptr, level, ladjfile
      integer :: mitot, mjtot, i1, i2
      integer :: i,j, ivar, z, loc,m
      integer :: allocsize, new_size
      character(len=*) :: adjfile
      logical foundFile, initial

      real(kind=8), allocatable, target, dimension(:) :: new_storage

      iadd(ivar,i,j)  = adjoints(k)%loc(mptr)
     .        + ivar - 1 + (adjoints(k)%meqn)*((j-1)*mitot+i-1)

c     Initializing all levels to zero
      adjoints(k)%gridlevel(:) = 0

c     ! Checking to see if fort.t file exists
      ladjfile = len(trim(adjfile))
      adjfile(ladjfile-4:ladjfile-4) = 't'
      write(6,*) 'Attempting to reload data '
      write(6,*) '  fort.t* file: ',trim(adjfile)
      inquire(file=trim(adjfile),exist=foundFile)
      if (.not. foundFile) then
        write(*,*)" Did not find fort.t* file!"
        stop
      endif
      open(9,file=trim(adjfile),status='old',form='formatted')
      rewind 9

c     ! Reading from fort.t file
      read(9, "(e18.8)") adjoints(k)%time
      read(9, "(i6)") adjoints(k)%meqn
      read(9, "(i6)") adjoints(k)%ngrids
      read(9, "(i6)") adjoints(k)%naux
      read(9, "(i6)") adjoints(k)%ndim
      read(9, "(i6)") adjoints(k)%nghost

      close(9)

c     ! Allocating memory for alloc array
      allocsize = 4000000
      allocate(adjoints(k)%alloc(allocsize))

c     ! Checking to see if fort.q file exists
      adjfile(ladjfile-4:ladjfile-4) = 'q'
      write(6,*) 'Attempting to reload data '
      write(6,*) '  fort.q* file: ',trim(adjfile)
      inquire(file=trim(adjfile),exist=foundFile)
      if (.not. foundFile) then
        write(*,*)" Did not find fort.q* file!"
        stop
      endif
      open(10,file=trim(adjfile),status='old',form='formatted')
      rewind 10

c     ! Checking to see if fort.b file exists
      adjfile(ladjfile-4:ladjfile-4) = 'b'
      write(6,*) 'Attempting to reload data '
      write(6,*) '  fort.b* file: ',trim(adjfile)
      inquire(file=trim(adjfile),exist=foundFile)
      if (.not. foundFile) then
         write(*,*)" Did not find fort.b* file!"
         stop
      endif
      open(20,file=trim(adjfile),status='unknown',access='stream')
      rewind 20

c     ! Reading from fort.q* file and fort.b* files
      loc = 1
      do z = 1, adjoints(k)%ngrids
        read(10,"(i6)") mptr
        adjoints(k)%gridpointer(z) = mptr

        read(10,"(i6)") level
        adjoints(k)%gridlevel(mptr) = level

        read(10,"(i6)") adjoints(k)%ncellsx(mptr)
        read(10,"(i6)") adjoints(k)%ncellsy(mptr)
        read(10,"(e18.8)") adjoints(k)%xlowvals(mptr)
        read(10,"(e18.8)") adjoints(k)%ylowvals(mptr)
        read(10,"(e18.8)") adjoints(k)%hxposs(level)
        read(10,"(e18.8)") adjoints(k)%hyposs(level)
        read(10,*)

        mitot = adjoints(k)%ncellsx(mptr) + 2*adjoints(k)%nghost
        mjtot = adjoints(k)%ncellsy(mptr) + 2*adjoints(k)%nghost

        adjoints(k)%loc(mptr) = loc
        loc = loc + mitot*mjtot*(adjoints(k)%meqn)

c       Checking to see if the alloc array is large enough
c       to hold the new grid
c       If not, making the alloc array larger
        if (allocsize .lt. loc) then
           new_size = 2*allocsize
           allocate(new_storage(new_size))

           new_storage(1:allocsize) = adjoints(k)%alloc
           call move_alloc(new_storage,adjoints(k)%alloc)
           allocsize = new_size
        endif

c       ! This is the bulk of the reading
        do j=1,mjtot
            do i=1,mitot
                read(20) adjoints(k)%alloc(iadd(1,i,j)),
     .        adjoints(k)%alloc(iadd(2,i,j)),
     .        adjoints(k)%alloc(iadd(3,i,j)),
     .        adjoints(k)%alloc(iadd(4,i,j))
            enddo
        enddo

      enddo

      close(10)
      close(20)

      adjoints(k)%lfine = maxval(adjoints(k)%gridlevel)

      return
      end
