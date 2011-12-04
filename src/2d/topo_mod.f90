! ============================================================================
!  Program:     /Users/mandli/src/claw44/branches/ktm-geoclaw-mod/2dxy/lib
!  File:        topo_mod
!  Created:     2010-04-22
!  Author:      Kyle Mandli
! ============================================================================
!      Copyright (C) 2010-04-22 Kyle Mandli <mandli@amath.washington.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================
!  Module for topography data
! ============================================================================
module topo_module

    implicit none

    ! Work array
    double precision, allocatable :: topowork(:)

    ! Topography file data
    character*150, allocatable :: topofname(:)
    integer :: mtopofiles,mtoposize
    double precision, allocatable :: xlowtopo(:), ylowtopo(:), tlowtopo(:)
    double precision, allocatable :: xhitopo(:), yhitopo(:), thitopo(:)
    double precision, allocatable :: dxtopo(:), dytopo(:)
    integer, allocatable ::  mxtopo(:), mytopo(:)

    integer, allocatable :: i0topo(:), mtopo(:), mtopoorder(:)
    integer, allocatable ::  minleveltopo(:), maxleveltopo(:), itopotype(:)

    ! Moving topography support
    integer :: imovetopo

contains

    ! ========================================================================
    ! Read topography files as specified in settopo.data
    !
    ! Each topography file has a type stored in topotype(i).
    !   topotype = 1:  standard GIS format: 3 columns: lon,lat,height(m)
    !   topotype = 2:  Header as in DEM file, height(m) one value per line
    !   topotype = 3:  Header as in DEM file, height(m) one row per line
    ! For other formats modify readtopo routine.
    !
    ! advancing northwest to northeast then from north to south. Values should
    ! be uniformly spaced.
    !
    ! Finest value of topography in a given region will be used for
    ! computation
    ! ========================================================================
    subroutine set_topo(fname)

        use geoclaw_module

        implicit none

        ! Input arguments
        character*25, intent(in), optional :: fname

        ! Locals
        integer, parameter :: iunit = 7
        integer :: i,j,itopo,finer_than,rank
        double precision :: area_i,area_j,x_junk,y_junk
        character*25 :: file_name
        logical :: found_file

        ! Open and begin parameter file output
        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'SETTOPO:'
        write(GEO_PARM_UNIT,*) '---------'


        if (present(fname)) then
            file_name = fname
        else
            file_name  = 'settopo.data'
        endif
        inquire(file=file_name,exist=found_file)
        if (.not. found_file) then
            print *, 'You must provide a file ', file_name
            stop
        endif

        call opendatafile(iunit, file_name)

        read(iunit,*) mtopofiles

        if (mtopofiles == 0) then
            write(GEO_PARM_UNIT,*) '   mtopofiles = 0'
            write(GEO_PARM_UNIT,*) '   No topo files specified, '
            write(GEO_PARM_UNIT,*) '          will set B(x,y) = 0 in setaux'
            return
        endif

        write(GEO_PARM_UNIT,*) '   mtopofiles = ',mtopofiles

        ! Read and allocate data parameters for each file
        allocate(mxtopo(mtopofiles),mytopo(mtopofiles))
        allocate(xlowtopo(mtopofiles),ylowtopo(mtopofiles))
        allocate(tlowtopo(mtopofiles),xhitopo(mtopofiles),yhitopo(mtopofiles))
        allocate(thitopo(mtopofiles),dxtopo(mtopofiles),dytopo(mtopofiles))
        allocate(topofname(mtopofiles),itopotype(mtopofiles))
        allocate(minleveltopo(mtopofiles),maxleveltopo(mtopofiles))
        allocate(i0topo(mtopofiles),mtopo(mtopofiles),mtopoorder(mtopofiles))

        do i=1,mtopofiles
            read(iunit,*) topofname(i)
            read(iunit,*) itopotype(i),minleveltopo(i), maxleveltopo(i), &
                tlowtopo(i),thitopo(i)

            write(GEO_PARM_UNIT,*) '   '
            write(GEO_PARM_UNIT,*) '   ',topofname(i)
            write(GEO_PARM_UNIT,*) '  itopotype = ', itopotype(i)
            write(GEO_PARM_UNIT,*) '  minlevel, maxlevel = ', &
                minleveltopo(i), maxleveltopo(i)
            write(GEO_PARM_UNIT,*) '  tlow, thi = ', tlowtopo(i),thitopo(i)
            if (abs(itopotype(i)) == 1) then
                print *, 'WARNING: topotype 1 has been deprecated'
                print *, 'converting to topotype > 1 is encouraged'
                print *, 'python tools for converting files are provided'
            endif
            call read_topo_header(topofname(i),itopotype(i),mxtopo(i), &
                mytopo(i),xlowtopo(i),ylowtopo(i),xhitopo(i),yhitopo(i), &
                dxtopo(i),dytopo(i))
            mtopo(i) = mxtopo(i)*mytopo(i)
        enddo

        ! Indexing into work array
        i0topo(1)=1
        if (mtopofiles > 1) then
            do i=2,mtopofiles
                i0topo(i)=i0topo(i-1) + mtopo(i-1)
            enddo
        endif

        ! Read and allocate topography for each file
        mtoposize = sum(mtopo)
        allocate(topowork(mtoposize))

        do i=1,mtopofiles
            call read_topo(mxtopo(i),mytopo(i),itopotype(i),topofname(i), &
                topowork(i0topo(i):i0topo(i)+mtopo(i)-1))
        enddo

        ! topography order...This determines which order to process topography
        !
        ! The finest topography will be given priority in any region
        ! mtopoorder(rank) = i means that i'th topography file has rank rank,
        ! where the file with rank=1 is the finest and considered first.

        do i=1,mtopofiles
            finer_than = 0
            do j=1,mtopofiles
                if (j /= i) then
                    area_i=dxtopo(i)*dytopo(i)
                    area_j=dxtopo(j)*dytopo(j)
                    if (area_i < area_j) finer_than = finer_than + 1
!                   if two files have the same resolution, order is arbitrarily chosen
                    if ((area_i == area_j).and.(j < i)) then
                        finer_than = finer_than + 1
                    endif
                endif
            enddo
!         # ifinerthan tells how many other files i is finer than
            rank = mtopofiles - finer_than
            mtopoorder(rank) = i
        enddo

        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '  Ranking of topography files', &
            '  finest to coarsest: ', &
            (mtopoorder(rank),rank=1,mtopofiles)
        write(GEO_PARM_UNIT,*) ' '

    end subroutine set_topo

    ! ========================================================================
    !  read_topo(mx,my,dx,dy,xlow,xhi,ylow,yhi,itopo,fname,topo_type)
    !
    !  Read topo file.
    !  New feature: topo_type < 0 means z values need to be negated.
    ! ========================================================================
    subroutine read_topo(mx,my,topo_type,fname,topo)

        use geoclaw_module

        implicit none

        ! Arguments
        integer, intent(in) :: mx,my,topo_type
        character*150, intent(in) :: fname
        double precision, intent(inout) :: topo(1:mx*my)

        ! Locals
        integer, parameter :: iunit = 19, miss_unit = 17
        double precision, parameter :: topo_missing = -150.d0
        logical, parameter :: maketype2 = .false.
        integer :: i,j,num_points,missing,status,topo_start
        double precision :: no_data_value,x,y,z

        print *, ' '
        print *, 'Reading topography file  ', fname

        open(unit=iunit, file=fname, status='unknown',form='formatted')

        select case(abs(topo_type))
            ! ASCII file with x,y,z values on each line.
            ! (progressing from upper left corner across rows, then down)
            ! Assumes a uniform rectangular grid of data values.
            case(1)
                i = 0
                status = 0
                do while (status == 0)
                    i = i + 1
                    read(iunit,fmt=*,iostat=status) x,y,topo(i)
                enddo

            ! ================================================================
            ! ASCII file with header followed by z data
            ! (progressing from upper left corner across rows, then down)
            ! one value per line if topo_type=2 or
            ! mx values per line if topo_type=3
            ! ================================================================
            case(2:3)
                ! Read header
                do i=1,5
                    read(iunit,*)
                enddo
                read(iunit,*) no_data_value

                ! Read in data
                missing = 0
                select case(abs(topo_type))
                    case(2)
                        do i=1,mx*my
                            read(iunit,*) topo(i)
                            if (topo(i) == no_data_value) then
                                missing = missing + 1
                                topo(i) = topo_missing
                            endif
                        enddo
                    case(3)
                        do j=1,my
                            read(iunit,*) (topo((j-1)*mx + i),i=1,mx)
                            do i=1,mx
                                if (topo((j-1)*mx + i) == no_data_value) then
                                    missing = missing + 1
                                    topo((j-1)*mx + i) = topo_missing
                                endif
                            enddo
                        enddo
                end select

                ! Write a warning if we found and missing values
                if (missing > 0)  then
                    print *, '   WARNING...some missing data values this file'
                    print *, '       ',missing,' missing data values'
                    print *, '              (see fort.missing)'
                    print *, '   These values have arbitrarily been set to ',&
                        topo_missing
                endif
        end select

        close(unit=iunit)

        ! Handle negative topo types
        if (topo_type < 0) then
            forall(i=1:mx*my)
                topo(i) = -topo(i)
            end forall
        endif

        ! ====================================================================
        ! when topo_type=1 and data has x,y,z columns,
        ! set maketype2 to true to create a file new.topo_type2 with only z
        ! values, so next time it will take less time to read in.
        ! only works if dx = dy.
!         if ((topo_type == 1).and.maketype2) then
!             open(unit=29,file='new.tt2',status='unknown',form='formatted')
!             write(29,*) mx, '       mx'
!             write(29,*) my, '       my'
!             write(29,*) xll, '     xllcorner'
!             write(29,*) yll, '     yllcorner'
!             write(29,*) dx, '     cellsize'
!             write(29,*) -9999, '     nodataval'
!             do i=topo_start,itopo
!                 write(29,'(d18.8)') topowork(i)
!             enddo
!             close(unit=29)
!             print *, 'New topo file new.tt2 created'
!             if (dx.ne.dy) then
!                 print *,  ' ** Warning, dx must equal dy for this'
!                 print *,  ' dx = ',dx,'  dy = ',dy
!             endif
!         endif
!         ! ====================================================================

    end subroutine read_topo

    ! ========================================================================
    ! subroutine read_topo_header(fname,topo_type,mx,my,xll,yll,xhi,yhi,dx,dy)
    ! ========================================================================
    !  Read topo file header to determine space needed in allocatable array
    !
    !  :Input:
    !   - fname - (char) Name of file
    !   - topo_type - (int) Type of topography file (-3 < topo_type < 3)
    !
    !  :Output:
    !   - mx,my - (int) Number of grid points
    !   - xll,yll,xhi,yhi - (float) Lower and upper coordinates for grid
    !   - dx,dy - (float) Spatial resolution of grid
    ! ========================================================================
    subroutine read_topo_header(fname,topo_type,mx,my,xll,yll,xhi,yhi,dx,dy)

        use geoclaw_module

        implicit none

        ! Input and Output
        character*150, intent(in) :: fname
        integer, intent(in) :: topo_type
        integer, intent(out) :: mx,my
        double precision, intent(out) :: xll,yll,xhi,yhi,dx,dy

        ! Local
        integer, parameter :: iunit = 19
        integer :: topo_size, status
        double precision :: x,y,z,nodata_value
        logical :: found_file

        inquire(file=fname,exist=found_file)
        if (.not. found_file) then
            print *, 'Missing topography file:'
            print *, '   ', fname
            stop
        endif

        open(unit=iunit, file=fname, status='unknown',form='formatted')

        select case(abs(topo_type))
            ! ASCII file with 3 columns
            ! determine data size
            case(1)
                ! Initial size variables
                topo_size = 0
                mx = 0

                ! Read in first values, determines xlow and yhi
                read(iunit,*) xll,yhi
                topo_size = topo_size + 1
                mx = mx + 1

                ! Go through first row figuring out mx, continue to count
                y = yhi
                do while (yhi == y)
                    read(iunit,*) x,y,z
                    topo_size = topo_size + 1
                    mx = mx + 1
                enddo
                mx = mx - 1
                ! Continue to count the rest of the lines
                status = 0
                do while (status == 0)
                    read(iunit,fmt=*,iostat=status) x,y,z
                    topo_size = topo_size + 1
                enddo
                if (status > 0) then
                    print *, "IO error occured in ",fname,", aborting!"
                    stop
                endif

                ! Calculate remaining values
                my = topo_size / mx
                xhi = x
                yll = y
                dx = (xhi-xll) / (mx-1)
                dy = (yhi-yll) / (my-1)

            ! ASCII file with header followed by z data
            case(2:3)
                read(iunit,*) mx
                read(iunit,*) my
                read(iunit,*) xll
                read(iunit,*) yll
                read(iunit,*) dx
                read(iunit,*) nodata_value
                dy = dx
                xhi = xll + (mx-1)*dx
                yhi = yll + (my-1)*dy

            case default
                print *, 'ERROR:  Unrecognized topo_type'
                print *, '    topo_type = ',topo_type
                print *, '  for topography file:'
                print *, '   ', fname
                stop
        end select

        close(iunit)

        write(GEO_PARM_UNIT,*) '  mx = ',mx,'  x = (',xll,',',xhi,')'
        write(GEO_PARM_UNIT,*) '  my = ',my,'  y = (',yll,',',yhi,')'
        write(GEO_PARM_UNIT,*) '  dx, dy (meters/degrees) = ', dx,dy

    end subroutine read_topo_header

end module topo_module
