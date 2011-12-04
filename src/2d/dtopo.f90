! ============================================================================
!  Program:     /Users/mandli/src/claw44/branches/ktm-geoclaw-mod/2dxy/lib
!  File:        dtopo_mod
!  Created:     2010-04-22
!  Author:      Kyle Mandli
! ============================================================================
!      Copyright (C) 2010-04-22 Kyle Mandli <mandli@amath.washington.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================
!  Module containing changing topography data
! ============================================================================
module dtopo_module

    implicit none

    ! Work array
    double precision, allocatable :: dtopowork(:)

    ! File data parameters
    character*150, allocatable :: dtopofname(:)
    double precision, allocatable :: xlowdtopo(:),ylowdtopo(:),xhidtopo(:)
    double precision, allocatable :: yhidtopo(:),t0dtopo(:),tfdtopo(:)
    double precision, allocatable :: dxdtopo(:),dydtopo(:),dtdtopo(:)
    integer, allocatable :: mxdtopo(:),mydtopo(:),mtdtopo(:),mdtopo(:)
    integer, allocatable :: minleveldtopo(:),maxleveldtopo(:),dtopotype(:)
    integer, allocatable :: i0dtopo(:)
    integer :: num_dtopo
    double precision dz
    logical, allocatable :: topoaltered(:)

contains
    ! ========================================================================
    !  set_dtopo(fname)
    ! ========================================================================
    ! Read moving topography info from setdtopo.data
    ! Time-dependend topography is used to initiate tsunami, for example.
    !
    ! If num_dtopo = 0, no movement of topography specified.
    !
    ! If num_dtopo = 1, the topo changes dynamically.
    ! the topofile can then be formated in the following ways
    ! dtopotype > 1: file contains a header, with a line for each record
    ! my; mx; mt; xlower; ylower; t0; dx; dy; dt;
    ! dtopotype = 1:
    ! Then the dtopofile should have 4 columns:
    ! time,longitude,latitude,vertical displacement(m) since t0
    !
    ! Longitude and latitude advance in the standard GIS way from
    ! upper left corner across in x and then down in y.
    ! Time column advances most slowly.
    ! ========================================================================
    subroutine set_dtopo(fname)

        use geoclaw_module
        use topo_module

        implicit none

        ! Input arguments
        character*25, optional, intent(in) :: fname

        ! Locals
        character*25 :: file_name
        integer, parameter :: iunit = 79
        logical :: found_file
        double precision :: xcell, xim, xip, ycell, yjm, yjp, ztopoij
        double precision :: capac_area, deg2rad
        integer :: i,m,ib,jb,ij,ijdtopo,jbr

        ! Function
        double precision :: topointegral

        ! Common block
        double precision :: tstart
        common /ctstart/ tstart

        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'SETDTOPO:'
        write(GEO_PARM_UNIT,*) '-------------'

        if (present(fname)) then
            file_name = fname
        else
            file_name  = 'setdtopo.data'
        endif
        inquire(file=file_name,exist=found_file)
        if (.not. found_file) then
            print *,'You must provide a file ', file_name
            stop
        endif

        call opendatafile(iunit, file_name)

        read(iunit,*) num_dtopo
        write(GEO_PARM_UNIT,*) '   num dtopo files = ',num_dtopo
        if (num_dtopo == 0) then
            return
        endif

        ! Allocate and read in dtopo info
        allocate(dtopofname(num_dtopo),minleveldtopo(num_dtopo))
        allocate(maxleveldtopo(num_dtopo),mxdtopo(num_dtopo))
        allocate(mydtopo(num_dtopo),mtdtopo(num_dtopo),mdtopo(num_dtopo))
        allocate(xlowdtopo(num_dtopo),ylowdtopo(num_dtopo),t0dtopo(num_dtopo))
        allocate(xhidtopo(num_dtopo),yhidtopo(num_dtopo),tfdtopo(num_dtopo))
        allocate(dtopotype(num_dtopo),i0dtopo(num_dtopo))
        allocate(dxdtopo(num_dtopo),dydtopo(num_dtopo),dtdtopo(num_dtopo))
        allocate(topoaltered(num_dtopo))

        do i=1,num_dtopo
            read(iunit,*) dtopofname(i)
            read(iunit,*) dtopotype(i),minleveldtopo(i), maxleveldtopo(i)
            write(GEO_PARM_UNIT,*) '   fname:',dtopofname(i)
            write(GEO_PARM_UNIT,*) '   topo type:',dtopotype(i)
            write(GEO_PARM_UNIT,*) '   minlevel, maxlevel:'
            write(GEO_PARM_UNIT,*)  minleveldtopo(i),maxleveldtopo(i)

            ! Read in header data
            call read_dtopo_header(dtopofname(i),dtopotype(i),mxdtopo(i), &
                mydtopo(i),mtdtopo(i),xlowdtopo(i),ylowdtopo(i),t0dtopo(i), &
                xhidtopo(i),yhidtopo(i),tfdtopo(i),dxdtopo(i),dydtopo(i), &
                dtdtopo(i))
            mdtopo(i) = mxdtopo(i) * mydtopo(i) * mtdtopo(i)
        enddo

        ! Indexing into work array
        i0dtopo(1) = 1
        if (num_dtopo > 1) then
            do i=2,num_dtopo
                i0dtopo(i) = i0dtopo(i-1) + mdtopo(i-1)
            enddo
        endif

        ! Allocate and read dtopo files
        allocate(dtopowork(sum(mdtopo)))

        do i=1,num_dtopo
            call read_dtopo(mxdtopo(i),mydtopo(i),mtdtopo(i),dtopotype(i), &
                dtopofname(i),dtopowork(i0dtopo(i):i0dtopo(i)+mdtopo(i)-1))
        enddo

        ! ==============Fix topography for a restart==========================
        ! topoaltered is ordinarily initialized as .false.
        ! this variable indicates when topo arrays have been altered
        ! based on a dtopo model.  For restarts after the end of motion,
        ! the topo arrays are altered here to match the final topo + dtopo.
        ! ====================================================================
        do i=1,num_dtopo
            if (tstart <= tfdtopo(i)) then
                topoaltered(i) = .false.
            elseif (tstart > tfdtopo(i)) then
                topoaltered(i) = .true.
                write(GEO_PARM_UNIT,*) '  Altering topo arrays at t=', tstart
                print *, 'SETDTOPO Resetting topo arrays at t=',tstart
                do m=1,mtopofiles
                    if ((xlowtopo(m) <= xhidtopo(i)).and. &
                            (xhitopo(m) >= xlowdtopo(i)).and. &
                            (ylowtopo(m) <= yhidtopo(i)).and. &
                            (yhitopo(m) >= ylowdtopo(i))) then
                        do ib=1,mxtopo(m)
                            xcell=xlowtopo(m) + (ib-0.5d0)*dxtopo(m)
                            xim=xlowtopo(m) + (ib-1.d0)*dxtopo(m)
                            xip=xlowtopo(m) + ib*dxtopo(m)

                            if ((xim >= xlowdtopo(i)).and. &
                                    (xip <= xhidtopo(i))) then
                                do jb=1,mytopo(m)
                                    ycell=ylowtopo(m) + (jb-.5d0)*dytopo(m)
                                    yjm=ylowtopo(m) + (jb-1.d0)*dytopo(m)
                                    yjp=ylowtopo(m) + jb*dytopo(m)

                                    if ((yjm >= ylowdtopo(i)).and. &
                                            (yjp <= yhidtopo(i))) then
                                        ztopoij = 0.d0
                                        ijdtopo = i0dtopo(i)+(mtdtopo(i)-1) &
                                            *mxdtopo(i)*mydtopo(i)

                                        ztopoij= topointegral(xim,xcell,xip, &
                                            yjm,ycell,yjp,xlowdtopo(i), &
                                            ylowdtopo(i),dxdtopo(i), &
                                            dydtopo(i),mxdtopo(i), &
                                            mydtopo(i),dtopowork(ijdtopo),1)

                                        ztopoij=ztopoij/((yjp-yjm)*(xip-xim))
                                        if (icoordsys == 2) then
                                            deg2rad = pi/180.d0
                                            capac_area = deg2rad*Rearth**2 &
                                                * (sin(yjp*deg2rad) &
                                                - sin(yjm*deg2rad))/(yjp-yjm)
                                            ztopoij = ztopoij / capac_area
                                        endif

                                        jbr=mytopo(m)-jb+1
                                        ij=i0topo(m) + (jbr-1)*mxtopo(m)+ib-1
                                        topowork(ij)=topowork(ij)+ztopoij

                                    endif
                                enddo
                            endif
                        enddo
                    endif
                enddo
            endif
        enddo

    end subroutine set_dtopo
    ! ========================================================================

    ! ========================================================================
    !  read_dtopo(fname)
    ! ========================================================================
    subroutine read_dtopo(mx,my,mt,dtopo_type,fname,dtopo)

      implicit none

      ! Arguments
      integer, intent(in) :: mx,my,mt,dtopo_type
      character*150, intent(in) :: fname
      double precision, intent(inout) :: dtopo(1:mx*my*mt)

      ! Local
      integer, parameter :: iunit = 29
      integer :: i,j,k,dtopo_size,status
      double precision :: t,x,y

      open(unit=iunit, file=fname, status = 'unknown',form='formatted')

      select case(abs(dtopo_type))
         case(1)
            ! ASCII file with 4 columns
            do i = 1,mx*my*mt
               read(iunit,fmt=*,iostat=status) t,x,y, dtopo(i)
            enddo

         case(2)
            ! read header
            do i = 1,9
               read(iunit,*)
            enddo
            ! read the data
            do i = 1,mx*my*mt
               read(iunit,*) dtopo(i)
            enddo
         case(3)
            ! read header
            do i = 1,9
               read(iunit,*)
            enddo
            do k = 1,mt
               do j = 1,my
                  read(iunit,*) (dtopo((k-1)*mx*my + (j-1)*mx + i) , i=1,mx)
               enddo
            enddo
      end select

    end subroutine read_dtopo

    ! ========================================================================
    !  subroutine read_dtopo_header(fname,topo_type,mx,my,mt,xlow,ylow,t0,xhi,
    !                               yhi,tf,dx,dy,dt)
    ! ========================================================================

    !  Read in dtopo file header and either read or calculate the grid info
    !
    !  :Input:
    !   - fname - (char) Name of the dtopo file
    !   - topo_type - (int) Topography file type (1-3 are valid)
    !
    !  :Output:
    !   - mx,my,mt - (int) Number of grid point in space (mx,my) and time (mt)
    !   - xlow,ylow - (dp) Lower corner spatial coordinate of grid
    !   - xhi,yhi - (dp) Upper corner spatial coodinate of grid
    !   - t0,tf - (dp) Beginning and end times for the file
    !   - dx,dy,dt - (dp) Distance between space (dx,dy) and time (dt) points
    ! ========================================================================
    subroutine read_dtopo_header(fname,topo_type,mx,my,mt,xlow,ylow,t0,xhi, &
        yhi,tf,dx,dy,dt)

        implicit none

        ! Input Arguments
        character*150, intent(in) :: fname
        integer, intent(in) :: topo_type

        ! Output Arguments
        integer, intent(out) :: mx,my,mt
        double precision, intent(out) :: xlow,ylow,t0,xhi,yhi,tf,dx,dy,dt

        ! Locals
        integer, parameter :: iunit = 7
        integer :: topo_size,status
        double precision :: x,y,t,y_old,t_old
        logical :: found_file

        ! Open file
        inquire(file=fname,exist=found_file)
        if (.not.found_file) then
            print *, 'Missing dtopo file:'
            print *, '    ', fname
            stop
        endif
        open(unit=iunit,file=fname,status='unknown',form='formatted')

        select case(topo_type)
            ! Old style ASCII dtopo files
            case(1)
                ! Initial size variables
                topo_size = 0
                my = 1
                mt = 1

                ! Read in first values, determines xlow, yhi and t0
                read(iunit,*) t0,xlow,yhi
                topo_size = topo_size + 1
                t = t0
                y_old = yhi
                ! Go through entire file figuring out my, mt and topo_size
                status = 0
                do while (status == 0.and. t.eq.t0)
                    read(iunit,fmt=*,iostat=status) t,x,y
                    topo_size = topo_size + 1
                    if (y /= y_old .and. t.eq.t0 ) then
                        my = my + 1
                        y_old = y
                    endif
                enddo
                mx = (topo_size-1)/my
                do while (status == 0)
                    read(iunit,fmt=*,iostat=status) t,x,y
                    topo_size = topo_size + 1
                enddo


                if (status > 0) then
                    print *, "IO error occured in ",fname,", aborting!"
                    stop
                endif

                ! Calculate remaining values
                mt = (topo_size-1)/ (my*mx)
                xhi = x
                ylow = y
                tf = t
                dx = (xhi-xlow) / (mx-1)
                dy = (yhi-ylow) / (my-1)
                dt = (tf - t0) / (mt-1)

            ! New ASCII headered dtopo files, similar to topography files type
            ! 2 and 3
            case(2:3)
                ! Read in header directly
                read(iunit,*) mx
                read(iunit,*) my
                read(iunit,*) mt
                read(iunit,*) xlow
                read(iunit,*) ylow
                read(iunit,*) t0
                read(iunit,*) dx
                read(iunit,*) dy
                read(iunit,*) dt

                xhi = xlow + dx*(mx-1)
                yhi = ylow + dy*(my-1)
                tf = t0 + dt*(mt-1)
            case default
                print *, 'ERROR:  Unrecognized topography type'
                print *, '    topo_type = ',topo_type
                print *, '  for dtopo file:'
                print *, '   ', fname
                stop
        end select

         close(iunit)
    end subroutine read_dtopo_header

end module dtopo_module