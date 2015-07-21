
subroutine fgmax_interpolate(mx,my,meqn,mbc,maux,q,aux,dx,dy, &
           xlower,ylower,ifg,level,fg_values,mask_fgrid,fg_npts)

    ! Given a grid patch, return fg_values containing values interpolated
    ! to the fixed grid fg => FG_fgrids(ifg).
    ! If there are aux arrays, also set fg%aux at any points where it has
    ! not yet been set on this level.

    use fgmax_module

    implicit none
    integer, intent(in) :: mx,my,meqn,mbc,maux,ifg,level,fg_npts
    real(kind=8), intent(in) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: dx,dy,xlower,ylower
    !real(kind=8), dimension(:,:), allocatable, intent(inout) :: fg_values
    !logical, dimension(:), allocatable, intent(inout) :: mask_fgrid
    real(kind=8), intent(inout) :: fg_values(FG_NUM_VAL,fg_npts)
    logical, intent(inout) :: mask_fgrid(fg_npts)

    type(fgrid), pointer :: fg
    integer :: i,j,k,mv,ma
    integer :: i1,i2,j1,j2
    logical, allocatable, dimension(:,:) :: mask_patch
    real(kind=8), allocatable, dimension(:,:,:) :: values
    real(kind=8), allocatable, dimension(:,:) :: a,b,c
    real(kind=8), allocatable, dimension(:) :: dxk, dyk
    integer, allocatable, dimension(:) :: ik, jk
    real(kind=8) :: x1,x2,y1,y2,x,y,xupper,yupper
    logical :: debug

    debug = FG_DEBUG
    if (debug) then
        write(61,*) '========================================'
        write(61,*) 'In fgmax_interpolate, Level = ',level
        endif

    fg => FG_fgrids(ifg)
    !fg_npts = fg%npts  ! now an input parameter
    !write(61,*) '++++ interpolate xNbb', fg%x1bb, fg%x2bb
    
    !write(6,61) fg%fgno, fg%npts
 61 format('Updating fgrid number ',i2,' with',i7,' points')


    ! Determine intersection of this patch with bounding box of fgrid:

    xupper = xlower+mx*dx
    yupper = ylower+my*dy
    x1 = max(xlower, fg%x1bb) !- dx
    x2 = min(xupper, fg%x2bb) !+ dx
    y1 = max(ylower, fg%y1bb) !- dy
    y2 = min(yupper, fg%y2bb) !+ dy
    if (debug) then
        write(61,*) 'xlower,xupper: ',xlower,xupper
        write(61,*) 'ylower,yupper: ',ylower,yupper
        write(61,*) 'x1bb, x2bb: ',fg%x1bb, fg%x2bb
        write(61,*) 'y1bb, y2bb: ',fg%y1bb, fg%y2bb
        write(61,*) 'Intersection x1, x2:',x1,x2
        write(61,*) 'Intersection y1, y2:',y1,y2
        endif

    ! If this grid does not intersect the bounding box of fg, do nothing:
    if ((x1 > x2) .or. (y1 > y2)) then
        if (debug) then
            write(61,*) '+++ No intersection found!'
            endif
        mask_fgrid = .false.
        return
        endif

    allocate(mask_patch(1-mbc:mx+mbc, 1-mbc:my+mbc))
    allocate(values(FG_NUM_VAL, 1-mbc:mx+mbc, 1-mbc:my+mbc))
    allocate(a(1-mbc:mx+mbc, 1-mbc:my+mbc))
    allocate(b(1-mbc:mx+mbc, 1-mbc:my+mbc))
    allocate(c(1-mbc:mx+mbc, 1-mbc:my+mbc))
    allocate(dxk(1:fg%npts))
    allocate(dyk(1:fg%npts))
    allocate(ik(1:fg%npts))
    allocate(jk(1:fg%npts))

    if (mbc < 1) then
        write(6,*) '*** mbc >= 1 required by fgmax_interpolate'
        stop
        endif

    ! Create a mask that is .true. only in part of patch intersecting fgrid:
    i1 = max(int((x1 - xlower + 0.5d0*dx) / dx), 0)
    i2 = min(int((x2 - xlower + 0.5d0*dx) / dx) + 1, mx+1)
    j1 = max(int((y1 - ylower + 0.5d0*dy) / dy), 0)
    j2 = min(int((y2 - ylower + 0.5d0*dy) / dy) + 1, my+1)

    if (debug) then
        write(61,*) 'patch intersecting fgrid: i1,i2: ',i1,i2
        write(61,*) 'patch intersecting fgrid: j1,j2: ',j1,j2
        endif
    forall (i=1-mbc:mx+mbc, j=1-mbc:my+mbc)
        mask_patch(i,j) = ((i >= i1) .and. (i <= i2) .and. &
                           (j >= j1) .and. (j <= j2))
        end forall

    ! Create a mask that is .true. only in part of fgrid intersecting patch:
    do k=1,fg%npts
        mask_fgrid(k) = ((fg%x(k) >= x1) .and. (fg%x(k) <= x2) .and. &
                         (fg%y(k) >= y1) .and. (fg%y(k) <= y2))
        enddo

    if (debug) then
        write(61,*) '+++ in fgmax_interpolate, count = ',count(mask_fgrid)
        write(61,*) '+++ mask_patch is true only at i,j,x,y: '
        do i=1-mbc,mx+mbc
            x = xlower + (i-0.5d0)*dx
            do j=1-mbc,my+mbc
                y = ylower + (j-0.5d0)*dy
                if (mask_patch(i,j)) then
                    write(61,63) i,j,x,y
 63                 format(2i4,2d16.6)
                    endif
                enddo
            enddo
     
        endif
 

    ! Set the array values to be the values that we want to update:
    ! Note that values will be set properly only where mask_patch == .true.
    values = FG_NOTSET
    call fgmax_values(mx,my,meqn,mbc,maux,q,aux,dx,dy, &
                   xlower,ylower,mask_patch,values)
    if (debug) then
        write(65,*) '+++ i,j,x,y,values(1,i,j) '
        write(65,*) '    at points where mask_patch(i,j) == .true.'
        do i=1-mbc,mx+mbc
            x = xlower + (i-0.5d0)*dx
            do j=1-mbc,my+mbc
                y = ylower + (j-0.5d0)*dy
                if (mask_patch(i,j)) then
                    write(65,65) i,j,x,y,values(1,i,j)
 65                 format(2i4,2d16.6,d20.9)
                    endif
                enddo
            enddo
        endif
        

    ! Determine indices of cell center at lower left of rectangle including
    ! fixed grid point (x(k),y(k)) that will be used for bilinear interp.
    ! Then determine dxk and dyk, distance from each fgrid point to lower
    ! left point.

    ! Note that these are vectorized over all k:
    !ik = int((fg%x - xlower + 0.5d0*dx) / dx)
    !jk = int((fg%y - ylower + 0.5d0*dy) / dy)
    !dxk = fg%x - (xlower + (ik-0.5d0)*dx)
    !dyk = fg%y - (ylower + (jk-0.5d0)*dy)

    do k=1,fg%npts 
        ik(k) = int((fg%x(k) - xlower + 0.5d0*dx) / dx)
        jk(k) = int((fg%y(k) - ylower + 0.5d0*dy) / dy)
        dxk(k) = fg%x(k) - (xlower + (ik(k)-0.5d0)*dx)
        dyk(k) = fg%y(k) - (ylower + (jk(k)-0.5d0)*dy)
        enddo

    if (debug) then
        write(61,*) '+++ mask_fgrid is true only at k,x,y,ik,jk,dxk,dyk: '
        do k=1,fg%npts
            if (mask_fgrid(k)) then
                write(61,62) k,fg%x(k),fg%y(k), ik(k),jk(k),dxk(k),dyk(k)
 62             format(i4,2d16.6,2i6,2d16.6)
                endif
            enddo
        endif

    do mv=1,FG_NUM_VAL
        ! loop over the different values we want to monitor

        ! Compute coefficients needed for bilinear interpolation at all patch
        ! points in the intersection with the fgrid bounding box:
        forall (i=0:mx, j=0:my, mask_patch(i,j))
            a(i,j) = (values(mv,i+1,j) - values(mv,i,j)) / dx
            b(i,j) = (values(mv,i,j+1) - values(mv,i,j)) / dy
            c(i,j) = (values(mv,i+1,j+1) + values(mv,i,j) &
                      - (values(mv,i+1,j) + values(mv,i,j+1))) / (dx*dy)
            end forall

        do k=1,fg%npts 
            if (mask_fgrid(k)) then
    
                i = ik(k)
                j = jk(k)
                if ((i==mx+1) .or. (j==my+1)) then
                    write(6,*) '**** Warning from fgmax_interpolate:'
                    write(6,*) '**** Expected i <= mx, j<= my'
                    write(6,*) 'i,j,mx,my: ',i,j,mx,my
                    endif
                if ((values(mv,i,j)==FG_NOTSET) .or. &
                        (values(mv,i+1,j)==FG_NOTSET) .or. &
                        (values(mv,i,j+1)==FG_NOTSET) .or. &
                        (values(mv,i+1,j+1)==FG_NOTSET)) then
                    fg_values(mv,k) = FG_NOTSET
                else
                    fg_values(mv,k) = values(mv,ik(k),jk(k)) &
                       + a(ik(k),jk(k))*dxk(k) &
                       + b(ik(k),jk(k))*dyk(k) &
                       + c(ik(k),jk(k))*dxk(k)*dyk(k)
                    endif

!               write(6,64) mv,values(mv,ik(k),jk(k)),a(ik(k),jk(k)), &
!                  b(ik(k),jk(k)),c(ik(k),jk(k))
!64             format('mv,v,a,b,c: ',i2,4d16.6)
                endif
            enddo

        enddo


    ! Set the fg%aux(level,:,:) if this hasn't yet been set.
    if ((maux > 0) .and. (.not. fg%auxdone(level))) then
        ! at least some points do not yet have fg%aux set on this level

        do ma=1,FG_NUM_AUX
    
            ! Compute coefficients needed for bilinear interpolation at all 
            ! patch points in the intersection with the fgrid bounding box:
            forall (i=0:mx, j=0:my, mask_patch(i,j))
                a(i,j) = (aux(ma,i+1,j) - aux(ma,i,j)) / dx
                b(i,j) = (aux(ma,i,j+1) - aux(ma,i,j)) / dy
                c(i,j) = (aux(ma,i+1,j+1) + aux(ma,i,j) &
                          - (aux(ma,i+1,j) + aux(ma,i,j+1))) / (dx*dy)
                end forall
    
            do k=1,fg%npts 
                !print *, '+++ loop on k, aux = ', fg%aux(level,ma,k)
                !print *, '+++ diff = ', fg%aux(level,ma,k) - FG_NOTSET
                if (mask_fgrid(k) .and. (fg%aux(level,ma,k) == FG_NOTSET)) then
                    fg%aux(level,ma,k) = aux(ma,ik(k),jk(k)) &
                           + a(ik(k),jk(k))*dxk(k) &
                           + b(ik(k),jk(k))*dyk(k) &
                           + c(ik(k),jk(k))*dxk(k)*dyk(k)
                    !print *, '+++ set aux to ', fg%aux(level,ma,k)
    !               write(6,64) ma,aux(ma,ik(k),jk(k)),a(ik(k),jk(k)), &
    !                  b(ik(k),jk(k)),c(ik(k),jk(k))
    !64             format('ma,aux,a,b,c: ',i2,4d16.6)
                    endif
                enddo
    
            enddo
        if (minval(fg%aux(level,1,:)) > FG_NOTSET) then
            ! Done with aux arrays at all fgrid points on this level
            !print *, '+++ level,fg%aux:',level,fg%aux(level,1,1)
            fg%auxdone(level) = .true.
            endif
        endif

    deallocate(a,b,c,dxk,dyk,ik,jk,mask_patch)

end subroutine fgmax_interpolate
