module fixedgrids_module

    implicit none
    save

    ! Parameters for fixed grid data arrays
    integer, parameter :: maxfgrids = 3
    integer, parameter :: maxfgridsize = 423803

    ! Fixed grid arrays and sizes
    integer :: mfgrids
    real(kind=8), dimension(maxfgridsize) :: fgridearly,fgridlate,fgridoften
    integer, dimension(maxfgrids) :: mfgridvars,mfgridvars2,mxfg,myfg

    ! Geometry
    real(kind=8), dimension(maxfgrids) :: xlowfg,xhifg,ylowfg,yhifg,dxfg,dyfg

    ! Work array indices
    integer, dimension(maxfgrids) :: i0fg,i0fg2
      
    ! Time tracking and output
    integer, dimension(maxfgrids) :: noutfg,ilastoutfg
    integer, dimension(maxfgrids) :: ioutarrivaltimes,ioutsurfacemax
    real(kind=8), dimension(maxfgrids) :: tlastoutfg,tstartfg,tendfg,dtfg,tcfmax

contains
    
    ! Setup routine that reads in the fixed grids data file and sets up the
    ! appropriate data structures
    subroutine set_fixed_grids(fname)

        use amr_module

        implicit none
        
        ! Subroutine arguments
        character(len=*), optional, intent(in) :: fname
        
        ! Local storage
        integer, parameter :: unit = 7
        integer :: i,mspace,mspace2


        write(parmunit,*) ' '
        write(parmunit,*) '--------------------------------------------'
        write(parmunit,*) 'SETFIXEDGRIDS:'
        write(parmunit,*) '-----------'

        ! Open data file
        if (present(fname)) then
            call opendatafile(unit,fname)
        else
            call opendatafile(unit,'setfixedgrids.data')
        endif

        ! Read in data
        read(unit,'(i2)') mfgrids
        write(parmunit,*) '  mfgrids = ',mfgrids
        if (mfgrids > maxfgrids) then
            print *,'SETFIXEDGRIDS: ERROR mfgrids > maxfgrids'
            print *,'Decrease the number of fixed grids or'
            print *,'Increase maxfgrids in fixedgrids.i'
            stop
        endif

        if (mfgrids == 0) then
            write(parmunit,*) '  No fixed grids specified for output'
            return
        endif
        
        ! Initialize work array indices
        i0fg(1) = 1
        i0fg2(1) = 1

        ! Read in data for each fixed grid
        do i=1,mfgrids
            ! Read in this grid's data
            read(unit,*) tstartfg(i),tendfg(i),noutfg(i),xlowfg(i),xhifg(i), &
                   ylowfg(i),yhifg(i),mxfg(i),myfg(i), &
                   ioutarrivaltimes(i), ioutsurfacemax(i)

           ! Setup data for this grid
           ! Set dtfg (the timestep length between outputs) for each grid
           if (tendfg(i) <= tstartfg(i)) then
               if (noutfg(i).gt.1) then 
                  print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
                  print *,'tstartfg=tendfg yet noutfg>1'
                  print *,'set tendfg > tstartfg or set noutfg = 1'
                  stop
               else
                   dtfg(i)=0.d0
               endif
           else
               if (noutfg(i) < 2) then
                   print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
                   print *,'tendfg>tstartfg, yet noutfg=1'
                   print *,'set noutfg > 2'
                   stop
               else
                   dtfg(i) = (tendfg(i) - tstartfg(i)) / (noutfg(i) - 1)
               endif
            endif

            ! Initialize tlastoutfg and ilastoutfg
            tlastoutfg(i) = tstartfg(i) - dtfg(i)
            ilastoutfg(i) = 0

            ! Set spatial intervals dx and dy on each grid
            if (mxfg(i) > 1) then
               dxfg(i) = (xhifg(i) - xlowfg(i)) / (mxfg(i) - 1)
            else if (mxfg(i) == 1) then
               dxfg(i) = 0.d0
            else
                 print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
                 print *,'x grid points mxfg<=0, set mxfg>= 1'
            endif

            if (myfg(i) > 1) then
                dyfg(i) = (yhifg(i) - ylowfg(i)) / (myfg(i) - 1)
            else if (myfg(i).eq.1) then
                dyfg(i) = 0.d0
            else
                print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
                print *,'y grid points myfg<=0, set myfg>= 1'
            endif 
       
            ! set the number of variables stored for each grid
            ! this should be (the number of variables you want to write out + 1)
            mfgridvars(i)  = 6  
            mfgridvars2(i) = 3*ioutsurfacemax(i) + ioutarrivaltimes(i)
            
            ! find entry point into work arrays for each fixed grid
            i0fg(i)= i0fg(i-1)   + mfgridvars(i-1)*mxfg(i-1)*myfg(i-1)
            i0fg2(i)= i0fg2(i-1) + mfgridvars2(i-1)*mxfg(i-1)*myfg(i-1)
       enddo
       close(unit)
       
       ! Make sure enough space has been alotted for fixed grids in memory
       mspace=i0fg(mfgrids) + mfgridvars(mfgrids)*mxfg(mfgrids)*myfg(mfgrids)
       mspace2=i0fg2(mfgrids) + mfgridvars2(mfgrids)*mxfg(mfgrids)*myfg(mfgrids)
       mspace=mspace+mspace2
       if (mspace > maxfgridsize) then
           print *,'SETFIXEDGRIDS: ERROR not enough memory allocated'
           print *,'Decrease the number and size of fixed grids or'
           print *,'set maxfgridsize in fixedgrids.i to:', mspace
           stop
       endif

       ! Initialize fixed grid work arrays to huge, this will prevent 
       ! non-filled values from being misinterpreted
       fgridearly = nan()
       fgridlate = nan()
       fgridoften = nan()
       
       tcfmax=-1.d16

    end subroutine set_fixed_grids
    
    !=====================FGRIDINTERP=======================================
    !         # This routine interpolates q and aux on a computational grid
    !         # to a fgrid not necessarily aligned with the computational grid
    !         # using bilinear interpolation defined on computation grid
    !=======================================================================
    subroutine fgrid_interp(fgrid,xlowfg,ylowfg, &
                            xhifg,yhifg,dxfg,dyfg,mxfg,myfg,t,mvarsfg,q,meqn, &
                            mxc,myc,mbc,dxc,dyc,nvar,xlowc,ylowc,maux,aux, &
                            ioutarrivaltimes,ioutsurfacemax,maxcheck)
    
        use geoclaw_module
    
        implicit none
    
        ! Subroutine arguments
        integer, intent(in) :: mxfg,myfg,mvarsfg,meqn,mxc,myc,mbc,nvar,maux
        integer, intent(in) :: ioutarrivaltimes,ioutsurfacemax,maxcheck
        real(kind=8), intent(in) :: xlowfg,ylowfg,xhifg,yhifg,dxfg,dyfg
        real(kind=8), intent(in) :: t,dxc,dyc,xlowc,ylowc
        real(kind=8), intent(inout) :: fgrid(1:mxfg,1:myfg,mvarsfg)
        real(kind=8), intent(in) :: q(1-mbc:mxc+mbc,1-mbc:myc+mbc, meqn)
        real(kind=8), intent(in) :: aux(1-mbc:mxc+mbc,1-mbc:myc+mbc, maux)
    
        ! Locals
        integer :: ifg,jfg,m,ic1,ic2,jc1,jc2
        integer :: indb,indeta,indarrive,indetamin,indetamax,indetanow
        real(kind=8) :: xfg,yfg,xc1,xc2,yc1,yc2,xhic,yhic,xterm,yterm,xyterm
        real(kind=8) :: tol,arrivaltol,totaldepth,depthindicator,check
        real(kind=8) :: z11,z12,z21,z22,z11w,z12w,z21w,z22w
        real(kind=8) :: a,b,c,d,h11,h12,h21,h22
            
    
        xhic=xlowc + dxc*mxc  
        yhic=ylowc + dyc*myc    
        
        tol=drytolerance
        arrivaltol=1.d-2
        
        indb=meqn+1
        indeta=meqn+2
    
        if (maxcheck.gt.0) then 
            indarrive=0
            indetamin=0
            indetamax=0
            indetanow=0
    
            if (ioutarrivaltimes.gt.0) then
                indarrive = 1
            endif
        
            if (ioutsurfacemax.gt.0) then
                indetanow = indarrive+1
                indetamin = indarrive+2
                indetamax = indarrive+3
            endif
        endif
    
        ! Primary interpolation loops 
        do ifg=1,mxfg
            xfg=xlowfg + (ifg-1)*dxfg
            do jfg=1,myfg
                yfg=ylowfg + (jfg-1)*dyfg
    
                if (.not.((xfg.lt.xlowc.or.xfg.gt.xhic).or.(yfg.lt.ylowc.or.yfg.gt.yhic))) then
    
                    ! find where xfg,yfg is in the computational grid
                    ic1 = int((xfg-(xlowc+0.5d0*dxc))/(dxc))+1
                    jc1 = int((yfg-(ylowc+0.5d0*dyc))/(dyc))+1
                    if (ic1.eq.mxc) ic1=mxc-1
                    if (jc1.eq.myc) jc1=myc-1 
                    ic2= ic1+1
                    jc2= jc1+1
                        
                    xc1= xlowc + dxc*(ic1-0.5d0)
                    yc1= ylowc + dyc*(jc1-0.5d0)
                    xc2= xlowc + dxc*(ic2-0.5d0)
                    yc2= ylowc + dyc*(jc2-0.5d0)
         
                    ! interpolate bilinear used to interpolate to xfg,yfg
                    ! define constant parts of bilinear
                    xterm=(xfg-xc1)/dxc
                    yterm=(yfg-yc1)/dyc
                    xyterm= xterm*yterm
        
                    if (maxcheck.eq.0) then 
                        do m=1,meqn
                            z11=q(ic1,jc1,m)
                            z21=q(ic2,jc1,m)
                            z12=q(ic1,jc2,m)
                            z22=q(ic2,jc2,m)
                            a=z21-z11
                            b=z12-z11
                            d=z11
                            c=z22-(a+b+d)
                            fgrid(ifg,jfg,m) = a*xterm + b*yterm + c*xyterm + d
                        enddo
                        z11=aux(ic1,jc1,1)
                        z21=aux(ic2,jc1,1)
                        z12=aux(ic1,jc2,1)
                        z22=aux(ic2,jc2,1) 
                        a=z21-z11
                        b=z12-z11
                        d=z11
                        c=z22-(a+b+d)
                        fgrid(ifg,jfg,indb) = a*xterm + b*yterm + c*xyterm + d
                    endif
                    ! This next output variable is the surface using bilinear interpolation,
                    ! using a surface that only uses the wet eta points near the shoreline
        
                    z11=aux(ic1,jc1,1)+q(ic1,jc1,1)
                    z21=aux(ic2,jc1,1)+q(ic2,jc1,1)
                    z12=aux(ic1,jc2,1)+q(ic1,jc2,1)
                    z22=aux(ic2,jc2,1)+q(ic2,jc2,1)
                        
                    h11=q(ic1,jc1,1)
                    h21=q(ic2,jc1,1)
                    h12=q(ic1,jc2,1)
                    h22=q(ic2,jc2,1)
                    depthindicator= min(h11,h12,h21,h22)
                    totaldepth= h11+h22+h21+h12
    
                    if (depthindicator.lt.tol.and.totaldepth.gt.4.d0*tol) then !near shoreline
                        if (h11.lt.tol) then
                            z11w=  (h12*z12 + h21*z21 + h22*z22)/totaldepth
                            z11=z11w
                        endif
                        if (h12.lt.tol) then
                            z12w=  (h11*z11 + h21*z21 + h22*z22)/totaldepth
                            z12=z12w
                        endif
                        if (h21.lt.tol) then
                            z21w=  (h11*z11 + h12*z12 + h22*z22)/totaldepth
                            z21=z21w
                        endif
                        if (h22.lt.tol) then
                            z22w=  (h12*z12 + h21*z21 + h11*z11)/totaldepth
                            z22=z22w
                        endif            
                    endif
                    if (totaldepth.le.4.d0*tol) then
                        z22=nan()
                    endif
    
                    a=z21-z11
                    b=z12-z11
                    d=z11
                    c=z22-(a+b+d)
    
                    ! If eta max/min are saved on this grid initialized if necessary
                    if (ioutsurfacemax.gt.0.and.maxcheck.eq.2) then 
                        if (.not.(fgrid(ifg,jfg,indetamin).eq.fgrid(ifg,jfg,indetamin))) fgrid(ifg,jfg,indetamin)=0.d0
                        if (.not.(fgrid(ifg,jfg,indetamax).eq.fgrid(ifg,jfg,indetamax))) fgrid(ifg,jfg,indetamax)=0.d0
                    endif
    
                ! check which task to perform
                    if (maxcheck.eq.0) then 
                        fgrid(ifg,jfg,indeta) = a*xterm + b*yterm + c*xyterm + d
                        fgrid(ifg,jfg,mvarsfg) = t
                    else if (maxcheck.eq.1.and.ioutsurfacemax.gt.0) then
                        fgrid(ifg,jfg,indetanow) = a*xterm + b*yterm + c*xyterm + d   
                    else if (maxcheck.eq.2.and.ioutsurfacemax.gt.0) then
                        fgrid(ifg,jfg,indetamin) = min(fgrid(ifg,jfg,indetamin),fgrid(ifg,jfg,indetanow))
                        fgrid(ifg,jfg,indetamax) = max(fgrid(ifg,jfg,indetamax),fgrid(ifg,jfg,indetanow))            
                    endif
    
                    ! If arrival times are saved on this grid
                    if (maxcheck.eq.1.and.ioutarrivaltimes.gt.0) then
                        check=fgrid(ifg,jfg,indarrive)
                        !# check=NaN: Waves haven't arrived previously
                        if (.not.(check == check)) then
                            if (abs(fgrid(ifg,jfg,indeta)).gt.arrivaltol) then
                                fgrid(ifg,jfg,indarrive)= t
                            endif
                        endif
                    endif
                endif
            enddo
        enddo
    
    end subroutine fgrid_interp

    !=====================FGRIDOUT==========================================
    ! This routine interpolates in time and then outputs a grid at
    ! time=toutfg
    !
    ! files have a header, followed by columns of data
    !=======================================================================
    subroutine fgrid_out(fgrid1,fgrid2,fgrid3,xlowfg,xhifg,ylowfg,yhifg, &
                        mxfg,myfg,mvarsfg,mvarsfg2,toutfg,ioutfg,ng, &
                        ioutarrival,ioutflag)

        implicit none
        
        ! Subroutine arguments
        integer, intent(in) :: mxfg,myfg,mvarsfg,mvarsfg2,ioutfg,ng
        integer, intent(in) :: ioutarrival,ioutflag
        real(kind=8), intent(in) :: xlowfg,xhifg,ylowfg,yhifg,toutfg
        real(kind=8), intent(inout) :: fgrid1(1:mxfg,1:myfg,mvarsfg)
        real(kind=8), intent(inout) :: fgrid2(1:mxfg,1:myfg,mvarsfg)
        real(kind=8), intent(inout) :: fgrid3(1:mxfg,1:myfg,mvarsfg2)
              
        ! Locals
        integer, parameter :: unit = 95
        integer :: ifg,jfg,iv
        integer :: ngridnumber,ipos,idigit,noutnumber,icolumns
        integer :: indetamin,indetamax
        real(kind=8) :: t0,tf,tau
        character(len=30) :: fgoutname

        character(len=*), parameter :: header_format = "(e18.8,'    time', /," // &
                                                        "i5,'    mx', /,"   // &
                                                        "i5,'    my', /,"   // &
                                                     "e18.8,'    xlow',/"   // &
                                                     "e18.8,'    ylow',/"   // &
                                                     "e18.8,'    xhi',/,"   // &
                                                     "e18.8,'    yhi',/,"   // &
                                                        "i5,'  columns',/)"
        character(len=*), parameter :: data_format = "(8e26.16)"
        character(len=*), parameter :: arrival_header_format = &
                                                       "(i5,'    mx', /," // &
                                                        "i5,'    my', /," // &
                                                     "e18.8,'    xlow',/" // &
                                                     "e18.8,'    ylow',/" // &
                                                     "e18.8,'    xhi',/," // &
                                                     "e18.8,'    yhi',/," // &
                                                        "i5,'  columns',/)"
    

        ! Make the file names and open output files
        fgoutname = 'fort.fgnn_xxxx'
        ngridnumber= ng
        do ipos = 9, 8, -1
            idigit= mod(ngridnumber,10)
            fgoutname(ipos:ipos) = char(ichar('0') + idigit)
            ngridnumber = ngridnumber/ 10
        enddo

        noutnumber=ioutfg
        do ipos = 14, 11, -1
            idigit = mod(noutnumber,10)
            fgoutname(ipos:ipos) = char(ichar('0') + idigit)
            noutnumber = noutnumber / 10
        enddo

        open(unit,file=fgoutname,status='unknown',form='formatted')

        icolumns = mvarsfg -1
        if (mvarsfg2.gt.1) then
           icolumns=icolumns + 2
        endif
        
        write(unit,header_format) toutfg,mxfg,myfg,xlowfg,ylowfg,xhifg,yhifg,icolumns

        indetamin = ioutarrival+2
        indetamax = ioutarrival+3
        
        ! interpolate the grid in time, to the output time, using 
        ! the solution in fgrid1 and fgrid2, which represent the 
        ! solution on the fixed grid at the two nearest computational times
        do jfg=1,myfg
            do ifg=1,mxfg
                t0=fgrid1(ifg,jfg,mvarsfg)
                tf=fgrid2(ifg,jfg,mvarsfg)
                tau=(toutfg-t0)/(tf-t0)
                
                do iv=1,mvarsfg-1
                    if (dabs(fgrid1(ifg,jfg,iv)) .lt. 1d-90) fgrid1(ifg,jfg,iv) = 0.d0
                    if (dabs(fgrid2(ifg,jfg,iv)) .lt. 1d-90) fgrid2(ifg,jfg,iv) = 0.d0
                enddo
                if (icolumns.eq.mvarsfg-1) then 
                    write(unit,data_format) ((1.d0 - tau)*fgrid1(ifg,jfg,iv)+tau*fgrid2(ifg,jfg,iv), iv=1,mvarsfg-1)
                else
                    if (abs(fgrid3(ifg,jfg,indetamin)) .lt. 1d-90) fgrid3(ifg,jfg,indetamin) = 0.d0
                    if (abs(fgrid3(ifg,jfg,indetamax)) .lt. 1d-90) fgrid3(ifg,jfg,indetamax) = 0.d0
                    write(unit,data_format) ((1.d0 - tau)*fgrid1(ifg,jfg,iv)+tau*fgrid2(ifg,jfg,iv), iv=1,mvarsfg-1), &
                                          fgrid3(ifg,jfg,indetamin), &
                                          fgrid3(ifg,jfg,indetamax)
                endif
            enddo
        enddo

        close(unit)
        print "(a,i2,a,i2,a,e18.8)",' FGRIDOUT: Fixed Grid  ',ng, '  frame ',ioutfg,' at time =',toutfg

        ! ==================== Output for arrival times============
        if (ioutflag.eq.1) then
            ! Make the file name and open output file for arrival times
            fgoutname = 'fort.fgnn_arrivaltimes'
            ngridnumber= ng
            do ipos = 9, 8, -1
                idigit= mod(ngridnumber,10)
                fgoutname(ipos:ipos) = char(ichar('0') + idigit)
                ngridnumber = ngridnumber/ 10
            enddo
            open(unit,file=fgoutname,status='unknown',form='formatted')

            write(95,arrival_header_format) mxfg,myfg,xlowfg,ylowfg,xhifg,yhifg

            do jfg=1,myfg
                do ifg=1,mxfg
                    write(unit,"(1e26.16)") fgrid3(ifg,jfg,1)
                enddo
            enddo
            close(unit)

            print "(a,i2,a)", ' FGRIDOUT: Fixed Grid  ', ng, '  arrival times output'
        endif
      
    end subroutine fgrid_out
    
    real(kind=8) function nan()
        real(kind=8) dnan
        integer inan(2)
        equivalence (dnan,inan)
        inan(1)=2147483647
        inan(2)=2147483647
        nan=dnan
    end function nan

end module fixedgrids_module

! module fixedgrids_module
! 
!     implicit none
!     save
! 
!     type fixed_grid_data_type
!         real(kind=8), pointer :: data(:,:,:)
!     end type fixed_grid_data_type
! 
!     ! Number of fixed grids
!     integer :: num_fixed_grids
!     
!     ! Primary data storage, these are each an array of pointers
!     integer, allocatable :: num_grid_vars(:,:)
!     type(fixed_grid_data_type), allocatable :: early_data_fg(:)
!     type(fixed_grid_data_type), allocatable :: often_data_fg(:)
!     type(fixed_grid_data_type), allocatable :: late_data_fg(:)
!     
!     ! Parameters for output of each fixed grid
!     integer, allocatable :: arrival_times_output(:),surface_max_output(:)
!     integer, allocatable :: num_output_fg(:)
!     real(kind=8), allocatable :: t_start_fg(:),t_end_fg(:)
! 
!     ! Output tracking
!     integer, allocatable :: t_last_output_index_fg(:)
!     real(kind=8), allocatable :: t_last_output_fg(:)
!     real(kind=8) :: max_fg_time
!      
!     ! Geometry and time step size
!     integer, allocatable :: mx_fg(:),my_fg(:)
!     real(kind=8), allocatable :: x_low_fg(:),x_hi_fg(:)
!     real(kind=8), allocatable :: y_low_fg(:),y_hi_fg(:)
!     real(kind=8), allocatable :: dt_fg(:),dx_fg(:),dy_fg(:)
!        
! contains
! 
! !     ! Read in fixed grid settings and setup data structures
! !     subroutine set_fixed_grids(fname)
! ! 
! !         use amr_module, only: parmunit
! ! 
! !         implicit none
! !         
! !         ! Subroutine arguments
! !         character(len=25), optional :: fname
! !         
! !         ! File opening
! !         integer, parameter :: unit = 7
! !         character(len=*), parameter :: fg_line_format = "(2d16.8,1i2,4d16.8,2i4,2d16.8)"
! !         
! !         ! Allocation pointer for fixed grid data
! !         real(kind=8), pointer :: temp_data(:,:,:)
! ! 
! !         ! Other locals
! !         integer :: i
! ! 
! !         write(parmunit,*) ' '
! !         write(parmunit,*) '--------------------------------------------'
! !         write(parmunit,*) 'SETFIXEDGRIDS:'
! !         write(parmunit,*) '-----------'
! ! 
! !         ! Open file
! !         if (present(fname)) then
! !             call opendatafile(unit,fname)
! !         else
! !             call opendatafile(unit,'setfixedgrids.data')
! !         endif
! ! 
! !         ! Read in data
! !         read(7,"(i2)") num_fixed_grids
! !         write(parmunit,*) '  num_fixed_grids = ',num_fixed_grids
! !         if (num_fixed_grids == 0) then
! !             write(parmunit,*) "  No fixed grids specified for output"
! !             return
! !         endif
! !         
! !         ! Allocate all fixed grid data and info arrays
! !         allocate(t_start_fg(num_fixed_grids),t_end_fg(num_fixed_grids))
! !         allocate(num_output_fg(num_fixed_grids))
! !         allocate(x_low_fg(num_fixed_grids),x_hi_fg(num_fixed_grids))
! !         allocate(y_low_fg(num_fixed_grids),y_hi_fg(num_fixed_grids))
! !         allocate(mx_fg(num_fixed_grids),my_fg(num_fixed_grids))
! !         allocate(dt_fg(num_fixed_grids))
! !         allocate(dx_fg(num_fixed_grids),dy_fg(num_fixed_grids))
! !         allocate(arrival_times_output(num_fixed_grids))
! !         allocate(surface_max_output(num_fixed_grids))
! !         allocate(num_grid_vars(num_fixed_grids,2))
! !         allocate(t_last_output_index_fg(num_fixed_grids))
! !         allocate(t_last_output_fg(num_fixed_grids))
! !         
! !         ! These are the data arrays themselves (only pointers)
! !         allocate(early_data_fg(num_fixed_grids))
! !         allocate(late_data_fg(num_fixed_grids))
! !         allocate(often_data_fg(num_fixed_grids))
! !         
! !         ! Read in parameters for each fixed grid
! !         do i=1,num_fixed_grids
! !             read(unit,*) t_start_fg(i),t_end_fg(i),num_output_fg(i), &
! !                                    x_low_fg(i),x_hi_fg(i),y_low_fg(i),y_hi_fg(i), &
! !                                    mx_fg(i),my_fg(i), &
! !                                    arrival_times_output(i),surface_max_output(i)
! !             write(parmunit,*) t_start_fg(i),t_end_fg(i),num_output_fg(i), &
! !                                    x_low_fg(i),x_hi_fg(i),y_low_fg(i),y_hi_fg(i), &
! !                                    mx_fg(i),my_fg(i), &
! !                                    arrival_times_output(i),surface_max_output(i)
! !         enddo
! !         close(unit)
! !        
! !         ! Set some parameters for each grid
! !         do i=1,num_fixed_grids
! !             ! Set time step length between outputs
! !             if (t_end_fg(i) <= t_start_fg(i)) then
! !                 if (num_output_fg(i) > 1) then
! !                     print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
! !                     print *,'tstartfg=tendfg yet noutfg>1'
! !                     print *,'set tendfg > tstartfg or set noutfg = 1'
! !                     stop
! !                 else
! !                     dt_fg(i) = 0.d0
! !                 endif
! !             else
! !                 if (num_output_fg(i) < 2) then
! !                     print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
! !                     print *,'tendfg>tstartfg, yet noutfg=1'
! !                     print *,'set noutfg > 2'
! !                     stop
! !                 else
! !                     dt_fg(i) = (t_end_fg(i) - t_start_fg(i)) / real(num_output_fg(i)-1,kind=8)
! !                 endif
! !             endif
! !             
! !             ! Counters for keeping track of output times
! !             t_last_output_fg(i) = t_start_fg(i) - dt_fg(i)
! !             t_last_output_index_fg(i) = 0
! ! 
! !             ! Set spatial intervals dx and dy on each grid
! !             if (mx_fg(i) > 1) then
! !                 dx_fg(i) = (x_hi_fg(i) - x_low_fg(i)) / real(mx_fg(i) - 1,kind=8)
! !             else if (mx_fg(i) == 1) then
! !                 dx_fg(i) = 0.d0
! !             else
! !                 print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
! !                 print *,'x grid points mxfg<=0, set mxfg>= 1'
! !             endif
! !             if (my_fg(i) > 1) then
! !                 dy_fg(i) = (y_hi_fg(i) - y_low_fg(i)) / real(my_fg(i) - 1,kind=8)
! !             else if (my_fg(i) == 1) then
! !                 dy_fg(i) = 0.d0
! !             else
! !                 print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
! !                 print *,'y grid points myfg<=0, set myfg>= 1'
! !             endif
! ! 
! !         enddo
! ! 
! !       ! Set the number of variables stored for each grid
! !       ! this should be (the number of variables you want to write out + 1)
! !       do i=1, num_fixed_grids 
! !          num_grid_vars(i,1)  = 6  
! !          num_grid_vars(i,2) = 3 * surface_max_output(i) + arrival_times_output(i)
! !       enddo
! ! 
! !       ! Allocate data arrays, here we use allocate data to the temprorary
! !       ! pointer and transfer the data to the array of pointers stored in
! !       ! each of the individual arrays of pointers.  Since we only do this once
! !       ! this seems like it should not be intolerable in terms of performance
! !       !
! !       ! This code used to fill grid data with NaNs to prevent non-filled values
! !       ! from being misinterpreted, this was not portable and so the arrays are
! !       ! filled with the intrinsic `huge` instead (largest representable number)
! !       do i=1,num_fixed_grids
! !           allocate(temp_data(num_grid_vars(i,1),mx_fg(i),my_fg(i)))
! !           temp_data = huge(1.d0)
! ! !           call move_alloc(temp_data,early_data_fg(i))
! !           early_data_fg(i)%data => temp_data
! !       enddo
! !       do i=1,num_fixed_grids
! !           allocate(temp_data(num_grid_vars(i,2),mx_fg(i),my_fg(i)))
! !           temp_data = huge(1.d0)
! ! !           call move_alloc(temp_data,late_data_fg(i))
! !           late_data_fg(i)%data => temp_data
! !       enddo
! !       do i=1,num_fixed_grids
! !           allocate(temp_data(num_grid_vars(i,2),mx_fg(i),my_fg(i)))
! !           temp_data = huge(1.d0)
! ! !           call move_alloc(temp_data,often_data_fg(i))
! !           often_data_fg(i)%data => temp_data
! !       enddo
! !       
! !       ! Set convenience maximum output for fixed grid time
! !       max_fg_time = 1d-16
! ! 
! !     end subroutine set_fixed_grids
! ! 
! ! 

! ! 
! !     ! Advance (output) all fgrids at all times that have not yet been output
! !     ! but that have been bracketed by computational times.
! !     subroutine fgrid_advance(t,dt)
! !         
! ! !         use amr_module
! ! !         use regions_module
! ! !         use qinit_module
! ! !         use gauges_module
! !         
! !         implicit none
! !     
! !         ! Subroutine arguments
! !         real(kind=8), intent(in) :: t,dt
! !         
! !         ! Local storage
! !         real(kind=8) :: tc_start,tc_final,out_time
! !         integer :: n,i,out_start_index,out_end_index,out_flag
! !         
! !         ! Store computational step times
! !         tc_start = t
! !         tc_final = tc_start + dt
! !         
! !         ! Check to see if any fixed grids should be outputed
! !         do n=1,num_fixed_grids
! !             if (tc_start > t_start_fg(n) .and. t_last_output_index_fg(n) < num_output_fg(n)) then
! !                 ! fgrid n may need to be written out
! !                 ! find the first output number that has not been written out and
! !                 ! find the first output number on a fixed grid that is >= tc0
! !                 ! which will not be written out
! !                 if (dt_fg(n) > 0.d0) then
! !                     out_end_index = 1 + max(0,int((tc_start - t_start_fg(n)) / dt_fg(n)))
! !                 else
! !                     out_end_index = 1
! !                 endif
! !                 out_end_index = min(out_end_index,num_output_fg(n))
! !                 out_start_index = t_last_output_index_fg(n)
! ! 
! !                 ! write-out fgrid times that are less than tc0, and have not been 
! !                 ! written yet these should be the most accurate values at any given 
! !                 ! point in the fgrid since tc0> output time
! !                 do i=out_start_index,out_end_index
! !                     out_time = t_start_fg(n) + (i-1) * dt_fg(n)
! !                     if (out_time < tc_start) then
! !                         ! Write out solution for fixed grid n
! !                         out_flag = arrival_times_output(n) * (num_output_fg(n) - t_last_output_index_fg(n))
! !                         
! !                         call fgrid_output(early_data_fg(n)%data, &
! !                                           late_data_fg(n)%data,  &
! !                                           often_data_fg(n)%data, &
! !                                           x_low_fg(n),x_hi_fg(n),y_low_fg(n),y_hi_fg(n), &
! !                                           mx_fg(n),my_fg(n),num_grid_vars(n,:), &
! !                                           out_time,i,n,arrival_times_output(n),out_flag)
! !                         
! !                         t_last_output_fg(n) = out_time
! !                         t_last_output_index_fg(n) = t_last_output_index_fg(n) + 1
! !                     endif
! !                 enddo
! !             endif
! !         enddo
! !         
! !     end subroutine fgrid_advance
! ! 
! !     subroutine fgrid_output(fgrid1,fgrid2,fgrid3, &
! !                         x_low,x_hi,y_low,y_hi,mx,my,num_vars, &
! !                         t_out,num_out,grid_number, &
! !                         out_arrival,out_flag)
! ! 
! !         use geoclaw_module
! ! 
! !         implicit none
! !     
! !         ! Subroutine arguments
! !         integer, intent(in) :: num_vars(2),mx,my
! !         integer, intent(in) :: num_out,grid_number,out_arrival,out_flag
! !         real(kind=8), intent(in) :: x_low,x_hi,y_low,y_hi,t_out
! !         real(kind=8), intent(inout), dimension(num_vars(1),1:mx,1:my) :: fgrid1,fgrid2,fgrid3
! ! 
! !         ! File handling
! !         character(len=30) :: fg_file_name
! !         integer, parameter :: unit = 90
! !         integer :: digit,nout,ng,pos,status
! !     
! !         ! Other locals
! !         integer :: i,j,m,deta_min,deta_max,num_columns
! !         real(kind=8) :: t0,tf,tau
! !     
! !         ! Format strings
! !         character(len=*), parameter :: header_format = "(e18.8,'    time', /," // &
! !                                                            "i5,'    mx', /,"   // &
! !                                                            "i5,'    my', /,"   // &
! !                                                         "e18.8,'    xlow',/"   // &
! !                                                         "e18.8,'    ylow',/"   // &
! !                                                         "e18.8,'    xhi',/,"   // &
! !                                                         "e18.8,'    yhi',/,"   // &
! !                                                            "i5,'  columns',/)"
! !         character(len=*), parameter :: data_format = "(8e26.16)"
! !         character(len=*), parameter :: arrival_header_format = &
! !                                                           "(i5,'    mx', /," // &
! !                                                            "i5,'    my', /," // &
! !                                                         "e18.8,'    xlow',/" // &
! !                                                         "e18.8,'    ylow',/" // &
! !                                                         "e18.8,'    xhi',/," // &
! !                                                         "e18.8,'    yhi',/," // &
! !                                                            "i5,'  columns',/)"
! !         character(len=*), parameter :: arrival_data_format = "(1e26.16)"
! ! 
! !         ! Open file for writing                               
! !         fg_file_name = 'fort.fgnn_xxxx'
! !         ng = grid_number
! !         do pos=9,8,-1
! !             digit = mod(ng,10)
! !             fg_file_name(pos:pos) = char(ichar('0') + digit)
! !             ng = ng / 10
! !         enddo
! !     
! !         nout = num_out
! !         do pos=14,11,-1
! !             digit = mod(nout,10)
! !             fg_file_name(pos:pos) = char(ichar('0') + digit)
! !             nout = nout / 10
! !         enddo
! !     
! !         open(unit=unit, file=fg_file_name, iostat=status, status="new", action="write")
! !         if ( status /= 0 ) then
! !             print *, "Error opening file fixed grid output file ",fg_file_name
! !             stop
! !         endif
! !     
! !         ! Determine the number of columns of the output file
! !         num_columns = num_vars(1) - 1
! !         if (num_vars(2) > 1) then
! !             num_columns = num_columns + 2
! !         endif
! !         deta_min = out_arrival + 2
! !         deta_max = out_arrival + 3
! !     
! !         ! Write out header
! !         write(unit,header_format) t_out,mx,my,x_low,y_low,x_hi,y_hi,num_columns
! ! 
! !         ! Interpolate the grid in time to hit the requested output time using the
! !         ! solution in fgrid1 and fgrid2 which represent the solution on the fixed
! !         ! grid at the two nearest computational times
! !         do j=1,my
! !             do i=1,mx
! !                 ! Figure out time interpolant
! !                 t0 = fgrid1(num_vars(1),i,j)
! !                 tf = fgrid2(num_vars(1),i,j)
! !                 print *,t0,tf
! !                 tau = (t_out - t0) / (tf - t0)
! !             
! !                 ! Zero out values that are too small
! !                 do m=1,num_vars(1)-1
! !                     if (abs(fgrid1(m,i,j)) < 1d-90) then
! !                         fgrid1(m,i,j) = 0.d0
! !                     endif
! !                     if (abs(fgrid2(m,i,j)) < 1d-90) then
! !                         fgrid2(m,i,j) = 0.d0
! !                     endif
! !                 enddo
! !             
! !                 ! Write out interpolants
! !                 if (num_columns == num_vars(1) - 1) then
! !                     write(unit,data_format) ((1.d0 - tau) * fgrid1(m,i,j) &
! !                                 + tau * fgrid2(m,i,j),m=1,num_vars(1)-1)
! !                 else
! !                     if (abs(fgrid3(deta_min,i,j)) < 1d-90) then
! !                         fgrid3(deta_min,i,j) = 0.d0
! !                     endif
! !                     if (abs(fgrid3(deta_max,i,j)) < 1d-90) then
! !                         fgrid3(deta_max,i,j) = 0.d0
! !                     endif
! !                 
! !                     write(unit,data_format) ((1.d0 - tau) * fgrid1(m,i,j) &
! !                             + tau * fgrid2(m,i,j), m=1,num_vars(1)-1), &
! !                             fgrid3(deta_min,i,j),fgrid3(deta_max,i,j)
! !                 endif
! !             enddo
! !         enddo
! !     
! !         close(unit)
! !         print "(' FGRIDOUT: Fixed Grid  ', i2, '  output at time =', e18.8)",grid_number,t_out
! !     
! !         ! Output arrival times
! !         if (out_flag == 1) then
! !             ! Open output file
! !             fg_file_name = 'fort.fgnn_arrivaltimes'
! !             ng = grid_number
! !             do pos=9,8,-1
! !                 digit = mod(ng,10)
! !                 fg_file_name(pos:pos) = char(ichar('0') + digit)
! !                 ng = ng / 10
! !             enddo
! !         
! !             open(unit=unit, file=fg_file_name, iostat=status, status="new", action="write")
! !             if ( status /= 0 ) then
! !                 print *,"Error opening file for arrival times file ",fg_file_name
! !                 stop
! !             endif
! !         
! !             write(unit,arrival_header_format) mx,my,x_low,y_low,x_hi,y_hi
! !             do j=1,my
! !                 do i=1,mx
! !                     write(unit,arrival_data_format) fgrid3(1,i,j)
! !                 enddo
! !             enddo
! !         
! !             close(unit)
! !             print "(' FGRIDOUT: Fixed Grid  ', i2, '  arrival times output')",grid_number
! !         endif
! !     
! !     end subroutine fgrid_output
! 
! end module fixedgrids_module
