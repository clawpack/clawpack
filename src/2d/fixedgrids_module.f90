module fixedgrids_module

    implicit none
    save

    ! Container for fixed grid data, geometry and output settings
    type fixedgrid_type
        ! Grid data
        real(kind=8), pointer :: early(:,:,:)
        real(kind=8), pointer :: late(:,:,:)
        real(kind=8), pointer :: often(:,:,:)
        
        ! Geometry
        integer :: num_vars(2),mx,my
        real(kind=8) :: dx,dy,x_low,x_hi,y_low,y_hi
        
        ! Time Tracking and output types
        integer :: num_output,last_output_index
        integer :: output_arrival_times,output_surface_max
        real(kind=8) :: last_output_time,start_time,end_time,dt
    end type fixedgrid_type    

    ! Fixed grid arrays and sizes
    integer :: num_fixed_grids
    type(fixedgrid_type), allocatable :: fgrids(:)
    real(kind=8) :: tcfmax

contains
    
    ! Setup routine that reads in the fixed grids data file and sets up the
    ! appropriate data structures
    subroutine set_fixed_grids(fname)

        use amr_module, only: parmunit

        implicit none
        
        ! Subroutine arguments
        character(len=*), optional, intent(in) :: fname
        
        ! Local storage
        integer, parameter :: unit = 7
        integer :: i

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
        read(unit,'(i2)') num_fixed_grids
        write(parmunit,*) '  mfgrids = ',num_fixed_grids
        if (num_fixed_grids == 0) then
            write(parmunit,*) '  No fixed grids specified for output'
            return
        endif
        
        ! Allocate fixed grids (not the data yet though)
        allocate(fgrids(num_fixed_grids))

        ! Read in data for each fixed grid
        do i=1,num_fixed_grids
            ! Read in this grid's data
            read(unit,*) fgrids(i)%start_time, &
                         fgrids(i)%end_time, &
                         fgrids(i)%num_output, &
                         fgrids(i)%x_low, &
                         fgrids(i)%x_hi , &
                         fgrids(i)%y_low, &
                         fgrids(i)%y_hi , &
                         fgrids(i)%mx  , &
                         fgrids(i)%my  , &
                         fgrids(i)%output_arrival_times, &
                         fgrids(i)%output_surface_max

           ! Setup data for this grid
           ! Set dtfg (the timestep length between outputs) for each grid
           if (fgrids(i)%end_time <= fgrids(i)%start_time) then
               if (fgrids(i)%num_output > 1) then 
                  print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
                  print *,'start_time <= end_time yet num_output > 1'
                  print *,'set end_time > start_time or set num_output = 1'
                  stop
               else
                   fgrids(i)%dt = 0.d0
               endif
           else
               if (fgrids(i)%num_output < 2) then
                   print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
                   print *,'end_time > start_time, yet num_output = 1'
                   print *,'set num_output > 2'
                   stop
               else
                   fgrids(i)%dt = (fgrids(i)%end_time - fgrids(i)%start_time) &
                                       / (fgrids(i)%num_output - 1)
               endif
            endif

            ! Initialize last_output_time and index
            fgrids(i)%last_output_time = fgrids(i)%start_time - fgrids(i)%dt
            fgrids(i)%last_output_index = 0

            ! Set spatial intervals dx and dy on each grid
            if (fgrids(i)%mx > 1) then
               fgrids(i)%dx = (fgrids(i)%x_hi - fgrids(i)%x_low) / (fgrids(i)%mx - 1)
            else if (fgrids(i)%mx == 1) then
               fgrids(i)%dx = 0.d0
            else
                 print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
                 print *,'x grid points mx <= 0, set mx >= 1'
            endif

            if (fgrids(i)%my > 1) then
                fgrids(i)%dy = (fgrids(i)%y_hi - fgrids(i)%y_low) / (fgrids(i)%my - 1)
            else if (fgrids(i)%my == 1) then
                fgrids(i)%dy = 0.d0
            else
                print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
                print *,'y grid points my <= 0, set my >= 1'
            endif 
       
            ! set the number of variables stored for each grid
            ! this should be (the number of variables you want to write out + 1)
            fgrids(i)%num_vars(1) = 6
            fgrids(i)%num_vars(2) = 3*fgrids(i)%output_surface_max &
                                          + fgrids(i)%output_arrival_times
            
            ! Allocate new fixed grid data array
            allocate(fgrids(i)%early(fgrids(i)%num_vars(1),fgrids(i)%mx,fgrids(i)%my))
            fgrids(i)%early = nan()
            allocate(fgrids(i)%late(fgrids(i)%num_vars(1),fgrids(i)%mx,fgrids(i)%my))
            fgrids(i)%late = nan()
            allocate(fgrids(i)%often(fgrids(i)%num_vars(2),fgrids(i)%mx,fgrids(i)%my))
            fgrids(i)%often = nan()
       enddo
       close(unit)
       
       tcfmax=-1.d16

    end subroutine set_fixed_grids
    
    
    !=====================FGRIDINTERP=======================================
    ! This routine interpolates q and aux on a computational grid
    ! to a fgrid not necessarily aligned with the computational grid
    ! using bilinear interpolation defined on computation grid
    !=======================================================================
    subroutine fgrid_interp(fgrid_type,fgrid, &
                            t,q,meqn,mxc,myc,mbc,dxc,dyc,xlowc,ylowc, &
                            maux,aux,maxcheck)
    
        use geoclaw_module, only: drytolerance
    
        implicit none
    
        ! Subroutine arguments
        integer, intent(in) :: fgrid_type
        type(fixedgrid_type), intent(inout) :: fgrid
        integer, intent(in) :: meqn,mxc,myc,mbc,maux,maxcheck
        real(kind=8), intent(in) :: t,dxc,dyc,xlowc,ylowc
        real(kind=8), intent(in) :: q(meqn,1-mbc:mxc+mbc,1-mbc:myc+mbc)
        real(kind=8), intent(in) :: aux(maux,1-mbc:mxc+mbc,1-mbc:myc+mbc)
    
        ! Indices
        integer :: ifg,jfg,m,ic1,ic2,jc1,jc2
        integer :: bathy_index,eta_index,arrival_index
        integer :: eta_min_index,eta_max_index,eta_now_index

        ! Tolerances
        real(kind=8), parameter :: arrival_tolerance = 1.d-2
        real(kind=8) :: total_depth,depth_indicator,nan_check

        ! Geometry
        real(kind=8) :: xfg,yfg,xc1,xc2,yc1,yc2,xhic,yhic
        real(kind=8) :: geometry(4)
        
        ! Work arrays for eta interpolation
        real(kind=8) :: eta(2,2),h(2,2)
        
        
        ! Alias to data in fixed grid
        integer :: num_vars
        real(kind=8), pointer :: fg_data(:,:,:)
        
        ! Setup aliases for specific fixed grid
        if (fgrid_type == 1) then
            num_vars = fgrid%num_vars(1)
            fg_data => fgrid%early
        else if (fgrid_type == 2) then
            num_vars = fgrid%num_vars(1)
            fg_data => fgrid%late
        else
            num_vars = fgrid%num_vars(2)
            fg_data => fgrid%often
        endif
            
        xhic = xlowc + dxc*mxc  
        yhic = ylowc + dyc*myc    
        
        ! Find indices of various quantities in the fgrid arrays
        bathy_index = meqn + 1
        eta_index = meqn + 2
    
        if (maxcheck > 0) then 
            arrival_index = 0
            eta_now_index = 0
            eta_min_index = 0
            eta_max_index = 0
    
            if (fgrid%output_arrival_times > 0) then
                arrival_index = 1
            endif
        
            if (fgrid%output_surface_max > 0) then
                eta_now_index = arrival_index + 1
                eta_min_index = arrival_index + 2
                eta_max_index = arrival_index + 3
            endif
        endif
    
        ! Primary interpolation loops 
        do ifg=1,fgrid%mx
            xfg=fgrid%x_low + (ifg-1)*fgrid%dx
            do jfg=1,fgrid%my
                yfg=fgrid%y_low + (jfg-1)*fgrid%dy
    
                ! Check to see if this coordinate is inside of this grid
                if (.not.((xfg < xlowc.or.xfg > xhic).or.(yfg < ylowc.or.yfg > yhic))) then
    
                    ! find where xfg,yfg is in the computational grid and compute the indices
                    ! and relevant coordinates of each corner
                    ic1 = int((xfg-(xlowc+0.5d0*dxc))/(dxc))+1
                    jc1 = int((yfg-(ylowc+0.5d0*dyc))/(dyc))+1
                    if (ic1.eq.mxc) ic1=mxc-1
                    if (jc1.eq.myc) jc1=myc-1 
                    ic2 = ic1 + 1
                    jc2 = jc1 + 1
                        
                    xc1 = xlowc + dxc * (ic1 - 0.5d0)
                    yc1 = ylowc + dyc * (jc1 - 0.5d0)
                    xc2 = xlowc + dxc * (ic2 - 0.5d0)
                    yc2 = ylowc + dyc * (jc2 - 0.5d0)
         
                    ! Calculate geometry of interpolant
                    ! interpolate bilinear used to interpolate to xfg,yfg
                    ! define constant parts of bilinear
                    geometry = [(xfg - xc1) / dxc, &
                                (yfg - yc1) / dyc, &
                                (xfg - xc1) * (yfg - yc1) / (dxc*dyc), &
                                1.d0]
        
                    ! Interpolate for all conserved quantities and bathymetry
                    if (maxcheck == 0) then 
                        forall (m=1:meqn)
                            fg_data(m,ifg,jfg) = interpolate([[q(m,ic1,jc1),q(m,ic1,jc2)], &
                                                              [q(m,ic2,jc1),q(m,ic2,jc2)]], geometry)
                        end forall
                        fg_data(bathy_index,ifg,jfg) = interpolate([[aux(1,ic1,jc1),aux(1,ic1,jc2)], &
                                                                    [aux(1,ic2,jc1),aux(1,ic2,jc2)]], geometry)
                    endif

                    ! If eta max/min are saved on this grid initialize if necessary
                    if (fgrid%output_surface_max > 0 .and. maxcheck == 2) then 
                        if (.not.(fg_data(eta_min_index,ifg,jfg) == fg_data(eta_min_index,ifg,jfg))) then
                            fg_data(eta_min_index,ifg,jfg) = 0.d0
                        endif
                        if (.not.(fg_data(eta_max_index,ifg,jfg) == fg_data(eta_max_index,ifg,jfg))) then
                            fg_data(eta_max_index,ifg,jfg) = 0.d0
                        endif
                    endif
                    
                    ! Interpolate surface eta, only use wet eta points near the shoreline
                    eta(1,:) = [aux(1,ic1,jc1) + q(1,ic1,jc1), aux(1,ic1,jc2) + q(1,ic1,jc2)]
                    eta(2,:) = [aux(1,ic2,jc1) + q(1,ic2,jc1), aux(1,ic2,jc2) + q(1,ic2,jc2)]
                    h(1,:) = [q(1,ic1,jc1),q(1,ic1,jc2)]
                    h(2,:) = [q(1,ic2,jc1),q(1,ic2,jc2)]
                         
                    depth_indicator= min(h(1,1),h(1,2),h(2,1),h(2,2))
                    total_depth = sum(h)
    
                    ! We are near shoreline
                    if (depth_indicator < drytolerance .and. total_depth > 4.d0*drytolerance) then
                        ! Check to see if each cell around fixed grid point is 
                        ! wet, if not re-balance
                        if (h(1,1) < drytolerance) then
                            eta(1,1) =  (h(1,2)*eta(1,2) &
                                       + h(2,1)*eta(2,1) &
                                       + h(2,2)*eta(2,2)) / total_depth
                        endif
                        if (h(1,2) < drytolerance) then
                            eta(1,2) =  (h(1,1)*eta(1,1) &
                                       + h(2,1)*eta(2,1) &
                                       + h(2,2)*eta(2,2)) / total_depth
                        endif
                        if (h(2,1) < drytolerance) then
                            eta(2,1) =  (h(1,1)*eta(1,1) &
                                       + h(1,2)*eta(1,2) &
                                       + h(2,2)*eta(2,2)) / total_depth
                        endif
                        if (h(2,2) < drytolerance) then
                            eta(2,2)=  (h(1,2)*eta(1,2) &
                                      + h(2,1)*eta(2,1) &
                                      + h(1,1)*eta(1,1)) / total_depth
                        endif            
                    endif
                    if (total_depth <= 4.d0*drytolerance) then
                        eta(2,2) = nan()
                    endif
    
                    ! Check which task to perform and evaluate the interpolant
                    ! or evaluate the eta min and max functions
                    if (maxcheck == 0) then 
                        fg_data(eta_index,ifg,jfg) = interpolate(eta,geometry)
                        fg_data(num_vars,ifg,jfg) = t
                    else if (maxcheck.eq.1.and.fgrid%output_surface_max.gt.0) then
                        fg_data(eta_now_index,ifg,jfg) = interpolate(eta,geometry)
                    else if (maxcheck.eq.2.and.fgrid%output_surface_max.gt.0) then
                        fg_data(eta_min_index,ifg,jfg) = min(fg_data(eta_min_index,ifg,jfg), &
                                                             fg_data(eta_now_index,ifg,jfg))
                        fg_data(eta_max_index,ifg,jfg) = max(fg_data(eta_max_index,ifg,jfg), &
                                                             fg_data(eta_now_index,ifg,jfg))            
                    endif
    
                    ! If arrival times are saved on this grid
                    if (maxcheck == 1 .and. fgrid%output_arrival_times > 0) then
                        ! TODO: It would be nice to replace this with an
                        ! intrinsic such as ieee_is_nan but this is not widely
                        ! implemented yet
                        nan_check = fg_data(arrival_index,ifg,jfg)
                        ! if nan_check = NaN: Waves haven't arrived previously
                        if (.not.(nan_check == nan_check)) then
                            if (abs(fg_data(eta_index,ifg,jfg)) > arrival_tolerance) then
                                fg_data(arrival_index,ifg,jfg)= t
                            endif
                        endif
                    endif
                    
                endif ! Enclosing if statement to see if fixed grid point is
                      ! in this computational grid
            enddo ! Fixed grid y-coordinate loop
        enddo ! Fixed grid x-coordinte loop
    
    end subroutine fgrid_interp
    
    ! Interpolation function
    ! Given 4 points (points) and geometry from x,y,and cross terms
    real(kind=8) pure function interpolate(points,geometry)
                            
        implicit none
                                
        ! Function signature
        real(kind=8), intent(in) :: points(2,2)
        real(kind=8), intent(in) :: geometry(4)
        
        ! This is set up as a dot product between the approrpriate terms in 
        ! the input data.  This routine could be vectorized or a BLAS routine
        ! used instead of the intrinsics to ensure that the fastest routine
        ! possible is being used
        interpolate = sum([points(2,1)-points(1,1), &
                           points(1,2)-points(1,1), &
                           points(1,1) + points(2,2) - (points(2,1) + points(1,2)), &
                           points(1,1)] * geometry)
                           
    end function interpolate
    

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
        real(kind=8), intent(inout) :: fgrid1(mvarsfg,1:mxfg,1:myfg)
        real(kind=8), intent(inout) :: fgrid2(mvarsfg,1:mxfg,1:myfg)
        real(kind=8), intent(inout) :: fgrid3(mvarsfg2,1:mxfg,1:myfg)
              
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
                t0=fgrid1(mvarsfg,ifg,jfg)
                tf=fgrid2(mvarsfg,ifg,jfg)
                tau=(toutfg-t0)/(tf-t0)
                
                do iv=1,mvarsfg-1
                    if (dabs(fgrid1(iv,ifg,jfg)) .lt. 1d-90) fgrid1(iv,ifg,jfg) = 0.d0
                    if (dabs(fgrid2(iv,ifg,jfg)) .lt. 1d-90) fgrid2(iv,ifg,jfg) = 0.d0
                enddo
                if (icolumns.eq.mvarsfg-1) then 
                    write(unit,data_format) ((1.d0 - tau)*fgrid1(iv,ifg,jfg)+tau*fgrid2(iv,ifg,jfg), iv=1,mvarsfg-1)
                else
                    if (abs(fgrid3(indetamin,ifg,jfg)) < 1d-90) fgrid3(indetamin,ifg,jfg) = 0.d0
                    if (abs(fgrid3(indetamax,ifg,jfg)) < 1d-90) fgrid3(indetamax,ifg,jfg) = 0.d0
                    write(unit,data_format) ((1.d0 - tau)*fgrid1(iv,ifg,jfg)+tau*fgrid2(iv,ifg,jfg), iv=1,mvarsfg-1), &
                                          fgrid3(indetamin,ifg,jfg), &
                                          fgrid3(indetamax,ifg,jfg)
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
                    write(unit,"(1e26.16)") fgrid3(1,ifg,jfg)
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
