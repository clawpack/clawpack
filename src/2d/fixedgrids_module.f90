module fixedgrids_module

    implicit none
    save

    ! Container for fixed grid data, geometry and output settings
    type fixedgrid_type
        ! Grid data
        real(kind=8), allocatable :: early(:,:,:)
        real(kind=8), allocatable :: late(:,:,:)
        real(kind=8), allocatable :: often(:,:,:)
        
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
        real(kind=8), intent(inout) :: fgrid(mvarsfg,1:mxfg,1:myfg)
        real(kind=8), intent(in) :: q(meqn,1-mbc:mxc+mbc,1-mbc:myc+mbc)
        real(kind=8), intent(in) :: aux(maux,1-mbc:mxc+mbc,1-mbc:myc+mbc)
    
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
                            z11=q(m,ic1,jc1)
                            z21=q(m,ic2,jc1)
                            z12=q(m,ic1,jc2)
                            z22=q(m,ic2,jc2)
                            a=z21-z11
                            b=z12-z11
                            d=z11
                            c=z22-(a+b+d)
                            fgrid(m,ifg,jfg) = a*xterm + b*yterm + c*xyterm + d
                        enddo
                        z11=aux(1,ic1,jc1)
                        z21=aux(1,ic2,jc1)
                        z12=aux(1,ic1,jc2)
                        z22=aux(1,ic2,jc2)
                        a=z21-z11
                        b=z12-z11
                        d=z11
                        c=z22-(a+b+d)
                        fgrid(indb,ifg,jfg) = a*xterm + b*yterm + c*xyterm + d
                    endif
                    ! This next output variable is the surface using bilinear interpolation,
                    ! using a surface that only uses the wet eta points near the shoreline
        
                    z11 = aux(1,ic1,jc1) + q(1,ic1,jc1)
                    z21 = aux(1,ic2,jc1) + q(1,ic2,jc1)
                    z12 = aux(1,ic1,jc2) + q(1,ic1,jc2)
                    z22 = aux(1,ic2,jc2) + q(1,ic2,jc2)
                        
                    h11 = q(1,ic1,jc1)
                    h21 = q(1,ic2,jc1)
                    h12 = q(1,ic1,jc2)
                    h22 = q(1,ic2,jc2)
                    depthindicator= min(h11,h12,h21,h22)
                    totaldepth= h11+h22+h21+h12
    
                    ! Near shoreline
                    if (depthindicator.lt.tol.and.totaldepth.gt.4.d0*tol) then
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
                    if (ioutsurfacemax > 0 .and. maxcheck == 2) then 
                        if (.not.(fgrid(indetamin,ifg,jfg) == fgrid(indetamin,ifg,jfg))) then
                            fgrid(indetamin,ifg,jfg) = 0.d0
                        endif
                        if (.not.(fgrid(indetamax,ifg,jfg) == fgrid(indetamax,ifg,jfg))) then
                            fgrid(indetamax,ifg,jfg) = 0.d0
                        endif
                    endif
    
                    ! check which task to perform
                    if (maxcheck == 0) then 
                        fgrid(indeta,ifg,jfg) = a*xterm + b*yterm + c*xyterm + d
                        fgrid(mvarsfg,ifg,jfg) = t
                    else if (maxcheck.eq.1.and.ioutsurfacemax.gt.0) then
                        fgrid(indetanow,ifg,jfg) = a*xterm + b*yterm + c*xyterm + d   
                    else if (maxcheck.eq.2.and.ioutsurfacemax.gt.0) then
                        fgrid(indetamin,ifg,jfg) = min(fgrid(indetamin,ifg,jfg),fgrid(indetanow,ifg,jfg))
                        fgrid(indetamax,ifg,jfg) = max(fgrid(indetamax,ifg,jfg),fgrid(indetanow,ifg,jfg))            
                    endif
    
                    ! If arrival times are saved on this grid
                    if (maxcheck == 1 .and. ioutarrivaltimes > 0) then
                        check=fgrid(indarrive,ifg,jfg)
                        !# check=NaN: Waves haven't arrived previously
                        if (.not.(check == check)) then
                            if (abs(fgrid(indeta,ifg,jfg)) > arrivaltol) then
                                fgrid(indarrive,ifg,jfg)= t
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
