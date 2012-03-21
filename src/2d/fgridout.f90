subroutine fgridout(fgrid1,fgrid2,fgrid3, &
                    x_low,x_hi,y_low,y_hi,mx,my,num_vars,num_vars_2, &
                    t_out,num_out,grid_number, &
                    out_arrival,out_flag)

    use geoclaw_module

    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: num_vars,num_vars_2,mx,my
    integer, intent(in) :: num_out,grid_number,out_arrival,out_flag
    real(kind=8), intent(in) :: x_low,x_hi,y_low,y_hi,t_out
    real(kind=8), intent(inout), dimension(num_vars,1:mx,1:my) :: fgrid1,fgrid2,fgrid3

    ! File handling
    character(len=30) :: fg_file_name
    integer, parameter :: unit = 90
    integer :: digit,nout,ng,pos,status
    
    ! Other locals
    integer :: i,j,m,deta_min,deta_max,num_columns
    real(kind=8) :: t0,tf,tau
    
    ! Format strings
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
    character(len=*), parameter :: arrival_data_format = "(1e26.16)"

    ! Open file for writing                               
    fg_file_name = 'fort.fgnn_xxxx'
    ng = grid_number
    do pos=9,8,-1
        digit = mod(ng,10)
        fg_file_name(pos:pos) = char(ichar('0') + digit)
        ng = ng / 10
    enddo
    
    nout = num_out
    do pos=14,11,-1
        digit = mod(nout,10)
        fg_file_name(pos:pos) = char(ichar('0') + digit)
        nout = nout / 10
    enddo
    
    open(unit=unit, file=fg_file_name, iostat=status, status="new", action="write")
    if ( status /= 0 ) then
        print *, "Error opening file fixed grid output file ",fg_file_name
        stop
    endif
    
    ! Determine the number of columns of the output file
    num_columns = num_vars - 1
    if (num_vars_2 > 1) then
        num_columns = num_columns + 2
    endif
    deta_min = out_arrival + 2
    deta_max = out_arrival + 3
    
    ! Write out header
    write(unit,header_format) t_out,mx,my,x_low,y_low,x_hi,y_hi,num_columns

    ! Interpolate the grid in time to hit the requested output time using the
    ! solution in fgrid1 and fgrid2 which represent the solution on the fixed
    ! grid at the two nearest computational times
    do j=1,my
        do i=1,mx
            ! Figure out time interpolant
            t0 = fgrid1(num_vars,i,j)
            tf = fgrid2(num_vars,i,j)
            tau = (t_out - t0) / (tf - t0)
            
            ! Zero out values that are too small
            do m=1,num_vars-1
                if (abs(fgrid1(m,i,j)) < 1d-90) then
                    fgrid1(m,i,j) = 0.d0
                endif
                if (abs(fgrid2(m,i,j)) < 1d-90) then
                    fgrid2(m,i,j) = 0.d0
                endif
            enddo
            
            ! Write out interpolants
            if (num_columns == num_vars - 1) then
                write(unit,data_format) ((1.d0 - tau) * fgrid1(m,i,j) &
                            + tau * fgrid2(m,i,j),m=1,num_vars-1)
            else
                if (abs(fgrid3(deta_min,i,j)) < 1d-90) then
                    fgrid3(deta_min,i,j) = 0.d0
                endif
                if (abs(fgrid3(deta_max,i,j)) < 1d-90) then
                    fgrid3(deta_max,i,j) = 0.d0
                endif
                
                write(unit,data_format) ((1.d0 - tau) * fgrid1(m,i,j) &
                        + tau * fgrid2(m,i,j), m=1,num_vars-1), &
                        fgrid3(deta_min,i,j),fgrid3(deta_max,i,j)
            endif
        enddo
    enddo
    
    close(unit)
    print "(' FGRIDOUT: Fixed Grid  ', i2, '  output at time =', e18.8)",grid_number,t_out
    
    ! Output arrival times
    if (out_flag == 1) then
        ! Open output file
        fg_file_name = 'fort.fgnn_arrivaltimes'
        ng = grid_number
        do pos=9,8,-1
            digit = mod(ng,10)
            fg_file_name(pos:pos) = char(ichar('0') + digit)
            ng = ng / 10
        enddo
        
        open(unit=unit, file=fg_file_name, iostat=status, status="new", action="write")
        if ( status /= 0 ) then
            print *,"Error opening file for arrival times file ",fg_file_name
            stop
        endif
        
        write(unit,arrival_header_format) mx,my,x_low,y_low,x_hi,y_hi
        do j=1,my
            do i=1,mx
                write(unit,arrival_data_format) fgrid3(1,i,j)
            enddo
        enddo
        
        close(unit)
        print "(' FGRIDOUT: Fixed Grid  ', i2, '  arrival times output')",grid_number
    endif
    
end subroutine fgridout