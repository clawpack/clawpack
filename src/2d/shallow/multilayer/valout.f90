!! Geoclaw specific output - adds eta to q array before writing out
!!
!! Write the results to the file fort.q<iframe>
!! Use format required by matlab script  plotclaw2.m or Python tools
!!
!! set outaux = .true. to also output the aux arrays to fort.a<iframe>
subroutine valout(level_begin, level_end, time, num_eqn, num_aux)

    use amr_module, only: alloc, t0, output_aux_onlyonce, output_aux_components
    use amr_module, only: frame => matlabu, num_ghost => nghost, lstart
    use amr_module, only: hxposs, hyposs, output_format, store1, storeaux
    use amr_module, only: node, rnode, ndilo, ndihi, ndjlo, ndjhi
    use amr_module, only: cornxlo, cornylo, levelptr, mxnest
    use amr_module, only: timeValout, timeValoutCPU, tvoll, tvollCPU, rvoll
    use amr_module, only: timeTick, tick_clock_start, t0

    use storm_module, only: storm_specification_type, output_storm_location
    use storm_module, only: output_storm_location
    use storm_module, only: landfall, display_landfall_time

    use geoclaw_module, only: rho
    use multilayer_module, only: num_layers

#ifdef HDF5
    use hdf5
#endif

    implicit none

    ! Input
    integer, intent(in) :: level_begin, level_end, num_eqn, num_aux
    real(kind=8), intent(in) :: time

    ! Locals
    logical :: timing_file_exists
    integer, parameter :: out_unit = 50
    integer :: i, j, k, m, level, output_aux_num, num_stop, digit
    integer :: index, grid_ptr, num_cells(2), num_grids, q_loc, aux_loc
    real(kind=8) :: lower_corner(2), delta(2)
    logical :: out_aux
    character(len=11) :: file_name(5)

    real(kind=8) :: h(num_layers), hu(num_layers), hv(num_layers)
    real(kind=8) :: eta(num_layers)
    real(kind=8), allocatable :: qeta(:)

#ifdef HDF5
    ! HDF writing
    integer :: hdf_error
    integer(hid_t) :: hdf_file, data_space, data_set
    integer(hsize_t) :: dims(2)
#endif

    ! Timing
    integer :: clock_start, clock_finish, clock_rate
    integer    tick_clock_finish, tick_clock_rate, timeTick_int
    real(kind=8) :: cpu_start, cpu_finish, t_CPU_overall, timeTick_overall
    character(len=128) :: console_format
    character(len=256) :: timing_line, timing_substr
    character(len=*), parameter :: timing_file_name = "timing.csv"

    character(len=*), parameter :: header_format =                             &
                                    "(i6,'                 grid_number',/," // &
                                     "i6,'                 AMR_level',/,"   // &
                                     "i6,'                 mx',/,"          // &
                                     "i6,'                 my',/"           // &
                                     "e26.16,'    xlow', /, "               // &
                                     "e26.16,'    ylow', /,"                // &
                                     "e26.16,'    dx', /,"                  // &
                                     "e26.16,'    dy',/)"
    character(len=*), parameter :: t_file_format = "(e18.8,'    time', /,"  // &
                                           "i6,'                 meqn'/,"   // &
                                           "i6,'                 ngrids'/," // &
                                           "i6,'                 naux'/,"   // &
                                           "i6,'                 ndim'/,"   // &
                                           "i6,'                 nghost'/,/)"

    ! Output timing
    call system_clock(clock_start,clock_rate)
    call cpu_time(cpu_start)

    ! Count how many aux components requested
    output_aux_num = 0
    do i=1, num_aux
        output_aux_num = output_aux_num + output_aux_components(i)
    end do

    ! Note:  Currently outputs all aux components if any are requested
    out_aux = ((output_aux_num > 0) .and.               &
              ((.not. output_aux_onlyonce) .or. (abs(time - t0) < 1d-90)))

    ! Output storm track if needed
    if (storm_specification_type /= 0) then
        call output_storm_location(time)
    end if

    ! Construct file names
    file_name(1) = 'fort.qxxxx'
    file_name(2) = 'fort.txxxx'
    file_name(3) = 'fort.axxxx'
    file_name(4) = 'fort.bxxxx'
    num_stop = frame
    do i = 10, 7, -1
        digit = mod(num_stop, 10)
        do j = 1, 4
            file_name(j)(i:i) = char(ichar('0') + digit)
        end do
        num_stop = num_stop / 10
    end do
    ! Slightly modified for HDF file output
    file_name(5) = 'clawxxxx.h5'
    num_stop = frame
    do i = 8, 5, -1
        digit = mod(num_stop, 10)
        file_name(5)(i:i) = char(ichar('0') + digit)
        num_stop = num_stop / 10
    end do

    ! ==========================================================================
    ! Write out fort.q file (and fort.bXXXX and clawxxxx.h5 files if necessary)
    ! Here we let fort.q be out_unit and the the other two be out_unit + 1
    open(unit=out_unit, file=file_name(1), status='unknown', form='formatted')
    if (output_format == 3) then
        open(unit=out_unit + 1, file=file_name(4), status="unknown",    &
             access='stream')
    else if (output_format == 4) then
#ifdef HDF5
        ! Note that we will use this file for both q and aux data
        call h5create_f(file_name(5), H5F_ACC_TRUNC_F, hdf_file, hdf_error)

        ! Create group for q
        call h5gcreate_f(hdf_file, "/q", q_group, hdf_error)
#endif
    end if
    num_grids = 0

    ! Loop over levels
    do level = level_begin, level_end
        grid_ptr = lstart(level)
        delta = [hxposs(level), hyposs(level)]

        ! Loop over grids on each level
        do while (grid_ptr /= 0)
            ! Extract grid data
            num_grids = num_grids + 1
            num_cells(1) = node(ndihi, grid_ptr) - node(ndilo, grid_ptr) + 1
            num_cells(2) = node(ndjhi, grid_ptr) - node(ndjlo, grid_ptr) + 1
            q_loc = node(store1, grid_ptr)
            aux_loc = node(storeaux, grid_ptr)
            lower_corner = [rnode(cornxlo, grid_ptr), rnode(cornylo, grid_ptr)]

            ! Write out header data
            write(out_unit, header_format) grid_ptr, level,             &
                                           num_cells(1),                &
                                           num_cells(2),                &
                                           lower_corner(1),             &
                                           lower_corner(2),             &
                                           delta(1), delta(2)

            ! Output grids
            select case(output_format)
                ! ASCII output
                case(1)
                    ! Round off if nearly zero
                    forall (m = 1:num_eqn,                              &
                            i=num_ghost + 1:num_cells(1) + num_ghost,   &
                            j=num_ghost + 1:num_cells(2) + num_ghost,   &
                            abs(alloc(iadd(m, i, j))) < 1d-90)

                        alloc(iadd(m, i, j)) = 0.d0
                    end forall

                    do j = num_ghost + 1, num_cells(2) + num_ghost
                        do i = num_ghost + 1, num_cells(1) + num_ghost

                            ! Extract depth and momenta
                            do k=1, num_layers
                                index = 3 * (k - 1)
                                h(k) = alloc(iadd(index + 1, i, j)) / rho(k)
                                hu(k) = alloc(iadd(index + 2, i, j)) / rho(k)
                                hv(k) = alloc(iadd(index + 3, i, j)) / rho(k)
                            end do

                            ! Calculate sufaces
                            eta(num_layers) = h(num_layers)         &
                                                + alloc(iaddaux(1, i, j))
                            do k=num_layers - 1, 1, -1
                                eta(k) = h(k) + eta(k + 1)
                                if (abs(eta(k)) < 1d-99) then
                                    eta(k) = 0.d0
                                end if
                            end do

                            write(out_unit, "(50e26.16)")                   &
                                (h(k), hu(k), hv(k), k=1, num_layers),      &
                                (eta(k), k=1, num_layers)
                        end do
                        write(out_unit, *) ' '
                    end do

                ! What is case 2?
                case(2)
                    stop "Unknown format."

                ! Binary output
                case(3)
                    ! Need to add eta to the output data
                    allocate(qeta((num_eqn + num_layers)                &
                             * (num_cells(1) + 2 * num_ghost)           &
                             * (num_cells(2) + 2 * num_ghost)))
                    do j = 1, num_cells(2) + 2 * num_ghost
                        do i = 1, num_cells(1) + 2 * num_ghost

                            ! Extract depth and momenta
                            do k=1,num_layers
                                 index = 3 * (k - 1)
                                 h(k) = alloc(iadd(index + 1,i,j)) / rho(k)
                                 hu(k)= alloc(iadd(index + 2,i,j)) / rho(k)
                                 hv(k) = alloc(iadd(index + 3,i,j)) / rho(k)
                            end do

                            ! Calculate surfaces
                            eta(num_layers) = h(num_layers)         &
                                                + alloc(iaddaux(1,i,j))
                            do k=num_layers-1,1,-1
                                eta(k) = h(k) + eta(k+1)
                            enddo

                            do m=1,num_layers
                                index = 3*(m - 1)
                                qeta(iaddqeta(index+1,i,j)) = h(m)
                                qeta(iaddqeta(index+2,i,j)) = hu(m)
                                qeta(iaddqeta(index+3,i,j)) = hv(m)
                                qeta(iaddqeta(3*num_layers+m,i,j)) = eta(m)
                            end do

                        end do
                    end do

                    ! Note: We are writing out ghost cell data also
                    write(out_unit + 1) qeta

                    deallocate(qeta)

                ! HDF5 output
                case(4)
#ifdef HDF5
                ! Create data space - handles dimensions of the corresponding 
                ! data set - annoyingling need to stick grid size into other
                ! data type
                dims = (/ num_eqn, num_cells(1) + 2 * num_ghost,               &
                                   num_cells(2) + 2 * num_ghost /)
                call h5screate_simple_f(2, dims, data_space, hdf_error)

                ! Create new dataset for this grid
                call h5dcreate_f(hdf_file, data_set, H5T_IEEE_F64LE,           &
                                 data_space, data_set, hdf_error)

                ! Write q into file
                i = (iadd(num_eqn, num_cells(1) + 2 * num_ghost,               &
                                   num_cells(2) + 2 * num_ghost))
                call h5dwrite_f(data_set, H5T_NATIVE_DOUBLE,                   &
                                alloc(iadd(1, 1, 1):i), hdf_error)
                call h5dclose_f(data_set, hdf_error)
                call h5sclose_f(data_space, hdf_error)
#endif
                case default
                    print *, "Unsupported output format", output_format,"."
                    stop 

            end select
            grid_ptr = node(levelptr, grid_ptr)
        end do
    end do
    close(out_unit)
#ifdef HDF5
    if (output_format == 4) then
        call h5gclose_f(q_group, hdf_error)
    end if
#endif

    ! ==========================================================================
    ! Write out fort.a file
    if (out_aux) then
        if (output_format == 1) then
            open(unit=out_unit, file=file_name(3), status='unknown',        &
                 form='formatted')
        else if (output_format == 3) then
            open(unit=out_unit, file=file_name(3), status='unknown',        &
                 access='stream')
        else if (output_format == 4) then
#ifdef HDF5
        ! Create group for aux
        call h5gcreate_f(hdf_file, "/aux", aux_group, hdf_error)
#endif            
        end if

        do level = level_begin, level_end
            grid_ptr = lstart(level)
            delta = [hxposs(level), hyposs(level)]

            ! Loop over grids on each level
            do while (grid_ptr /= 0)
                ! Extract grid data
                num_cells(1) = node(ndihi, grid_ptr) - node(ndilo, grid_ptr) + 1
                num_cells(2) = node(ndjhi, grid_ptr) - node(ndjlo, grid_ptr) + 1
                aux_loc = node(storeaux, grid_ptr)
                lower_corner = [rnode(cornxlo, grid_ptr),           &
                                rnode(cornylo, grid_ptr)]

                ! Output grids
                select case(output_format)
                    ! ASCII output
                    case(1)

                        ! We only output header info for aux data if writing 
                        ! ASCII data
                        write(out_unit, header_format) grid_ptr, level, &
                                                       num_cells(1),    &
                                                       num_cells(2),    &
                                                       lower_corner(1), &
                                                       lower_corner(2), &
                                                       delta(1), delta(2)

                        ! Round off if nearly zero
                        forall (m = 1:num_aux,                              &
                                i=num_ghost + 1:num_cells(1) + num_ghost,   &
                                j=num_ghost + 1:num_cells(2) + num_ghost,   &
                                abs(alloc(iaddaux(m, i, j))) < 1d-90)

                            alloc(iaddaux(m, i, j)) = 0.d0
                        end forall

                        do j = num_ghost + 1, num_cells(2) + num_ghost
                            do i = num_ghost + 1, num_cells(1) + num_ghost
                                write(out_unit, "(50e26.16)")                   &
                                         (alloc(iaddaux(m, i, j)), m=1, num_aux)
                            end do
                            write(out_unit, *) ' '
                        end do

                    ! What is case 2?
                    case(2)
                        stop "Unknown format."

                    ! Binary output
                    case(3)
                        ! Note: We are writing out ghost cell data also
                        i = (iaddaux(num_aux, num_cells(1) + 2 * num_ghost, &
                                              num_cells(2) + 2 * num_ghost))
                        write(out_unit) alloc(iaddaux(1, 1, 1):i)

                    ! HDF5 output
                    case(4)
#ifdef HDF5
                ! Create data space - handles dimensions of the corresponding 
                ! data set - annoyingling need to stick grid size into other
                ! data type
                dims = (/ num_aux, num_cells(1) + 2 * num_ghost,               &
                                   num_cells(2) + 2 * num_ghost /)
                call h5screate_simple_f(2, dims, data_space, hdf_error)

                ! Create new dataset for this grid
                call h5dcreate_f(hdf_file, data_set, H5T_IEEE_F64LE,           &
                                 data_space, data_set, hdf_error)

                call h5dcreate_f(hdf_file, data_set, H5T_IEEE_F64LE,           &
                                 data_space, data_set, hdf_error)

                ! Write q into file
                i = (iadd_aux(num_aux, num_cells(1) + 2 * num_ghost,           &
                                       num_cells(2) + 2 * num_ghost))
                call h5dwrite_f(data_set, H5T_NATIVE_DOUBLE,                   &
                                alloc(iadd_aux(1, 1, 1):i), hdf_error)
                call h5dclose_f(data_set, hdf_error)
                call h5sclose_f(data_space, hdf_error)
#endif
                    case default
                        print *, "Unsupported output format", output_format,"."
                        stop 

                end select
                grid_ptr = node(levelptr, grid_ptr)
            end do
        end do
    end if
#ifdef HDF5
    if (out_aux) then
        call h5gclose_f(aux_group, hdf_error)
    end if
    call h5fclose_f(hdf_file, hdf_error)
#endif


    ! ==========================================================================
    ! Write fort.t file
    open(unit=out_unit, file=file_name(2), status='unknown', form='formatted')

    ! Note:  We need to print out num_ghost too in order to strip ghost cells
    !        from q array when reading in pyclaw.io.binary
    write(out_unit, t_file_format) time, num_eqn + num_layers, num_grids, &
                                   num_aux, 2, num_ghost
    close(out_unit)

    ! ==========================================================================
    ! Write out timing stats
    open(unit=out_unit, file=timing_file_name, form='formatted',         &
             status='old', action='write', position='append')
    
    timing_line = "(e16.6, ', ', e16.6, ', ', e16.6,"
    do level=1, mxnest
        timing_substr = "', ', e16.6, ', ', e16.6, ', ', e16.6"
        timing_line = trim(timing_line) // timing_substr
    end do
    timing_line = trim(timing_line) // ")"

    if (time == t0) then
        t_CPU_overall = 0.d0
        timeTick_overall = 0.d0
      else
        call cpu_time(t_CPU_overall)
        call system_clock(tick_clock_finish,tick_clock_rate)
        timeTick_int = timeTick + tick_clock_finish - tick_clock_start
        timeTick_overall = real(timeTick_int, kind=8)/real(clock_rate,kind=8)
      endif

    write(out_unit, timing_line) time, timeTick_overall, t_CPU_overall, &
        (real(tvoll(i), kind=8) / real(clock_rate, kind=8), &
         tvollCPU(i), rvoll(i), i=1,mxnest)
    
    close(out_unit)

    ! ==========================================================================
    ! Print output info
    if (display_landfall_time) then
        ! Convert time to days relative to landfall
        console_format = "('AMRCLAW: Frame ',i4,' output files done at " // &
                         "time t = ', f5.2,/)"
        print console_format, frame, time / (3.3d3 * 24.d0)
    else
        console_format = "('AMRCLAW: Frame ',i4,' output files done at " // &
                         "time t = ', d13.6,/)"
        print console_format, frame, time
    end if

    ! Increment frame counter
    frame = frame + 1

    ! Ouptut timing
    call system_clock(clock_finish,clock_rate)
    call cpu_time(cpu_finish)
    timeValout = timeValout + clock_finish - clock_start
    timeValoutCPU = timeValoutCPU + cpu_finish - cpu_start
    

contains

    ! Index into q array
    pure integer function iadd(m, i, j)
        implicit none
        integer, intent(in) :: m, i, j
        iadd = q_loc + m - 1 + num_eqn * ((j - 1) * (num_cells(1) + 2 * num_ghost) + i - 1)
    end function iadd

    ! Index into aux array
    pure integer function iaddaux(m, i, j)
        implicit none
        integer, intent(in) :: m, i, j
        iaddaux = aux_loc + m - 1 + num_aux * (i - 1) + num_aux * (num_cells(1) + 2 * num_ghost) * (j - 1)
    end function iaddaux

    ! Index into qeta for binary output
    ! Note that this implicitly assumes that we are outputting only h, hu, hv
    ! and will not output more (change num_eqn parameter above)
    pure integer function iaddqeta(m, i, j)
        implicit none
        integer, intent(in) :: m, i, j
        iaddqeta = 1 + m - 1 + (num_eqn + num_layers) * ((j - 1) * (num_cells(1) + 2 * num_ghost) + i - 1) 
    end function iaddqeta

end subroutine valout
