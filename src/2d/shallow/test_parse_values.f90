

program test_parse_values

    ! To test: select a test below and then:
    !    gfortran utility_module.f90 test_parse_values.f90 
    !    ./a.out

    use utility_module, only: parse_values
    use topo_module, only: read_topo_header

    implicit none
    character(len=150) str,fname
    integer :: n,i,topo_type,mx,my
    real(kind=8) :: values(10), xll,yll,xhi,yhi,dx,dy

    if (.true.) then
        write(6,*) 'input line with mix of character strings and numbers...'
        read(5,'(a)') str
        write(6,*) 'read in str  = ',str

        call parse_values(trim(str), n, values)
        write(6,*) 'n = ',n
        write(6,*) 'values = ',(values(i), i=1,n)
        write(6,*) 'integer values = ',(nint(values(i)), i=1,n)
        endif

    if (.false.) then
        ! requires also compiling amr_module, topo_module, and setting path:
        fname = '/Users/rjl/topo/etopo/etopo4min100E69W70S70N.asc'
        !fname = '/Users/rjl/topo/etopo/etopo4min180E65W65S0N.asc'
        topo_type = 3
        call read_topo_header(fname,topo_type,mx,my,xll,yll,xhi,yhi,dx,dy)
        write(6,*) 'mx = ',mx
        write(6,*) 'my = ',my
        write(6,*) 'xll = ',xll
        write(6,*) 'yll = ',yll
        write(6,*) 'dx = ',dx
        write(6,*) 'dy = ',dy
        endif

end program test_parse_values
