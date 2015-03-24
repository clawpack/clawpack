! Contains a number of useful functions
module utility_module

    implicit none
    
Contains

    ! Returns number of arguments in list, assumes arbitrary number of white 
    ! space as delimiter
    integer function get_value_count(line) result(count)

        implicit none

        character(len=*), intent(in) :: line

        character(len=200) :: search_buffer
        integer :: i
        logical :: found

        count = 1

        search_buffer = trim(line(1:index(line,"=:") - 2))

        do i=1,len_trim(search_buffer)
            if (search_buffer(i:i) == " ") then
                if (.not. found) then
                    found = .true.
                    count = count + 1
                endif
            else
                found = .false.
            endif
        enddo

    end function get_value_count

    ! Converts seconds to days truncating at the second decimal place
    real(kind=8) pure function convert2days(seconds) result(days)
        
        implicit none
        real(kind=8), intent(in) :: seconds

        days = real(int(seconds * 1.d2 / 8.64d4) / 1.d2, kind=8)

    end function convert2days


    !======================================
    subroutine parse_values(str, n, values)
    !======================================

    ! Take the input string `str` and parse it to extract any numerical values,
    ! ignoring any character strings. 
    ! Returns `n` values in array `values`. 
    ! Assumes n <= 10.
    ! If you expect value(i) to be an integer, evaluate as nint(value(i)).

    implicit none

    character(len=*), intent(in) :: str
    integer, intent(out) :: n
    real(kind=8), intent(out) :: values(10)

    integer :: pos2,nw,i,e
    character(len=80) :: word(10), str2
    real(kind=8) :: x

    ! First break into words / tokens based on white space.  
    ! Each might be character or numerical:

    nw = 0
    str2 = trim(adjustl(str))
    do while (len(trim(adjustl(str2))) > 0) 
        pos2 = index(str2, " ")
        
        if (pos2 == 0) then
           nw = nw + 1
           word(nw) = trim(adjustl(str2))
           exit
           endif

        nw = nw + 1
        word(nw) = trim(adjustl(str2(1:pos2-1)))
        str2 = trim(adjustl(str2(pos2+1:)))
        enddo

    ! now extract numerical values:
    n = 0
    do i=1,nw
        read(word(i),*,IOSTAT=e) x
        if (e == 0) then
            ! this token is numerical
            n = n+1
            values(n) = x
            endif
        enddo


    end subroutine parse_values

end module utility_module
