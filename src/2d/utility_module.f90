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

end module utility_module