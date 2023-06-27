!
! Module for working with string arrays of key/value pairs
!
module keysvalues

    public
    integer, parameter :: MAXSTR=512

contains

    function getStrParam(keys,values,k)
        character(*), dimension(:) :: keys
        character(*), dimension(:) :: values
        character(*) :: k

        character(len=MAXSTR) :: getStrParam

        integer :: i

        do i=1,size(keys)
            if (trim(keys(i))==trim(k)) then
                 getStrParam=trim(values(i))
                 return
            end if
        end do

        getStrParam=""

    end function


    function getRealParam(keys,values,k)
        character(*), dimension(:) :: keys
        character(*), dimension(:) :: values
        character(*) :: k

        real(8) :: getRealParam

        character(MAXSTR) :: buf

        getRealParam = 0

        buf = getStrParam(keys,values,k)
        read (buf,*) getRealParam

    end function


    function getRealArrParam(keys,values,k,n)
        character(*), dimension(:) :: keys
        character(*), dimension(:) :: values
        character(*) :: k
        integer :: n

        real(8), dimension(n) :: getRealArrParam

        character(MAXSTR) :: buf

        getRealArrParam = 0

        buf = getStrParam(keys,values,k)
        read (buf,*) getRealArrParam

    end function


    function getIntParam(keys,values,k)
        character(*), dimension(:) :: keys
        character(*), dimension(:) :: values
        character(*) :: k

        integer :: getIntParam

        character(MAXSTR) :: buf

        getIntParam = 0

        buf = getStrParam(keys,values,k)
        read (buf,*) getIntParam

    end function

    function hasParam(keys,values,k)
        character(*), dimension(:) :: keys
        character(*), dimension(:) :: values
        character(*) :: k

        logical :: hasParam

        integer :: i

        do i=1,size(keys)
            if (trim(keys(i))==trim(k)) then
                hasParam=.true.
                return
            end if
        end do

        hasParam=.false.

    end function

    subroutine setRealParam(keys,values,k,v)
        character(*), dimension(:) :: keys
        character(*), dimension(:) :: values
        character(*) :: k
        real(8) :: v

        character(MAXSTR) :: buf

        integer :: i

        write (buf,*) v

        do i=1,size(keys)
            if (trim(keys(i))==trim(k)) then
                 values(i)=trim(buf)
                 return
            end if
        end do

    end subroutine


end module
