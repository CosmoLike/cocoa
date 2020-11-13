    module ArrayTests
    use ArrayUtils
    implicit none
    contains


    function RunArrayTests() result (fails)
    integer fails
    character(LEN=:), allocatable :: S, S2
    real, allocatable :: arr(:)
    integer, allocatable :: arrI(:)

    fails = 0

    allocate(arr(3), source = [1.,2.,3.])
    call reallocate(arr,2)
    call reallocate(arr,5)
    call reallocate(arrI,5)
    arrI(3)=4
    call reallocate(arrI,8)
    if (arr(2)/=2. .or. arrI(3)/=4) then
        fails = fails+1
        print *, 'Error in realloc'
    end if
    call reallocate(arrI,0)
    if (fails==0) print *, 'OK array tests'


    end function

    end module
