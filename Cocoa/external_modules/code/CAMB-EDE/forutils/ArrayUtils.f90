    module ArrayUtils
    use MiscUtils
    use MpiUtils
    implicit none

    INTERFACE Reallocate
    module procedure Realloc_R, Realloc_D, Realloc_I
    END INTERFACE Reallocate

    contains

    subroutine realloc_R(arr, new_size, keep)
    real, allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: new_size
    logical, intent(in), optional :: keep
    integer sz
    real, allocatable :: tmp(:)

    if (.not. allocated(arr)) then
        allocate(arr(new_size))
        return
    end if
    if (LBOUND(arr,1)/=1) call MpiStop('Realloc only works on arrays starting at 1')
    if (new_size/=size(arr)) then
        allocate(tmp(new_size))
        if (DefaultTrue(keep)) then
            sz = min(new_size, size(arr))
            tmp(:sz) = arr(:sz)
        end if
        call move_alloc(tmp, arr)
    end if

    end subroutine realloc_R

    subroutine realloc_D(arr, new_size, keep)
    double precision, allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: new_size
    logical, intent(in), optional :: keep
    integer sz
    double precision, allocatable :: tmp(:)

    if (.not. allocated(arr)) then
        allocate(arr(new_size))
        return
    end if
    if (LBOUND(arr,1)/=1) call MpiStop('Realloc only works on arrays starting at 1')
    if (new_size/=size(arr)) then
        allocate(tmp(new_size))
        if (DefaultTrue(keep)) then
            sz = min(new_size, size(arr))
            tmp(:sz) = arr(:sz)
        end if
        call move_alloc(tmp, arr)
    end if

    end subroutine realloc_D

    subroutine realloc_I(arr, new_size, keep)
    integer, allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: new_size
    logical, intent(in), optional :: keep
    integer sz
    integer, allocatable :: tmp(:)

    if (.not. allocated(arr)) then
        allocate(arr(new_size))
        return
    end if
    if (LBOUND(arr,1)/=1) call MpiStop('Realloc only works on arrays starting at 1')
    if (new_size/=size(arr)) then
        allocate(tmp(new_size))
        if (DefaultTrue(keep)) then
            sz = min(new_size, size(arr))
            tmp(:sz) = arr(:sz)
        end if
        call move_alloc(tmp, arr)
    end if

    end subroutine realloc_I

    function IndexOf(aval,arr, n)
    integer, intent(in) :: n, arr(n), aval
    integer IndexOf, i

    do i=1,n
        if (arr(i)==aval) then
            IndexOf= i
            return
        end if
    end do
    IndexOf = 0

    end function IndexOf

    function MaxIndex(arr, n)
    integer, intent(in) :: n
    real, intent(in) :: arr(n)
    integer locs(1:1), MaxIndex

    locs = maxloc(arr(1:n))
    MaxIndex = locs(1)

    end function MaxIndex


    function MinIndex(arr, n)
    integer, intent(in) :: n
    real, intent(in) :: arr(n)
    integer locs(1:1), MinIndex

    locs = minloc(arr(1:n))
    MinIndex = locs(1)

    end function MinIndex

    end module ArrayUtils
