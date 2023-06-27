    module ListTests
    use ObjectLists
    use FileUtils
    use FileTests
    implicit none
    character(LEN=0), target :: Empty_String = ''

    contains

    function RunListTests() result (fails)
    integer fails
    fails = 0
    call test_TRealArrayList(fails)
    call test_TRealList(fails)
    call test_TStringLIst(fails)
    call test_TObjectList(fails)
    end function RunListTests

    subroutine test_TObjectList(fails)
    !This will stop if internal error
    integer fails
    Type(TObjectList) :: L
    Type(TBinaryFile) F
    double precision arr(9)
    real arr2(9)
    integer arr3(3)
    logical arr4(1)

    arr=1
    call L%Add(arr)
    arr2=2
    call L%Add(arr2)
    arr3=3
    call L%Add(arr3)
    arr4 = .true.
    call L%Add(arr4)
    call checkL()

    call F%CreateFile(temp_file)
    call L%SaveBinary(F%unit)
    call F%Close()
    call L%Clear()
    call F%Open(temp_file)
    call L%ReadBinary(F%unit)
    call F%Close()
    call F%CreateFile(temp_file)
    call L%SaveBinary(F%unit)
    call F%Close(del = .true.)
    call checkL()
    contains

    subroutine checkL()
    class(*), pointer :: PP(:)

    PP => L%ArrayItem(1)
    select type(P=> PP)
    type is (double precision)
        if (all(P==1)) then
            print *, 'OK list double'
        else
            fails = fails+1
            print *, 'Error double value'
        end if
        class default
        fails = fails+1
        print *,'Error getting double real'
    end select
    PP => L%ArrayItem(2)
    select type(P=>PP)
    type is (real)
        if (all(P==2)) then
            print *, 'OK list real'
        else
            fails = fails+1
            print *, 'Error real value'
        end if
        class default
        fails = fails+1
        print *,'Error getting list real'
    end select
    PP => L%ArrayItem(3)
    select type(P=>PP)
    type is (integer)
        if (all(P==3)) then
            print *, 'OK list int'
        else
            fails = fails+1
            print *, 'Error int value'
        end if
        class default
        fails = fails+1
        print *,'Error getting int real'
    end select
    PP=>L%ArrayItem(4)
    select type(P=>PP)
    type is (logical)
        if (all(P .eqv. .true.)) then
            print *, 'OK list logical'
        else
            fails = fails+1
            print *, 'Error logical value'
        end if
        class default
        fails = fails+1
        print *,'Error getting list logical'
    end select
    select type(P=>L%ArrayItemIndex(2,2))
    type is (real)
        if (P/=2.0) then
            fails = fails +1
            print *, 'Error with ArrayItemIndex'
        end if
        class default
        fails = fails+1
        print *, 'Error getting ArrayItemIndex'
    end select
    end subroutine checkL

    end  subroutine test_TObjectList

    subroutine test_TRealLIst(fails)
    Type(TRealList) T
    integer i, fails

    print *, 'test_TRealLIst'
    call T%Add(0.5d0)
    call T%Add(-3.d0)
    call T%Add(12.d0)
    do i=1, 40
        call T%Add(i*1.d0)
    end do
    call T%Add(-5.d0)

    call T%Sort()
    if (T%Item(1) == -5.d0 .and. T%item(3) ==0.5d0) then
        print *,'test_TRealLIst OK'
    else
        fails= fails + 1
        print *,'error'
        print *, T%AsArray()
    end if

    end subroutine

    subroutine test_TRealArrayList(fails)
    Type(TRealArrayList) T
    integer fails

    print *, 'test_TRealArrayList'
    call T%Add([0.5d0])
    call T%Add([-3.d0, 1.d0])
    if (all(T%Item(1) == [0.5d0]) .and. T%item(2,2) ==1.d0) then
        print *,'test_TRealArrayList OK'
    else
        fails= fails+1
        print *,'error'
    end if

    end subroutine test_TRealArrayList


    subroutine ReadWrite(T)
    class(TObjectList) T
    integer unit

    open (newunit=unit,file=temp_file,status='replace', form='unformatted')
    call T%SaveBinary(unit)
    close(unit)
    open (newunit=unit,file=temp_file,status='old', form='unformatted')
    call T%Clear()
    call T%ReadBinary(unit)
    close(unit, status = 'DELETE')

    end subroutine

    subroutine test_TStringLIst(fails)
    Type(TStringList) :: T
    integer i, fails

    print *, 'test_TStringLIst'
    call T%Add('here')
    call T%Add('there')
    call T%Add('alpha')
    call T%sort()

    do i=1,2
        if (T%Item(1) /= 'alpha' .or. T%Item(3) /='there') then
            print *, 'Error'
            fails = fails + 1
            call T%WriteItems()
        else
            print *, 'test_TStringLIst OK'
        end if
        if (i==2) exit
        call ReadWrite(T)
    end do
    call T%Clear()
    call T%Add('here','there')
    call T%Add('hot','alpha')
    if (T%ValueOf('here')/='there') then
        print *,'value Of error'
        fails = fails + 1
    else
        print *,'Value of OK'
    end if

    end subroutine test_TStringLIst

    end module ListTests
