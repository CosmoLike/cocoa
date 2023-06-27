    module FileTests
    use FileUtils
    implicit none

    character(LEN=*), parameter :: temp_file='test.tmp'

    !Some things are already tested in ObjectLists_tests

    contains

    function RunFileTests() result (fails)
    integer fails
    character(LEN=7400) :: St
    Type(TTextFile) :: T
    character(LEN=:), allocatable :: InLine, comment
    integer line
    real testval

    fails = 0

    St = 'test string'
    St(len(St):len(St)) = 'c'

    call T%CreateFile(temp_file)

    call T%Write('#comment header')
    call T%Write(' ')
    call T%Write(St)
    call T%Write('name = ',3.04)
    call T%Write(' tester values')
    call T%Close()

    line = 0
    do while (T%ReadNextContentLine(temp_file,InLine))
        line = line+1
        if (line==1) then
            if (Inline/=St) then
                fails = fails+1
                print *, 'Error reading string'
            end if
        else if (line==2) then
            read(InLine(7:),*) testval
            if (testval /= 3.04) then
                fails = fails+1
                print *, 'Error reading real'
            end if
        end if
    end do

    call T%Open(temp_file)
    if (T%ReadLineSkipEmptyAndComments(InLine, Comment)) then
        if (comment/='comment header') then
            fails = fails+1
            print *, 'Error reading header'
        end if
    else
        fails= fails+1
        print *, 'Error reading line'
    end if
    call T%Close(del=.true.)
    if (fails==0) print *,'OK file tests'

    end function RunFileTests


    end module FileTests