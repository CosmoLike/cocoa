    module IniTests
    use IniObjects
    use StringUtils
    use FileUtils
    use FileTests
    implicit none
    character(LEN=0), target :: Empty_String = ''

    contains

    function RunIniTests() result (fails)
    integer fails
    fails = 0

    call test_readini(fails)
    call test_Read_env(fails)

    end function RunIniTests

    subroutine test_readini(fails)
    integer fails
    Type(TIniFile) :: Ini
    character (LEN=200), dimension(2) :: Lines
    real x
    double precision y
    Type(TTextFile) T

    call T%CreateFile(temp_file)
    Lines(1) = 'parameter = ../tester'
    Lines(2) = 'x = 3.04'
    call T%Write(Lines(1))
    call T%Write(Lines(2))
    call T%Close()
    print *, 'Wrote ini OK'
    call Ini%Open(temp_file)
    print *, 'Read ini OK'
    call File%Delete(temp_file)
    if (Ini%Read_String('parameter')/='../tester') then
        print *, 'error reading '//trim(Lines(1))
        fails = fails + 1
    else
        print *, 'Ini OK read'
    end if
    if (Ini%Read_String('parameter', NotFoundFail=.true.) /='../tester') then
        print *, 'error reading '//trim(Lines(1))
        fails = fails + 1
    end if
    call Ini%Read('x', x)
    call Ini%Read('x', y)
    if (x/=3.04 .or. y/=3.04d0) then
        print *, 'Ini error reading real'
        fails = fails + 1
    else
        print *, 'Ini OK read real'
    end if
    call Ini%Close()
    call Ini%Open_FromLines(Lines,2)
    if (Ini%Read_String('parameter')/='../tester') then
        print *, 'error reading lines '//trim(Lines(1))
        fails = fails + 1
    else
        print *, 'Ini OK read lines'
    end if

    call Ini%close()

    end subroutine test_readini

    subroutine test_read_env(fails)
    integer fails
    Type(TIniFile) :: Ini
    character(LEN=:), allocatable :: S
    character (LEN=200), dimension(1) :: Lines

    Lines(1) = 'parameter = test$(PATH)/mypath$(PATH)'

    call Ini%Open_FromLines(Lines,1)
    S = Ini%Read_String('parameter')
    if (S /= 'test'//GetEnvironmentVariable('PATH')//'/mypath'//GetEnvironmentVariable('PATH')) then
        fails = fails + 1
        print *, 'error'
    else
        print *, 'OK Ini path'
    end if

    end subroutine test_read_env

    end module IniTests
