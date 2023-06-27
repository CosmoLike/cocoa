    module StringTests
    use StringUtils
    implicit none
    contains


    function RunStringTests() result (fails)
    integer fails
    character(LEN=:), allocatable :: S, S2
    fails = 0

    if (UpperCase('Low@er')/='LOW@ER' .or. LowerCase('1L oW@eR')/='1l ow@er') then
        print *, 'Error in uppercase/lowercase'
        fails = fails+1
    end if

    S = FormatString('Test %d for %f %s', 91, 3.04, 'result')
    if (S/='Test 91 for 3.0400 result') then
        fails = fails+1
        print *, 'Error in FormatString: ', S
    end if
    
    S2 = Join(',', S,'XP')
    if (S2 /= S//',XP') then
        fails = fails+1
        print *, 'Error in Join: ', S2
    end if
    
    if (fails==0) print *, 'OK string tests'
    

    end function

    end module
