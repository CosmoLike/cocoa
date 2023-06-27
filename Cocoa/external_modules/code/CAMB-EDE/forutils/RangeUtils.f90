    !A collection of ranges, consisting of sections of minimum step size
    !Useful for getting piecewise uniform point samplings for numerical integration
    !(with possible mix of log and linear spacing)
    !Antony Lewis, http://cosmologist.info/

    module RangeUtils
    use MiscUtils
    use MpiUtils
    use ArrayUtils
    implicit none
    private

    type, public :: TRange
        integer :: start_index
        integer :: steps
        logical :: IsLog
        double precision :: Low, High
        double precision :: delta
        double precision :: delta_max, delta_min !for log spacing, the non-log max and min step size
    contains

    end type TRange

    type, public :: TRanges
        integer :: count = 0
        integer :: npoints = 0
        double precision :: Lowest, Highest
        type(TRange), allocatable :: R(:)
        logical :: has_dpoints = .false.
        double precision, dimension(:), allocatable :: points, dpoints
        !dpoints is (points(i+1)-points(i-1))/2
        double precision :: RangeTol = 0.1d0
        !fraction of bin width we are prepared for merged bin widths to increase by
        logical, private :: changed = .true.
    contains
    procedure :: Init => TRanges_Free
    procedure :: Free => TRanges_Free
    procedure :: IndexOf => TRanges_IndexOf
    procedure :: Array => TRanges_Array
    procedure :: dArray => TRanges_dArray
    procedure :: GetArray => TRanges_GetArray
    procedure :: Getdpoints => TRanges_Getdpoints
    procedure :: Add_delta => TRanges_Add_delta
    procedure :: Add => TRanges_Add
    procedure :: Write => TRanges_Write
    end type TRanges

    contains

#ifdef __GFORTRAN__
    impure elemental subroutine TRanges_Free(this)
#else
    subroutine TRanges_Free(this)
#endif
    class(TRanges), intent(inout) :: this

    if (allocated(this%R)) deallocate(this%R)
    if (allocated(this%points)) deallocate(this%points)
    if (allocated(this%dpoints)) deallocate(this%dpoints)
    this%count = 0
    this%npoints = 0
    this%has_dpoints = .false.
    this%changed = .true.
    end subroutine TRanges_Free


    function TRanges_IndexOf(this, tau) result(pointstep)
    class(TRanges), intent(in) :: this
    double precision, intent(in) :: tau
    integer :: pointstep, i

    pointstep=0
    do i=1, this%count
        associate(AReg => this%R(i))
            if (tau < AReg%High .and. tau >= AReg%Low) then
                if (AReg%IsLog) then
                    pointstep = AReg%start_index + int(log(tau / AReg%Low) / AReg%delta)
                else
                    pointstep = AReg%start_index + int((tau - AReg%Low) / AReg%delta)
                end if
                return
            end if
        end associate
    end do

    if (tau >= this%Highest) then
        pointstep = this%npoints
    else
        print *, "tau=", tau, ",this%Highest=", this%Highest
        call MpiStop('TRanges_IndexOf: value out of range')
    end if
    end function TRanges_IndexOf

    function TRanges_Array(this)
    class(TRanges), intent(inout), target :: this
    double precision, pointer :: TRanges_Array(:)

    if (this%changed) call this%GetArray(this%has_dpoints)
    TRanges_Array => this%points

    end function TRanges_Array

    function TRanges_dArray(this)
    class(TRanges), intent(inout), target :: this
    double precision, pointer :: TRanges_dArray(:)

    if (this%changed) then
        call this%GetArray(.true.)
    else if (.not. this%has_dpoints) then
        call this%Getdpoints()
        this%has_dpoints = .true.
    end if
    TRanges_dArray => this%dpoints

    end function TRanges_dArray

    subroutine TRanges_GetArray(this, want_dpoints)
    class(TRanges), intent(inout) :: this
    logical, intent(in), optional :: want_dpoints
    integer :: i,j,ix

    this%has_dpoints = DefaultTrue(want_dpoints)

    call reallocate(this%points, this%npoints)

    ix=0
    do i=1, this%count
        associate (AReg => this%R(i))
            do j = 0, AReg%steps-1
                ix=ix+1
                if (AReg%IsLog) then
                    this%points(ix) = AReg%Low*exp(j*AReg%delta)
                else
                    this%points(ix) = AReg%Low + AReg%delta*j
                end if
            end do
        end associate
    end do
    ix =ix+1
    this%points(ix) = this%Highest
    if (ix /= this%npoints) call MpiStop('TRanges_GetArray: ERROR')
    this%changed= .false.
    if (this%has_dpoints) call this%Getdpoints()
    end subroutine TRanges_GetArray


    subroutine TRanges_Getdpoints(this, half_ends)
    class(TRanges), intent(inout) :: this
    logical, intent(in), optional :: half_ends
    integer :: i
    logical :: halfs

    halfs = DefaultTrue(half_ends)

    call reallocate(this%dpoints, this%npoints)
    do i=2, this%npoints-1
        this%dpoints(i) = (this%points(i+1) - this%points(i-1))/2
    end do
    if (halfs) then
        this%dpoints(1) = (this%points(2) - this%points(1))/2
        this%dpoints(this%npoints) = (this%points(this%npoints) - this%points(this%npoints-1))/2
    else
        this%dpoints(1) = (this%points(2) - this%points(1))
        this%dpoints(this%npoints) = (this%points(this%npoints) - this%points(this%npoints-1))
    end if
    end subroutine TRanges_Getdpoints


    subroutine TRanges_Add_delta(this, t_start, t_end, t_approx_delta, IsLog)
    class(TRanges), intent(inout) :: this
    logical, intent(in), optional :: IsLog
    double precision, intent(in) :: t_start, t_end, t_approx_delta
    integer :: n
    logical :: WantLog

    WantLog = DefaultFalse(IsLog)

    if (t_end <= t_start) &
        call MpiStop('TRanges_Add_delta: end must be larger than start')
    if (t_approx_delta <=0) call MpiStop('TRanges_Add_delta: delta must be > 0')

    if (WantLog) then
        n  = max(1,int(log(t_end/t_start)/t_approx_delta + 1.d0 - this%RangeTol))
    else
        n  = max(1,int((t_end-t_start)/t_approx_delta + 1.d0 - this%RangeTol))
    end if
    call this%Add(t_start, t_end, n, WantLog)
    end subroutine TRanges_Add_delta


    subroutine TRanges_Add(this, t_start, t_end, nstep, IsLog)
    class(TRanges), intent(inout) :: this
    logical, intent(in), optional :: IsLog
    double precision, intent(in) :: t_start, t_end
    integer, intent(in) :: nstep
    type(TRange), allocatable, target :: NewRanges(:)
    double precision, allocatable :: EndPoints(:), RequestDelta(:)
    integer :: ixin, nreg, ix, i,j, nsteps
    double precision :: min_request, max_request, min_log_step, max_log_step
    double precision :: delta, diff, max_delta, splitpt
    logical :: WantLog
    type(TRange) :: NewR, TmpR

    WantLog = DefaultFalse(IsLog)

    if (WantLog) then
        delta = log(t_end / t_start) / nstep
    else
        delta = (t_end - t_start) / nstep
    end if

    if (t_end <= t_start) call MpiStop('TRanges_Add: end must be larger than start')
    if (nstep <= 0) call MpiStop('TRanges_Add: nstep must be > 0')

    nreg = 0
    allocate(NewRanges(this%count+1))
    NewR%Low   = t_start
    NewR%High  = t_end
    NewR%delta = delta
    NewR%steps = nstep
    NewR%IsLog = WantLog
    call Setmins(NewR)

    !Add existing ranges adjusting existing bounds if overlapping non-trivially
    do i=1, this%count
        associate(R => this%R(i))
            if (NewR%Low >= R%Low .and. NewR%High <= R%High .and. &
                (NewR%Delta_min >= R%Delta_max .or. (NewR%IsLog .eqv. R%IsLog)  .and. NewR%delta >= R%delta)) then
                !New range wider than existing
                return
            else if (NewR%Low <= R%Low .and. NewR%High >= R%High .and. &
                (NewR%Delta_min <= R%Delta_min .or. (NewR%IsLog .eqv. R%IsLog)  .and. NewR%delta <= R%delta)) then
                !R is completely inside and wider than new one
                cycle
            else if (NewR%Delta_min < R%Delta_max .and. NewR%High > R%Low .and. NewR%Low < R%High) then
                if (R%isLog .and. .not. NewR%isLog) then
                    !Find point where log and linear spacing become equally fine
                    nsteps = int(log(NewR%delta/(exp(R%delta)-1)/R%Low)/R%delta+1)
                    splitpt = R%Low*exp(R%delta*nsteps)
                    if (splitpt*(1-exp(-R%delta)) > NewR%delta*(1+this%RangeTol)) then
                        splitpt = NewR%delta/(1-exp(-R%delta))
                    end if

                    if (splitpt < R%High) then
                        if (R%High > NewR%High) then
                            TmpR%Low = NewR%High
                            TmpR%High = R%High
                            TmpR%IsLog = R%IsLog
                            call SetDelta(TmpR, R%Delta)
                            call AddRange(TmpR)
                        end if
                        R%High = max(NewR%Low, splitpt)
                        call SetDelta(R, R%delta)
                    end if
                else if (.not. R%IsLog .and. NewR%IsLog) then
                    nsteps = int(log(R%delta/(exp(NewR%delta)-1)/NewR%Low)/NewR%delta+1)
                    splitpt = NewR%Low*exp(NewR%delta*nsteps)
                    if (splitpt*(1-exp(-NewR%delta)) > R%delta*(1+this%RangeTol)) then
                        splitpt = R%delta/(1-exp(-NewR%delta))
                    end if
                    if (splitpt < R%High .and. splitpt > R%Low) then
                        if (R%Low < NewR%Low) then
                            TmpR%Low = R%Low
                            TmpR%High = NewR%Low
                            TmpR%IsLog = R%IsLog
                            call SetDelta(TmpR, R%delta)
                            call AddRange(TmpR)
                        end if
                        R%Low = min(NewR%High, splitpt)
                        call SetDelta(R, R%delta)
                    end if
                end if
            end if
            call AddRange(R)
        end associate
    end do
    this%changed = .true.
    if (allocated(this%R)) deallocate(this%R)
    if (nreg>0) allocate(this%R, source = NewRanges(1:nreg))
    this%count = nreg
    call AddRange(NewR)

    allocate(EndPoints(0:nreg * 2))

    !Get end points in order
    ix = 0
    do i=1, nreg
        associate (AReg => NewRanges(i))
            if (ix==0) then
                EndPoints(1) = AReg%Low
                EndPoints(2) = AReg%High
                ix = 2
            else
                ixin = ix
                do j=1,ixin
                    if (AReg%Low < EndPoints(j)) then
                        EndPoints(j+1:ix+1) = EndPoints(j:ix)
                        EndPoints(j) = AReg%Low
                        ix=ix+1
                        exit
                    end if
                end do
                if (ixin == ix) then
                    ix = ix+1
                    EndPoints(ix) = AReg%Low
                    ix = ix+1
                    EndPoints(ix) = AReg%High
                else
                    ixin = ix
                    do j=1,ixin
                        if (AReg%High < EndPoints(j)) then
                            EndPoints(j+1:ix+1) = EndPoints(j:ix)
                            EndPoints(j) = AReg%High
                            ix=ix+1
                            exit
                        end if
                    end do
                    if (ixin == ix) then
                        ix = ix+1
                        EndPoints(ix) = AReg%High
                    end if
                end if
            end if
        end associate
    end do

    !remove duplicate points
    ixin = ix
    ix = 1
    do i=2, ixin
        if (EndPoints(i) /= EndPoints(ix)) then
            ix=ix+1
            EndPoints(ix) = EndPoints(i)
        end if
    end do

    !ix is the number of end points
    this%Lowest = EndPoints(1)
    this%Highest = EndPoints(ix)
    this%count = 0

    max_delta = this%Highest - this%Lowest

    if (.not. allocated(this%R) .or. size(this%R) < ix) then
        if (allocated (this%R)) deallocate (this%R)
        allocate (this%R(ix))
    end if

    allocate(RequestDelta(ix))

    do i=1, ix - 1
        associate (AReg => this%R(i))
            AReg%Low = EndPoints(i)
            AReg%High = EndPoints(i+1)

            delta = max_delta  ! max_delta = EndPoints(i+1) - EndPoints(i)
            AReg%IsLog = .false.
            !Loop over requested ranges to find step size "delta" required to satisfy requirements for this AReg
            do j=1, nreg
                if (AReg%Low >= NewRanges(j)%Low .and. AReg%Low < NewRanges(j)%High) then
                    if (NewRanges(j)%IsLog) then
                        if (AReg%IsLog) then
                            delta = min(delta,NewRanges(j)%delta)
                        else
                            min_log_step = AReg%Low*(exp(NewRanges(j)%delta)-1)
                            if (min_log_step < delta) then
                                max_log_step = AReg%High*(1-exp(-NewRanges(j)%delta))
                                AReg%IsLog = .true.
                                if (max_log_step <= delta) then
                                    delta = NewRanges(j)%delta
                                else
                                    delta = -log(1-delta/AReg%High)
                                end if
                            end if
                        end if
                    else !New Range is not log
                        if (AReg%IsLog) then
                            max_log_step = AReg%High*(1-exp(-delta))
                            if (NewRanges(j)%delta < max_log_step) then
                                min_log_step = AReg%Low*(exp(delta)-1)
                                if (min_log_step <  NewRanges(j)%delta) then
                                    !This is obviously not as good as splitting up the range
                                    delta = -log(1-NewRanges(j)%delta/AReg%High)
                                else
                                    AReg%IsLog = .false.
                                    delta = NewRanges(j)%delta
                                end if
                            end if
                        else
                            delta = min(delta, NewRanges(j)%delta)
                        end if
                    end if
                end if
            end do
            call SetDelta(AReg, delta)
            this%count = this%count + 1
            RequestDelta(this%count) = delta
        end associate
    end do


    !Get rid of tiny TRanges
    ix = this%count
    do i=ix, 1, -1
        associate (AReg => this%R(i))
            if (AReg%steps ==1) then
                Diff = AReg%High - AReg%Low
                if (AReg%IsLog) then
                    min_request = AReg%Low*(exp(RequestDelta(i))-1)
                    max_request = AReg%High*(1-exp(-RequestDelta(i)))
                else
                    min_request = RequestDelta(i)
                    max_request = min_request
                end if
                if (i/= this%count) then
                    associate (LastReg => this%R(i+1))
                        if (RequestDelta(i) >= AReg%delta .and. Diff <= LastReg%Delta_min &
                            .and. LastReg%Delta_min <= max_request .and. &
                            (.not. LastReg%IsLog .or. AReg%Low > LastReg%Low*exp(-LastReg%delta))) then

                            LastReg%Low = AReg%Low
                            if (Diff > LastReg%Delta_min*this%RangeTol) then
                                LastReg%steps =  LastReg%steps + 1
                            end if
                            if (LastReg%IsLog) then
                                LastReg%delta = log(LastReg%High/LastReg%Low) / LastReg%steps
                            else
                                LastReg%delta = (LastReg%High -LastReg%Low) / LastReg%steps
                            end if
                            this%R(i:this%Count-1) = this%R(i+1:this%Count)
                            this%Count = this%Count -1
                            cycle
                        end if
                    end associate
                end if
                if (i/=1) then
                    associate (LastReg => this%R(i-1))
                        if (RequestDelta(i) >= AReg%delta .and. Diff <= LastReg%Delta_max &
                            .and. LastReg%Delta_max <= min_request) then
                            LastReg%High = AReg%High
                            if (Diff > LastReg%Delta_max*this%RangeTol) then
                                LastReg%steps =  LastReg%steps + 1
                            end if
                            if (LastReg%IsLog) then
                                LastReg%delta = log(LastReg%High/LastReg%Low) / LastReg%steps
                            else
                                LastReg%delta = (LastReg%High -LastReg%Low) / LastReg%steps
                            end if
                            this%R(i:this%Count-1) = this%R(i+1:this%Count)
                            this%Count = this%Count -1
                        end if
                    end associate
                end if
            end if
        end associate
    end do


    !Set up start indices and get total number of steps
    nsteps = 1
    do i = 1, this%Count
        associate (AReg => this%R(i))
            AReg%Start_index = nsteps
            nsteps = nsteps + AReg%steps
            if (AReg%IsLog) then
                if (AReg%steps ==1) then
                    AReg%Delta_min = AReg%High - AReg%Low
                    AReg%Delta_max = AReg%Delta_min
                else
                    AReg%Delta_min = AReg%Low*(exp(AReg%delta)-1)
                    AReg%Delta_max = AReg%High*(1-exp(-AReg%delta))
                end if
            else
                AReg%Delta_max = AReg%delta
                AReg%Delta_min = AReg%delta
            end if
        end associate
    end do

    this%npoints = nsteps

    contains

    subroutine Setmins(AReg)
    type(TRange), intent(inout) :: Areg

    if (AReg%IsLog) then
        if (AReg%steps ==1) then
            AReg%Delta_min = AReg%High - AReg%Low
            AReg%Delta_max = AReg%Delta_min
        else
            AReg%Delta_min = AReg%Low*(exp(AReg%delta)-1)
            AReg%Delta_max = AReg%High*(1-exp(-AReg%delta))
        end if
    else
        AReg%Delta_max = AReg%delta
        AReg%Delta_min = AReg%delta
    end if
    end subroutine

    subroutine SetDelta(AReg, delta)
    Type(TRange) :: AReg
    double precision, intent(in) :: delta
    double precision Diff

    if (AReg%IsLog) then
        Diff = log(AReg%High/AReg%Low)
    else
        Diff = AReg%High - AReg%Low
    endif
    if (delta >= Diff) then
        AReg%delta = Diff
        AReg%steps = 1
    else
        AReg%steps  = max(1,int(Diff/delta + 1.d0 - this%RangeTol))
        AReg%delta = Diff / AReg%steps
    end if
    call Setmins(AReg)

    end subroutine SetDelta


    subroutine AddRange(R)
    Type(TRange), intent(in) :: R
    Type(TRange), allocatable :: Tmp(:)

    nreg = nreg +1
    if (nreg> size(NewRanges)) then
        allocate(tmp(nreg+3))
        tmp(:nreg-1) = NewRanges(:nreg-1)
        call move_alloc(tmp, NewRanges)
    end if
    NewRanges(nreg) = R
    end subroutine

    end subroutine TRanges_Add


    subroutine TRanges_Write(this)
    class(TRanges), intent(in), target :: this
    integer :: i

    do i=1,this%count
        associate (AReg => this%R(i))
            if (AReg%IsLog) then
                Write (*,'("Range ",I3,":", 3E14.4," log")') i, AReg%Low, AReg%High, AReg%delta
            else
                Write (*,'("Range ",I3,":", 3E14.4," linear")') i, AReg%Low, AReg%High, AReg%delta
            end if
        end associate
    end do
    end subroutine TRanges_Write


    end module RangeUtils
