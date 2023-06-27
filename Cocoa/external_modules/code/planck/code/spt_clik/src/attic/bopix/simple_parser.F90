MODULE simple_parser

  IMPLICIT NONE

  PRIVATE

  PUBLIC &
       &   parse      &
       & , parse_file

  INTERFACE parse
     MODULE PROCEDURE &
          &           parse_single  &
          &         , parse_double  &
          &         , parse_i4      &
          &         , parse_i8      &
          &         , parse_logical &
          &         , parse_string
  END INTERFACE

  CHARACTER(len=256) :: line,name,value
  INTEGER(4)         :: eqpos,ios
  LOGICAL            :: lfound

CONTAINS

!!$-------------------------------------------------------------------

  SUBROUTINE parse_file(filename,default,found)

    CHARACTER(len=*),            INTENT(out) :: filename
    CHARACTER(len=*),  OPTIONAL, INTENT(in)  :: default
    LOGICAL,           OPTIONAL, INTENT(out) :: found


    INTEGER :: iarg

    INTERFACE
       FUNCTION iargc()
         INTEGER iargc
       END FUNCTION iargc

       SUBROUTINE getarg (num, res)
         INTEGER, INTENT(in) :: num
         CHARACTER(len=*), INTENT(out) :: res
       END SUBROUTINE getarg
    END INTERFACE


    lfound = .FALSE.

    !DO iarg=1,iargc()
    !   CALL getarg(iarg,line)
!
!    !   eqpos = SCAN(line,"=")
!    !   name  = TRIM(ADJUSTL(line(:eqpos-1)))
!    !   value = TRIM(ADJUSTL(line(eqpos+1:)))
!
!    !   IF (TRIM(name) .EQ. "--parfile") THEN
!    !      IF (LEN_TRIM(value) .NE. 0) THEN
!    !         filename = TRIM(ADJUSTL(value))
!    !         lfound   = .TRUE.
!    !         EXIT
!    !      END IF
!    !   END IF
!
!    !END DO
!
    IF (.NOT. lfound) THEN

       IF (PRESENT(default)) THEN
          filename=TRIM(ADJUSTL(default))
       ELSE
          WRITE(unit=*,fmt="(a)",advance="no")&
               &"Input the name of the def. file: "
          READ(unit=*,fmt="(a)")filename
       END IF

    END IF

    INQUIRE(file=TRIM(ADJUSTL(filename)), exist=lfound)

    IF (PRESENT(found)) found = lfound


!!$    IF (.NOT. found ) THEN
!!$       WRITE(*,*) TRIM(filename)//" does not exist"
!!$       STOP
!!$    END IF

  END SUBROUTINE parse_file

!!$-------------------------------------------------------------------

  SUBROUTINE parse_single(filename,string,variable,default,found)

    CHARACTER(len=*),            INTENT(in)  :: filename
    CHARACTER(len=*),            INTENT(in)  :: string
    REAL(4),                     INTENT(out) :: variable
    REAL(4),           OPTIONAL, INTENT(in)  :: default
    LOGICAL,           OPTIONAL, INTENT(out) :: found

    lfound   = .FALSE.
    variable = 0

    OPEN (unit=111,file=trim(filename),status="OLD",action="READ")

    DO
       READ (unit=111,fmt="(a)",iostat=ios) line

       IF (ios .NE. 0) EXIT

       eqpos=SCAN(line,"=")

       IF ( (eqpos .EQ. 0) .OR. (line(1:1) .EQ. "#") ) CYCLE

       name  = TRIM(ADJUSTL(line(:eqpos-1)))
       value = TRIM(ADJUSTL(line(eqpos+1:)))

       IF ( TRIM(ADJUSTL(name)) .EQ. TRIM(ADJUSTL(string)) ) THEN
          READ(unit=value,fmt=*) variable
          lfound= .TRUE.
          EXIT
       END IF

    END DO

    CLOSE(unit=111)

    IF ( (.NOT. lfound) .AND. (PRESENT(default)) )  THEN
       variable = default
       lfound   = .TRUE.
    END IF

    IF (PRESENT(found)) found=lfound

  END SUBROUTINE parse_single

!!$-------------------------------------------------------------------

  SUBROUTINE parse_double(filename,string,variable,default,found)

    CHARACTER(len=*),            INTENT(in)  :: filename
    CHARACTER(len=*),            INTENT(in)  :: string
    REAL(8),                     INTENT(out) :: variable
    REAL(8),           OPTIONAL, INTENT(in)  :: default
    LOGICAL,           OPTIONAL, INTENT(out) :: found

    lfound   = .FALSE.
    variable = 0

    OPEN (unit=111,file=trim(filename),status="OLD",action="READ")

    DO
       READ (unit=111,fmt="(a)",iostat=ios) line

       IF(ios .NE. 0) EXIT

       eqpos=SCAN(line,"=")

       IF ( (eqpos .EQ. 0) .OR. (line(1:1) .EQ. "#") ) CYCLE

       name  = TRIM(ADJUSTL(line(:eqpos-1)))
       value = TRIM(ADJUSTL(line(eqpos+1:)))

       IF ( TRIM(ADJUSTL(name)) .EQ. TRIM(ADJUSTL(string)) ) THEN
          READ(unit=value,fmt=*) variable
          lfound= .TRUE.
          EXIT
       END IF

    END DO

    CLOSE(unit=111)

    IF ( (.NOT. lfound) .AND. (PRESENT(default)) )  THEN
       variable = default
       lfound   = .TRUE.
    END IF

    IF (PRESENT(found)) found=lfound

  END SUBROUTINE parse_double

!!$-------------------------------------------------------------------

  SUBROUTINE parse_i4(filename,string,variable,default,found)

    CHARACTER(len=*),            INTENT(in)  :: filename
    CHARACTER(len=*),            INTENT(in)  :: string
    INTEGER(4),           INTENT(out) :: variable
    INTEGER(4), OPTIONAL, INTENT(in)  :: default
    LOGICAL,           OPTIONAL, INTENT(out) :: found

    lfound   = .FALSE.
    variable = 0

    OPEN (unit=111,file=trim(filename),status="OLD",action="READ")

    DO
       READ (unit=111,fmt="(a)",iostat=ios) line

       IF(ios .NE. 0) EXIT

       eqpos=SCAN(line,"=")

       IF ( (eqpos .EQ. 0) .OR. (line(1:1) .EQ. "#") ) CYCLE

       name  = TRIM(ADJUSTL(line(:eqpos-1)))
       value = TRIM(ADJUSTL(line(eqpos+1:)))

       IF ( TRIM(ADJUSTL(name)) .EQ. TRIM(ADJUSTL(string)) ) THEN
          READ(unit=value,fmt=*) variable
          lfound= .TRUE.
          EXIT
       END IF

    END DO

    CLOSE(unit=111)

    IF ( (.NOT. lfound) .AND. (PRESENT(default)) )  THEN
       variable = default
       lfound   = .TRUE.
    END IF

    IF (PRESENT(found)) found=lfound

  END SUBROUTINE parse_i4

!!$-------------------------------------------------------------------

  SUBROUTINE parse_i8(filename,string,variable,default,found)

    CHARACTER(len=*),            INTENT(in)  :: filename
    CHARACTER(len=*),            INTENT(in)  :: string
    INTEGER(8),           INTENT(out) :: variable
    INTEGER(8), OPTIONAL, INTENT(in)  :: default
    LOGICAL,           OPTIONAL, INTENT(out) :: found

    lfound   = .FALSE.
    variable = 0

    OPEN (unit=111,file=trim(filename),status="OLD",action="READ")

    DO
       READ (unit=111,fmt="(a)",iostat=ios) line

       IF(ios .NE. 0) EXIT

       eqpos=SCAN(line,"=")

       IF ( (eqpos .EQ. 0) .OR. (line(1:1) .EQ. "#") ) CYCLE

       name  = TRIM(ADJUSTL(line(:eqpos-1)))
       value = TRIM(ADJUSTL(line(eqpos+1:)))

       IF ( TRIM(ADJUSTL(name)) .EQ. TRIM(ADJUSTL(string)) ) THEN
          READ(unit=value,fmt=*) variable
          lfound= .TRUE.
          EXIT
       END IF

    END DO

    CLOSE(unit=111)

    IF ( (.NOT. lfound) .AND. (PRESENT(default)) )  THEN
       variable = default
       lfound   = .TRUE.
    END IF

    IF (PRESENT(found)) found=lfound

  END SUBROUTINE parse_i8

!!$-------------------------------------------------------------------

  SUBROUTINE parse_logical(filename,string,variable,default,found)

    CHARACTER(len=*),            INTENT(in)  :: filename
    CHARACTER(len=*),            INTENT(in)  :: string
    LOGICAL,                     INTENT(out) :: variable
    LOGICAL,           OPTIONAL, INTENT(in)  :: default
    LOGICAL,           OPTIONAL, INTENT(out) :: found

    lfound   = .FALSE.
    variable = .FALSE.

    OPEN (unit=111,file=trim(filename),status="OLD",action="READ")

    DO
       READ (unit=111,fmt="(a)",iostat=ios) line

       IF(ios .NE. 0) EXIT

       eqpos=SCAN(line,"=")

       IF ( (eqpos .EQ. 0) .OR. (line(1:1) .EQ. "#") ) CYCLE

       name  = TRIM(ADJUSTL(line(:eqpos-1)))
       value = TRIM(ADJUSTL(line(eqpos+1:)))

       IF ( TRIM(ADJUSTL(name)) .EQ. TRIM(ADJUSTL(string)) ) THEN
          READ(unit=value,fmt=*) variable
          lfound= .TRUE.
          EXIT
       END IF

    END DO

    CLOSE(unit=111)

    IF ( (.NOT. lfound) .AND. (PRESENT(default)) )  THEN
       variable = default
       lfound   = .TRUE.
    END IF

    IF (PRESENT(found)) found=lfound

  END SUBROUTINE parse_logical

!!$-------------------------------------------------------------------

  SUBROUTINE parse_string(filename,string,variable,default,found)

    CHARACTER(len=*),            INTENT(in)  :: filename
    CHARACTER(len=*),            INTENT(in)  :: string
    CHARACTER(len=*),            INTENT(out) :: variable
    CHARACTER(len=*),  OPTIONAL, INTENT(in)  :: default
    LOGICAL,           OPTIONAL, INTENT(out) :: found

    lfound   = .FALSE.
    variable = ""

    OPEN (unit=111,file=trim(filename),status="OLD",action="READ")

    DO
       READ (unit=111,fmt="(a)",iostat=ios) line

       IF(ios .NE. 0) EXIT

       eqpos=SCAN(line,"=")

       IF ( (eqpos .EQ. 0) .OR. (line(1:1) .EQ. "#") ) CYCLE

       name  = TRIM(ADJUSTL(line(:eqpos-1)))
       value = TRIM(ADJUSTL(line(eqpos+1:)))

       IF ( TRIM(ADJUSTL(name)) .EQ. TRIM(ADJUSTL(string)) ) THEN
          variable = TRIM(ADJUSTL(value))
          lfound   = .TRUE.
          EXIT
       END IF

    END DO

    CLOSE(unit=111)

    IF ( (.NOT. lfound) .AND. (PRESENT(default)) )  THEN
       variable = default
       lfound   = .TRUE.
    END IF

    IF (PRESENT(found)) found=lfound

  END SUBROUTINE parse_string

!!$-------------------------------------------------------------------

END MODULE simple_parser
