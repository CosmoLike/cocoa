   SUBROUTINE READ_PARAMETER()

      USE PARAMETER_MODULE
      USE SIMPLE_PARSER
      USE HEALPIX_TYPES
      USE PIX_TOOLS
      USE ALM_TOOLS
      USE FITSTOOLS, ONLY: INPUT_MAP

     IMPLICIT NONE

    defs_file='input_file.txt'
!=================================================================================================================================
!  key="NSIDEOUT"
!  CALL parse(defs_file,key,NSIDEOUT,found=found)
!  IF (.NOT. FOUND) THEN
!       WRITE(6,*)'NSIDEOUT NOT FOUND'
!       STOP
!   END IF

!  key="BOPIX_LMAX_COV"
!  CALL parse(defs_file,key,BOPIX_LMAX_COV,found=found)
!  IF (.NOT. FOUND) THEN
!       WRITE(6,*)'BOPIX_LMAX_COV NOT FOUND'
!       STOP
!   END IF

  key="MENO_LOGLIK_FACTOR"
  CALL parse(defs_file,key, MENO_LOGLIK_FACTOR,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'MENO_LOGLIK_FACTOR NOT FOUND'
       STOP
   END IF

 key="FWHM_DEG"
  CALL parse(defs_file,key, FWHM_DEG,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'FWHM_DEG NOT FOUND'
       STOP
   END IF

!  key="BOPIX_NTHREADS"
!  CALL parse(defs_file,key,BOPIX_NTHREADS,found=found)
!  IF (.NOT. FOUND) THEN
!       WRITE(6,*)'BOPIX_NTHREADS NOT FOUND'
!       STOP
!   END IF

  key="BOPIX_DATA_DIR"
  CALL parse(defs_file,key,BOPIX_DATA_DIR,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'BOPIX_DATA_DIR NOT FOUND'
       STOP
   END IF

!================================================================================================================================
  key="NOISE"
  CALL parse(defs_file,key,NOISE,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'NOISE NOT FOUND'
       STOP
   END IF

  key="ORDERING_COV_SIGNAL"
  CALL parse(defs_file,key,ORDERING_COV_SIGNAL,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'ORDERING_COV_SIGNAL NOT FOUND'
       STOP
   END IF

  key="COV_FILE"
  CALL parse(defs_file,key,COV_FILE,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'COV_FILE NOT FOUND'
       STOP
   END IF
   COV_FILE=TRIM(BOPIX_DATA_DIR)//COV_FILE

  key="TDFACT_COVNOISE"
  CALL parse(defs_file,key,TDFACT_COVNOISE,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'TDFACT_COVNOISE NOT FOUND'
       STOP
   END IF

  key="UNIT_FACT_COV"
  CALL parse(defs_file,key,UNIT_FACT_COV,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'UNIT_FACT_COV NOT FOUND'
       STOP
   END IF

  key="TD_FACT"
  CALL parse(defs_file,key,TD_FACT,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'TD_FACT NOT FOUND'
       STOP
   END IF

  key="UNIT_FACT"
  CALL parse(defs_file,key,UNIT_FACT,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'UNIT_FACT NOT FOUND'
       STOP
   END IF

  key="WINDOW_FILE"
  CALL parse(defs_file,key,WINDOW_FILE,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'WINDOW_FILE NOT FOUND'
       STOP
   END IF
   WINDOW_FILE=TRIM(BOPIX_DATA_DIR)//WINDOW_FILE

!================================================================================================================================

  key="CONVERT_ORD_MAP"
  CALL parse(defs_file,key,CONVERT_ORD_MAP,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'CONVERT_ORD_MAP NOT FOUND'
       STOP
   END IF

  key="MAP_TYPE"
  CALL parse(defs_file,key,MAP_TYPE,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'MAP_TYPE NOT FOUND'
       STOP
   END IF

   SELECT CASE (MAP_TYPE)
       CASE (1)
           key="MAP_FILE"
           CALL parse(defs_file,key,MAP_FILE,found=found)
           IF (.NOT. FOUND) THEN
              WRITE(6,*)'MAP_FILE NOT FOUND'
              STOP
           END IF
         MAP_FILE=TRIM(BOPIX_DATA_DIR)//MAP_FILE
        CASE (3) 
          key="MASKED_MAP_T"
          CALL parse(defs_file,key,MASKED_MAP_T,found=found)
          IF (.NOT. FOUND) THEN
             WRITE(6,*)'MASKED_MAP_T NOT FOUND'
             STOP
          END IF

          key="MASKED_MAP_P"
          CALL parse(defs_file,key,MASKED_MAP_P,found=found)
          IF (.NOT. FOUND) THEN
             WRITE(6,*)'MASKED_MAP_P NOT FOUND'
             STOP
          END IF
            MASKED_MAP_T=TRIM(BOPIX_DATA_DIR)//MASKED_MAP_T
            MASKED_MAP_P=TRIM(BOPIX_DATA_DIR)//MASKED_MAP_P

    END SELECT           

!=================================================================================================================================
  key="CONVERT_ORD_MASK"
  CALL parse(defs_file,key,CONVERT_ORD_MASK,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'CONVERT_ORD_MASK NOT FOUND'
       STOP
   END IF

  key="NFILES_MASK"
  CALL parse(defs_file,key,NFILES_MASK,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'NFILES_MASK NOT FOUND'
       STOP
   END IF
  
   SELECT CASE (NFILES_MASK)
    CASE (1)
        key="MASKFILE"
        CALL parse(defs_file,key,MASKFILE,found=found)
         IF (.NOT. FOUND) THEN
           WRITE(6,*)'MASKFILE NOT FOUND'
           STOP
         END IF
        MASKFILE=TRIM(BOPIX_DATA_DIR)//MASKFILE

    CASE (2)
        key="MASKFILE_T"
        CALL parse(defs_file,key,MASKFILE_T,found=found)
          IF (.NOT. FOUND) THEN
            WRITE(6,*)'MASKFILE_T NOT FOUND'
            STOP
          END IF
        MASKFILE_T=TRIM(BOPIX_DATA_DIR)//MASKFILE_T

         key="MASKFILE_P"
         CALL parse(defs_file,key,MASKFILE_P,found=found)
          IF (.NOT. FOUND) THEN
            WRITE(6,*)'MASKFILE_P NOT FOUND'
            STOP
          END IF
        MASKFILE_P=TRIM(BOPIX_DATA_DIR)//MASKFILE_P
     END SELECT

!=================================================================================================================================
  key="BOPIX_CL_FILE"
  CALL parse(defs_file,key,BOPIX_CL_FILE,found=found)
  IF (.NOT. FOUND) THEN
       WRITE(6,*)'BOPIX_CL_FILE NOT FOUND'
       STOP
   END IF
   
   SELECT CASE (BOPIX_CL_FILE)
   
    CASE (1)
    key="BOPIX_CL_LMAX"
    CALL parse(defs_file,key,BOPIX_CL_LMAX,found=found)
    IF (.NOT. FOUND) THEN
        WRITE(6,*)'BOPIX_CL_LMAX NOT FOUND'
        STOP
     END IF

    key="BOPIX_CL_FILENAME"
    CALL parse(defs_file,key,BOPIX_CL_FILENAME,found=found)
    IF (.NOT. FOUND) THEN
       WRITE(6,*)'BOPIX_CL_FILENAME NOT FOUND'
       STOP
    END IF
    BOPIX_CL_FILENAME=TRIM(BOPIX_DATA_DIR)//BOPIX_CL_FILENAME      
   
   END SELECT

   END SUBROUTINE READ_PARAMETER

