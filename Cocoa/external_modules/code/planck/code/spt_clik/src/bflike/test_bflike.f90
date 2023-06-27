program test_bflike
  use bflike
  use healpix_types  

  implicit none

  real(dp) :: loglike
  integer :: ell,l
  real(dp) :: clin(2:64,4),cltest(2:64,6)
  character(len=256) :: clikfile
  character(len=*),parameter :: testfile='/global/scratch2/sd/gerbino/commander/commander/src/comm_process_resfiles/camb_tau0.06_r0.00_Aprior.dat'
  integer :: t1,t2,t3,t4, clock_rate, clock_max
  real(dp),parameter :: expvalue=-9356.5402665494148096

  open(unit=2010,file=trim(testfile),status='old')
  do l=2,64
      read(2010,*)ell,clin(l,:)
  enddo
  close(2010)

  cltest=0.d0
  cltest(:,1:4)=clin

  !if (nArguments() == 1) call getArgument(1, clikfile) 
  clikfile='bflike.clik'
  if (iargc() .eq. 1) call getarg(1, clikfile)

  call system_clock ( t1, clock_rate, clock_max )
  call init_pix_like(trim(clikfile))  
  call system_clock ( t2, clock_rate, clock_max )
  call get_pix_loglike(cltest,loglike)
  call system_clock ( t3, clock_rate, clock_max )
  call clean_pix_like()
  call system_clock ( t4, clock_rate, clock_max )

  write (*,*) 'Init time = ', real ( t2 - t1 ) / real ( clock_rate )
  write (*,*) 'Like time = ', real ( t3 - t2 ) / real ( clock_rate )
  write (*,*) 'Finalize time = ', real ( t4 - t3 ) / real ( clock_rate )

  write(*,'(a,f,a,f,a,f)') 'Like: ',loglike,'. Expected: ',expvalue, '. Difference: ',expvalue-loglike



end program test_bflike
