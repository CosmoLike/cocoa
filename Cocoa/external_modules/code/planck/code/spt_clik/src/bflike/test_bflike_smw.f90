program test_bflike

  use healpix_types
  use bflike_smw

  implicit none

  integer(i4b) ,parameter :: myTT=1 ,myEE=2 ,myBB=3 ,myTE=4 ,myTB=5 ,myEB=6
  integer(i4b) ,parameter :: maxl = 64 ,npt = 1000

  integer(i4b)            :: j ,l 
  real(dp) ,parameter     :: reflike = -5639.55369105528_dp
  real(dp)                :: like ,cls(2:maxl,6) 
  real(sp)                :: t2,t1

!reads parameter file for bflike
  call init_pix_like_smw('./')
  print *,"kidding"
  cls = 0.d0

!fiducial cls
  open(10,file='./camb_tau0.06_r0.00_Aprior.dat',status='old',action='read')
  do j=2,maxl
     read(10,*) l,cls(l,myTT),cls(l,myEE),cls(l,myBB),cls(l,myTE)
  end do
  close(10)
  
  call cpu_time(t1)
  !input to get_pix_loglike are l(l+1)C_l/2pi in a (2:lmax,1:6) array with
  ! TT EE BB TE TB EB
  print *,"here we go"
  call get_pix_loglike_smw(cls,like)
  call cpu_time(t2)

  write(*,*) 'you have like = ',like
  write(*,*) 'it should be  = ',reflike
  write(*,*) 'relative difference = ',(like-reflike)/reflike

  call clean_pix_like_smw()

end program test_bflike
