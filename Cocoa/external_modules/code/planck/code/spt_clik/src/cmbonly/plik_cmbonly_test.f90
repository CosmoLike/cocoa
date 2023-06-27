 !===========================================================
 program test

 !===========================================================
 !E. Calabrese, J. Dunkley 2014. 
 !===========================================================

 use Plik_CMBonly 

 implicit none
 real(8), dimension(:), allocatable :: cell_tt,cell_ee,cell_te
 character(LEN=128) :: filename
 real(8)            :: plike,cp
 integer            :: lun, il, dummy,i,j

 !---------------------------------------------------

 print *,""
 print *,"Planck CMB-only likelihood test"
 print *,"==================================="

 data_dir = "/Users/benabed/Boulot/clik-hg/plik_lite18/like_cmbonly_plikv18/data/"
 use_ee = .true.
 use_te = .true.
 
 call like_init_cmbonly
 !---------------------------------------------------
 ! read in test Cls
 !---------------------------------------------------
 filename = trim(data_dir)//'base_plikHM_TT_lowTEB.dat'
 write(*,*)"Reading in Cls from: ",trim(filename)
 lun = 1000
 allocate(cell_tt(2:tt_lmax),cell_ee(2:tt_lmax),cell_te(2:tt_lmax)) 
 cell_tt(2:tt_lmax)=0.d0
 cell_ee(2:tt_lmax)=0.d0
 cell_te(2:tt_lmax)=0.d0

 open(unit=lun,file=filename,action='read',status='old')
 do il=2,tt_lmax
    read(lun,*)dummy,cell_tt(il),cell_te(il),cell_ee(il)
 enddo
 close(lun)
cp=1
 call calc_like_cmbonly(plike,cell_tt,cell_te,cell_ee,cp)

 write(*,*) '-------------------------------------'
 write(*,*) 'Planck chi2 = ', plike*2.
 write(*,*) '-------------------------------------'

 end program test
