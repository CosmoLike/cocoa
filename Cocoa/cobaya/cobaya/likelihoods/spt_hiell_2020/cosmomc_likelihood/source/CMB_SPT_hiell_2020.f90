!Likelihood code used in 
  !SPTpol+SZ l=2000-11000 power spectrum
  !For questions, please contact Christian Reichardt
  module CMB_SPT_hiell_2020
    use Likelihood_Cosmology
    use CosmoTheory
    use CosmologyTypes
    use FileUtils
    use MatrixUtils
    use foregrounds

    implicit none

    ! The bandpowers, in order 90x90, 90x150, 90x220, 150x150, 150x220, 220x220
    integer :: nfreq, nband, nall,maxnbin
    integer,dimension(:),allocatable :: nbins,offsets
    ! cov == The cholesky factored bandpower covariance  
    double precision, dimension(:,:),allocatable :: cov, windows, beam_err, cov_w_beam,cov_tmp
    integer :: spt_windows_lmin,spt_windows_lmax
    double precision, dimension(:), allocatable :: cl_to_dl_conversion,spec,ells
    double precision, dimension(:,:), allocatable :: spt_eff_fr
    double precision, dimension(:), allocatable :: spt_prefactor
    double precision, dimension(5) ::  spt_norm_fr
    double precision, dimension(3,3) :: cal_cov,tmp
    double precision :: CalibLnL0
    integer, dimension(:,:),allocatable :: indices
    
    logical :: SuccessfulSPTInitialization
    logical :: normalizeSZ_143GHz
    logical :: CallFGPrior
    logical :: ApplyFTSPrior
    logical :: use_dZ

    logical :: SuccessfulSPTHiellInitialization
    logical :: printDlSPT, printDlSPTComponents
    logical :: plankish
    
    logical :: binaryCov, binaryWindows, binaryBeamErr
    
    !real(mcp) :: meanBeam

    Type, extends(TCMBLikelihood) :: TSPTHiEllLike
  contains
    procedure :: ReadIni => SPT_HiEll_ReadIni
    procedure :: InitSPTHiEllData
    procedure :: LogLike => SPTHiEllLnLike
 end Type TSPTHiEllLike

contains

  subroutine setSPTUninitialized()
    SuccessfulSPTHiellInitialization=.false.
  end subroutine setSPTUninitialized


 subroutine SPT_HiEll_ReadIni(this, Ini)
   use IniFile
   use IniObjects
   implicit none
   class(TSPTHiEllLike) :: this
   class(TSettingIni) :: Ini
   character (LEN=Ini_max_string_len) :: desc_file
   character (LEN=Ini_max_string_len) :: bp_file, param_file
   character (LEN=Ini_max_string_len) :: cov_file,beamerr_file
   character (LEN=Ini_max_string_len) :: window_folder


   
   call InitFGModel(Ini)
   

   param_file = Ini%Read_String_Default( 'spt_hiell_params_file','')
   call this%loadParamNames(param_file)

   

   desc_file = Ini%Read_String_Default('spt_hiell_description','')
   bp_file = Ini%Read_String_Default('spt_hiell_bandpowers','')
   
   
   !default to binary since faster i/o
   binaryCov=.True.
   binaryWindows=.True.
   binaryBeamErr=.True.
   cov_file =  Ini%Read_String_Default('spt_hiell_binary_covariance','')
   if (cov_file == '') then 
      cov_file = Ini%Read_String_Default('spt_hiell_covariance','')
      binaryCov=.False.
   endif
   window_folder = Ini%Read_String_Default('spt_hiell_binary_window','')
   if (window_folder == '') then 
      window_folder = Ini%Read_String_Default('spt_hiell_windir','')
      binaryWindows=.False.
   endif
   beamerr_file = Ini%Read_String_Default('spt_hiell_binary_beamerr','')
   if (beamerr_file == '') then 
      beamerr_file = Ini%Read_String_Default('spt_hiell_beamerr','')
      binaryBeamErr=.False.
   endif
      
   normalizeSZ_143GHz = .true. !Ini%Read_Logical('normalizeSZ_143ghz',.false.)
   !do we want extra debug prints
   CallFGPrior   = Ini%Read_Logical('apply_fg_prior',.true.)
   ApplyFTSPrior = Ini%Read_Logical('apply_prior_FTS_bands',.true.)
   
   printDlSPT = Ini%Read_Logical('print_spectrum',.false.)
   printDlSPTComponents = Ini%Read_Logical('print_spectrum_components',.false.)
   use_dZ = Ini%Read_Logical('interpret_ksz2_as_dz',.false.)
   if ((printDlSPT .or. printDlSPTComponents) .and. MPIRank /= 0) then
      call MPIStop('Warning - print_spectrum/print_spectrum_components is not MPI thread-safe, quitting...')
   endif
   
   if (bp_file=='' .or. desc_file=='' .or. window_folder=='' .or. cov_file=='' .or. beamerr_file == '') then
      print*,'Missing required spt hiell key: received: ',bp_file,desc_file,window_folder,cov_file,beamerr_file
      stop
   endif
   
   call this%InitSPTHiEllData(desc_file, bp_file, cov_file, beamerr_file, window_folder)
 end subroutine SPT_HiEll_ReadIni
 
 subroutine InitSPTHiEllData(this, desc_file, bp_file, cov_file, beamerr_file, window_folder)
   use IniFile
   implicit none
   class(TSPTHiEllLike) :: this
   character(LEN=Ini_max_string_len) :: desc_file, bp_file, cov_file,beamerr_file
   character(LEN=Ini_max_string_len) :: window_folder
   integer, parameter :: tmp_file_unit=82
   integer i,j,k,l,dum,lwin
   integer*4 :: neff
   double precision, allocatable, dimension(:) :: locwin
   integer*8 :: offset,delta
   integer*4 :: efflmin,efflmax,j0,j1
   real*4 :: arr(2)
   real*8 rtmp
   Type(TTextFile) :: F
   integer*4 :: errcode
   logical wexist
   double precision, dimension(3):: zeros

   zeros=0
   cal_cov(1,1) = 1.1105131e-05
   cal_cov(1,2) = 3.5551351e-06
   cal_cov(2,1) = cal_cov(1,2)
   cal_cov(1,3) = 1.1602891e-06
   cal_cov(3,1) = cal_cov(1,3)
   cal_cov(2,2) = 3.4153547e-06
   cal_cov(2,3) = 2.1348018e-06
   cal_cov(3,2) = cal_cov(2,3)
   cal_cov(3,3) = 1.7536000e-05
   tmp=cal_cov
   CalibLnL0 = Matrix_GaussianLogLikeDouble(tmp,zeros)
   !Obtain necessary info from the desc_file pertaining
   !to which freqs we're using, ell range, and number of bins per spectrum.
   inquire(FILE=trim(desc_file),EXIST=wexist)
   if (.not. wexist) then
      print*,'SPT hiell 2020, missing desc file:', trim(desc_file)
      call mpistop()
   endif
   
   call F%Open(desc_file)
   read(F%unit,*) nall,nfreq !number of bandpowers, number of freqs
   if (nfreq > MaxNFreq) &
        call MpiStop('spt initialized with more than allowed Nfrequencies')
   nband = (nfreq)*(nfreq+1)/2
   allocate(nbins(nband))
   allocate(spt_eff_fr(5,nfreq),spt_prefactor(nfreq),indices(2,nband),offsets(nband))
   do i=1,nband
      read (F%unit,*) j
      nbins(i)=j
   end do
   if (nall .ne. sum(nbins)) call MpiStop('mismatched number of bandpowers')
   maxnbin=maxval(nbins(:))
   do i=1,5
      read (F%unit,*) rtmp
      spt_norm_fr(i) = rtmp
   enddo
   if (normalizeSZ_143GHz) then 
      spt_norm_fr(5) = 143.
      print*,'using 143 as tSZ center freq'
   endif
   read (F%unit,*) spt_windows_lmin, spt_windows_lmax !Min and Max ell in window file
   
   do j=1,nfreq
       do i=1,5
          read (F%unit,*) rtmp
          spt_eff_fr(i,j) = rtmp
       enddo
    enddo
    do j=1,nfreq
       read (F%unit,*) spt_prefactor(j)
    enddo
   call F%Close()

   if (feedback > 1) then 
      print *, 'nall: ', nall
      print *, 'nfreq: ', nfreq
      print *, 'spt_windows_lmin: ', spt_windows_lmin
      print *, 'spt_windows_lmax: ', spt_windows_lmax
      print *, 'window_folder: ', trim(window_folder)
   endif

   if (spt_windows_lmax .gt. ReportFGLmax()) then
      print*,'Hard-wired lmax in foregrounds.f90 is too low for CMB_SPT_hiell_2020.f90'
      call mpistop()
   endif


   allocate(this%cl_lmax(CL_T,CL_T), source=0)
   this%cl_lmax(CL_T,CL_T) = spt_windows_lmax

   if (spt_windows_lmin < 2 .or. spt_windows_lmin >= spt_windows_lmax) then
      call mpistop('Invalid lranges for sptpol')
   end if

   !ells vector is 2 ell longer in order to do cl derivatives.
   !As a result, so is cl_to_dl_conversion
   allocate( ells(spt_windows_lmin:spt_windows_lmax), &
             cl_to_dl_conversion(spt_windows_lmin:spt_windows_lmax) )

   allocate(windows(spt_windows_lmin:spt_windows_lmax,nall), &
        spec(nall))

   allocate(cov(nall,nall), beam_err(nall,nall),cov_w_beam(nall,nall),cov_tmp(nall,nall))

   !Define an array with the l*(l+1)/2pi factor to convert to Dl from Cl.
   do j=spt_windows_lmin,spt_windows_lmax
      ells(j) = j
   enddo
   cl_to_dl_conversion(:) = (ells*(ells+1d0))/TWOPI

   ! Read in bandpowers
   !Should be 90x90, 90x`150, 90x220,150x150, 150x220, 220x220 in that order
   inquire(FILE=trim(bp_file),EXIST=wexist)
   if (.not. wexist) then
      print*,'SPT hiell 2020, missing bp file:', trim(bp_file)
      call mpistop()
   endif
   call F%Open(bp_file)
   do i=1,nall
      read (F%unit,*) dum,spec(i)
   end do
   call F%close()
   
   
   
   ! Read in covariance
   inquire(FILE=trim(cov_file),EXIST=wexist)
   if (.not. wexist) then
      print*,'SPT hiell 2020, missing cov file:', trim(cov_file)
      call mpistop()
   endif
   if (binaryCov) then 
      call OpenReadBinaryFile(cov_file,tmp_file_unit,nall*8_8)
      do i=1,nall
         read(tmp_file_unit,rec=i)cov(:,i)
      enddo
      close (tmp_file_unit)
   else
      call F%open(cov_file)
      do i=1,nall
         do j=1,nall
            read (F%unit,*) cov(j,i)
         end do
      end do
      call F%close()
   endif
   if (feedback > 1) print *, 'First entry of covariance matrix: ', cov(1,1)
   inquire(FILE=trim(beamerr_file),EXIST=wexist)
   if (.not. wexist) then
      print*,'SPT hiell 2020, missing beamerr file:', trim(beamerr_file)
      call mpistop()
   endif
   if (binaryBeamErr) then 
      call OpenReadBinaryFile(beamerr_file,tmp_file_unit,nall*8_8)
      do i=1,nall
         read(tmp_file_unit,rec=i)beam_err(:,i)
      enddo
      close (tmp_file_unit)
   else
      call F%open(beamerr_file)
      do i=1,nall
         do j=1,nall
            read (F%unit,*) beam_err(j,i)
         end do
      end do
      call F%close()
   endif
   if (feedback > 1) print *, 'First entry of beam correlation matrix: ', beam_err(1,1)
   
   
   
   ! Read in windows
   if (binaryWindows) then
      inquire(FILE=trim(window_folder),EXIST=wexist)
      if (.not. wexist) then
         print*,'SPT hiell 2020, missing window file:', trim(window_folder)
         call mpistop()
      endif
      call OpenReadBinaryStreamFile(trim(window_folder),tmp_file_unit)
      read(tmp_file_unit,pos=1)efflmin,efflmax
      allocate(locwin(efflmin:efflmax))
      if ((efflmax .lt. spt_windows_lmin) .or. (efflmin .gt. spt_windows_lmax)) &
           call MpiStop('unallowed l-ranges for binary window functions')
      j0=efflmin
      if (spt_windows_lmin > j0) j0=spt_windows_lmin
      j1=efflmax
      if (spt_windows_lmax < j1) j1=spt_windows_lmax
      if (j1 < j0) &
           call MpiStop('unallowed l-ranges for binary window functions - no allowed ells')
      delta=(efflmax-efflmin+1)*8_8
      offset=2 * 4_8+1
      do i=1,nall
         read(tmp_file_unit,pos=((i-1)*delta + offset)) locwin
         windows(j0:j1,i)=locwin(j0:j1)
      end do
      close(tmp_file_unit)
      deallocate(locwin)
   else
      do i=1,nall
         inquire(FILE=trim(window_folder)//trim(numcat('window_',i)),EXIST=wexist)
         if (.not. wexist) then
            print*,'SPTpol, missing window file:', trim(window_folder)//trim(numcat('window_',i))
            call mpistop()
         endif
         call F%Open(trim(window_folder)//trim(numcat('window_',i)))
         do j=spt_windows_lmin,spt_windows_lmax
            read (F%unit,*) dum, windows(j,i)
         end do
         call F%Close()
      end do
   end if

   i=0
   do j=1,nfreq
      do k=j,nfreq
         i=i+1
         indices(1,i)=j
         indices(2,i)=k
      end do
   end do
   offsets(1)=1
   do i=2,nband
      offsets(i)=offsets(i-1)+nbins(i-1)
   enddo
   

   SuccessfulSPTHiellInitialization = .true.

   if (feedback > 1) then
      print *, 'Successfully initialized SPT_HIELL data...'
   endif

 end subroutine InitSPTHiEllData



 function SPTHiEllLnLike(this, CMB, Theory, DataParams) 
   use MpiUtils
   implicit none
   
   class(TSPTHiEllLike) :: this
   Class(CMBParams) :: CMB
   Class(TCosmoTheoryPredictions), target :: Theory
   double precision :: DataParams(:) 
   double precision, dimension(spt_windows_lmax) :: dl_cmb
   double precision, dimension(spt_windows_lmax,7) :: component_spectra
   double precision :: PriorLnLike
   double precision :: dum
   double precision :: SPTHiEllLnLike,CalibLnLike, FTSLnLike, FTSfactor
   double precision :: FGPriorLnLike, NoCalLnLike
   double precision, parameter :: d3000 = 3000*3001/TWOPI
   double precision, parameter :: beta = 0.0012309
   double precision, parameter :: dipole_cosine = -0.4033
   double precision, dimension(1:nall) :: deltacb,cbs
   double precision, dimension(1:maxnbin) :: tmpcb
   double precision, dimension(1) :: junk, detcov
   double precision, dimension(2) :: PoissonLevels 
   double precision, dimension(2) :: ADust
   double precision, dimension(2) :: alphaDust
   double precision, dimension(3) ::  CalFactors,delta_calib !90, 150, 220
   real*4, dimension(10) ::  comp_arr
   type(foreground_params) :: foregroundParams
   integer :: i,j,k, l,kk, thisoffset,thisnbin
   double precision :: norm
   integer fid
   real*4, dimension(2) :: arr
   real*4, dimension(7) :: arr7
   integer*4 :: errcode
   integer, dimension(1)::ivec
   double precision :: kszfac, tszfac

   double precision, dimension(spt_windows_lmin:spt_windows_lmax) :: dl_fgs
   integer, parameter :: iFG = 5

   ! get CMB spectrum
   call Theory%ClArray(dl_cmb(:),CL_T,CL_T)      

   if (HaveForegroundsBeenInitialized() .eq. .false.) then
      write(*,*)'trying to call SPT likelihood w/o initializing foregrounds'
      call mpistop()
   endif
   if (iFG .gt. 100) call mpistop() !haven't done it yet
   CalFactors = DataParams(1:3)
   FTSFactor  = DataParams(4)
   foregroundParams = GetForegroundParamsFromArray(DataParams(iFG:iFG+nForegroundParams))

   tszfac = cosmo_scale_tsz(CMB%H0,Theory%sigma_8,CMB%omb)
   foregroundParams.czero_tsz = foregroundParams.czero_tsz * tszfac
   kszfac = cosmo_scale_ksz(CMB%H0,Theory%sigma_8,CMB%omb,CMB%omc+CMB%omb+CMB%omnu,CMB%InitPower(ns_index),CMB%tau)
   foregroundParams.czero_ksz = foregroundParams.czero_ksz * kszfac

   if (use_dZ) then
      kszfac = pkSZ(foregroundParams.czero_ksz2)
      foregroundParams.czero_ksz2 = kszfac
   endif


   !add this to use negative for correlation

!   if (printDlSPT) then
!      call printForegrounds(foregroundParams)
!   endif


   !$OMP PARALLEL DO  DEFAULT(NONE), &
   !$OMP  SHARED(cbs,indices,foregroundParams,spt_eff_fr,spt_norm_fr,cl_to_dl_conversion,nbins,nfreq,dl_cmb,spt_windows_lmax,spt_windows_lmin,windows,spt_prefactor,CalFactors,deltacb,offsets,spec,nband,printDlSPT,printDlSPTComponents, FTSfactor), &
   !$OMP  private(i,j,k,dl_fgs,tmpcb,thisnbin,l,thisoffset,fid,arr,component_spectra,comp_arr), &
   !$OMP SCHEDULE(STATIC)
   do i=1,nband
      j=indices(1,i)
      k=indices(2,i)
      thisoffset=offsets(i)
      thisnbin=nbins(i)
      tmpcb(:)=0

      !first get theory spectra
      if (printDlSPTComponents) then
         dl_fgs(spt_windows_lmin:spt_windows_lmax) = dl_foreground(foregroundParams,j,k,nfreq,spt_eff_fr+FTSfactor, &
              spt_norm_fr,spt_windows_lmin,spt_windows_lmax,component_spectra) 
      else
         dl_fgs(spt_windows_lmin:spt_windows_lmax) = dl_foreground(foregroundParams,j,k,nfreq,spt_eff_fr+FTSfactor, &
              spt_norm_fr,spt_windows_lmin,spt_windows_lmax) 
      endif
      !add CMB
      dl_fgs(spt_windows_lmin:spt_windows_lmax)=dl_fgs(spt_windows_lmin:spt_windows_lmax)+dl_cmb(spt_windows_lmin:spt_windows_lmax)
      


      if (printDlSPT) then
         fid=33+GetMpiRank()+i
         print*,'printing Dl spt'
         call OpenWriteBinaryFile(trim(numcat('suxp_spt_',i)),fid,4_8 * 2)
         do l=spt_windows_lmin,spt_windows_lmax
            arr(1)=l
            arr(2)=dl_fgs(l)
            write(fid,rec=l-spt_windows_lmin+1) arr(1:2)
         enddo
         close(fid)
      endif
      if (printDlSPTComponents) then
         fid=33+GetMpiRank()+i
         call OpenWriteBinaryFile(trim(numcat('suxp_spt_components_',i)),fid,4_8 * 10)
         print*,spt_windows_lmin,spt_windows_lmax
         do l=spt_windows_lmin,spt_windows_lmax
            comp_arr(1)=l
            comp_arr(2)=dl_fgs(l)
            comp_arr(3) = dl_cmb(l)
            comp_arr(4) = component_spectra(l,1) 
            comp_arr(5) = component_spectra(l,2) 
            comp_arr(6) = component_spectra(l,3) 
            comp_arr(7) = component_spectra(l,4) 
            comp_arr(8) = component_spectra(l,5) 
            comp_arr(9) = component_spectra(l,6) 
            comp_arr(10) = component_spectra(l,7)
            write(fid,rec=l-spt_windows_lmin+1) comp_arr(1:10)
            if ((l/1000)*1000 .eq. l) then 
               print*,'Test: ',comp_arr(:)
            endif
         enddo
         close(fid)
      endif


      !now bin with window functions
      call dgemv('T',spt_windows_lmax-spt_windows_lmin+1,thisnbin,1.0d0,&
           windows(:,thisoffset:thisoffset+thisnbin-1),spt_windows_lmax-spt_windows_lmin+1,&
           dl_fgs,1,0.0d0,tmpcb,1)

      !apply prefactors
      tmpcb = tmpcb * spt_prefactor(k)*spt_prefactor(j)*CalFactors(j)*CalFactors(k)
      if (printDlSPT) then
         open(fid,file=trim(numcat('est_bandpowers_spt_',i)))
         write(fid,*)'# ',j,k
         do l=1,thisnbin
            write(fid,*)l,tmpcb(l),spec(thisoffset+l-1)
         enddo
         close(fid)
      endif

      cbs(thisoffset:thisoffset+thisnbin-1) = tmpcb(1:thisnbin) 
   enddo

   deltacb = cbs - spec

   do i=1,nall
      do j=1,nall
         cov_w_beam(i,j) = beam_err(i,j)*cbs(i)*cbs(j)
      enddo
   enddo

   cov_w_beam = cov + cov_w_beam
   if (feedback > 1)    cov_tmp=cov_w_beam
   SPTHiEllLnLike =  Matrix_GaussianLogLikeDouble(cov_w_beam, deltacb)
   NoCalLnLike=SPTHiEllLnLike   
   if (CallFGPrior) then
      FGPriorLnLike = getForegroundPriorLnL(foregroundParams)
      SPTHiEllLnLike = SPTHiEllLnLike + FGPriorLnLike
   endif

   delta_calib = log(CalFactors)
   !can take off cov term since cov is fixed, and constant for all points
   tmp=cal_cov
   CalibLnLike = Matrix_GaussianLogLikeDouble(tmp,delta_calib) - CalibLnL0
   

   SPTHiEllLnLike = SPTHiEllLnLike + CalibLnLike

   FTSLnLike = 0
   if (ApplyFTSPrior) &
        FTSLnLike = 0.5* (FTSfactor / 0.3)**2
   SPTHiEllLnLike = SPTHiEllLnLike + FTSLnLike
   

   if (feedback > 1)  then
      print *, 'SPTHiEllLnLike lnlike = ', SPTHiEllLnLike
      print*, 'Calibration chisq (mistakenly double subtracting)',2*(CalibLnLike-CalibLnL0)
   print*, 'Calibration chisq',2*(CalibLnLike)
      detcov = Matrix_GaussianLogLikeDouble(cov_tmp, deltacb*0)
      print*,'lnLcov term',detcov
      print*,'chisq for cov only:',   2*(NoCalLnLike-detcov)
      print*,'chisq for fts bit:',2*FTSLnLike
      print*,'chisq for FG prior:',2*FGPriorLnLike
      print *, 'SPTHiEllLike chisq (after priors) = ', 2*(SPTHiEllLnLike-detcov)
   endif
 end function SPTHiEllLnLike

end module CMB_SPT_hiell_2020
