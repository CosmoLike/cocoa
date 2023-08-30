! SPT-3G 2018 TTTEEE likelihood
! Written by Lennart Balkenhol 2022
! Main module

!---------------------------------------------------------!
!          MODULE DEFINITION, IMPORTS, CONSTANTS          !
!---------------------------------------------------------!

! Imports
module CMB_SPT3G_2018_TTTEEE
#ifdef _STANDALONE_
use SPT3G_utils
#else
use Likelihood_Cosmology
use CosmoTheory
use CosmologyTypes
use FileUtils
use MatrixUtils
use settings
#endif

use CMB_SPT3G_2018_TTTEEE_foregrounds


! Parameters that decide how the covariance is manipulated to be positive definite
real(mcp) :: cov_eval_cut_threshold, cov_eval_large_number_replacement

! Level of feedback
logical :: print_spectrum

! Only print out the chisq
logical :: print_chisq

! Overwrite CAMB's theory spectrum?
character(LEN=:), allocatable :: fixed_theory_input_file

! Data set parameters
integer :: SPT3G_windows_lmin, SPT3G_windows_lmax ! Will be read in from param file
integer, parameter :: N_freq_0 = 3 ! Number of frequency bands
integer, parameter :: N_b_0_TT = 44 ! Maximum (uncropped) number of bandpowers
integer, parameter :: N_b_0_TE = 44 ! Maximum (uncropped) number of bandpowers
integer, parameter :: N_b_0_EE = 44 ! Maximum (uncropped) number of bandpowers
integer, parameter :: N_s_0 = 18 ! Maximum (uncropped) number of spectra that can be used for parameter estimation
integer, parameter :: N_b_0_total = 6*N_b_0_TT+6*N_b_0_TE+6*N_b_0_EE ! Maximum number of bins

! Variables that help juggling the indexing
Type(TStringList) :: spectra_to_fit_list ! List specifying which spectra to use for parameter estimation
integer, allocatable :: spectra_start_ix(:), spectra_stop_ix(:) ! List of start and stop ix of the spectra in the order they are used for fitting

! Late cropping of data - only use this to test non-continuous band powers vectors (i.e. skipping a random data point in the middle)
logical :: do_late_crop
integer, allocatable :: late_crop_msk(:) ! Optional list that can be used to crop down the final band powers in the logl call to allow for non-continuous bdp vectors (0 = use, 1 = dont use)

Type(TStringList) :: full_bandpower_list ! List specifying the order of bandpowers in the bandpower file
integer, allocatable :: full_bandpower_start_ix(:), full_bandpower_stop_ix(:) ! List specifying the start and stop indices for the spectra in the band power file

Type(TStringList) :: full_window_list ! List specifying the order of bandpowers in the window files

Type(TStringList) :: full_covariance_list ! List specifying the order of spectra in the covariance matrix file
integer, allocatable :: full_covariance_start_ix(:), full_covariance_stop_ix(:) ! List specifying the start and stop indices for the spectra in the covariance file

Type(TStringList) :: fiducial_covariance_list ! List specifying the order of spectra in the fiducial covariance matrix file
integer, allocatable :: fiducial_covariance_start_ix(:), fiducial_covariance_stop_ix(:) ! List specifying the start and stop indices for the spectra in the fiducial covariance file

Type(TStringList) :: full_beam_covariance_list ! List specifying the order of spectra in the beam covariance matrix file
integer, allocatable :: full_beam_covariance_start_ix(:), full_beam_covariance_stop_ix(:) ! List specifying the start and stop indices for the spectra in the beam covariance file

integer, allocatable :: calibration_ix(:,:) ! List specifying the indices of calibration paramters for the spectra requested in the fit
integer, allocatable :: spectra_to_fit_Poisson_ix(:) ! List specifying the Poisson power indices of the spectra requested in the fit
integer, allocatable :: spectra_to_fit_CIB_cl_decorr_ix(:,:) ! Array holding the indices of CIB clustering decorrelation parameters

integer, dimension(6) :: cal_rows_use ! Convenience array indicating whether cal 90/150/220 is used

! Effective band centres
Type(TStringList) :: nu_eff_list ! List specifying the order of spectra in the central frequency file
real(mcp), allocatable :: nu_eff_DSFG(:,:)
real(mcp), allocatable :: nu_eff_tSZ(:,:)
real(mcp), allocatable :: nu_eff_pol_gal_dust(:,:)
real(mcp), allocatable :: nu_eff_gal_cirrus(:,:)

! Indices of nuisance parameters sampled by CosmoMC
integer, parameter :: iKappa = 1 ! Super sample lensing
integer, parameter :: iTcal90=2, iTcal150=3, iTcal220=4, iEcal90=5, iEcal150=6, iEcal220=7 ! Calibration parameters
integer, parameter :: iEE_Poisson_90x90=8, iEE_Poisson_90x150=9, iEE_Poisson_90x220=10, iEE_Poisson_150x150=11, iEE_Poisson_150x220=12, iEE_Poisson_220x220=13 ! EE Poisson power
integer, parameter :: iEE_GalDust_Amp=14, iEE_GalDust_Alpha=15, iEE_GalDust_Beta=16 ! EE Polarised Galactic Dust
integer, parameter :: iTE_GalDust_Amp=17, iTE_GalDust_Alpha=18, iTE_GalDust_Beta=19 ! TE Polarised Galactic Dust
integer, parameter :: iTT_Poisson_90x90=20, iTT_Poisson_90x150=21, iTT_Poisson_90x220=22, iTT_Poisson_150x150=23, iTT_Poisson_150x220=24, iTT_Poisson_220x220=25 ! TT Poisson power
integer, parameter :: iTT_GalDust_Amp=26, iTT_GalDust_Alpha=27, iTT_GalDust_Beta=28 ! TT Galactic Dust (Cirrus)
integer, parameter :: iTT_CIB_Clustering_Amp=29, iTT_CIB_Clustering_Alpha=30, iTT_CIB_Clustering_Beta=31 ! TT CIB clustering
integer, parameter :: iTT_CIB_Clustering_decorr_90=32, iTT_CIB_Clustering_decorr_150=33, iTT_CIB_Clustering_decorr_220=34 ! TT CIB clustering decorrealtion
integer, parameter :: iTT_tSZ_Amp=35 ! TT tSZ Amplitude
integer, parameter :: iTT_tSZ_CIB_corr=36 ! TT tSZ-CIB correlation
integer, parameter :: iTT_kSZ_Amp=37 ! TT kSZ Amplitude
! Indices of cosmo paramters
integer, parameter :: id_H0       = 1
integer, parameter :: id_omb      = 2
integer, parameter :: id_sigma_8  = 3
integer, parameter :: id_omm      = 4
integer, parameter :: id_ns       = 5
integer, parameter :: id_tau      = 6

! Indices specifying order of fg components in output array
integer, parameter :: N_fg_max = 9
integer, parameter :: iOutCal = 1, iOutAberration = 2, iOutSSL = 3, iOutGalDust = 4, iOutPoisson = 5
integer, parameter :: iOutCIBClustering = 6, iOuttSZ = 7, iOuttSZCIB = 8, iOutkSZ = 9 ! TT exclusive foregrounds

!---------------------------------------------------------!
!                      MAIN FUNCTIONS                     !
!---------------------------------------------------------!

! Module definition
#ifdef _STANDALONE_
Type :: TSPT3G_2018_TTTEEE_Likelihood
#else
Type, extends(TCMBLikelihood) :: TSPT3G_2018_TTTEEE_Likelihood
#endif
  integer :: N_b_total ! Total number of bins to use
  integer, allocatable :: spec_bin_min(:), spec_bin_max(:) ! Bin selection
  integer :: N_freq ! Number of frequencies that get used
  integer :: N_s, N_s_TT, N_s_TE, N_s_EE ! Number of spectra in total and broken up by TT, TE, EE
  real(mcp) :: aberration_coefficient
  real(mcp), allocatable :: bandpowers(:)
  real(mcp), allocatable :: bdp_covariance(:,:)
  real(mcp), allocatable :: bdp_covariance_pos_def(:,:) ! Band power covariance after ensuring it has no negative eigenvalues
  real(mcp), allocatable :: beam_covariance(:,:)
  real(mcp), allocatable :: full_windows(:,:,:)
  real(mcp), allocatable :: inv_cal_covariance(:,:)
  real(mcp), allocatable :: theory_spectrum_overwrite(:,:)
  real(mcp) :: beam_cov_scale ! Scales the beam covariance artificially
  integer :: include_logdet
contains
#ifndef _STANDALONE_
  procedure :: ReadIni => SPT3G_2018_TTTEEE_ReadIni
  procedure :: LogLike => SPT3G_2018_TTTEEE_LogLike
  procedure :: SPT3G_2018_TTTEEE_InitData
  procedure :: MakeCovariancePositiveDefinite
#endif
  procedure :: GenerateStartStopIx
  procedure :: GenerateStartStopIxForFit
  procedure :: BinPowerSpectrum
  procedure :: AddBeamCovariance
  procedure :: AddEEForegrounds
  procedure :: AddTEForegrounds
  procedure :: AddTTForegrounds
  procedure :: GrabTheorySpectrumOverwrite
  procedure :: CheckTTBinSelection
  procedure :: SPT3G_2018_TTTEEE_LogLike_external
end Type TSPT3G_2018_TTTEEE_Likelihood

!#ifdef _STANDALONE_
!class(TSPT3G_2018_TTTEEE_Likelihood) :: single_Lkl
!#endif

contains


#ifndef _STANDALONE_
! Read ini (gets called at the start by CosmoMC)
! Read in the .ini file here and grab all the names of files to use (band powers, covariance, calibration covariance, etc.)
! Hand off to data initialisation function at the end
subroutine SPT3G_2018_TTTEEE_ReadIni(this, Ini)
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  class(TSettingIni) :: Ini

  ! Define variables
  character(LEN=:), allocatable :: parameter_filename

  character(LEN=:), allocatable :: bandpower_filename
  character(LEN=:), allocatable :: full_bandpower_list_string

  character(LEN=:), allocatable :: covariance_filename
  character(LEN=:), allocatable :: full_covariance_list_string

  character(LEN=:), allocatable :: fiducial_covariance_filename
  character(LEN=:), allocatable :: fiducial_covariance_list_string

  character(LEN=:), allocatable :: beam_covariance_filename
  character(LEN=:), allocatable :: full_beam_covariance_list_string

  character(LEN=:), allocatable :: cal_covariance_filename

  character(LEN=:), allocatable :: window_folder
  character(LEN=:), allocatable :: full_window_list_string

  character(LEN=:), allocatable :: nu_eff_filename
  character(LEN=:), allocatable :: nu_eff_list_string

  character(LEN=:), allocatable :: spectra_to_fit_list_string

  Type(TStringList) :: spec_bin_min_list
  character(LEN=:), allocatable :: spec_bin_min_list_string
  Type(TStringList) :: spec_bin_max_list
  character(LEN=:), allocatable :: spec_bin_max_list_string

  character(LEN=:), allocatable :: late_crop_msk_string
  Type(TStringList) :: late_crop_msk_list

  integer i, stat

  print *, "SPT-3G 2018 TTTEEE: Requested Likelihood. Beginning initialisation."

  !KARIM: Fortran untils here that read from the ini file

  SPT3G_windows_lmin = Ini%Read_Int("SPT3G_2018_TTTEEE_window_l_min")
  SPT3G_windows_lmax = Ini%Read_Int("SPT3G_2018_TTTEEE_window_l_max")

  ! Specify to what ell_max we need the CAMB spectra
  ! Increased from what is needed for the window functions to match CosmoPower training set
  ! cosmomc specific, don't need
  allocate(this%cl_lmax(CL_E,CL_E), source=0)
  this%cl_lmax(CL_T,CL_T) = 5000!SPT3G_windows_lmax+1
  this%cl_lmax(CL_E,CL_T) = 5000!SPT3G_windows_lmax+1
  this%cl_lmax(CL_E,CL_E) = 5000!SPT3G_windows_lmax+1
  !---

  ! Read in spectrum selection
  spectra_to_fit_list_string = Ini%Read_String("SPT3G_2018_TTTEEE_spectra_to_fit", .true.)
  call spectra_to_fit_list%SetFromString(spectra_to_fit_list_string)
  this%N_s = spectra_to_fit_list%Count

  ! Read in bin selection for each spectrum and the total number of bins we are fitting
  allocate(this%spec_bin_min(this%N_s))
  this%spec_bin_min = -1
  allocate(this%spec_bin_max(this%N_s))
  this%spec_bin_max = 45
  spec_bin_min_list_string = Ini%Read_String("SPT3G_2018_TTTEEE_spectra_to_fit_bin_min", .true.)
  call spec_bin_min_list%SetFromString(spec_bin_min_list_string)
  spec_bin_max_list_string = Ini%Read_String("SPT3G_2018_TTTEEE_spectra_to_fit_bin_max", .true.)
  call spec_bin_max_list%SetFromString(spec_bin_max_list_string)

  this%N_b_total = 0
  do i=1, this%N_s
    call str2int(spec_bin_min_list%Item(i), this%spec_bin_min(i), stat)
    call str2int(spec_bin_max_list%Item(i), this%spec_bin_max(i), stat)
    if (((this%spec_bin_min(i) .le. 0) .or. (this%spec_bin_max(i) .ge. 45)) .or. (this%spec_bin_max(i) .lt. this%spec_bin_min(i))) then
      call MPIStop("SPT-3G 2018 TTTEEE: bad ell range selection for spectrum: "//spectra_to_fit_list%Item(i))
    end if
    this%N_b_total = this%N_b_total + this%spec_bin_max(i) - this%spec_bin_min(i) + 1
  end do

  ! Check if a late crop is requested and read in the mask if necessary
  do_late_crop = Ini%Read_Logical("SPT3G_2018_TTTEEE_late_crop", .false.)
  if (do_late_crop .eqv. .true.) then
    late_crop_msk_string = Ini%Read_String("SPT3G_2018_TTTEEE_late_crop_msk", .true.)
    call late_crop_msk_list%SetFromString(late_crop_msk_string)
    allocate(late_crop_msk(this%N_b_total))
    do i=1, this%N_b_total
      call str2int(late_crop_msk_list%Item(i), late_crop_msk(i), stat)
    end do
  end if

  ! Determine how many spectra are TT vs TE vs EE and the total number of bins we are fitting
  call DetermineNumberOfSpectraByField(spectra_to_fit_list, this%N_s_TT, this%N_s_TE, this%N_s_EE)

  ! Determine how many different frequencies get used
  call DetermineNumberofFrequencies(spectra_to_fit_list, this%N_freq)


  ! cosmomc specific, don't need
  ! Read in nuisance parameter file
  parameter_filename = Ini%Read_String_Default("SPT3G_2018_TTTEEE_params_file","")
  call this%loadParamNames(parameter_filename)
  !---

  ! Read in band power file name
  bandpower_filename = Ini%ReadFileName("SPT3G_2018_TTTEEE_bandpower_file",relative = .true., NotFoundFail=.true.)
  full_bandpower_list_string = Ini%Read_String("SPT3G_2018_TTTEEE_bandpower_file_order", .true.)
  call full_bandpower_list%SetFromString(full_bandpower_list_string)

  ! Read in band power covariance matrix file name
  covariance_filename = Ini%ReadFileName("SPT3G_2018_TTTEEE_covariance_matrix",relative = .true., NotFoundFail=.true.)
  full_covariance_list_string = Ini%Read_String("SPT3G_2018_TTTEEE_covariance_matrix_order", .true.)
  call full_covariance_list%SetFromString(full_covariance_list_string)

  ! Read in fiducial band power covariance matrix file name
  fiducial_covariance_filename = Ini%ReadFileName("SPT3G_2018_TTTEEE_fiducial_covariance_matrix",relative = .true., NotFoundFail=.false.)
  fiducial_covariance_list_string = Ini%Read_String("SPT3G_2018_TTTEEE_fiducial_covariance_matrix_order", .true.)
  call fiducial_covariance_list%SetFromString(fiducial_covariance_list_string)

  cov_eval_cut_threshold = Ini%Read_Real("SPT3G_2018_TTTEEE_cov_eval_cut_threshold", 0.2)
  cov_eval_large_number_replacement = Ini%Read_Real("SPT3G_2018_TTTEEE_cov_eval_replacement", 1e3)

  ! Read in beam covariance matrix file name and scale factor
  beam_covariance_filename = Ini%ReadFileName("SPT3G_2018_TTTEEE_beam_covariance_matrix",relative = .true., NotFoundFail=.true.)
  full_beam_covariance_list_string = Ini%Read_String("SPT3G_2018_TTTEEE_beam_covariance_matrix_order", .true.)
  call full_beam_covariance_list%SetFromString(full_beam_covariance_list_string)

  this%beam_cov_scale = Ini%Read_Real("SPT3G_2018_TTTEEE_beam_covariance_scale", 1.0)

  ! Read in calibration covariance file name
  cal_covariance_filename = Ini%ReadFileName("SPT3G_2018_TTTEEE_cal_covariance_matrix",relative = .true., NotFoundFail=.true.)

  ! Read in winfow functions file name
  window_folder = Ini%Read_String("SPT3G_2018_TTTEEE_window_folder", NotFoundFail=.true.)
  full_window_list_string = Ini%Read_String("SPT3G_2018_TTTEEE_window_folder_order", .true.)
  call full_window_list%SetFromString(full_window_list_string)

  ! Read in central frequencies file name
  nu_eff_filename = Ini%ReadFileName("SPT3G_2018_TTTEEE_central_frequency_file", relative = .true., NotFoundFail=.true.)
  nu_eff_list_string = Ini%Read_String("SPT3G_2018_TTTEEE_central_frequency_file_order", .true.)
  call nu_eff_list%SetFromString(nu_eff_list_string)

  ! Read in aberration correction coefficient
  this%aberration_coefficient = Ini%Read_Real("SPT3G_2018_TTTEEE_aberration_coefficient", 0.0)

  ! Check if only the chisq is requested
  print_chisq = Ini%Read_Logical("SPT3G_2018_TTTEEE_print_chisq", .false.)

  ! I don't need that for clik
  ! Check if we want to overwrite the theory spectrum
  fixed_theory_input_file = Ini%Read_String_Default("SPT3G_2018_TTTEEE_theory_spectrum", "USE_CAMB")
  if (fixed_theory_input_file /= "USE_CAMB") then
    allocate(this%theory_spectrum_overwrite(4,1+SPT3G_windows_lmax-SPT3G_windows_lmin))
    call File%ReadTextMatrix(fixed_theory_input_file, this%theory_spectrum_overwrite) ! Read in the theory spectrum
    print *, "SPT-3G 2018 TTTEEE: using theory spectrum at ", fixed_theory_input_file
  end if

  ! Check level of feedback
  print_spectrum = Ini%Read_Logical("SPT3G_2018_TTTEEE_print_spectrum", .false.)
  if (print_spectrum .and. MPIRank /= 0) then
    call MPIStop("SPT-3G 2018 TTTEEE: print_spectrum is not MPI thread-safe!")
  endif
  ! ---
  
  ! Initialise foreground model
  call SPT3G_2018_TTTEEE_ReadIni_Foregrounds(Ini)

  ! Hand off to SPT3G_2018_TTTEEE_InitData
  call this%SPT3G_2018_TTTEEE_InitData(bandpower_filename, covariance_filename, fiducial_covariance_filename, beam_covariance_filename, cal_covariance_filename, window_folder, nu_eff_filename)

end subroutine SPT3G_2018_TTTEEE_ReadIni

! Initialise the data
! This reads in all the files and crops them down to the requested spectra and ell range
! Manipulates data in any other way that only needs to be done once
subroutine SPT3G_2018_TTTEEE_InitData(this, bandpower_filename, covariance_filename, fiducial_covariance_filename, beam_covariance_filename, cal_covariance_filename, window_folder, nu_eff_filename)
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  character(LEN=:), allocatable :: bandpower_filename
  real(mcp), allocatable :: full_bandpowers(:,:)

  character(LEN=:), allocatable :: covariance_filename
  real(mcp), allocatable :: full_covariance_matrix(:,:)
  character(LEN=:), allocatable :: fiducial_covariance_filename

  character(LEN=:), allocatable :: beam_covariance_filename
  real(mcp), allocatable :: full_beam_covariance_matrix(:,:)

  character(LEN=:), allocatable :: cal_covariance_filename
  real(mcp), allocatable :: full_cal_covariance_matrix(:,:)

  character(LEN=:), allocatable :: window_folder
  real(mcp), allocatable :: current_window(:,:)

  character(LEN=:), allocatable :: nu_eff_filename
  real(mcp), allocatable :: nu_eff_matrix(:,:)

  integer i, j, ii, jj, i_spec, j_spec, ix, spec_ok

  real(mcp), allocatable :: test_array(:)
  real(mcp), allocatable :: test_array_binned(:)

  print *, "SPT-3G 2018 TTTEEE: Initialising data. STANDALONE"

  !KARIM: A bunch of calls to the CosmomMC helper function ReadTextMatrix

  ! Grab indices of spectra in fit order
  call this%GenerateStartStopIxForFit(spectra_to_fit_list, spectra_start_ix, spectra_stop_ix)

  ! Band Powers
  allocate(full_bandpowers(1,N_b_0_total))
  call File%ReadTextMatrix(bandpower_filename, full_bandpowers, N_b_0_total,1) ! Read in full band power vector - this is actually a vector, but ReadTextMatrix is nicer than ReadText
  call this%GenerateStartStopIx(spectra_to_fit_list, full_bandpower_list, full_bandpower_start_ix, full_bandpower_stop_ix) ! Get the indices of the requested bins in the fitting order
  allocate(this%bandpowers(this%N_b_total))
  do i_spec=1, spectra_to_fit_list%Count ! Loop over full spectra and populate band power vector that will be used
    this%bandpowers(spectra_start_ix(i_spec):spectra_stop_ix(i_spec)) = full_bandpowers(1,full_bandpower_start_ix(i_spec):full_bandpower_stop_ix(i_spec))
  end do

  ! Covariance Matrix
  allocate(full_covariance_matrix(N_b_0_total,N_b_0_total))
  call File%ReadTextMatrix(covariance_filename, full_covariance_matrix) ! Read in the full covariance matrix
  call this%GenerateStartStopIx(spectra_to_fit_list, full_covariance_list, full_covariance_start_ix, full_covariance_stop_ix)  ! Get the indices of the requested bins in the fitting order
  allocate(this%bdp_covariance(this%N_b_total,this%N_b_total))
  do i_spec=1, spectra_to_fit_list%Count ! Loop over full covariance block and populate cov that will be used
    do j_spec=1, spectra_to_fit_list%Count
      this%bdp_covariance(spectra_start_ix(i_spec):spectra_stop_ix(i_spec),spectra_start_ix(j_spec):spectra_stop_ix(j_spec)) = full_covariance_matrix(full_covariance_start_ix(i_spec):full_covariance_stop_ix(i_spec),full_covariance_start_ix(j_spec):full_covariance_stop_ix(j_spec))
    end do
  end do
  allocate(this%bdp_covariance_pos_def(this%N_b_total,this%N_b_total))
  call this%MakeCovariancePositiveDefinite(this%bdp_covariance, fiducial_covariance_filename, fiducial_covariance_list, this%bdp_covariance_pos_def) ! Ensure covariance is positive definite

  ! Beam Covariance Matrix
  allocate(full_beam_covariance_matrix(N_b_0_total,N_b_0_total))
  call File%ReadTextMatrix(beam_covariance_filename, full_beam_covariance_matrix) ! Read in the full beam covariance matrix
  call this%GenerateStartStopIx(spectra_to_fit_list, full_beam_covariance_list, full_beam_covariance_start_ix, full_beam_covariance_stop_ix)  ! Get the indices of the requested bins in the fitting order
  allocate(this%beam_covariance(this%N_b_total,this%N_b_total))
  do i_spec=1, spectra_to_fit_list%Count ! Loop over full covariance block and populate cov that will be used
    do j_spec=1, spectra_to_fit_list%Count
      this%beam_covariance(spectra_start_ix(i_spec):spectra_stop_ix(i_spec),spectra_start_ix(j_spec):spectra_stop_ix(j_spec)) = full_beam_covariance_matrix(full_beam_covariance_start_ix(i_spec):full_beam_covariance_stop_ix(i_spec),full_beam_covariance_start_ix(j_spec):full_beam_covariance_stop_ix(j_spec))
    end do
  end do
  this%beam_covariance = this%beam_covariance * this%beam_cov_scale

  ! Window Functions
  ! These are a bit trickier to handle due to the independent cuts possible for TT/TE/EE
  ! The windows for low ell TT spectra exist in the files so that we can read these in in a nice array
  ! Re-order/crop later when the binning is performed
  allocate(this%full_windows(N_b_0_EE,1+SPT3G_windows_lmax-SPT3G_windows_lmin,N_s_0)) ! only referencing EE here because it has 44 bins
  allocate(current_window(1+SPT3G_windows_lmax-SPT3G_windows_lmin,N_s_0+1)) ! +1 more specs than needed because of ells in the file
  do i=1,N_b_0_EE ! only referencing EE here because it has 44 bins
    call file%ReadTextMatrix(trim(window_folder)//trim(numcat("window_",i))//trim(".txt"), current_window)
    this%full_windows(i,:,:) = current_window(:,2:)
  end do

  ! Calibration Covariance
  ! This is typically a 6x6 matrix, but here we crop out any unnecessary elements
  allocate(full_cal_covariance_matrix(N_freq_0*2,N_freq_0*2))
  call File%ReadTextMatrix(cal_covariance_filename, full_cal_covariance_matrix) ! Read in the full covariance matrix
  cal_rows_use = 0
  call MaskCalibrationCovarianceMatrix(spectra_to_fit_list, cal_rows_use) ! Find out which elements we need to keep
  allocate(this%inv_cal_covariance(sum(cal_rows_use),sum(cal_rows_use)))
  call CropAndInvertCalibrationCovarianceMatrix(full_cal_covariance_matrix, cal_rows_use, this%inv_cal_covariance) ! Crop down and invert

  ! Effective band centres
  allocate(nu_eff_matrix(5,N_freq_0)) ! pol gal dust, DSFG, Radio galaxies, tSZ
  call File%ReadTextMatrix(nu_eff_filename, nu_eff_matrix) ! Read in the effective band centres
  allocate(nu_eff_gal_cirrus(this%N_s,2))
  allocate(nu_eff_pol_gal_dust(this%N_s,2))
  allocate(nu_eff_DSFG(this%N_s,2))
  allocate(nu_eff_tSZ(this%N_s,2))
  call GrabEffectiveBandCentres(spectra_to_fit_list, nu_eff_matrix, nu_eff_gal_cirrus, nu_eff_pol_gal_dust, nu_eff_DSFG, nu_eff_tSZ) ! Sort the band centres from the big matrix into arrays that can be indexed by spectrum requested in the fit

  ! Prepare mapping from spectra to nuisance parameter indices
  call GetCalibrationParamIx(spectra_to_fit_list, calibration_ix) ! Calibration indices
  call GetPoissonParamIx(spectra_to_fit_list, spectra_to_fit_Poisson_ix) ! Poisson power indices
  call GrabCIBClusteringDecorrelationParamIx(spectra_to_fit_list, spectra_to_fit_CIB_cl_decorr_ix) ! CIB clustering decorrelation indices

  this%include_logdet = 1

  ! Give feedback
  print *, "SPT-3G 2018 TTTEEE: Likelihood successfully initialised!"
  print *, "Fitting spectra: (bins)"
  do i=1,this%N_s
    ! Check that no low-ell TT data is requested
    call this%CheckTTBinSelection(spectra_to_fit_list%Item(i), this%spec_bin_min(i), spec_ok)
    print *, spectra_to_fit_list%Item(i), "(", this%spec_bin_min(i), "-", this%spec_bin_max(i), ")"
    if (spec_ok .eq. 0) then
      print *, "The above spectrum/bin selection is invalid."
      call MPIStop("SPT-3G 2018 TTTEEE: invalid spectrum/bin selection!")
    end if
  end do
  if (do_late_crop .eqv. .true.) then
    print *, "(with some last-second omissions)"
  end if
  print *, "--------------------"

end subroutine SPT3G_2018_TTTEEE_InitData

# else

subroutine SPT3G_2018_TTTEEE_Ini_external(this,eSPT3G_windows_lmin, eSPT3G_windows_lmax, full_bandpower, full_bandpower_list_string, full_covariance_matrix, &
                                          full_covariance_list_string, full_beam_covariance_matrix,  full_beam_covariance_list_string ,      &
                                          full_cal_covariance_matrix, full_windows, full_window_list_string,     &
                                          nu_eff_matrix,nu_eff_list_string, spectra_to_fit_list_string,      &
                                          spec_bin_min_list_string,  spec_bin_max_list_string,late_crop_msk_string, &
                                          ecov_eval_cut_threshold,ecov_eval_large_number_replacement,beam_cov_scale, &
                                          aberration_coefficient, enu_0_galdust, eT_galdust, enu_0_CIB, eT_CIB, enu_0_tSZ, etSZCosmologyScalingEnabled,full_tSZ_template,ekSZCosmologyScalingEnabled,full_kSZ_template,einclude_logdet)
  class(TSPT3G_2018_TTTEEE_Likelihood),intent(inout) :: this
  !issue here....

  integer,intent(in) :: eSPT3G_windows_lmin, eSPT3G_windows_lmax
  real(mcp),intent(in) :: full_bandpower(N_b_0_total)
  character(LEN=*),intent(in) :: full_bandpower_list_string

  real(mcp),intent(in) :: full_covariance_matrix(N_b_0_total,N_b_0_total)
  character(LEN=*),intent(in) :: full_covariance_list_string

  real(mcp),intent(in) :: full_beam_covariance_matrix(N_b_0_total,N_b_0_total)
  character(LEN=*),intent(in) :: full_beam_covariance_list_string

  real(mcp),intent(in)::full_cal_covariance_matrix(N_freq_0*2,N_freq_0*2)
  

  real(mcp),intent(in):: full_windows(N_b_0_EE,1+eSPT3G_windows_lmax-eSPT3G_windows_lmin,N_s_0)
  character(LEN=*),intent(in) :: full_window_list_string

  real(mcp),intent(in) :: nu_eff_matrix(5,N_freq_0)
  character(LEN=*),intent(in) :: nu_eff_list_string

  character(LEN=*),intent(in) :: spectra_to_fit_list_string

  character(LEN=*),intent(in) :: spec_bin_min_list_string
  character(LEN=*),intent(in) :: spec_bin_max_list_string 

  character(LEN=*),intent(in) :: late_crop_msk_string
  real(mcp),intent(in):: ecov_eval_cut_threshold,ecov_eval_large_number_replacement,beam_cov_scale,aberration_coefficient
  real(mcp),intent(in):: enu_0_galdust, eT_galdust, enu_0_CIB, eT_CIB, enu_0_tSZ
  integer,intent(in) :: etSZCosmologyScalingEnabled,ekSZCosmologyScalingEnabled,einclude_logdet

  real(mcp),intent(in) :: full_tSZ_template(1+eSPT3G_windows_lmax-eSPT3G_windows_lmin)
  real(mcp),intent(in) :: full_kSZ_template(1+eSPT3G_windows_lmax-eSPT3G_windows_lmin)

  Type(TStringList) :: late_crop_msk_list
  Type(TStringList) :: spec_bin_min_list
  Type(TStringList) :: spec_bin_max_list
  
  integer i, stat,edo_late_crop
    
  SPT3G_windows_lmin = eSPT3G_windows_lmin
  SPT3G_windows_lmax = eSPT3G_windows_lmax

  ! Read in spectrum selection
  call spectra_to_fit_list%SetFromString(spectra_to_fit_list_string)
  this%N_s = spectra_to_fit_list%Count
  
  ! Read in bin selection for each spectrum and the total number of bins we are fitting
  allocate(this%spec_bin_min(this%N_s))
  this%spec_bin_min = -1
  allocate(this%spec_bin_max(this%N_s))
  this%spec_bin_max = 45
  call spec_bin_min_list%SetFromString(spec_bin_min_list_string)
  call spec_bin_max_list%SetFromString(spec_bin_max_list_string)
  

  this%N_b_total = 0
  do i=1, this%N_s
    call str2int(spec_bin_min_list%Item(i), this%spec_bin_min(i), stat)
    call str2int(spec_bin_max_list%Item(i), this%spec_bin_max(i), stat)
    if (((this%spec_bin_min(i) .le. 0) .or. (this%spec_bin_max(i) .ge. 45)) .or. (this%spec_bin_max(i) .lt. this%spec_bin_min(i))) then
      write(*,*) "SPT-3G 2018 TTTEEE: bad ell range selection for spectrum: "//spectra_to_fit_list%Item(i)
      STOP
    end if
    this%N_b_total = this%N_b_total + this%spec_bin_max(i) - this%spec_bin_min(i) + 1
  end do

  ! Check if a late crop is requested and read in the mask if necessary
  edo_late_crop = len(trim(late_crop_msk_string))
  if (edo_late_crop .gt. 0) then
    do_late_crop = .true.
    call late_crop_msk_list%SetFromString(late_crop_msk_string)
    allocate(late_crop_msk(this%N_b_total))
    do i=1, this%N_b_total
      call str2int(late_crop_msk_list%Item(i), late_crop_msk(i), stat)
    end do
  end if

  ! Determine how many spectra are TT vs TE vs EE and the total number of bins we are fitting
  call DetermineNumberOfSpectraByField(spectra_to_fit_list, this%N_s_TT, this%N_s_TE, this%N_s_EE)

  ! Determine how many different frequencies get used
  call DetermineNumberofFrequencies(spectra_to_fit_list, this%N_freq)

  ! Read in band power file name
  call full_bandpower_list%SetFromString(full_bandpower_list_string)

  ! Read in band power covariance matrix file name
  call full_covariance_list%SetFromString(full_covariance_list_string)

  
  cov_eval_cut_threshold = ecov_eval_cut_threshold
  cov_eval_large_number_replacement = ecov_eval_large_number_replacement

  call full_beam_covariance_list%SetFromString(full_beam_covariance_list_string)

  this%beam_cov_scale = beam_cov_scale

  call full_window_list%SetFromString(full_window_list_string)

  call nu_eff_list%SetFromString(nu_eff_list_string)

  this%aberration_coefficient = aberration_coefficient

  print_chisq = .false.

  print_spectrum = .False.
  
  ! Grab indices of spectra in fit order
  call this%GenerateStartStopIxForFit(spectra_to_fit_list, spectra_start_ix, spectra_stop_ix)

  ! Band Powers
  call this%GenerateStartStopIx(spectra_to_fit_list, full_bandpower_list, full_bandpower_start_ix, full_bandpower_stop_ix) ! Get the indices of the requested bins in the fitting order
  allocate(this%bandpowers(this%N_b_total))
  do i_spec=1, spectra_to_fit_list%Count ! Loop over full spectra and populate band power vector that will be used
    this%bandpowers(spectra_start_ix(i_spec):spectra_stop_ix(i_spec)) = full_bandpower(full_bandpower_start_ix(i_spec):full_bandpower_stop_ix(i_spec))
  end do

    ! Covariance Matrix
  call this%GenerateStartStopIx(spectra_to_fit_list, full_covariance_list, full_covariance_start_ix, full_covariance_stop_ix)  ! Get the indices of the requested bins in the fitting order
  allocate(this%bdp_covariance(this%N_b_total,this%N_b_total))
  do i_spec=1, spectra_to_fit_list%Count ! Loop over full covariance block and populate cov that will be used
    do j_spec=1, spectra_to_fit_list%Count
      this%bdp_covariance(spectra_start_ix(i_spec):spectra_stop_ix(i_spec),spectra_start_ix(j_spec):spectra_stop_ix(j_spec)) = full_covariance_matrix(full_covariance_start_ix(i_spec):full_covariance_stop_ix(i_spec),full_covariance_start_ix(j_spec):full_covariance_stop_ix(j_spec))
    end do
  end do
  allocate(this%bdp_covariance_pos_def(this%N_b_total,this%N_b_total))
  !!!call this%MakeCovariancePositiveDefinite(this%bdp_covariance, fiducial_covariance_filename, fiducial_covariance_list, this%bdp_covariance_pos_def) ! Ensure covariance is positive definite
  this%bdp_covariance_pos_def(:,:) = this%bdp_covariance(:,:)

  ! Beam Covariance Matrix
  call this%GenerateStartStopIx(spectra_to_fit_list, full_beam_covariance_list, full_beam_covariance_start_ix, full_beam_covariance_stop_ix)  ! Get the indices of the requested bins in the fitting order
  allocate(this%beam_covariance(this%N_b_total,this%N_b_total))
  do i_spec=1, spectra_to_fit_list%Count ! Loop over full covariance block and populate cov that will be used
    do j_spec=1, spectra_to_fit_list%Count
      this%beam_covariance(spectra_start_ix(i_spec):spectra_stop_ix(i_spec),spectra_start_ix(j_spec):spectra_stop_ix(j_spec)) = full_beam_covariance_matrix(full_beam_covariance_start_ix(i_spec):full_beam_covariance_stop_ix(i_spec),full_beam_covariance_start_ix(j_spec):full_beam_covariance_stop_ix(j_spec))
    end do
  end do
  this%beam_covariance = this%beam_covariance * this%beam_cov_scale

  ! Window Functions
  ! These are a bit trickier to handle due to the independent cuts possible for TT/TE/EE
  ! The windows for low ell TT spectra exist in the files so that we can read these in in a nice array
  ! Re-order/crop later when the binning is performed
  allocate(this%full_windows(N_b_0_EE,1+SPT3G_windows_lmax-SPT3G_windows_lmin,N_s_0)) ! only referencing EE here because it has 44 bins
  this%full_windows(:,:,:)  = full_windows(:,:,:) 
  
  ! Calibration Covariance
  ! This is typically a 6x6 matrix, but here we crop out any unnecessary elements
  cal_rows_use = 0
  call MaskCalibrationCovarianceMatrix(spectra_to_fit_list, cal_rows_use) ! Find out which elements we need to keep
  allocate(this%inv_cal_covariance(sum(cal_rows_use),sum(cal_rows_use)))
  call CropAndInvertCalibrationCovarianceMatrix(full_cal_covariance_matrix, cal_rows_use, this%inv_cal_covariance) ! Crop down and invert

  ! Effective band centres
  allocate(nu_eff_gal_cirrus(this%N_s,2))
  allocate(nu_eff_pol_gal_dust(this%N_s,2))
  allocate(nu_eff_DSFG(this%N_s,2))
  allocate(nu_eff_tSZ(this%N_s,2))
  call GrabEffectiveBandCentres(spectra_to_fit_list, nu_eff_matrix, nu_eff_gal_cirrus, nu_eff_pol_gal_dust, nu_eff_DSFG, nu_eff_tSZ) ! Sort the band centres from the big matrix into arrays that can be indexed by spectrum requested in the fit

! Prepare mapping from spectra to nuisance parameter indices
  call GetCalibrationParamIx(spectra_to_fit_list, calibration_ix) ! Calibration indices
  call GetPoissonParamIx(spectra_to_fit_list, spectra_to_fit_Poisson_ix) ! Poisson power indices
  call GrabCIBClusteringDecorrelationParamIx(spectra_to_fit_list, spectra_to_fit_CIB_cl_decorr_ix) ! CIB clustering decorrelation indices

  ! NEED TO DO THE FG PART AS WELL

  call SPT3G_2018_TTTEEE_Ini_Foregrounds(eSPT3G_windows_lmin,eSPT3G_windows_lmax, enu_0_galdust, eT_galdust, enu_0_CIB, eT_CIB, enu_0_tSZ, etSZCosmologyScalingEnabled,full_tSZ_template,ekSZCosmologyScalingEnabled,full_kSZ_template)

  if (einclude_logdet.NE.0) then
    this%include_logdet = 1
  else
    this%include_logdet = 0
  endif

end subroutine SPT3G_2018_TTTEEE_Ini_external
# endif



! Likelihood function
! Loops over the spectra, adding foregrounds etc., bins them
! Gets the final covariance and returns logL

#ifndef _STANDALONE_
function SPT3G_2018_TTTEEE_LogLike(this, CMB, Theory, DataParams) result(SPT_LogLike)
  implicit none
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  class(CMBParams) :: CMB
  class(TCosmoTheoryPredictions), target :: Theory
  real(mcp) :: DataParams(:)
  real(mcp) :: SPT_LogLike, detcov, chisq

  ! Helpers for shuffling the spectra around
  integer i, j, ii, jj, i_spec, ix, ix_1, ix_2, spec_length
  real(mcp), dimension((SPT3G_windows_lmax+1)*3) :: Theory_Cl
  real(mcp),dimension(6) :: CMBParams 

  call Theory%ClArray(Theory_Cl(:SPT3G_windows_lmax+1),CL_T,CL_T)
  call Theory%ClArray(Theory_Cl(SPT3G_windows_lmax+1:(SPT3G_windows_lmax+1)*2),CL_E,CL_E)
  call Theory%ClArray(Theory_Cl((SPT3G_windows_lmax+1)*2:(SPT3G_windows_lmax+1)*3),CL_T,CL_E)
  
  CMBParams( id_H0       ) = CMB%H0
  CMBParams( id_omb      ) = CMB%omb
  CMBParams( id_sigma_8  ) = Theory%sigma_8
  CMBParams( id_omm      ) = CMB%omc + CMB%omnu + CMB%omb
  CMBParams( id_ns       ) = CMB%InitPower(ns_index)
  CMBParams( id_tau      ) = CMB%tau

  SPT_LogLike = this%SPT3G_2018_TTTEEE_LogLike_external(Theory_Cl,CMBParams,DataParams)
end function
#endif


function SPT3G_2018_TTTEEE_LogLike_external(this, Theory_Cl,CMBParams,DataParams) result(SPT_LogLike)
  implicit none
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  real(mcp) :: DataParams(:),CMBParams(:)
  real(mcp), dimension((SPT3G_windows_lmax+1)*3) :: Theory_Cl
  real(mcp) :: SPT_LogLike, detcov, chisq

  ! Helpers for shuffling the spectra around
  integer i, j, ii, jj, i_spec, ix, ix_1, ix_2, spec_length
  character(LEN=:), allocatable :: current_field, current_freq_1, current_freq_2

  ! Unbinned spectra and the various foreground components/corrections
  real(mcp), dimension(SPT3G_windows_lmax+1) :: current_Dl_theory_unbinned_CMB_only
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: current_Dl_theory_unbinned
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_foregrounds ! Holds copies of foregrounds for printing

  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: ssl_correction
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: aberration_correction

  ! Binned spectra
  real(mcp), allocatable :: current_Dl_theory_binned(:)
  real(mcp), dimension(this%N_b_total) :: Dl_theory_binned
  real(mcp), dimension(this%N_b_total) :: Delta_data_model,tmp
  real(mcp), allocatable :: Delta_data_model_final(:)

  ! Complete covariance (band powers and beam)
  real(mcp), dimension(this%N_b_total,this%N_b_total) :: cov_for_logl !, cov_for_logl_copy ! Copy needed for chisq calculation if requested
  real(mcp), allocatable :: cov_for_logl_final(:,:), cov_for_logl_copy_final(:,:) ! Needed for late crop to ignore band powers in the middle of the vector

  ! Vector of calibration parameters
  real(mcp), dimension(sum(cal_rows_use)) :: cal_vec, cal_vec_buf
  integer :: cal_ix_arr(2*N_freq_0) = (/iTcal90, iTcal150, iTcal220, iEcal90, iEcal150, iEcal220/)

  ! Number of bdps cropped out in the last step
  integer num_cropped,n,info
  character(len=20)::chumber

  !!!open(122,file="cl_retest.txt")
  !!!write(122,*) Theory_Cl
  !!!write(122,*) DataParams
  !!!write(122,*) CMBParams
  !!!close(122)

  ! Loop over requested spectra
  do i_spec=1,this%N_s

    ! Identify spectrum (get field and frequencies)
    call IdSpec(spectra_to_fit_list%Item(i_spec), current_field, current_freq_1, current_freq_2)

    ! Grab theory Dls
    ! these start at ell = 1 (which is always zero) and go to ell = 3201 (as requested)
#ifndef _STANDALONE_
    if (fixed_theory_input_file .eq. "USE_CAMB") then
#endif
      
      if (current_field .eq. "TT") then
        current_Dl_theory_unbinned_CMB_only = Theory_Cl(:SPT3G_windows_lmax+1)
      else if (current_field .eq. "TE") then
        current_Dl_theory_unbinned_CMB_only = Theory_Cl((SPT3G_windows_lmax+1)*2:(SPT3G_windows_lmax+1)*3)
      else if (current_field .eq. "EE") then
        current_Dl_theory_unbinned_CMB_only = Theory_Cl((SPT3G_windows_lmax+1):(SPT3G_windows_lmax+1)*2)
      end if

      ! Crop down to spectrum we will modify (1-3200 expected by window functions)
      current_Dl_theory_unbinned(SPT3G_windows_lmin:SPT3G_windows_lmax) = current_Dl_theory_unbinned_CMB_only(SPT3G_windows_lmin:SPT3G_windows_lmax)
#ifndef _STANDALONE_
    else

      ! Use supplied spectrum
      call this%GrabTheorySpectrumOverwrite(current_field, current_Dl_theory_unbinned)
      current_Dl_theory_unbinned_CMB_only(SPT3G_windows_lmin:SPT3G_windows_lmax) = current_Dl_theory_unbinned

    end if
#endif

    ! Add foregrounds
    if (current_field .eq. "TT") then
      !KARIM: Passing the CMB and Theory objects to the TT foregrounds function because tSZ and kSZ have a weak cosmology scaling (i.e. need the current values of sigma_8, H0, omega_m, and some other parameters)
      call this%AddTTForegrounds(i_spec, CMBParams, DataParams, current_Dl_theory_unbinned, Dl_foregrounds)
    else if (current_field .eq. "TE") then
      call this%AddTEForegrounds(i_spec, DataParams, current_Dl_theory_unbinned, Dl_foregrounds)
    else if (current_field .eq. "EE") then
      call this%AddEEForegrounds(i_spec, DataParams, current_Dl_theory_unbinned, Dl_foregrounds)
    end if

    !!write(chumber,*) i_spec
    !!chumber=adjustl(chumber)
    !!open(unit=334,file="fg_"//trim(chumber)//".dat")
    !!write(*,*) "fg_"//trim(chumber)//".dat"
    !!write(334,*) Dl_foregrounds
    !!close(334)

#ifndef _STANDALONE_
    ! Write smooth theory out with all foreground components
    if (print_spectrum) then
      call WriteOutTheorySpectrum(spectra_to_fit_list%Item(i_spec), 33+i_spec, current_Dl_theory_unbinned_CMB_only, current_Dl_theory_unbinned, Dl_foregrounds)
    end if
#endif

    ! Bin
    allocate(current_Dl_theory_binned(spectra_stop_ix(i_spec)-spectra_start_ix(i_spec)+1))
    call this%BinPowerSpectrum(current_Dl_theory_unbinned, i_spec, current_Dl_theory_binned)

    ! Slot into long vector of binned theory
    Dl_theory_binned(spectra_start_ix(i_spec):spectra_stop_ix(i_spec)) = current_Dl_theory_binned

    ! Tidy up
    deallocate(current_Dl_theory_binned) ! Different spectra have different numbers of bins, so this needs to be flexible in length
    Dl_foregrounds = 0

  end do

  ! Calculate difference of theory and data
  Delta_data_model = this%bandpowers - Dl_theory_binned

  ! Add the beam coariance to the band power covariance
  cov_for_logl = this%bdp_covariance_pos_def
  call this%AddBeamCovariance(Dl_theory_binned, cov_for_logl)
  !cov_for_logl_copy = cov_for_logl

  ! Final crop to ignore select band powers
  if (do_late_crop .eqv. .true.) then
    num_cropped = sum(late_crop_msk)
    allocate(Delta_data_model_final(this%N_b_total-num_cropped))
    allocate(cov_for_logl_final(this%N_b_total-num_cropped,this%N_b_total-num_cropped))

    ii = 1
    do i=1, this%N_b_total
      if (late_crop_msk(i) .eq. 0) then
        Delta_data_model_final(ii) = Delta_data_model(i)
        jj = 1
        do j=1, this%N_b_total
          if (late_crop_msk(j) .eq. 0) then
            cov_for_logl_final(ii,jj) = cov_for_logl(i,j)
            jj = jj + 1
          end if
        end do
        ii = ii+1
      end if
    end do
    !allocate(cov_for_logl_copy_final(this%N_b_total-num_cropped,this%N_b_total-num_cropped))
    !cov_for_logl_copy_final = cov_for_logl_final
  else
    allocate(Delta_data_model_final(this%N_b_total))
    Delta_data_model_final = Delta_data_model
    allocate(cov_for_logl_final(this%N_b_total,this%N_b_total))
    cov_for_logl_final = cov_for_logl
    !allocate(cov_for_logl_copy_final(this%N_b_total,this%N_b_total))
    !cov_for_logl_copy_final = cov_for_logl_copy
  end if


  
# ifndef _STANDALONE_
  if (print_chisq) then
    allocate(cov_for_logl_copy_final, source=cov_for_logl_final)
    !cov_for_logl_copy_final = cov_for_logl_final
  end if
  SPT_LogLike = Matrix_GaussianLogLikeDouble(cov_for_logl_final, Delta_data_model_final)
  ! Print chisq
  if (print_chisq) then
     ! Get the contribution from the covariance determinant part to the logL
     !KARIM: This is a CosmoMC wrapper for some Lapack routines
     detcov = Matrix_GaussianLogLikeDouble(cov_for_logl_copy_final, Delta_data_model_final*0)
     chisq = 2*(SPT_LogLike - detcov)
     print *, "SPT-3G 2018 TTTEEE: chi square = ", chisq
     call MPIStop("SPT-3G 2018 TTTEEE: Completed chi square calculation.")
     deallocate(cov_for_logl_copy_final)
  end if
# else
    n = size(Delta_data_model_final)
    call dpotrf ('L', n, cov_for_logl_final, n, info)
    if (info/=0) then
      write(*,*) "total covariance not positive definite"
      SPT_LogLike = -1e30
      stop
      return
    endif
    SPT_LogLike = 0
    if (this%include_logdet.NE.0) then
      !Log Det term:
      do i=1, n
          SPT_LogLike = SPT_LogLike  + log(cov_for_logl_final(i,i))
      end do
    endif
    tmp = Delta_data_model_final
    call DPOTRS('L', N, 1, cov_for_logl_final, n, tmp, n, INFO )
    if (info/=0) then
      write(*,*) "total covariance solve failed"
      SPT_LogLike = -1e30

      return
    endif

    !Add together
    !print *,SPT_LogLike
    SPT_LogLike = SPT_LogLike + dot_product(tmp,Delta_data_model_final)/2.
    !print *,SPT_LogLike
#endif
  ! Grab vector of calibration parameters used and offset from unity
  ix = 1
  do i=1, N_freq_0*2
    if (cal_rows_use(i) .eq. 1) then
      cal_vec(ix) = DataParams(cal_ix_arr(i))
      ix = ix + 1
    end if
  end do
  cal_vec = cal_vec-1.0_mcp

  ! Apply calibration prior
#ifndef _STANDALONE_
  !KARIM: This is a CosmoMC wrapper for some Lapack routines
  SPT_LogLike = SPT_LogLike + 0.5d0 * Matrix_QuadForm(this%inv_cal_covariance, cal_vec)

#else
  cal_vec_buf = MatMul(this%inv_cal_covariance,cal_vec)
  SPT_LogLike = SPT_LogLike + 0.5d0 * dot_product(cal_vec,cal_vec_buf)
#endif
  deallocate(Delta_data_model_final)
  deallocate(cov_for_logl_final)
end function SPT3G_2018_TTTEEE_LogLike_external


!function SPT3G_TTTEEE_LogLike
! This is where the magic happens
! Get CAMB spectra
! Modify as needed (foregrounds, binning etc.)
! Calculate and return logL

!---------------------------------------------------------!
!             FOREGROUND FUNCTIONS AND BINNING            !
!---------------------------------------------------------!

! Start and stop ell expected by binning function is 1-3200 (matching what's in the window files)
! i_spec is the index of the spectrum in the list of specs to fit
subroutine BinPowerSpectrum(this, power_spectrum, i_spec_to_bin, power_spectrum_binned)
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: power_spectrum
  real(mcp) :: power_spectrum_binned(*)
  integer, intent(in) :: i_spec_to_bin
  character(LEN=:), allocatable :: field, freq_1, freq_2
  integer i, i_spec, j_spec, ix, ix_1, ix_2, spec_length, first_bin, last_bin
  real(mcp), allocatable :: window(:,:)

  ! Identify the correct window function to use for this spectrum
  allocate(window(1+SPT3G_windows_lmax-SPT3G_windows_lmin,N_s_0))
  do i_spec=1, N_s_0 ! This loops over the spectra given in the window file
    if (full_window_list%Item(i_spec) .eq. spectra_to_fit_list%Item(i_spec_to_bin)) then
      window = this%full_windows(:,:,i_spec)
    end if
  end do

  first_bin = this%spec_bin_min(i_spec_to_bin)
  last_bin = this%spec_bin_max(i_spec_to_bin)
  spec_length = last_bin-first_bin+1

  ! Bin power spectrum
  call dgemv('N',spec_length,SPT3G_windows_lmax-SPT3G_windows_lmin+1,1.0d0,&
    window(first_bin:last_bin,:),&
    spec_length,power_spectrum,1,0d0,power_spectrum_binned,1)

end subroutine

! Master method for TT foregrounds
! Slightly different calling signature to others, since it includes cosmology-dependent foregrounds
subroutine AddTTForegrounds(this, i_spec, CMBParams, DataParams, Dl_theory, Dl_foregrounds)
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  integer, intent(in) :: i_spec
  real(mcp), intent(in) :: DataParams(:),CMBParams(:)
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_theory
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_foregrounds

  call ApplySuperSampleLensing(DataParams(iKappa), Dl_theory, Dl_foregrounds)

  call ApplyAberrationCorrection(this%aberration_coefficient, Dl_theory, Dl_foregrounds)

  call AddPoissonPower(DataParams(spectra_to_fit_Poisson_ix(i_spec)), Dl_theory, Dl_foregrounds)

  call AddGalacticDust(DataParams(iTT_GalDust_Amp), DataParams(iTT_GalDust_Alpha), DataParams(iTT_GalDust_Beta), nu_eff_gal_cirrus(i_spec,1), nu_eff_gal_cirrus(i_spec,2), Dl_theory, Dl_foregrounds)

  call AddCIBClustering(DataParams(iTT_CIB_Clustering_Amp), DataParams(iTT_CIB_Clustering_Alpha), DataParams(iTT_CIB_Clustering_Beta), nu_eff_DSFG(i_spec,1), nu_eff_DSFG(i_spec,2), DataParams(spectra_to_fit_CIB_cl_decorr_ix(i_spec,1)), DataParams(spectra_to_fit_CIB_cl_decorr_ix(i_spec,2)), Dl_theory, Dl_foregrounds)

  call AddtSZ(DataParams(iTT_tSZ_Amp), nu_eff_tSZ(i_spec,1), nu_eff_tSZ(i_spec,2), CMBParams(id_H0), CMBParams(id_sigma_8), CMBParams(id_omb), Dl_theory, Dl_foregrounds)

  call AddtSZCIBCorrelation(DataParams(iTT_tSZ_CIB_corr),&
                            DataParams(iTT_tSZ_Amp),&
                            DataParams(iTT_CIB_Clustering_Amp),&
                            DataParams(iTT_CIB_Clustering_Alpha),&
                            DataParams(iTT_CIB_Clustering_Beta),&
                            DataParams(spectra_to_fit_CIB_cl_decorr_ix(i_spec,1)),&
                            DataParams(spectra_to_fit_CIB_cl_decorr_ix(i_spec,2)),&
                            nu_eff_DSFG(i_spec,1),&
                            nu_eff_DSFG(i_spec,2),&
                            nu_eff_tSZ(i_spec,1),&
                            nu_eff_tSZ(i_spec,2),&
                            CMBParams(id_H0),&
                            CMBParams(id_sigma_8),&
                            CMBParams(id_omb),&
                            Dl_theory,&
                            Dl_foregrounds)
  call addkSZ(DataParams(iTT_kSZ_Amp), CMBParams(id_H0), CMBParams(id_sigma_8), CMBParams(id_omb), CMBParams(id_omm),CMBParams(id_ns_index),CMBParams(id_tau), Dl_theory, Dl_foregrounds)

  call ApplyCalibration(DataParams(calibration_ix(i_spec,1)), DataParams(calibration_ix(i_spec,2)), DataParams(calibration_ix(i_spec,3)), DataParams(calibration_ix(i_spec,4)), Dl_theory, Dl_foregrounds)

end subroutine AddTTForegrounds

! Master method for TE foregrounds
subroutine AddTEForegrounds(this, i_spec, DataParams, Dl_theory, Dl_foregrounds)
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  integer, intent(in) :: i_spec
  real(mcp), intent(in) :: DataParams(:)
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_theory
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_foregrounds

  call ApplySuperSampleLensing(DataParams(iKappa), Dl_theory, Dl_foregrounds)

  call ApplyAberrationCorrection(this%aberration_coefficient, Dl_theory, Dl_foregrounds)

  call AddGalacticDust(DataParams(iTE_GalDust_Amp), DataParams(iTE_GalDust_Alpha), DataParams(iTE_GalDust_Beta), nu_eff_pol_gal_dust(i_spec,1), nu_eff_pol_gal_dust(i_spec,2), Dl_theory, Dl_foregrounds)

  call ApplyCalibration(DataParams(calibration_ix(i_spec,1)), DataParams(calibration_ix(i_spec,2)), DataParams(calibration_ix(i_spec,3)), DataParams(calibration_ix(i_spec,4)), Dl_theory, Dl_foregrounds)

end subroutine AddTEForegrounds

! Master method for EE foregrounds
subroutine AddEEForegrounds(this, i_spec, DataParams, Dl_theory, Dl_foregrounds)
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  integer, intent(in) :: i_spec
  real(mcp), intent(in) :: DataParams(:)
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_theory
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_foregrounds

  call ApplySuperSampleLensing(DataParams(iKappa), Dl_theory, Dl_foregrounds)

  call ApplyAberrationCorrection(this%aberration_coefficient, Dl_theory, Dl_foregrounds)

  call AddPoissonPower(DataParams(spectra_to_fit_Poisson_ix(i_spec)), Dl_theory, Dl_foregrounds)

  call AddGalacticDust(DataParams(iEE_GalDust_Amp), DataParams(iEE_GalDust_Alpha), DataParams(iEE_GalDust_Beta), nu_eff_pol_gal_dust(i_spec,1), nu_eff_pol_gal_dust(i_spec,2), Dl_theory, Dl_foregrounds)

  call ApplyCalibration(DataParams(calibration_ix(i_spec,1)), DataParams(calibration_ix(i_spec,2)), DataParams(calibration_ix(i_spec,3)), DataParams(calibration_ix(i_spec,4)), Dl_theory, Dl_foregrounds)

end subroutine AddEEForegrounds

!---------------------------------------------------------!
!               COVARIANCE HELPER FUNCTIONS               !
!---------------------------------------------------------!

! Master routine to make the covariance positive definite
! Checks if the matrix already is positive definite, otherwise:
! Change basis, rescale, fix eigenvalues, undo rescaling, undo change of basis
# ifndef _STANDALONE_
subroutine MakeCovariancePositiveDefinite(this, cov_input, fiducial_covariance_filename, fiducial_covariance_list, cov_output)
  implicit none
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  real(mcp), intent(in) :: cov_input(:,:)
  character(LEN=:), allocatable :: fiducial_covariance_filename
  Type(TStringList), intent(in) :: fiducial_covariance_list
  real(mcp), intent(out) :: cov_output(:,:)

  integer N
  logical is_pos_def

  real(mcp), allocatable :: full_fiducial_covariance_matrix(:,:)
  real(mcp), allocatable :: fiducial_covariance_matrix(:,:)

  real(mcp), allocatable :: tmp(:)
  integer i, j, i_spec, j_spec

  real(mcp), allocatable :: cov_in_new_basis(:,:), cov_in_new_basis_scaled(:,:), cov_in_new_basis_scaled_fixed(:,:), scaling_mat(:,:)

  !print *, "before pos def check"
  ! Check if the matrix is already positive definite
  call IsMatrixPosDef(cov_input, is_pos_def)
  if (is_pos_def .eqv. .true.) then
    cov_output = cov_input
    return
  end if
  
  print *, "SPT-3G 2018 TTTEEE: attempting to make covariance positive definite."

  ! Get the size of the matrix and allocate all the intermediate arrays needed
  N = size(cov_input, DIM=1)
  allocate(cov_in_new_basis(N,N))
  allocate(cov_in_new_basis_scaled(N,N))
  allocate(cov_in_new_basis_scaled_fixed(N,N))
  allocate(scaling_mat(N,N))

  ! Read in the fiducial covariance and crop it the same way the data covariance has been cropped
  allocate(full_fiducial_covariance_matrix(N_b_0_total,N_b_0_total))
  call File%ReadTextMatrix(fiducial_covariance_filename, full_fiducial_covariance_matrix) ! Read in the full covariance matrix
  call this%GenerateStartStopIx(spectra_to_fit_list, fiducial_covariance_list, fiducial_covariance_start_ix, fiducial_covariance_stop_ix)  ! Get the indices of the requested bins in the fitting order
  allocate(fiducial_covariance_matrix(this%N_b_total,this%N_b_total))
  do i_spec=1, spectra_to_fit_list%Count ! Loop over full covariance block and populate cov that will be used
    do j_spec=1, spectra_to_fit_list%Count
      fiducial_covariance_matrix(spectra_start_ix(i_spec):spectra_stop_ix(i_spec),spectra_start_ix(j_spec):spectra_stop_ix(j_spec)) = full_fiducial_covariance_matrix(fiducial_covariance_start_ix(i_spec):fiducial_covariance_stop_ix(i_spec),fiducial_covariance_start_ix(j_spec):fiducial_covariance_stop_ix(j_spec))
    end do
  end do

  ! Change basis!
  call ChangeBasisToFiducial(cov_input, fiducial_covariance_matrix, .true., cov_in_new_basis)

  ! Rescale
  do i=1, N
    do j=1, N
      scaling_mat(i,j) = SQRT(cov_in_new_basis(i,i)*cov_in_new_basis(j,j))
    end do
  end do
  cov_in_new_basis_scaled = cov_in_new_basis/scaling_mat

  ! Fix the negative eigenvalues
  call FixEvalsOfMatrix(cov_in_new_basis_scaled, cov_in_new_basis_scaled_fixed)

  ! Undo rescaling
  cov_in_new_basis_scaled_fixed = cov_in_new_basis_scaled_fixed * scaling_mat

  ! Undo change of basis
  call ChangeBasisToFiducial(cov_in_new_basis_scaled_fixed, fiducial_covariance_matrix, .false., cov_output)

end subroutine MakeCovariancePositiveDefinite

! Checks if a matrix is positive definite
subroutine IsMatrixPosDef(cov_input, is_pos_def)
  implicit none
  real(mcp), intent(in) :: cov_input(:,:)
  logical, intent(out) :: is_pos_def
  integer N

  real(mcp), allocatable :: tmp(:)
  integer i, ierr, tmpsize
  integer, allocatable :: ipiv(:)

  real(mcp), allocatable :: cov_eigenvalues(:)
  real(mcp), allocatable :: cov_eigenvectors(:,:)

  ! Get the size of the matrix and allocate all the intermediate arrays needed
  N = size(cov_input, DIM=1)
  allocate(ipiv(N))
  allocate(cov_eigenvalues(N))
  allocate(cov_eigenvectors(N,N))

  cov_eigenvectors = cov_input
  tmpsize =  max( (ILAENV_wrap(1,'DSYTRD','U',N,N,N,N)+2)*N,max(1,3*N-1)) ! Taken from Matrix_utils
  allocate(tmp(tmpsize))
  call DSYEV('V','U',N,cov_eigenvectors,N,cov_eigenvalues,tmp,tmpsize,ierr)

  ! Check if any eigenvalues are negative
  is_pos_def = .true.
  do i=1,N
    if (cov_eigenvalues(i) .LT. 0) then
      is_pos_def = .false.
    end if
  end do

end subroutine IsMatrixPosDef

! Changes the basis of the matrix using the eigenvectors of another
! Goes either fowards or backwards
subroutine ChangeBasisToFiducial(cov_input, cov_fiducial, forward, cov_output)
  implicit none
  real(mcp), intent(in) :: cov_input(:,:)
  real(mcp), intent(in) :: cov_fiducial(:,:)
  logical, intent(in) :: forward
  real(mcp), intent(out) :: cov_output(:,:)

  integer N
  real(mcp), allocatable :: tmp(:)
  integer i, ierr, tmpsize
  integer, allocatable :: ipiv(:)

  real(mcp), allocatable :: cov_eigenvalues(:)
  real(mcp), allocatable :: cov_eigenvectors(:,:)
  real(mcp), allocatable :: cov_eigenvectors_inv(:,:)
  real(mcp), allocatable :: cov_intermediate_matmul(:,:)

  ! Get the size of the matrix and allocate all the intermediate arrays needed
  N = size(cov_input, DIM=1)
  allocate(ipiv(N))
  allocate(cov_eigenvalues(N))
  allocate(cov_eigenvectors(N,N))
  allocate(cov_eigenvectors_inv(N,N))
  allocate(cov_intermediate_matmul(N,N))

  ! Calculate the eigenvectors
  cov_eigenvectors = cov_fiducial
  tmpsize =  max( (ILAENV_wrap(1,'DSYTRD','U',N,N,N,N)+2)*N,max(1,3*N-1)) ! Taken from Matrix_utils
  allocate(tmp(tmpsize))
  call DSYEV('V','U',N,cov_eigenvectors,N,cov_eigenvalues,tmp,tmpsize,ierr)

  ! Get the inverse of the eigenvectors matrix
  cov_eigenvectors_inv = cov_eigenvectors
  call DGETRF(N,N,cov_eigenvectors_inv,N,ipiv,ierr)
  call DGETRI(N,cov_eigenvectors_inv,N,ipiv,tmp,tmpsize,ierr)

  ! Change basis (either to fiducial or back)
  !KARIM: This is a CosmoMC wrapper for some Lapack routines
  if (forward .eqv. .true.) then

    call Matrix_Mult(cov_input,cov_eigenvectors,cov_intermediate_matmul,1.0d0,0d0)
    call Matrix_Mult(cov_eigenvectors_inv,cov_intermediate_matmul,cov_output,1.0d0,0d0)

  else

    call Matrix_Mult(cov_input,cov_eigenvectors_inv,cov_intermediate_matmul,1.0d0,0d0)
    call Matrix_Mult(cov_eigenvectors,cov_intermediate_matmul,cov_output,1.0d0,0d0)
  end if

end subroutine ChangeBasisToFiducial

! Fix negative eigenvalues of a matrix by setting anything below a threshold to some fixed value
subroutine FixEvalsOfMatrix(cov_input, cov_output)
  implicit none
  real(mcp), intent(in) :: cov_input(:,:)
  real(mcp), intent(out) :: cov_output(:,:)
  integer N

  real(mcp), allocatable :: tmp(:)
  integer i, ierr, tmpsize
  integer, allocatable :: ipiv(:)

  real(mcp), allocatable :: cov_eigenvalues(:)
  real(mcp), allocatable :: cov_eigenvectors(:,:), cov_eigenvectors_inv(:,:), cov_eigenvalues_2d(:,:), cov_intermediate_matmul(:,:)

  ! Get the size of the matrix and allocate all the intermediate arrays needed
  N = size(cov_input, DIM=1)
  allocate(ipiv(N))
  allocate(cov_eigenvalues(N))
  allocate(cov_eigenvectors(N,N))
  allocate(cov_eigenvectors_inv(N,N))
  allocate(cov_eigenvalues_2d(N,N))
  allocate(cov_intermediate_matmul(N,N))

  ! Get the eigenvalues and eigenvectors of the covariance

  ! Alt way of calculating eigenvectors!
  cov_eigenvectors = cov_input
  tmpsize =  max( (ILAENV_wrap(1,'DSYTRD','U',N,N,N,N)+2)*N,max(1,3*N-1)) ! Taken from Matrix_utils
  allocate(tmp(tmpsize))
  call DSYEV('V','U',N,cov_eigenvectors,N,cov_eigenvalues,tmp,tmpsize,ierr)

  ! Set negative eigenvalues to a large positive number
  cov_eigenvalues_2d = 0
  do i=1,N
    cov_eigenvalues_2d(i,i) = cov_eigenvalues(i)
    if (cov_eigenvalues(i) .LT. cov_eval_cut_threshold) then
      cov_eigenvalues_2d(i,i) = cov_eval_large_number_replacement
    end if
  end do

  ! Get the inverse of the eigenvectors matrix
  cov_eigenvectors_inv = cov_eigenvectors
  call DGETRF(N,N,cov_eigenvectors_inv,N,ipiv,ierr)
  call DGETRI(N,cov_eigenvectors_inv,N,ipiv,tmp,tmpsize,ierr)

  ! Cast back into a covariance matrix
  !KARIM: This is a CosmoMC wrapper for some Lapack routines
  call Matrix_Mult(cov_eigenvalues_2d,cov_eigenvectors_inv,cov_intermediate_matmul,1.0d0,0d0)
  call Matrix_Mult(cov_eigenvectors,cov_intermediate_matmul,cov_output,1.0d0,0d0)

end subroutine FixEvalsOfMatrix
#endif

! Add the beam coariance to the band power covariance
subroutine AddBeamCovariance(this, binned_spec, bdp_cov)
  implicit none
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  real(mcp), dimension(this%N_b_total), intent(in) :: binned_spec
  real(mcp), dimension(this%N_b_total,this%N_b_total), intent(out) :: bdp_cov
  integer i, j

  do i=1,this%N_b_total
    do j=1,this%N_b_total
      bdp_cov(i,j) = bdp_cov(i,j) + binned_spec(i)*binned_spec(j)*this%beam_covariance(i,j)
    end do
  end do

end subroutine AddBeamCovariance

!---------------------------------------------------------!
!                      MISC HELPERS                       !
!---------------------------------------------------------!

! Check the TT bin selection
subroutine CheckTTBinSelection(this, check_spec, check_bin_min, spec_ok)
  implicit none
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  character(*), intent(in) :: check_spec
  integer, intent(in) :: check_bin_min
  integer, intent(inout) :: spec_ok
  integer this_min_req
  character(LEN=:), allocatable :: current_field, current_freq_1, current_freq_2

  call IdSpec(check_spec, current_field, current_freq_1, current_freq_2)

  spec_ok = 1
  this_min_req = 1
  if (current_field .eq. "TT") then
    if ((current_freq_1 .eq. "90") .and. (current_freq_2 .eq. "90")) then
      this_min_req = 10
    else if ((current_freq_1 .eq. "90") .and. (current_freq_2 .eq. "150")) then
      this_min_req = 10
    else if ((current_freq_1 .eq. "90") .and. (current_freq_2 .eq. "220")) then
      this_min_req = 10
    else if ((current_freq_1 .eq. "150") .and. (current_freq_2 .eq. "150")) then
      this_min_req = 10
    else if ((current_freq_1 .eq. "150") .and. (current_freq_2 .eq. "220")) then
      this_min_req = 15
    else if ((current_freq_1 .eq. "220") .and. (current_freq_2 .eq. "220")) then
      this_min_req = 15
    end if
  end if

  if (check_bin_min .lt. this_min_req) then
    spec_ok = 0
  end if

end subroutine

! Grab the right component of the overwrite theory spectrum for a given field
subroutine GrabTheorySpectrumOverwrite(this, current_field, theory_spec)
  implicit none
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  character(LEN=2), intent(in) :: current_field
  real(mcp), dimension(SPT3G_windows_lmax+1), intent(out) :: theory_spec

  if (current_field .eq. "TT") then
    theory_spec = this%theory_spectrum_overwrite(2,:)
  else if (current_field .eq. "TE") then
    theory_spec = this%theory_spectrum_overwrite(3,:)
  else if (current_field .eq. "EE") then
    theory_spec = this%theory_spectrum_overwrite(4,:)
  end if
end subroutine GrabTheorySpectrumOverwrite

! Generates the start and stop indices given the requested spectra to fit for a specific file
! (this is to deal with potentially different ordering of requested spectra and how they appear in a file and make cropping easier)
subroutine GenerateStartStopIx(this, spec_to_fit, spec_in_file, start_ix_fit_order, stop_ix_fit_order)
  implicit none
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  Type(TStringList), intent(in) :: spec_to_fit, spec_in_file ! List of spectra to fit, List of spectra in the file
  integer, allocatable, intent(out) :: start_ix_fit_order(:), stop_ix_fit_order(:) ! Output indices translating between fit list and file order

  integer, allocatable :: start_ix(:), stop_ix(:) ! Intermediate indices denoting file ordering
  integer i_spec, j_spec, ix, ix_1, ix_2, spec_length, start_offset, stop_offset
  character(LEN=:), allocatable :: current_field, current_freq_1, current_freq_2

  ! First generate the list of start and stop indices in the file order
  allocate(start_ix(spec_in_file%Count))
  allocate(stop_ix(spec_in_file%Count))
  allocate(start_ix_fit_order(spec_to_fit%Count))
  allocate(stop_ix_fit_order(spec_to_fit%Count))

  start_ix = 0

  do i_spec=2, spec_in_file%Count ! Can skip over first element because it always starts at zero

    spec_length = 44

    ! Previous start ix plus length of size of last spec
    start_ix(i_spec) = start_ix(i_spec-1) + spec_length

  end do
  start_ix = start_ix + 1 ! Fortran indexing starts at 1

  ! Stop indices are given by start indixes with offset
  stop_ix(:spec_in_file%Count-1) = start_ix(2:) - 1
  stop_ix(spec_in_file%Count) = N_b_0_total

  ! Now reshuffle to get the (fit) order requested
  ! Also adjust for ell bin cropping here
  do i_spec=1, this%N_s ! This loops over the spectra requested in the fit
    do j_spec=1, N_s_0 ! This loops over the spectra given in the file
      if (spec_to_fit%Item(i_spec) .eq. spec_in_file%Item(j_spec)) then

        start_offset = this%spec_bin_min(i_spec)-1
        stop_offset = N_b_0_TT-this%spec_bin_max(i_spec)

        start_ix_fit_order(i_spec) = start_ix(j_spec)+start_offset
        stop_ix_fit_order(i_spec) = stop_ix(j_spec)-stop_offset

      end if
    end do
  end do

end subroutine GenerateStartStopIx

! Same thing, but simplified for only the fitting order
! Returns start and stop indices for the spectra
subroutine GenerateStartStopIxForFit(this, spec_to_fit, start_ix, stop_ix)
  implicit none
  class(TSPT3G_2018_TTTEEE_Likelihood) :: this
  Type(TStringList), intent(in) :: spec_to_fit ! List of spectra to fit, List of spectra in the file
  integer, allocatable, intent(out) :: start_ix(:), stop_ix(:) ! Output indices

  integer i_spec, j_spec, ix, ix_1, ix_2, spec_length, start_offset, stop_offset
  character(LEN=:), allocatable :: current_field, current_freq_1, current_freq_2

  ! Generate a list of start and stop ix for the spectra in the file ordering
  allocate(start_ix(spec_to_fit%Count))
  allocate(stop_ix(spec_to_fit%Count))
  start_ix = 0

  ! Catch if only one spectrum is requested
  if (spec_to_fit%Count .eq. 1) then
    start_ix = 1
    stop_ix = this%N_b_total
    return
  end if

  do i_spec=2, spec_to_fit%Count ! Can skip over first element because it always starts at zero

    spec_length = this%spec_bin_max(i_spec-1) - this%spec_bin_min(i_spec-1) + 1

    ! Previous start ix plus length of size of last spec
    start_ix(i_spec) = start_ix(i_spec-1) + spec_length

  end do
  start_ix = start_ix + 1 ! Fortran indexing starts at 1

  ! Stop indices are given by start indixes with offset
  stop_ix(:spec_to_fit%Count-1) = start_ix(2:) - 1
  stop_ix(spec_to_fit%Count) = this%N_b_total

end subroutine GenerateStartStopIxForFit

! Determine number of frequencies requested in the fit
subroutine DetermineNumberofFrequencies(spec_to_fit, N_freq_req)
  implicit none
  Type(TStringList), intent(in) :: spec_to_fit ! List of spectra to fit
  integer, intent(out) :: N_freq_req
  integer i_spec, ix, ix_1, ix_2, used90, used150, used220
  character(LEN=:), allocatable :: current_field, current_freq_1, current_freq_2

  ! Start count at zero
  N_freq_req = 0
  used90 = 0
  used150 = 0
  used220 = 0

  ! Loop over all the spectra requested and determine their type
  do i_spec=1, spec_to_fit%Count
    ! Identify spectrum
    call IdSpec(spectra_to_fit_list%Item(i_spec), current_field, current_freq_1, current_freq_2)

    if ((current_freq_1 .eq. "90") .or. (current_freq_2 .eq. "90")) then
      used90 = 1
    end if

    if ((current_freq_1 .eq. "150") .or. (current_freq_2 .eq. "150")) then
      used150 = 1
    end if

    if ((current_freq_1 .eq. "220") .or. (current_freq_2 .eq. "220")) then
      used220 = 1
    end if

  end do

  N_freq_req = used90 + used150 + used220

end subroutine DetermineNumberofFrequencies

! Determine how many TT, TE, and EE spectra were requested in the fit
subroutine DetermineNumberOfSpectraByField(spec_to_fit, NTT, NTE, NEE)
  implicit none
  Type(TStringList), intent(in) :: spec_to_fit ! List of spectra to fit
  integer, intent(out) :: NTT, NTE, NEE
  integer i_spec, ix, ix_1, ix_2
  character(LEN=:), allocatable :: current_field, current_freq_1, current_freq_2

  ! Start count at zero
  NTT = 0
  NTE = 0
  NEE = 0

  ! Loop over all the spectra requested and determine their type
  do i_spec=1, spec_to_fit%Count

    ! Identify spectrum
    call IdSpec(spectra_to_fit_list%Item(i_spec), current_field, current_freq_1, current_freq_2)

    ! Count spectra
    if (current_field .eq. "TT") then
      NTT = NTT + 1
    else if (current_field .eq. "TE") then
      NTE = NTE + 1
    else if (current_field .eq. "EE") then
      NEE = NEE + 1
    end if

  end do

end subroutine DetermineNumberOfSpectraByField

! Get the indices of the calibration parameters for each spectrum
subroutine GetCalibrationParamIx(spec_to_fit, spectra_to_fit_cal_indices)
  implicit none
  Type(TStringList), intent(in) :: spec_to_fit ! List of spectra to fit
  integer, allocatable :: spectra_to_fit_cal_indices(:,:)
  integer i_spec, ix, ix_1, ix_2
  character(LEN=:), allocatable :: current_field, current_freq_1, current_freq_2

  allocate(spectra_to_fit_cal_indices(spec_to_fit%Count, 4))
  spectra_to_fit_cal_indices = 0

  ! Loop over all the spectra requested and determine their type
  do i_spec=1, spec_to_fit%Count
    ! Identify spectrum
    call IdSpec(spectra_to_fit_list%Item(i_spec), current_field, current_freq_1, current_freq_2)

    ! Get the corresponding indices of the calibration parameters
    ! Do the first map
    if (current_freq_1 .eq. "90") then
      if (current_field .eq. "TT") then
        spectra_to_fit_cal_indices(i_spec,1) = iTcal90
        spectra_to_fit_cal_indices(i_spec,3) = iTcal90
      else if (current_field .eq. "TE") then
        spectra_to_fit_cal_indices(i_spec,1) = iTcal90
        spectra_to_fit_cal_indices(i_spec,4) = iEcal90
      else if (current_field .eq. "EE") then
        spectra_to_fit_cal_indices(i_spec,1) = iEcal90
        spectra_to_fit_cal_indices(i_spec,3) = iEcal90
      end if
    end if
    if (current_freq_1 .eq. "150") then
      if (current_field .eq. "TT") then
        spectra_to_fit_cal_indices(i_spec,1) = iTcal150
        spectra_to_fit_cal_indices(i_spec,3) = iTcal150
      else if (current_field .eq. "TE") then
        spectra_to_fit_cal_indices(i_spec,1) = iTcal150
        spectra_to_fit_cal_indices(i_spec,4) = iEcal150
      else if (current_field .eq. "EE") then
        spectra_to_fit_cal_indices(i_spec,1) = iEcal150
        spectra_to_fit_cal_indices(i_spec,3) = iEcal150
      end if
    end if
    if (current_freq_1 .eq. "220") then
      if (current_field .eq. "TT") then
        spectra_to_fit_cal_indices(i_spec,1) = iTcal220
        spectra_to_fit_cal_indices(i_spec,3) = iTcal220
      else if (current_field .eq. "TE") then
        spectra_to_fit_cal_indices(i_spec,1) = iTcal220
        spectra_to_fit_cal_indices(i_spec,4) = iEcal220
      else if (current_field .eq. "EE") then
        spectra_to_fit_cal_indices(i_spec,1) = iEcal220
        spectra_to_fit_cal_indices(i_spec,3) = iEcal220
      end if
    end if

    ! Now do the second map
    if (current_freq_2 .eq. "90") then
      if (current_field .eq. "TT") then
        spectra_to_fit_cal_indices(i_spec,2) = iTcal90
        spectra_to_fit_cal_indices(i_spec,4) = iTcal90
      else if (current_field .eq. "TE") then
        spectra_to_fit_cal_indices(i_spec,2) = iEcal90
        spectra_to_fit_cal_indices(i_spec,3) = iTcal90
      else if (current_field .eq. "EE") then
        spectra_to_fit_cal_indices(i_spec,2) = iEcal90
        spectra_to_fit_cal_indices(i_spec,4) = iEcal90
      end if
    end if
    if (current_freq_2 .eq. "150") then
      if (current_field .eq. "TT") then
        spectra_to_fit_cal_indices(i_spec,2) = iTcal150
        spectra_to_fit_cal_indices(i_spec,4) = iTcal150
      else if (current_field .eq. "TE") then
        spectra_to_fit_cal_indices(i_spec,2) = iEcal150
        spectra_to_fit_cal_indices(i_spec,3) = iTcal150
      else if (current_field .eq. "EE") then
        spectra_to_fit_cal_indices(i_spec,2) = iEcal150
        spectra_to_fit_cal_indices(i_spec,4) = iEcal150
      end if
    end if
    if (current_freq_2 .eq. "220") then
      if (current_field .eq. "TT") then
        spectra_to_fit_cal_indices(i_spec,2) = iTcal220
        spectra_to_fit_cal_indices(i_spec,4) = iTcal220
      else if (current_field .eq. "TE") then
        spectra_to_fit_cal_indices(i_spec,2) = iEcal220
        spectra_to_fit_cal_indices(i_spec,3) = iTcal220
      else if (current_field .eq. "EE") then
        spectra_to_fit_cal_indices(i_spec,2) = iEcal220
        spectra_to_fit_cal_indices(i_spec,4) = iEcal220
      end if
    end if

  end do

end subroutine GetCalibrationParamIx

! Get the Poisson power nuisance parameter indices for each spectrum
! TE spectra get 0
subroutine GetPoissonParamIx(spec_to_fit, spectra_to_fit_pow_indices)
  implicit none
  Type(TStringList), intent(in) :: spec_to_fit ! List of spectra to fit
  integer, allocatable :: spectra_to_fit_pow_indices(:)
  integer :: Poisson_ix_arr(12) ! helper that holds the Poisson ix nuisance params TT first, then EE
  integer i, i_spec, ix, ix_1, ix_2, EE_mod, TE_mod
  character(LEN=:), allocatable :: current_field, current_freq_1, current_freq_2

  ! Assemble helper
  Poisson_ix_arr = (/iTT_Poisson_90x90,iTT_Poisson_90x150,iTT_Poisson_90x220,iTT_Poisson_150x150,iTT_Poisson_150x220,iTT_Poisson_220x220,iEE_Poisson_90x90,iEE_Poisson_90x150,iEE_Poisson_90x220,iEE_Poisson_150x150,iEE_Poisson_150x220,iEE_Poisson_220x220/)

  allocate(spectra_to_fit_pow_indices(spec_to_fit%Count))
  spectra_to_fit_pow_indices = 0

  ! Loop over all the spectra requested and determine their type
  do i_spec=1, spec_to_fit%Count

    ! Identify spectrum
    call IdSpec(spectra_to_fit_list%Item(i_spec), current_field, current_freq_1, current_freq_2)

    ! No Poisson power in TE
    TE_mod = 1
    if (current_field .eq. "TE") then
      TE_mod = 0
    end if

    ! EE indices are offset by 6 wrt TT
    EE_mod = 0
    if (current_field .eq. "EE") then
      EE_mod = 6
    end if

    if ((current_freq_1 .eq. "90") .and. (current_freq_2 .eq. "90")) then
      i = 1
    else if ((current_freq_1 .eq. "90") .and. (current_freq_2 .eq. "150")) then
      i = 2
    else if ((current_freq_1 .eq. "90") .and. (current_freq_2 .eq. "220")) then
      i = 3
    else if ((current_freq_1 .eq. "150") .and. (current_freq_2 .eq. "150")) then
      i = 4
    else if ((current_freq_1 .eq. "150") .and. (current_freq_2 .eq. "220")) then
      i = 5
    else if ((current_freq_1 .eq. "220") .and. (current_freq_2 .eq. "220")) then
      i = 6
    end if

    spectra_to_fit_pow_indices(i_spec) = Poisson_ix_arr((i+EE_mod))*TE_mod

  end do

end subroutine GetPoissonParamIx

! Help to crop the caibration covariance matrix
! Hardcoded to format of the calibration covariance (T90, T150, T220, E90, E150, E220) - sorry!
! There's probably a more elegenat way to code this, but this only gets called once
subroutine MaskCalibrationCovarianceMatrix(spec_to_fit, cal_row_mask)
  implicit none
  Type(TStringList), intent(in) :: spec_to_fit ! List of spectra to fit
  integer, dimension(6), intent(out) :: cal_row_mask
  integer i_spec, ix, ix_1, ix_2
  character(LEN=:), allocatable :: current_field, current_freq_1, current_freq_2

  cal_row_mask = 0

  ! Loop over all the spectra requested and determine their type
  do i_spec=1, spec_to_fit%Count
    ! Identify spectrum
    call IdSpec(spectra_to_fit_list%Item(i_spec), current_field, current_freq_1, current_freq_2)

    ! Get the corresponding indices in the cal covariance
    if ((current_freq_1 .eq. "90") .or. (current_freq_2 .eq. "90")) then
      if (current_field .eq. "TT") then
        cal_row_mask(1) = 1
      else if (current_field .eq. "TE") then
        cal_row_mask(1) = 1
        cal_row_mask(4) = 1
      else if (current_field .eq. "EE") then
        cal_row_mask(4) = 1
      end if
    end if
    if ((current_freq_1 .eq. "150") .or. (current_freq_2 .eq. "150")) then
      if (current_field .eq. "TT") then
        cal_row_mask(2) = 1
      else if (current_field .eq. "TE") then
        cal_row_mask(2) = 1
        cal_row_mask(5) = 1
      else if (current_field .eq. "EE") then
        cal_row_mask(5) = 1
      end if
    end if
    if ((current_freq_1 .eq. "220") .or. (current_freq_2 .eq. "220")) then
      if (current_field .eq. "TT") then
        cal_row_mask(3) = 1
      else if (current_field .eq. "TE") then
        cal_row_mask(3) = 1
        cal_row_mask(6) = 1
      else if (current_field .eq. "EE") then
        cal_row_mask(6) = 1
      end if
    end if

  end do

end subroutine MaskCalibrationCovarianceMatrix

! Crops and Inverts, like advertised, to be used together with masking function
subroutine CropAndInvertCalibrationCovarianceMatrix(full_cal_covariance_matrix, cal_rows_use, inv_cal_covariance)
  implicit none
  integer, dimension(6), intent(in) :: cal_rows_use
  real(mcp), intent(in) :: full_cal_covariance_matrix(:,:)
  real(mcp), intent(out) :: inv_cal_covariance(:,:)
  integer i, j, ii, jj,n,info

  ii = 1
  do i=1, N_freq_0*2 ! Loop over full matrix and copy desired elements into smaller matrix
    if (IntToLogical(cal_rows_use(i))) then
      jj = 1
      do j=1, N_freq_0*2
        if (IntToLogical(cal_rows_use(j))) then
          inv_cal_covariance(ii,jj) = full_cal_covariance_matrix(i,j)
          jj = jj + 1
        end if
      end do
      ii = ii + 1
    end if
  end do
#ifndef _STANDALONE_
  call Matrix_Inverse(inv_cal_covariance) ! Invert (only needs to be done once at the start)
#else
  n = ii-1
  call dpotrf ('L', n, inv_cal_covariance, n, info)
  if (info/=0) then
    write(*,*) "calibration covariance not positive definite"
    stop
  endif
  call dpotri ('L', n, inv_cal_covariance, n, info)
  if (info/=0) then
    write(*,*) "calibration covariance invert failed"
    stop
  endif
  ! symmetrize
  do i=1,n
        do j=1,i-1
            inv_cal_covariance(j,i) = inv_cal_covariance(i,j)
        end do
    end do
#endif

end subroutine CropAndInvertCalibrationCovarianceMatrix

! Extract the effective band centres from the file
! Relies on specific ordering (90, 150, 220) (pol gal dust, DSFG, raio galaxies, tSZ)
subroutine GrabEffectiveBandCentres(spec_to_fit, nu_matrix, nu_cirrus, nu_pol_gal, nu_DSFG, nu_tSZ)
  implicit none
  Type(TStringList), intent(in) :: spec_to_fit ! List of spectra to fit
  real(mcp), intent(in) :: nu_matrix(:,:)
  real(mcp), intent(out) :: nu_cirrus(:,:), nu_pol_gal(:,:), nu_DSFG(:,:), nu_tSZ(:,:)
  integer i_spec, ix, ix_1, ix_2
  character(LEN=:), allocatable :: current_field, current_freq_1, current_freq_2

  ! Loop over all the spectra requested and determine their type
  do i_spec=1, spec_to_fit%Count

    ! Identify spectrum
    call IdSpec(spectra_to_fit_list%Item(i_spec), current_field, current_freq_1, current_freq_2)

    ! Deal with the first frequency
    if (current_freq_1 .eq. "90") then
      nu_cirrus(i_spec,1) = nu_matrix(1,1)
      nu_pol_gal(i_spec,1) = nu_matrix(2,1)
      nu_DSFG(i_spec,1) = nu_matrix(3,1)
      nu_tSZ(i_spec,1) = nu_matrix(5,1)
    else if (current_freq_1 .eq. "150") then
      nu_cirrus(i_spec,1) = nu_matrix(1,2)
      nu_pol_gal(i_spec,1) = nu_matrix(2,2)
      nu_DSFG(i_spec,1) = nu_matrix(3,2)
      nu_tSZ(i_spec,1) = nu_matrix(5,2)
    else if (current_freq_1 .eq. "220") then
      nu_cirrus(i_spec,1) = nu_matrix(1,3)
      nu_pol_gal(i_spec,1) = nu_matrix(2,3)
      nu_DSFG(i_spec,1) = nu_matrix(3,3)
      nu_tSZ(i_spec,1) = nu_matrix(5,3)
    end if

    ! Now the second frequency
    if (current_freq_2 .eq. "90") then
      nu_cirrus(i_spec,2) = nu_matrix(1,1)
      nu_pol_gal(i_spec,2) = nu_matrix(2,1)
      nu_DSFG(i_spec,2) = nu_matrix(3,1)
      nu_tSZ(i_spec,2) = nu_matrix(5,1)
    else if (current_freq_2 .eq. "150") then
      nu_cirrus(i_spec,2) = nu_matrix(1,2)
      nu_pol_gal(i_spec,2) = nu_matrix(2,2)
      nu_DSFG(i_spec,2) = nu_matrix(3,2)
      nu_tSZ(i_spec,2) = nu_matrix(5,2)
    else if (current_freq_2 .eq. "220") then
      nu_cirrus(i_spec,2) = nu_matrix(1,3)
      nu_pol_gal(i_spec,2) = nu_matrix(2,3)
      nu_DSFG(i_spec,2) = nu_matrix(3,3)
      nu_tSZ(i_spec,2) = nu_matrix(5,3)
    end if

  end do

end subroutine GrabEffectiveBandCentres

! Grab indices of CIB clustering decorrelation parameters
subroutine GrabCIBClusteringDecorrelationParamIx(spec_to_fit, cib_cl_decorr_ix)
  implicit none
  Type(TStringList), intent(in) :: spec_to_fit ! List of spectra to fit
  integer, allocatable :: cib_cl_decorr_ix(:,:)
  integer i_spec, ix, ix_1, ix_2
  character(LEN=:), allocatable :: current_field, current_freq_1, current_freq_2

  allocate(cib_cl_decorr_ix(spec_to_fit%Count,2))

  ! Loop over all the spectra requested and determine their type
  do i_spec=1, spec_to_fit%Count

    ! Identify spectrum
    call IdSpec(spectra_to_fit_list%Item(i_spec), current_field, current_freq_1, current_freq_2)

    ! Deal with the first frequency
    if (current_freq_1 .eq. "90") then
      cib_cl_decorr_ix(i_spec,1) = iTT_CIB_Clustering_decorr_90
    else if (current_freq_1 .eq. "150") then
      cib_cl_decorr_ix(i_spec,1) = iTT_CIB_Clustering_decorr_150
    else if (current_freq_1 .eq. "220") then
      cib_cl_decorr_ix(i_spec,1) = iTT_CIB_Clustering_decorr_220
    end if

    ! Now the second frequency
    if (current_freq_2 .eq. "90") then
      cib_cl_decorr_ix(i_spec,2) = iTT_CIB_Clustering_decorr_90
    else if (current_freq_2 .eq. "150") then
      cib_cl_decorr_ix(i_spec,2) = iTT_CIB_Clustering_decorr_150
    else if (current_freq_2 .eq. "220") then
      cib_cl_decorr_ix(i_spec,2) = iTT_CIB_Clustering_decorr_220
    end if

  end do

end subroutine GrabCIBClusteringDecorrelationParamIx

! Helper that returns the map types and frequencies in a certain spectrum
subroutine IdSpec(spec_to_id, current_field, current_freq_1, current_freq_2)
  implicit none
  character(*), intent(in) :: spec_to_id
  character(LEN=:), allocatable, intent(out) :: current_field, current_freq_1, current_freq_2
  character(LEN=:), allocatable :: current_map_1_str, current_map_2_str
  integer ix, ix_1, ix_2

  ! Identify spectrum
  ix = scan(spec_to_id, "x")
  current_map_1_str = spec_to_id(1:ix-1)
  current_map_2_str = spec_to_id(ix+1:)
  ix_1 = scan(current_map_1_str, "_")
  ix_2 = scan(current_map_2_str, "_")
  current_field = trim(current_map_1_str(ix_1+1:)) // trim(current_map_2_str(ix_2+1:))
  current_freq_1 = current_map_1_str(1:ix_1-1)
  current_freq_2 = current_map_2_str(1:ix_2-1)

end subroutine IdSpec

#ifndef _STANDALONE_
! Write out smooth model spectrum with all the foreground components
subroutine WriteOutTheorySpectrum(spec_id, write_channel, current_Dl_theory_unbinned_CMB_only, current_Dl_theory_unbinned, Dl_foregrounds)
  implicit none
  character(*), intent(in) :: spec_id
  integer, intent(in) :: write_channel
  real(mcp), dimension(SPT3G_windows_lmax+1), intent(in) :: current_Dl_theory_unbinned_CMB_only
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax), intent(in) :: current_Dl_theory_unbinned
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax), intent(in) :: Dl_foregrounds
  real*4, dimension(N_fg_max+3) :: save_arr
  integer l

  save_arr=0
  call OpenWriteBinaryFile("likelihood_tests/unbinned_theory_"//trim(spec_id)//".dat",write_channel,12_8)
  do l=SPT3G_windows_lmin,SPT3G_windows_lmax
    save_arr(1)=l
    save_arr(2)=current_Dl_theory_unbinned_CMB_only(l)
    save_arr(3)=current_Dl_theory_unbinned(l)
    save_arr(4)=Dl_foregrounds(iOutCal,l)
    save_arr(5)=Dl_foregrounds(iOutAberration,l)
    save_arr(6)=Dl_foregrounds(iOutSSL,l)
    save_arr(7)=Dl_foregrounds(iOutGalDust,l)
    save_arr(8)=Dl_foregrounds(iOutPoisson,l)
    save_arr(9)=Dl_foregrounds(iOutCIBClustering,l)
    save_arr(10)=Dl_foregrounds(iOuttSZ,l)
    save_arr(11)=Dl_foregrounds(iOuttSZCIB,l)
    save_arr(12)=Dl_foregrounds(iOutkSZ,l)
    write(write_channel,rec=l) save_arr(1:12)
  enddo
  close(write_channel)

end subroutine WriteOutTheorySpectrum



! Used for writing out spectra
subroutine OpenWriteBinaryFile(aname,aunit,record_length)
  character(LEN=*), intent(IN) :: aname
  integer, intent(in) :: aunit
  integer*8,intent(in) :: record_length
  open(unit=aunit,file=aname,form='unformatted',status='replace',access='direct',recl=record_length, err=500)
  return
500 call MpiStop('File not able to be written to: '//trim(aname))
end subroutine OpenWriteBinaryFile
#endif

! Convert int to boolean
function IntToLogical(int_var) result(log_var)
  integer, intent(in) :: int_var
  logical :: log_var
  log_var = .false.
  if (int_var .eq. 1) then
    log_var = .true.
  end if
end function IntToLogical
subroutine str2int(str,int,stat)
  character(LEN=*),intent(IN) :: str
  integer,intent(out) :: int
  integer,intent(out) :: stat
  !print *,"inst"
  !print *,str
  read(str,*,iostat=stat)  int
end subroutine str2int

end module CMB_SPT3G_2018_TTTEEE
