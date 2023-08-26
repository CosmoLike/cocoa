! SPT-3G 2018 TTTEEE likelihood
! Written by Lennart Balkenhol 2022
! Foregrounds module

!---------------------------------------------------------!
!          MODULE DEFINITION, IMPORTS, CONSTANTS          !
!---------------------------------------------------------!

! Imports
module CMB_SPT3G_2018_TTTEEE_foregrounds
#ifndef _STANDALONE_
use settings
use FileUtils
use CosmologyTypes
#else
use SPT3G_utils
#endif
implicit none

private

public :: ApplyCalibration,&
          AddPoissonPower,&
          AddGalacticDust,&
          AddCIBClustering,&
          AddtSZ,&
          AddtSZCIBCorrelation,&
          AddkSZ,&
          ApplySuperSampleLensing,&
          ApplyAberrationCorrection, &
#ifdef _STANDALONE_
          SPT3G_2018_TTTEEE_Ini_Foregrounds, tSZCosmologyScalingEnabled, kSZCosmologyScalingEnabled 
#else
          SPT3G_2018_TTTEEE_ReadIni_Foregrounds
#endif

! Physical constants
real(mcp), parameter :: T_CMB = 2.72548_mcp    ! CMB temperature
real(mcp), parameter :: h = 6.62606957e-34_mcp ! Planck's constant
real(mcp), parameter :: kB = 1.3806488e-23_mcp ! Boltzmann constant
real(mcp), parameter :: Ghz_Kelvin = h/kB*1e9_mcp

! Indices specifying order of fg components in output array
integer, parameter :: N_fg_max = 9
integer, parameter :: iOutCal = 1, iOutAberration = 2, iOutSSL = 3, iOutGalDust = 4, iOutPoisson = 5
integer, parameter :: iOutCIBClustering = 6, iOuttSZ = 7, iOuttSZCIB = 8, iOutkSZ = 9 ! TT exclusive foregrounds

! Default ell range matching window files, but can be adjusted
integer :: SPT3G_windows_lmin = 1, SPT3G_windows_lmax = 3200

! Cosmology scaling of foregrounds
logical :: tSZCosmologyScalingEnabled
logical :: kSZCosmologyScalingEnabled

! Reference parameters for foreground model
real(mcp) :: nu_0_galdust
real(mcp) :: T_galdust
real(mcp) :: nu_0_CIB
real(mcp) :: T_CIB
real(mcp) :: nu_0_tSZ

! Foreground templates
real(mcp), allocatable :: tSZ_template(:) ! Normalised tSZ template
real(mcp), allocatable :: kSZ_template(:) ! Normalised kSZ template

contains

!---------------------------------------------------------!
!                  INITIALISATION METHOD                  !
!---------------------------------------------------------!

! Initialise the foreground data
#ifndef _STANDALONE_
subroutine SPT3G_2018_TTTEEE_ReadIni_Foregrounds(Ini)
  class(TSettingIni) :: Ini

  character(LEN=:), allocatable :: tSZ_template_name
  real(mcp), allocatable :: full_tSZ_template(:,:)

  character(LEN=:), allocatable :: kSZ_template_name
  real(mcp), allocatable :: full_kSZ_template(:,:)

  !KARIM: Fortran untils here that read from the ini file or read in a txt file via ReadTextMatrix

  ! Read in ell range
  SPT3G_windows_lmin = Ini%Read_Int("SPT3G_2018_TTTEEE_window_l_min")
  SPT3G_windows_lmax = Ini%Read_Int("SPT3G_2018_TTTEEE_window_l_max")

  ! Read in reference parameters for foregrounds
  nu_0_galdust = Ini%Read_Real("SPT3G_2018_TTTEEE_galdust_nu0", 150.0)
  T_galdust = Ini%Read_Real("SPT3G_2018_TTTEEE_galdust_T", 19.6)
  nu_0_CIB = Ini%Read_Real("SPT3G_2018_TTTEEE_CIB_nu0", 150.0)
  T_CIB = Ini%Read_Real("SPT3G_2018_TTTEEE_CIB_T", 25.0)
  nu_0_tSZ = Ini%Read_Real("SPT3G_2018_TTTEEE_tSZ_nu0", 143.0)

  ! Read in tSZ template and normalise
  ! Cosmology scaling is not supported!
  tSZ_template_name = Ini%ReadFileName("SPT3G_2018_TTTEEE_tSZ_template_file",relative = .true., NotFoundFail=.true.)
  tSZCosmologyScalingEnabled = Ini%Read_Logical("SPT3G_2018_TTTEEE_tSZ_cosmology_scaling", .false.)! not supported!

  allocate(full_tSZ_template(1+SPT3G_windows_lmax-SPT3G_windows_lmin,2)) ! File contains column of ell
  call File%ReadTextMatrix(tSZ_template_name, full_tSZ_template)
  allocate(tSZ_template(1+SPT3G_windows_lmax-SPT3G_windows_lmin))
  tSZ_template = full_tSZ_template(:,2) / full_tSZ_template(3000,2) ! Ensure normalisation

  ! Read in kSZ template and normalise
  ! Cosmology scaling is not supported!
  kSZ_template_name = Ini%ReadFileName("SPT3G_2018_TTTEEE_kSZ_template_file",relative = .true., NotFoundFail=.true.)
  kSZCosmologyScalingEnabled = Ini%Read_Logical("SPT3G_2018_TTTEEE_kSZ_cosmology_scaling", .false.)! not supported!

  allocate(full_kSZ_template(1+SPT3G_windows_lmax-SPT3G_windows_lmin,2)) ! File contains column of ell
  call File%ReadTextMatrix(kSZ_template_name, full_kSZ_template)
  allocate(kSZ_template(1+SPT3G_windows_lmax-SPT3G_windows_lmin))
  kSZ_template = full_kSZ_template(:,2) / full_kSZ_template(3000,2) ! Ensure normalisation

end subroutine SPT3G_2018_TTTEEE_ReadIni_Foregrounds
#else
subroutine SPT3G_2018_TTTEEE_Ini_Foregrounds(eSPT3G_windows_lmin,eSPT3G_windows_lmax, enu_0_galdust, eT_galdust, enu_0_CIB, eT_CIB, enu_0_tSZ, etSZCosmologyScalingEnabled,full_tSZ_template,ekSZCosmologyScalingEnabled,full_kSZ_template)
  real(mcp)::enu_0_galdust, eT_galdust, enu_0_CIB, eT_CIB, enu_0_tSZ
  integer :: eSPT3G_windows_lmin,eSPT3G_windows_lmax, etSZCosmologyScalingEnabled,ekSZCosmologyScalingEnabled

  real(mcp) :: full_tSZ_template(1+eSPT3G_windows_lmax-eSPT3G_windows_lmin)
  real(mcp) :: full_kSZ_template(1+eSPT3G_windows_lmax-eSPT3G_windows_lmin)


  ! Read in ell range
  SPT3G_windows_lmin = eSPT3G_windows_lmin
  SPT3G_windows_lmax = eSPT3G_windows_lmax

  ! Read in reference parameters for foregrounds
  nu_0_galdust = enu_0_galdust
  T_galdust    = eT_galdust
  nu_0_CIB     = enu_0_CIB
  T_CIB        = eT_CIB
  nu_0_tSZ     = enu_0_tSZ

  ! Read in tSZ template and normalise
  ! Cosmology scaling not supported
  tSZCosmologyScalingEnabled = .false.
  if (etSZCosmologyScalingEnabled .eq. 1) then
    tSZCosmologyScalingEnabled = .false.
  endif

  allocate(tSZ_template(1+SPT3G_windows_lmax-SPT3G_windows_lmin))
  tSZ_template = full_tSZ_template(:) / full_tSZ_template(3000) ! Ensure normalisation

  ! Read in kSZ template and normalise
  ! Cosmology scaling not supported
  kSZCosmologyScalingEnabled = .false.
  if (ekSZCosmologyScalingEnabled .eq. 1) then
    kSZCosmologyScalingEnabled = .false.
  endif


  allocate(kSZ_template(1+SPT3G_windows_lmax-SPT3G_windows_lmin))
  kSZ_template = full_kSZ_template(:) / full_kSZ_template(3000 ) ! Ensure normalisation

end subroutine SPT3G_2018_TTTEEE_Ini_Foregrounds
#endif

!---------------------------------------------------------!
!                  FOREGROUND FUNCTIONS                   !
!---------------------------------------------------------!

! Calibration
! Data is scaled as: TT: T1*T2, TE: 0.5*(T1*E2+T2*E1), EE: E1*E2
! Theory is scaled by the inverse
! In function this is calculated as  0.5*(cal1*cal2+cal3*cal4)
subroutine ApplyCalibration(cal1, cal2, cal3, cal4, Dl_theory, Dl_foregrounds)
  real(mcp), intent(in) :: cal1, cal2, cal3, cal4
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_theory
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_foregrounds
  real(mcp) calibration

  ! This is how the data spectra are calibrated
  calibration = 0.5*(cal1*cal2+cal3*cal4)

  ! So theory gets the inverse of this
  Dl_theory = Dl_theory/calibration
  Dl_foregrounds(iOutCal,:) = 1/calibration

end subroutine ApplyCalibration

! Add Poisson power, referenced at ell=3000
subroutine AddPoissonPower(pow_at_3000, Dl_theory, Dl_foregrounds)
  real(mcp), intent(in) :: pow_at_3000
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_theory
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_foregrounds
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_Poisson
  real(mcp), dimension(SPT3G_windows_lmax) :: ells
  integer i

  ! Grab ells helper (1-3200)
  do i=1,SPT3G_windows_lmax
   ells(i) = i
  end do

  ! Calculate and add Poisson power
  Dl_Poisson = ells*ells*pow_at_3000/(3000*3000)
  Dl_theory = Dl_theory + Dl_Poisson
  Dl_foregrounds(iOutPoisson,:) = Dl_Poisson

end subroutine AddPoissonPower


! Add galactic dust (intensity and polarisation)
! Referenced at ell=80, with power law dependence (alpha+2)
! At effective frequencies nu1 and nu2, with spectral index beta
subroutine AddGalacticDust(pow_at_80, alpha, beta, nu1, nu2, Dl_theory, Dl_foregrounds)
  real(mcp), intent(in) :: pow_at_80, alpha, beta, nu1, nu2
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_theory
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_foregrounds
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_galdust
  real(mcp), dimension(SPT3G_windows_lmax) :: ells
  integer i

  ! Grab ells helper (1-3200)
  do i=1,SPT3G_windows_lmax
   ells(i) = i
  end do

  ! Calculate and add galactic dust power
  Dl_galdust = pow_at_80 * (ells/80d0)**(alpha + 2.0d0)
  Dl_galdust = Dl_galdust * DustFreqScaling(beta, T_galdust, nu_0_galdust, nu1) * DustFreqScaling(beta, T_galdust, nu_0_galdust, nu2)
  Dl_theory = Dl_theory + Dl_galdust
  Dl_foregrounds(iOutGalDust,:) = Dl_galdust

end subroutine AddGalacticDust


! Add CIB clustering
! Referenced at ell=300, with power law dependence (alpha)
! At effective frequencies nu1 and nu2, with spectral index beta
! Decorrelation parameters zeta1 and zeta2
subroutine AddCIBClustering(pow_at_3000, alpha, beta, nu1, nu2, zeta1, zeta2, Dl_theory, Dl_foregrounds)
  real(mcp), intent(in) :: pow_at_3000, alpha, beta, nu1, nu2, zeta1, zeta2
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_theory
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_foregrounds
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_cib_clustering
  real(mcp), dimension(SPT3G_windows_lmax) :: ells
  integer i

  ! Grab ells helper (1-3200)
  do i=1,SPT3G_windows_lmax
   ells(i) = i
  end do

  ! Calculate and add polarised galactic dust power
  Dl_cib_clustering = pow_at_3000 * (ells/3000d0)**(alpha)
  Dl_cib_clustering = Dl_cib_clustering * DustFreqScaling(beta, T_CIB, nu_0_CIB, nu1) * DustFreqScaling(beta, T_CIB, nu_0_CIB, nu2)
  Dl_cib_clustering = Dl_cib_clustering * SQRT(zeta1*zeta2)
  Dl_theory = Dl_theory + Dl_cib_clustering
  Dl_foregrounds(iOutCIBClustering,:) = Dl_cib_clustering

end subroutine AddCIBClustering


! Add tSZ contribution
! Template normalised at ell=3000
! Cosmology scaling not supported
subroutine AddtSZ(pow_at_3000, nu1, nu2, H0, sigma_8, omb, Dl_theory, Dl_foregrounds)
  real(mcp), intent(in) :: pow_at_3000, nu1, nu2, H0, sigma_8, omb
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_theory
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_foregrounds
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_tSZ

  ! Calculate tSZ power
  Dl_tSZ = pow_at_3000 * tSZ_template ! Template
  Dl_tSZ = Dl_tSZ * tSZFrequencyScaling(nu1, nu_0_tSZ, T_CMB) * tSZFrequencyScaling(nu2, nu_0_tSZ, T_CMB) ! Frequency scaling

  ! Cosmology scaling
  !if (tSZCosmologyScalingEnabled) then
  !  Dl_tSZ = Dl_tSZ * tSZCosmologyScaling(H0, sigma_8, omb)
  !end if

  ! Add to model
  Dl_theory = Dl_theory + Dl_tSZ
  Dl_foregrounds(iOuttSZ,:) = Dl_tSZ

end subroutine AddtSZ


! Correlation between tSZ and CIB
! Only use clustered CIB component here and simplified model
! Sorry for the horrible call signature!
subroutine AddtSZCIBCorrelation(xi_tsz_CIB, tsz_pow_at_3000, CIB_pow_at_3000, alpha, beta, zeta1, zeta2, CIB_nu1, CIB_nu2, tSZ_nu1, tSZ_nu2, H0, sigma_8, omb, Dl_theory, Dl_foregrounds)
  real(mcp), intent(in) :: xi_tsz_CIB, tsz_pow_at_3000, CIB_pow_at_3000, alpha, beta, zeta1, zeta2, CIB_nu1, CIB_nu2, tSZ_nu1, tSZ_nu2, H0, sigma_8, omb
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_theory
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_foregrounds
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_foregrounds_dummy ! Dummy array to simplify calls
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_tSZ_11, Dl_tSZ_22 ! Different frequency combinations
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_cib_clustering_11, Dl_cib_clustering_22 ! Different frequency combinations
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_tSZ_CIB_corr

  ! Start at zero
  Dl_cib_clustering_11 = 0
  Dl_cib_clustering_22 = 0
  Dl_tSZ_11 = 0
  Dl_tSZ_22 = 0

  ! Calculate CIB components
  call AddCIBClustering(CIB_pow_at_3000, alpha, beta, CIB_nu1, CIB_nu1, zeta1, zeta1, Dl_cib_clustering_11, Dl_foregrounds_dummy)
  call AddCIBClustering(CIB_pow_at_3000, alpha, beta, CIB_nu2, CIB_nu2, zeta2, zeta2, Dl_cib_clustering_22, Dl_foregrounds_dummy)

  ! Calculate the tSZ components
  call AddtSZ(tsz_pow_at_3000, tSZ_nu1, tSZ_nu1, H0, sigma_8, omb, Dl_tSZ_11, Dl_foregrounds_dummy)
  call AddtSZ(tsz_pow_at_3000, tSZ_nu2, tSZ_nu2, H0, sigma_8, omb, Dl_tSZ_22, Dl_foregrounds_dummy)

  ! Calculate tSZ-CIB correlation
  ! Sign defined such that a positive xi corresponds to a reduction at 150GHz
  Dl_tSZ_CIB_corr = -1.0_mcp * xi_tsz_CIB * (SQRT(Dl_tSZ_11*Dl_cib_clustering_22) + SQRT(Dl_tSZ_22*Dl_cib_clustering_11))
  Dl_theory = Dl_theory + Dl_tSZ_CIB_corr
  Dl_foregrounds(iOuttSZCIB,:) = Dl_tSZ_CIB_corr

end subroutine AddtSZCIBCorrelation


! Add kSZ contribution
! Template normalised at ell=3000
! Cosmology scaling not supported
subroutine AddkSZ(pow_at_3000, H0, sigma_8, omb, omegam, ns, tau, Dl_theory, Dl_foregrounds)
  real(mcp), intent(in) :: pow_at_3000, H0, sigma_8, omb, omegam, ns, tau
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_theory
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_foregrounds
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: Dl_kSZ

  ! Calculate kSZ power
  Dl_kSZ = pow_at_3000 * kSZ_template ! Template

  ! Cosmology scaling
  !if (kSZCosmologyScalingEnabled) then
  !  Dl_kSZ = Dl_kSZ * kSZCosmologyScaling(H0, sigma_8, omb, omegam, ns, tau)
  !end if

  ! Add to model
  Dl_theory = Dl_theory + Dl_kSZ
  Dl_foregrounds(iOutkSZ,:) = Dl_kSZ

end subroutine AddkSZ

! Super sample lensing
! Based on Manzotti et al. 2014 (https://arxiv.org/pdf/1401.7992.pdf) Eq. 32
! Applies correction to the spectrum and returns the correction slotted into the fg array
subroutine ApplySuperSampleLensing(kappa, Dl_theory, Dl_foregrounds)
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: Cl_derivative
  real(mcp), intent(in) :: kappa
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_theory
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_foregrounds
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: ssl_correction
  real(mcp), dimension(SPT3G_windows_lmax) :: ells
  integer i

  ! Grab ells helper (1-3200)
  do i=1,SPT3G_windows_lmax
   ells(i) = i
  end do

  ! Grab Cl derivative
  call GetClDerivative(Dl_theory, Cl_derivative)

  ! Calculate super sample lensing correction
  ! (In Cl space) SSL = -k/l^2 d/dln(l) (l^2Cl) = -k(l*dCl/dl + 2Cl)
  ssl_correction = ells*Cl_derivative ! l*dCl/dl
  ssl_correction = ssl_correction * ells * (ells+1)/(2*PI) ! Convert this part to Dl space already
  ssl_correction = ssl_correction + 2*Dl_theory ! 2Cl - but already converted to Dl
  ssl_correction = -1 * ssl_correction * kappa ! -kappa

  ! Apply the correction
  Dl_theory = Dl_theory + ssl_correction
  Dl_foregrounds(iOutSSL,:) = ssl_correction

end subroutine ApplySuperSampleLensing

! Aberration Correction
! Based on Jeong et al. 2013 (https://arxiv.org/pdf/1309.2285.pdf) Eq. 23
! Applies correction to the spectrum and returns the correction by itself
subroutine ApplyAberrationCorrection(ab_coeff, Dl_theory, Dl_foregrounds)
  real(mcp), intent(in) :: ab_coeff
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: Cl_derivative
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_theory
  real(mcp), dimension(N_fg_max,SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Dl_foregrounds
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: aberration_correction
  real(mcp), dimension(SPT3G_windows_lmax) :: ells
  integer i

  ! AC = beta*l(l+1)dCl/dln(l)/(2pi)

  ! Grab ells helper (1-3200)
  do i=1,SPT3G_windows_lmax
   ells(i) = i
  end do

  ! Grab Cl derivative
  call GetClDerivative(Dl_theory, Cl_derivative)

  ! Calculate aberration correction
  ! (In Cl space) AC = -coeff*dCl/dln(l) = -coeff*l*dCl/dl
  ! where coeff contains the boost amplitude and direction (beta*<cos(theta)> in Jeong+ 13)
  aberration_correction = -1*ab_coeff*Cl_derivative*ells
  aberration_correction = aberration_correction*ells*(ells+1)/(2*PI) ! Convert to Dl

  ! Apply correction
  Dl_theory = Dl_theory + aberration_correction
  Dl_foregrounds(iOutAberration,:) = aberration_correction

end subroutine ApplyAberrationCorrection

! Helper to get the derivative of the spectrum
! Takes Dl in, but returns Cl derivative!
! Handles end points approximately
! Smoothes any spike at the ell_max
subroutine GetClDerivative(Dl_theory, Cl_derivative)
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax), intent(in) :: Dl_theory
  real(mcp), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax), intent(out) :: Cl_derivative
  real(mcp), dimension(SPT3G_windows_lmax) :: ells
  integer i

  ! Grab ells helper (1-3200)
  do i=1,SPT3G_windows_lmax
   ells(i) = i
  end do

  ! Calculate derivative
  Cl_derivative = Dl_theory*2*PI/(ells*(ells+1)) ! Convert to Cl
  Cl_derivative(2:SPT3G_windows_lmax-1) = 0.5*( Cl_derivative(3:SPT3G_windows_lmax) - Cl_derivative(1:SPT3G_windows_lmax-2) ) ! Find gradient
  Cl_derivative(1) = Cl_derivative(2) ! Handle start approximately
  Cl_derivative(SPT3G_windows_lmax) = Cl_derivative(SPT3G_windows_lmax-1) ! Handle end approximately

#ifndef _STANDALONE_
  ! Smooth over spike at lmax
  ! Transition point between Boltzmann solver Cl and where the spectrum comes from a lookup table/interpolation can cause a spike in derivative
  if (CosmoSettings%lmax_computed_cl .LT. SPT3G_windows_lmax-1) then
    Cl_derivative(CosmoSettings%lmax_computed_cl) = 0.75*Cl_derivative(CosmoSettings%lmax_computed_cl-1)&
      + 0.25*Cl_derivative(CosmoSettings%lmax_computed_cl+2)
    Cl_derivative(CosmoSettings%lmax_computed_cl+1) = 0.75*Cl_derivative(CosmoSettings%lmax_computed_cl+2)&
      + 0.25*Cl_derivative(CosmoSettings%lmax_computed_cl-1)
  end if
#endif

end subroutine GetClDerivative

!---------------------------------------------------------!
!                     SCALING HELPERS                     !
!---------------------------------------------------------!

! Galactic Dust Frequency Scaling
function DustFreqScaling(beta,Tdust,nu0,nu_eff) result(fdust)
  real(mcp), intent(in) :: beta
  real(mcp), intent(in) :: Tdust
  real(mcp), intent(in) :: nu_eff
  real(mcp), intent(in) :: nu0 ! Pivot frequency
  real(mcp) :: fdust

  fdust = (nu_eff/nu0)**beta
  fdust = fdust*Bnu(nu_eff,nu0,Tdust)/dBdT(nu_eff,nu0,T_CMB)

end function DustFreqScaling

! Planck function normalised to 1 at nu0
function Bnu(nu,nu0,T)
  real(mcp) Bnu,nu,nu0,T

  Bnu = (nu/nu0)**3
  Bnu = Bnu * (exp( Ghz_Kelvin*nu0/T)-1d0) / (exp( Ghz_Kelvin*nu/T)-1d0)

end function Bnu

! Derivative of Planck function normalised to 1 at nu0
function dBdT(nu,nu0,T)
  real(mcp) dBdT,dBdT0,nu,nu0,T,x,x0

  x0 = Ghz_Kelvin*nu0/T
  x = Ghz_Kelvin*nu/T

  dBdT0 = x0**4 * exp(x0) / (exp(x0)-1)**2
  dBdT =  x**4 * exp(x) / (exp(x)-1)**2

  dBdT = dBdT/dBdT0

end function dBdT

! tSZ Frequency Scaling
! Gives conversion factor for frequency nu from reference nu0
function tSZFrequencyScaling(nu, nu0, T) result(tSZfac)
  real(mcp) tSZfac,tSZfac0,nu,nu0,T,x,x0

  x0 = Ghz_Kelvin*nu0/T
  x = Ghz_Kelvin*nu/T

  tSZfac0 = x0*(exp(x0)+1)/(exp(x0)-1) - 4
  tSZfac = x*(exp(x)+1)/(exp(x)-1) - 4

  tSZfac = tSZfac/tSZfac0

end function tSZFrequencyScaling

! Taken from Reichardt et al. 2020 likelihood
! NOT SUPPORTED - DO NOT USE
function tSZCosmologyScaling(H0, sigma8, omegab) result(tSZfac)
  double precision :: H0, sigma8, omegab
  double precision :: tSZfac

  tSZfac = ((H0/71.0)**1.73 ) * ( (sigma8/.8)**8.34 ) * ( (omegab/.044)**2.81 )

end function tSZCosmologyScaling

! Taken from Reichardt et al. 2020 likelihood
! NOT SUPPORTED - DO NOT USE
function kSZCosmologyScaling(H0, sigma8, omegab, omegam, ns, tau) result(kSZfac)
  double precision :: H0, sigma8, omegab, omegam, ns, tau
  double precision :: kSZfac

  kSZfac = ((H0/71.0)**1.7 ) * ( (sigma8/.8)**4.7 ) * ( (omegab/.044)**2.1 ) * ( (omegam/.264)**(-0.44) ) * ( (ns/.96)**(-0.19) )

end function kSZCosmologyScaling

end module CMB_SPT3G_2018_TTTEEE_foregrounds
