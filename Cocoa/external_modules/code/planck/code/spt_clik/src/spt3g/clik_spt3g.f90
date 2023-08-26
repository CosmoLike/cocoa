module CMB_SPT3G_EETE_2020_clik

implicit None

real(8), parameter :: T_CMB = 2.72548    ! CMB temperature
real(8), parameter :: h = 6.62606957e-34 ! Planck's constant
real(8), parameter :: kB = 1.3806488e-23 ! Boltzmann constant
real(8), parameter :: Ghz_Kelvin = h/kB*1e9
real(8),                                parameter     :: twoPI = 6.283185307179586476925286766559005768394

integer :: iKappa=1,&
  iPoisson_90x90=2, iPoisson_90x150=3, iPoisson_90x220=4, iPoisson_150x150=5, iPoisson_150x220=6, iPoisson_220x220=7,&
  iTDust=8, iADust_TE_150=9, iBetaDust_TE=10, iAlphaDust_TE=11, iADust_EE_150=12,&
  iBetaDust_EE=13, iAlphaDust_EE=14,&
  iMapTcal90=15, iMapTcal150=16, iMapTcal220=17,&
  iMapEcal90=18, iMapEcal150=19, iMapEcal220=20
! Indices from old foreground model
! Remaining adjustments are made when paramfilename is read in according to model specified
integer :: iDl_Radio_TE_150=2, iAlphaRadio_TE=3, iDl_Radio_EE_150=4, iAlphaRadio_EE=5,&
  iTDSFG=6, iDl_DSFG_TE_150=7, iBetaDSFG_TE=8, iDl_DSFG_EE_150=9, iBetaDSFG_EE=10

integer, parameter :: cl_te_kind = 3, cl_ee_kind = 1
INTEGER:: BOK = 0
INTEGER:: CLIK_LMAX,CLIK_LMIN,big_clik_lmax
real(8), dimension(:), allocatable :: cl_clik_TE,cl_clik_EE,cl_clik_param
  
! Declare global variables
integer :: SPT3G_windows_lmin, SPT3G_windows_lmax
real(8) :: aberration_coefficient
integer :: poisson_switch, ssl_switch, radio_gal_switch, dsfg_switch, dust_switch
integer :: N_s, N_b, N_freq, bin_min, bin_max
  
real(8) :: nu_0_radio ! Unresolved radio galaxies
real(8) :: nu_0_dust ! Polarised galactic dust
real(8) :: nu_0_dsfg ! Dusty star-forming galaxies


integer, parameter :: N_freq_0 = 3 ! Number of frequency bands
integer, parameter :: N_b_0 = 44 ! Maximum (uncropped) number of bandpowers
integer, parameter :: N_s_0 = 12 ! Maximum (uncropped) number of spectra that can be used for parameter estimation

! Nuisance parameter priors
integer :: Cal_prior_switch,Kappa_prior_switch
integer :: AlphaDustEE_prior_switch, AlphaDustTE_prior_switch,BetaDustEE_prior_switch, BetaDustTE_prior_switch
real(8) :: Kappa_prior_mean, Kappa_prior_sigma
real(8) :: AlphaDustEE_prior_mean, AlphaDustEE_prior_sigma
real(8) :: AlphaDustTE_prior_mean, AlphaDustTE_prior_sigma
real(8) :: BetaDustEE_prior_mean, BetaDustEE_prior_sigma
real(8) :: BetaDustTE_prior_mean, BetaDustTE_prior_sigma
real(8) :: BeamCovScaling

integer,allocatable :: spectra_to_fit_bandpower_indices(:)

real(8), allocatable :: full_bandpowers(:,:)
real(8), allocatable :: bdp_covariance(:,:)
real(8), allocatable :: beam_covariance(:,:)
real(8), allocatable :: cal_inv_covariance(:,:)
real(8), allocatable :: full_windows(:,:,:)
real(8), allocatable :: spectra_to_fit_nu_eff(:,:)

integer, allocatable :: spectra_to_fit_tcal_indices(:,:) ! Gives the indices in the DataParams object for the TT spec calibration
integer, allocatable :: spectra_to_fit_ecal_indices(:,:) ! Gives the indices in the DataParams object for the EE spec calibration
integer, allocatable :: spectra_to_fit_kind(:)
integer, dimension(6) :: cal_row_use

logical :: use_simple_poisson_foregrounds
integer, allocatable :: spectra_to_fit_poisson_indices(:) ! Gives the indices in the DataParams object for the poisson amplitude

contains
subroutine save_1D(arr,name,id)
  real(8),intent(in),dimension(:)::arr
  character(len=*),intent(in)::name
  integer,intent(in)::id
  character(len=1024) :: filename, did
  integer::i

  write (did,"(I4,A4)") id,".txt"
  filename = trim(name)//trim(adjustl(did))
  open(unit=44,file=trim(filename))
  do i=1,size(arr)
    write(44,*) i,arr(i)
  end do
  close(44)
end subroutine save_1D

subroutine save_2D(arr,name,id)
  real(8),intent(in),dimension(:,:)::arr
  character(len=*),intent(in)::name
  integer,intent(in)::id
  character(len=1024) :: filename, did
  integer::i,j
  integer, dimension(2) :: ShapeArray

  write (did,"(I4,A4)") id,".txt"
  filename = trim(name)//trim(adjustl(did))
  open(unit=44,file=trim(filename))
  ShapeArray = shape(arr)
  do i=1,ShapeArray(1)
    write(44,*) (arr(i,j) ,j=1,ShapeArray(2))
  end do
  close(44)
end subroutine save_2D

function SPT3G_EETE_LogLike(ClTE, ClEE, DataParams) result(SPT_LogLike)
  implicit none

  real(8), intent(in),dimension(0:SPT3G_windows_lmax+1) :: ClTE,ClEE
  real(8), intent(in) :: DataParams(20)
  real(8) :: SPT_LogLike, SPT_PriorLogLike, detcov

  integer i, j, i_spec, ix, ix_1, ix_2,l
  character(LEN=:), allocatable :: current_spec_str, current_map_1_str, current_map_2_str, current_field, current_freq_1, current_freq_2
  
  real(8), dimension(SPT3G_windows_lmax+1) :: current_Dl_theory_unbinned_CMB_only ! Length is a buffer for derivatives
  real(8), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: current_Dl_theory_unbinned
  real(8), dimension(SPT3G_windows_lmax+1) :: Cl_derivative 
  real(8), dimension(N_b) :: current_Dl_theory_binned
  real(8), dimension(N_b*N_s) :: full_Dl_data_theory_binned, tmp
  real(8), dimension(N_b*N_s) :: Dl_data_theory_difference
  real(8), dimension(N_b*N_s,N_b*N_s) :: covariance, cov_copy

  real(8), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: super_sample_lensing
  real(8), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: aberration_correction
  real(8), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: simple_poisson
  real(8), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: radio_galaxies
  real(8), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: dsfg
  real(8), dimension(SPT3G_windows_lmin:SPT3G_windows_lmax) :: dust
  real(8) :: calibration
  
  real(8), dimension(sum(cal_row_use)) :: cal_vec,tmpcal

  integer te_poisson_switch
  integer current_iRadio, current_iAlphaRadio
  integer current_iDSFG, current_iBetaDSFG
  integer current_iDust, current_iAlphaDust, current_iBetaDust

  real(8), dimension(SPT3G_windows_lmax+1) :: ells

  integer :: info

  !print *, "SPT-3G Y1 EE/TE: Calculating likelihood.",SPT3G_windows_lmax
  !call save_1D(DataParams,"Params_",1000)

  ! Grab ells helper
  do i=1,SPT3G_windows_lmax+1
   ells(i) = i
  end do

  ! Loop over all spectra and get frequency scaling
  ! Only do this here if the scaling depends on the nuisance parameters
  ! Calculate the scaling for each map
  ! Add an array in the ini that reduces the requested spectra to what maps we have to do
  ! Loop over those maps and calculate the scaling
  ! Save in a way easily accessible in the spectra loop below (maybe just 2 N_s long freq_1/2 arrays)

  ! Loop over all spectra and add foregrounds
  ! Everything is in Dl space unless otherwise explicitly mentioned
  do i_spec=1,N_s
    !print *, i_spec,spectra_to_fit_kind(i_spec)
    ! Grab field specifics
    ! Theory CMB and foregrounds that are handled differently
    if (spectra_to_fit_kind(i_spec) .eq. cl_ee_kind) then
      current_Dl_theory_unbinned_CMB_only = ClEE(1:) /twoPI*(ells*(ells+1))
      
      te_poisson_switch = 1

      current_iRadio = iDl_Radio_EE_150 ! If simple foregrounds are used these point to the wrong params, but they never get used
      current_iAlphaRadio = iAlphaRadio_EE
      current_iDSFG = iDl_DSFG_EE_150
      current_iBetaDSFG = iBetaDSFG_EE
      
      current_iDust = iADust_EE_150
      current_iAlphaDust = iAlphaDust_EE
      current_iBetaDust = iBetaDust_EE

      ! Calibration for EE: 1/(Ecal_1*Ecal_2) since we matched the EE spectrum to Planck's
      calibration = 1/( DataParams(spectra_to_fit_ecal_indices(i_spec,1))*DataParams(spectra_to_fit_ecal_indices(i_spec,2)) )
      !print *,"go",calibration
    else if (spectra_to_fit_kind(i_spec) .eq. cl_te_kind) then
      current_Dl_theory_unbinned_CMB_only = ClTE(1:) /twoPI*(ells*(ells+1))
      
      
      te_poisson_switch = 0

      current_iRadio = iDl_Radio_TE_150
      current_iAlphaRadio = iAlphaRadio_TE
      current_iDSFG = iDl_DSFG_TE_150
      current_iBetaDSFG = iBetaDSFG_TE

      current_iDust = iADust_TE_150
      current_iAlphaDust = iAlphaDust_TE
      current_iBetaDust = iBetaDust_TE
    
      ! Calibration for TE: 0.5*(1/(Tcal_1*Ecal_2) + 1/(Tcal_2*Ecal_1))
      calibration = 1 /( DataParams(spectra_to_fit_tcal_indices(i_spec,1))*DataParams(spectra_to_fit_ecal_indices(i_spec,2)) )+&
        1 /( DataParams(spectra_to_fit_tcal_indices(i_spec,2))*DataParams(spectra_to_fit_ecal_indices(i_spec,1)) )
      calibration = calibration * 0.5

    else
      !print *, "Cannot recognise field in spectrum: ", spectra_to_fit_kind(i_spec) 
      stop
    end if
    !print *,"derivative"
    ! Grab Cl derivative for convenience
    Cl_derivative = current_Dl_theory_unbinned_CMB_only*twoPI/(ells*(ells+1)) ! Convert to Cl
    Cl_derivative(2:SPT3G_windows_lmax) = 0.5*( Cl_derivative(3:SPT3G_windows_lmax+1)&
      - Cl_derivative(1:SPT3G_windows_lmax-1) ) ! Take derivative
    Cl_derivative(1) = Cl_derivative(2) ! Handle endpoints approximately
    Cl_derivative(SPT3G_windows_lmax+1) = Cl_derivative(SPT3G_windows_lmax)
    !!if (CosmoSettings%lmax_computed_cl .LT. SPT3G_windows_lmax) then
    !!  Cl_derivative(CosmoSettings%lmax_computed_cl) = 0.75*Cl_derivative(CosmoSettings%lmax_computed_cl-1)&
    !!    + 0.25*Cl_derivative(CosmoSettings%lmax_computed_cl+2) ! Smooth spike from Boltzmann solver to look-up table in CAMB
    !!  Cl_derivative(CosmoSettings%lmax_computed_cl+1) = 0.75*Cl_derivative(CosmoSettings%lmax_computed_cl+2)&
    !!    + 0.25*Cl_derivative(CosmoSettings%lmax_computed_cl-1)
    !!end if

    ! Add CMB
    current_Dl_theory_unbinned = current_Dl_theory_unbinned_CMB_only(SPT3G_windows_lmin:SPT3G_windows_lmax)
    
    !print *,"ssl"
    ! Add super sample lensing
    ! (In Cl space) SSL = -k/l^2 d/dln(l) (l^2Cl) = -k(l*dCl/dl + 2Cl)
    super_sample_lensing = ells(SPT3G_windows_lmin:SPT3G_windows_lmax)*Cl_derivative(SPT3G_windows_lmin:SPT3G_windows_lmax) ! l*dCl/dl
    super_sample_lensing = super_sample_lensing * ells(SPT3G_windows_lmin:SPT3G_windows_lmax)&
      *(ells(SPT3G_windows_lmin:SPT3G_windows_lmax)+1)/(twoPI) ! Convert this part to Dl space already
    super_sample_lensing = super_sample_lensing + 2*current_Dl_theory_unbinned_CMB_only(SPT3G_windows_lmin:SPT3G_windows_lmax) ! 2Dl
    super_sample_lensing = -1 * super_sample_lensing * DataParams(iKappa) ! -kappa
    super_sample_lensing = super_sample_lensing * ssl_switch ! Switch on/off
    current_Dl_theory_unbinned = current_Dl_theory_unbinned + super_sample_lensing
  
    ! Aberration correction
    ! AC = beta*l(l+1)dCl/dln(l)/(2pi)
    ! Note that the CosmoMC internal aberration correction and the SPTpol Henning likelihood differ
    ! CosmoMC uses dCl/dl, Henning et al dDl/dl
    ! In fact, CosmoMC is correct: https://journals-aps-org.eu1.proxy.openathens.net/prd/pdf/10.1103/PhysRevD.89.023003
    !print *,"aberration"
    aberration_correction = -1*aberration_coefficient*Cl_derivative(SPT3G_windows_lmin:SPT3G_windows_lmax)&
      *ells(SPT3G_windows_lmin:SPT3G_windows_lmax)&
      *ells(SPT3G_windows_lmin:SPT3G_windows_lmax)&
      *(ells(SPT3G_windows_lmin:SPT3G_windows_lmax)+1)/(twoPI)
    current_Dl_theory_unbinned = current_Dl_theory_unbinned + aberration_correction
  
    ! Check the poisson foreground model
    if (use_simple_poisson_foregrounds) then
      !print *,"poisson"
        
        ! Simple poisson foregrounds
        ! This is any poisson power. Meant to describe both radio galaxies and DSFG. By giving each frequency combination an amplitude
        ! to play with this gives complete freedom to the data
        simple_poisson = te_poisson_switch*ells(SPT3G_windows_lmin:SPT3G_windows_lmax)*(ells(SPT3G_windows_lmin:SPT3G_windows_lmax)+1)&
          *DataParams(spectra_to_fit_poisson_indices(i_spec))/(3000*3001)
        simple_poisson = simple_poisson * poisson_switch ! Switch on/off
        current_Dl_theory_unbinned = current_Dl_theory_unbinned + simple_poisson

    else

        ! Unresovled radio galaxies
        ! Amplitude is defined at l=3000 in Dl space, so we canvert to Cl power first and then to Dl
        radio_galaxies = ells(SPT3G_windows_lmin:SPT3G_windows_lmax)*(ells(SPT3G_windows_lmin:SPT3G_windows_lmax)+1)*DataParams(current_iRadio)/(3000*3001)
        radio_galaxies = radio_galaxies * RadioGalFreqScaling(DataParams(current_iAlphaRadio),nu_0_radio,spectra_to_fit_nu_eff(i_spec,1))&
          *RadioGalFreqScaling(DataParams(current_iAlphaRadio),nu_0_radio,spectra_to_fit_nu_eff(i_spec,2))
        radio_galaxies = radio_galaxies * poisson_switch ! Switch on/off
        current_Dl_theory_unbinned = current_Dl_theory_unbinned + radio_galaxies

        ! Dusty star-forming galaxies
        ! Same thing as unresolved radio galaxies, but have a different frequency scaling
        dsfg = ells(SPT3G_windows_lmin:SPT3G_windows_lmax)*(ells(SPT3G_windows_lmin:SPT3G_windows_lmax)+1)*DataParams(current_iDSFG)/(3000*3001)
        dsfg = dsfg * DSFGFreqScaling(DataParams(current_iBetaDSFG),DataParams(iTDSFG),nu_0_dsfg,spectra_to_fit_nu_eff(i_spec,1))&
          *DSFGFreqScaling(DataParams(current_iBetaDSFG),DataParams(iTDSFG),nu_0_dsfg,spectra_to_fit_nu_eff(i_spec,2))
        dsfg = dsfg * poisson_switch ! Switch on/off
        current_Dl_theory_unbinned = current_Dl_theory_unbinned + dsfg

    end if
      
      !print *,"dust"

    ! Polarised galactic dust
    dust = DataParams(current_iDust)*(ells(SPT3G_windows_lmin:SPT3G_windows_lmax)/80d0)**(DataParams(current_iAlphaDust)+2.0d0)
    !print *,spectra_to_fit_nu_eff(i_spec,1),spectra_to_fit_nu_eff(i_spec,2),i_spec
    dust = dust * DustFreqScaling(DataParams(current_iBetaDust),DataParams(iTDust),nu_0_dust,spectra_to_fit_nu_eff(i_spec,1))&
      *DustFreqScaling(DataParams(current_iBetaDust),DataParams(iTDust),nu_0_dust,spectra_to_fit_nu_eff(i_spec,2))
    dust = dust * dust_switch ! Switch on/off
      
    current_Dl_theory_unbinned = current_Dl_theory_unbinned + dust
  
    ! Scale by calibration
    current_Dl_theory_unbinned = current_Dl_theory_unbinned * calibration
    
    ! Bin
    ! Select the right bins from the window function in this step
    ! (Difference to SPTpol likelihood due to the changed window function format (already transposed)

    !print *,N_b,SPT3G_windows_lmax-SPT3G_windows_lmin+1,1.0d0,spectra_to_fit_bandpower_indices(i_spec)
    !call save_1D(current_Dl_theory_unbinned,"unbinned",i_spec)
    call dgemv('N',N_b,SPT3G_windows_lmax-SPT3G_windows_lmin+1,1.0d0,&
      full_windows(bin_min:bin_max,:,spectra_to_fit_bandpower_indices(i_spec)),&
      N_b,current_Dl_theory_unbinned,1,0d0,current_Dl_theory_binned,1)
    !call save_1D(current_Dl_theory_binned,"binned",i_spec)
    
    full_Dl_data_theory_binned((i_spec-1)*N_b+1:i_spec*N_b) = current_Dl_theory_binned

    !!  
    ! Write smooth and binned theory out
    !!if (print_spectrum) then
    !!  save_arr=0
    !!  call OpenWriteBinaryFile(trim("likelihood_tests/binned_theory_")//trim(current_spec_str),fid,4_8 * 2)
    !!  do l=1,N_b
    !!    save_arr(1)=l
    !!    save_arr(2)=current_Dl_theory_binned(l)
    !!    write(fid,rec=l) save_arr(1:2)
    !!  enddo
    !!  close(fid)
    !!end if

    !!if (print_spectrum) then
    !!  save_arr=0
    !!  call OpenWriteBinaryFile(trim("likelihood_tests/unbinned_theory_")//trim(current_spec_str),fid,4_8 * 2)
    !!  do l=SPT3G_windows_lmin,SPT3G_windows_lmax
    !!    save_arr(1)=l
    !!    save_arr(2)=current_Dl_theory_unbinned_CMB_only(l)
    !!    save_arr(3)=aberration_correction(l)
    !!    save_arr(4)=super_sample_lensing(l)
    !!    if (use_simple_poisson_foregrounds) then
    !!      save_arr(6) = simple_poisson(l)
    !!    else
    !!      save_arr(5)=radio_galaxies(l)
    !!      save_arr(6)=dsfg(l)
    !!    end if
    !!    save_arr(7)=dust(l)
    !!    write(fid,rec=l-SPT3G_windows_lmin+1) save_arr(1:7)
    !!  enddo
    !!  close(fid)
    !!end if

    ! Take the difference to the measured bandpowers and preapre long (data-model) vector
    ! Crop data bandpowers to the right bins in this step too
    Dl_data_theory_difference((i_spec-1)*N_b+1:i_spec*N_b) = &
      current_Dl_theory_binned - full_bandpowers(bin_min:bin_max,spectra_to_fit_bandpower_indices(i_spec))


    ! Write binned model-data out
    !!if (print_spectrum) then
    !!  save_arr=0
    !!  call OpenWriteBinaryFile(trim("likelihood_tests/delta_data_model_")//trim(current_spec_str),fid,4_8 * 2)
    !!  do l=this%bin_min,this%bin_max
    !!    save_arr(1)=l
    !!    save_arr(2)=current_Dl_theory_binned(l-this%bin_min+1)
    !!    save_arr(3)=this%full_bandpowers(l,spectra_to_fit_bandpower_indices(i_spec))
    !!    save_arr(4)=current_Dl_theory_binned(l-this%bin_min+1) - this%full_bandpowers(l,spectra_to_fit_bandpower_indices(i_spec))
    !!    save_arr(5)=Dl_data_theory_difference((i_spec-1)*this%N_b+1+l-this%bin_min)
    !!    write(fid,rec=l-this%bin_min+1) save_arr(1:5)
    !!  enddo
    !!  close(fid)
    !!end if

  end do

  ! Construct the full covariance matrix
  covariance=0
  do i=1,N_b*N_s
    do j=1,N_b*N_s
      covariance(i,j)=beam_covariance(i,j)*full_Dl_data_theory_binned(i)*full_Dl_data_theory_binned(j)
    end do
  end do
  
  covariance = covariance + bdp_covariance
  cov_copy = covariance

  ! Calculate data likelihood
  
  !call save_2D(covariance,"covpre",122)
  call dpotrf ('L', N_b*N_s, covariance, N_b*N_s, info)
  if (info/=0) then 
    print *,"argl, dpotrf failed with",info
    stop
  endif
  !call save_2D(covariance,"covpost",122)
  SPT_LogLike = 0
  do i=1, N_b*N_s
    SPT_LogLike = SPT_LogLike  + log(covariance(i,i))
  end do
  detcov = SPT_LogLike
  tmp = Dl_data_theory_difference
  !call save_1D(Dl_data_theory_difference,"diff_data",122)

  call DPOTRS('L', N_b*N_s, 1, covariance, N_b*N_s, tmp, N_b*N_s, info )
  if (info/=0) then 
    print *,"argl dpotrs failed with",info
    stop
  endif
  !call save_1D(tmp,"tmp",122)
  
    
  SPT_LogLike = SPT_LogLike + dot_product(tmp,Dl_data_theory_difference)/2.0
  

  !!SPT_LogLike = Matrix_GaussianLogLikeDouble(covariance, Dl_data_theory_difference)

  ! Get the contribution from the determinant part
  !detcov = Matrix_GaussianLogLikeDouble(cov_copy, Dl_data_theory_difference*0)

  ! Get cal vector
  ix=1
  do i=1,N_freq_0*2
    if (cal_row_use(i) .eq. 1) then
      cal_vec(ix) = DataParams(iMapTcal90+i-1)! Cheeky, relies on the order of calibration parameters in DataParams
      ix = ix + 1
    end if
  end do
  cal_vec = log(cal_vec)

  ! Calculate prior likelihood
  SPT_PriorLogLike = 0.0d0
  !SPT_PriorLogLike = SPT_PriorLogLike + Cal_prior_switch*0.5d0*Matrix_QuadForm(cal_inv_covariance,cal_vec)
  SPT_PriorLogLike = SPT_PriorLogLike + Cal_prior_switch*0.5d0*dot_product(cal_vec,MatMul(cal_inv_covariance,cal_vec))
  SPT_PriorLogLike = SPT_PriorLogLike + Kappa_prior_switch*0.5d0*((DataParams(iKappa) - Kappa_prior_mean)/Kappa_prior_sigma)**2
  SPT_PriorLogLike = SPT_PriorLogLike + AlphaDustEE_prior_switch*0.5d0*((DataParams(iAlphaDust_EE) - AlphaDustEE_prior_mean)/AlphaDustEE_prior_sigma)**2
  SPT_PriorLogLike = SPT_PriorLogLike + AlphaDustTE_prior_switch*0.5d0*((DataParams(iAlphaDust_TE) - AlphaDustTE_prior_mean)/AlphaDustTE_prior_sigma)**2
  !! Print out likelihood contributions
  !if (feedback > 1) then
  !  print *, "SPT-3G Y1 EE/TE: Data LogLike: ", SPT_LogLike
  !  print *, "SPT-3G Y1 EE/TE: Data LogLike Detcov: ", detcov
  !  print *, "SPT-3G Y1 EE/TE: Prior LogLike: ", SPT_PriorLogLike
  !end if

  ! Add data and prior likelihoods and return
  SPT_LogLike = SPT_LogLike + SPT_PriorLogLike 

end function SPT3G_EETE_LogLike

! Dust Frequency Scaling
function DustFreqScaling(beta,Tdust,nu0,nu_eff) result(fdust)
  real(8), intent(in) :: beta
  real(8), intent(in) :: Tdust
  real(8), intent(in) :: nu_eff
  real(8), intent(in) :: nu0 ! Pivot frequency
  real(8) :: fdust

  fdust = (nu_eff/nu0)**beta
  fdust = fdust*Bnu(nu_eff,nu0,Tdust)/dBdT(nu_eff,nu0,T_CMB)

end function DustFreqScaling

! Poisson Radio Galaxies Frequency scaling
function RadioGalFreqScaling(alpha,nu0,nu_eff) result(fradiogal)
  real(8), intent(in) :: nu_eff
  real(8), intent(in) :: alpha
  real(8), intent(in) :: nu0 ! Pivot frequency
  real(8) :: fradiogal
  
  fradiogal = (nu_eff/nu0)**alpha
  fradiogal = fradiogal/dBdT(nu_eff,nu0,T_CMB)

end function RadioGalFreqScaling


! Dusty Star-Forming Galaxies Frequency scaling
function DSFGFreqScaling(beta,Tdust,nu0,nu_eff) result(fdsfg)
  real(8), intent(in) :: nu_eff
  real(8), intent(in) :: beta
  real(8), intent(in) :: Tdust
  real(8), intent(in) :: nu0 ! Pivot frequency
  real(8)  :: fdsfg

  fdsfg = (nu_eff/nu0)**beta
  fdsfg = fdsfg*Bnu(nu_eff,nu0,Tdust)/dBdT(nu_eff,nu0,T_CMB)


end function DSFGFreqScaling


! Planck function normalised to 1 at nu0
function Bnu(nu,nu0,T)
  real(8) Bnu,nu,nu0,T

  Bnu = (nu/nu0)**3
  Bnu = Bnu * (exp( Ghz_Kelvin*nu0/T)-1d0) / (exp( Ghz_Kelvin*nu/T)-1d0)

end function Bnu


! Derivative of Planck function normalised to 1 at nu0
function dBdT(nu,nu0,T)
  real(8) dBdT,dBdT0,nu,nu0,T,x,x0

  x0 = Ghz_Kelvin*nu0/T
  x = Ghz_Kelvin*nu/T

  dBdT0 = x0**4 * exp(x0) / (exp(x0)-1)**2
  dBdT =  x**4 * exp(x) / (exp(x)-1)**2

  dBdT = dBdT/dBdT0

end function dBdT


end module

SUBROUTINE spt3g_parameter_init(windows_lmin,windows_lmax,iaberration_coefficient, &
                                issl_switch, ipoisson_switch, idust_switch,inu_0_radio, inu_0_dsfg, inu_0_dust, &
                                ibin_min,ibin_max,iCal_prior_switch,iKappa_prior_switch, &
                                iKappa_prior_mean, iKappa_prior_sigma, &
                                iAlphaDustEE_prior_switch,iAlphaDustEE_prior_mean,iAlphaDustEE_prior_sigma, &
                                iAlphaDustTE_prior_switch,iAlphaDustTE_prior_mean,iAlphaDustTE_prior_sigma,iBeamCovScaling, &
                                ifull_bandpowers,iN_s, ispectra_to_fit_bandpower_indices,&
                                ibdp_covariance,ibeam_covariance, &
                                ispectra_to_fit_tcal_indices,ispectra_to_fit_ecal_indices, &
                                ical_row_use,ical_inv_covariance,ifull_windows,ispectra_to_fit_nu_eff, &
                                ispectra_to_fit_poisson_indices,ispectra_to_fit_kind)
  USE CMB_SPT3G_EETE_2020_clik

  implicit none

  integer,intent(in):: windows_lmin,windows_lmax,ibin_min,ibin_max
  real(8),intent(in):: iaberration_coefficient,inu_0_radio, inu_0_dsfg, inu_0_dust, iKappa_prior_mean 
  real(8),intent(in):: iKappa_prior_sigma, iAlphaDustEE_prior_mean, iAlphaDustEE_prior_sigma,iAlphaDustTE_prior_mean, iAlphaDustTE_prior_sigma, iBeamCovScaling
  integer,intent(in) :: ipoisson_switch, issl_switch,  idust_switch,iCal_prior_switch
  integer,intent(in) :: iKappa_prior_switch ,iN_s
  integer,intent(in) :: iAlphaDustEE_prior_switch, iAlphaDustTE_prior_switch
  real(8),intent(in),dimension(N_b_0*N_s_0) :: ifull_bandpowers
  integer,intent(in),dimension(iN_s) :: ispectra_to_fit_bandpower_indices
  real(8),intent(in),dimension((bin_max-bin_min+1)*iN_s*(bin_max-bin_min+1)*iN_s) :: ibdp_covariance,ibeam_covariance
  integer,intent(in),dimension(iN_s*2) ::ispectra_to_fit_tcal_indices,ispectra_to_fit_ecal_indices
  integer,intent(in),dimension(6) :: ical_row_use
  real(8),intent(in),dimension(sum(ical_row_use),sum(ical_row_use)) :: ical_inv_covariance
  real(8),intent(in),dimension(N_b_0*(1+windows_lmin-windows_lmax)*N_s_0) :: ifull_windows
  
  integer,intent(in),dimension(iN_s) :: ispectra_to_fit_poisson_indices
  integer,intent(in),dimension(iN_s) :: ispectra_to_fit_kind

  real(8),intent(in),dimension(N_s*2) :: ispectra_to_fit_nu_eff
  integer::i,j,k

  SPT3G_windows_lmin = windows_lmin
  SPT3G_windows_lmax = windows_lmax
  !print *,windows_lmin,windows_lmin,SPT3G_windows_lmin,SPT3G_windows_lmax
  allocate(cl_clik_EE((SPT3G_windows_lmax+1)))
  allocate(cl_clik_TE((SPT3G_windows_lmax+1)))
  allocate(cl_clik_Param(20))

  aberration_coefficient = iaberration_coefficient
  ssl_switch = issl_switch
  poisson_switch = ipoisson_switch
  dust_switch = idust_switch
  nu_0_radio = inu_0_radio 
  nu_0_dsfg = inu_0_dsfg 
  nu_0_dust = inu_0_dust 

  use_simple_poisson_foregrounds = .true.

  !call spectra_to_fit_list%SetFromString(spectra_to_fit_list_str)
  !N_s = spectra_to_fit_list%Count

  bin_min = ibin_min
  bin_max = ibin_max
  N_b = bin_max-bin_min+1

  Cal_prior_switch = iCal_prior_switch
  
  Kappa_prior_switch = iKappa_prior_switch
  Kappa_prior_mean = iKappa_prior_mean
  Kappa_prior_sigma = iKappa_prior_sigma

  AlphaDustEE_prior_switch = iAlphaDustEE_prior_switch
  AlphaDustEE_prior_mean = iAlphaDustEE_prior_mean
  AlphaDustEE_prior_sigma = iAlphaDustEE_prior_sigma


  AlphaDustTE_prior_switch = iAlphaDustTE_prior_switch
  AlphaDustTE_prior_mean = iAlphaDustTE_prior_mean
  AlphaDustTE_prior_sigma = iAlphaDustTE_prior_sigma

  BeamCovScaling = iBeamCovScaling

  allocate(full_bandpowers(N_b_0,N_s_0))
  do i=1,N_b_0
    do j=1,N_s_0
      full_bandpowers(i,j) = ifull_bandpowers((i-1)*N_s_0+j)
    end do
  end do

  N_s = iN_s
  allocate(spectra_to_fit_bandpower_indices(N_s))
  do i=1,N_s
    spectra_to_fit_bandpower_indices(i) = ispectra_to_fit_bandpower_indices(i)
    !print *,spectra_to_fit_bandpower_indices(i)
  enddo

  allocate(bdp_covariance(N_b*N_s,N_b*N_s))
  allocate(beam_covariance(N_b*N_s,N_b*N_s))
  do i=1,N_b*N_s
    do j=1,N_b*N_s
      bdp_covariance(i,j) = ibdp_covariance((i-1)*N_b*N_s+j)
      beam_covariance(i,j) = ibeam_covariance((i-1)*N_b*N_s+j)
    end do
  end do


  allocate(spectra_to_fit_tcal_indices(N_s,2))
  allocate(spectra_to_fit_ecal_indices(N_s,2))
  do i=1,N_s
  spectra_to_fit_tcal_indices(i,1) = ispectra_to_fit_tcal_indices((i-1)*2+1)
  spectra_to_fit_tcal_indices(i,2) = ispectra_to_fit_tcal_indices((i-1)*2+2)  
  spectra_to_fit_ecal_indices(i,1) = ispectra_to_fit_ecal_indices((i-1)*2+1)
  spectra_to_fit_ecal_indices(i,2) = ispectra_to_fit_ecal_indices((i-1)*2+2)
  end do

  do i=1,6
    cal_row_use(i) = ical_row_use(i)
  end do

  allocate(cal_inv_covariance(sum(cal_row_use),sum(cal_row_use)))
  do i=1,sum(cal_row_use)
    do j=1,sum(cal_row_use)
    cal_inv_covariance(i,j) = ical_inv_covariance(i,j)
    end do
  end do

  allocate(full_windows(N_b_0,1+SPT3G_windows_lmax-SPT3G_windows_lmin,N_s_0))
  do i=1,N_b_0
    do j=1,1+SPT3G_windows_lmax-SPT3G_windows_lmin
      do k=1,N_s_0
        full_windows(i,j,k) = ifull_windows((i-1)*N_s_0*(1+SPT3G_windows_lmax-SPT3G_windows_lmin)+(j-1)*N_s_0+k)
      end do
    end do
  end do

  allocate(spectra_to_fit_nu_eff(N_s,2))
  do i=1,N_s
    spectra_to_fit_nu_eff(i,1) = ispectra_to_fit_nu_eff((i-1)*2+1)
    spectra_to_fit_nu_eff(i,2) = ispectra_to_fit_nu_eff((i-1)*2+2)
  enddo

  allocate(spectra_to_fit_poisson_indices(N_s))
  do i=1,N_s
    spectra_to_fit_poisson_indices(i) = ispectra_to_fit_poisson_indices(i)
  end do

  allocate(spectra_to_fit_kind(N_s))
  do i=1,N_s
    spectra_to_fit_kind(i) = ispectra_to_fit_kind(i)
  end do
end subroutine spt3g_parameter_init

subroutine spt3g_lkl(LKL,CL)
  USE CMB_SPT3G_EETE_2020_clik
  REAL(8),INTENT(OUT)::LKL
  REAL(8),INTENT(IN),DIMENSION((SPT3G_windows_lmax+2)*2+20)::CL
  integer::i
  lkl = SPT3G_EETE_LogLike(CL(SPT3G_windows_lmax+2+1:(SPT3G_windows_lmax+2)*2),CL(:SPT3G_windows_lmax+2),CL((SPT3G_windows_lmax+2)*2+1:)) 
end subroutine spt3g_lkl

SUBROUTINE SPT3G_ONLY_ONE(MOK)
  USE CMB_SPT3G_EETE_2020_clik
  INTEGER,INTENT(OUT)::MOK
  MOK = BOK
  BOK = 1
END SUBROUTINE  SPT3G_ONLY_ONE

SUBROUTINE SPT3G_FREE()
  USE CMB_SPT3G_EETE_2020_clik
  BOK =0
  !deallocate(cltt)
END SUBROUTINE  SPT3G_FREE


Module CMB_SPT3G_TTEEE_2018_clik
use CMB_SPT3G_2018_TTTEEE
use SPT3G_utils
implicit None



class(TSPT3G_2018_TTTEEE_Likelihood), allocatable :: single_lkl

integer::bok=0
real(8), dimension(:), allocatable :: cl_clik,Dataparam,CMBparam 
integer::clik_lmax



end module 

SUBROUTINE spt3g_ttteee2018_parameter_init(eSPT3G_windows_lmin, eSPT3G_windows_lmax, full_bandpower, full_bandpower_list_string, l_full_bandpower_list_string, full_covariance_matrix, &
                                          full_covariance_list_string, l_full_covariance_list_string, full_beam_covariance_matrix,  full_beam_covariance_list_string ,l_full_beam_covariance_list_string ,      &
                                          full_cal_covariance_matrix, full_windows, full_window_list_string, l_full_window_list_string,     &
                                          nu_eff_matrix,nu_eff_list_string, l_nu_eff_list_string, spectra_to_fit_list_string, l_spectra_to_fit_list_string,      &
                                          spec_bin_min_list_string, l_spec_bin_min_list_string,  spec_bin_max_list_string, l_spec_bin_max_list_string,late_crop_msk_string, l_late_crop_msk_string,&
                                          ecov_eval_cut_threshold,ecov_eval_large_number_replacement,beam_cov_scale, &
                                          aberration_coefficient, enu_0_galdust, eT_galdust, enu_0_CIB, eT_CIB, enu_0_tSZ, etSZCosmologyScalingEnabled,full_tSZ_template,ekSZCosmologyScalingEnabled,full_kSZ_template,einclude_logdet)

  USE CMB_SPT3G_TTEEE_2018_clik
  use CMB_SPT3G_2018_TTTEEE
  use SPT3G_utils
  implicit none

  integer,intent(in) :: eSPT3G_windows_lmin, eSPT3G_windows_lmax
  real(mcp),intent(in) :: full_bandpower(N_b_0_total)
  integer,intent(in)::l_full_bandpower_list_string,l_full_covariance_list_string, &
                      l_full_beam_covariance_list_string,l_full_window_list_string, &
                      l_nu_eff_list_string, l_spectra_to_fit_list_string, l_spec_bin_min_list_string, &
                      l_spec_bin_max_list_string, l_late_crop_msk_string

  character(LEN=l_full_bandpower_list_string),intent(in) :: full_bandpower_list_string

  real(mcp),intent(in) :: full_covariance_matrix(N_b_0_total,N_b_0_total)
  character(LEN=l_full_covariance_list_string),intent(in) :: full_covariance_list_string

  real(mcp),intent(in) :: full_beam_covariance_matrix(N_b_0_total,N_b_0_total)
  character(LEN=l_full_beam_covariance_list_string),intent(in) :: full_beam_covariance_list_string

  real(mcp),intent(in)::full_cal_covariance_matrix(N_freq_0*2,N_freq_0*2)
  

  real(mcp),intent(in):: full_windows(N_b_0_EE,1+eSPT3G_windows_lmax-eSPT3G_windows_lmin,N_s_0)
  character(LEN=l_full_window_list_string),intent(in) :: full_window_list_string

  real(mcp),intent(in) :: nu_eff_matrix(5,N_freq_0)
  character(LEN=l_nu_eff_list_string),intent(in) :: nu_eff_list_string

  character(LEN=l_spectra_to_fit_list_string),intent(in) :: spectra_to_fit_list_string

  character(LEN=l_spec_bin_min_list_string),intent(in) :: spec_bin_min_list_string
  character(LEN=l_spec_bin_max_list_string),intent(in) :: spec_bin_max_list_string 

  character(LEN=l_late_crop_msk_string),intent(in) :: late_crop_msk_string
  real(mcp),intent(in):: ecov_eval_cut_threshold,ecov_eval_large_number_replacement,beam_cov_scale,aberration_coefficient
  real(mcp),intent(in):: enu_0_galdust, eT_galdust, enu_0_CIB, eT_CIB, enu_0_tSZ
  integer,intent(in) :: etSZCosmologyScalingEnabled,ekSZCosmologyScalingEnabled,einclude_logdet

  real(mcp),intent(in) :: full_tSZ_template(1+eSPT3G_windows_lmax-eSPT3G_windows_lmin)
  real(mcp),intent(in) :: full_kSZ_template(1+eSPT3G_windows_lmax-eSPT3G_windows_lmin)


  allocate(single_lkl)
  call SPT3G_2018_TTTEEE_Ini_external(single_lkl,eSPT3G_windows_lmin, eSPT3G_windows_lmax, full_bandpower, full_bandpower_list_string, full_covariance_matrix, &
                                          full_covariance_list_string, full_beam_covariance_matrix,  full_beam_covariance_list_string ,      &
                                          full_cal_covariance_matrix, full_windows, full_window_list_string,     &
                                          nu_eff_matrix,nu_eff_list_string, spectra_to_fit_list_string,      &
                                          spec_bin_min_list_string,  spec_bin_max_list_string,late_crop_msk_string, &
                                          ecov_eval_cut_threshold,ecov_eval_large_number_replacement,beam_cov_scale, &
                                          aberration_coefficient, enu_0_galdust, eT_galdust, enu_0_CIB, eT_CIB, enu_0_tSZ, etSZCosmologyScalingEnabled,full_tSZ_template,ekSZCosmologyScalingEnabled,full_kSZ_template,einclude_logdet)
  
  clik_lmax = eSPT3G_windows_lmax

  allocate(cl_clik(3*(clik_lmax+1)))
  allocate(Dataparam(37))
  allocate(CMBparam(6))

end subroutine spt3g_ttteee2018_parameter_init

subroutine spt3g_ttteee2018_lkl(LKL,CL)
  USE CMB_SPT3G_TTEEE_2018_clik
  use CMB_SPT3G_2018_TTTEEE
  use CMB_SPT3G_2018_TTTEEE_foregrounds
  use SPT3G_utils
  REAL(8),INTENT(OUT)::LKL
  REAL(8),INTENT(IN),DIMENSION((SPT3G_windows_lmax+1)*3+37+7)::CL
  integer::i,ell

  do ell=1,SPT3G_windows_lmax
    cl_clik(ell) = CL(ell+1)*ell*(ell+1.)/2./pi
    cl_clik(ell+SPT3G_windows_lmax+1) = CL(ell+1+SPT3G_windows_lmax+1)*ell*(ell+1.)/2./pi
    cl_clik(ell+(SPT3G_windows_lmax+1)*2) = CL(ell+1+(SPT3G_windows_lmax+1)*2)*ell*(ell+1.)/2./pi
  end do
  Dataparam(:) = CL((SPT3G_windows_lmax+1)*3+1:(SPT3G_windows_lmax+1)*3+37)
  if (tSZCosmologyScalingEnabled.or.kSZCosmologyScalingEnabled) then
    CMBparam = CL((SPT3G_windows_lmax+1)*3+37+1:(SPT3G_windows_lmax+1)*3+37+6)
  else
    CMBparam = (/0,0,0,0,0,0/)
  endif

  lkl = SPT3G_2018_TTTEEE_LogLike_external(single_lkl, cl_clik,CMBparam,Dataparam)
end subroutine spt3g_ttteee2018_lkl

SUBROUTINE SPT3G_ttteee2018_ONLY_ONE(MOK)
  USE CMB_SPT3G_TTEEE_2018_clik
  INTEGER,INTENT(OUT)::MOK
  MOK = BOK
  BOK = 1
END SUBROUTINE  SPT3G_ttteee2018_ONLY_ONE

SUBROUTINE SPT3G_ttteee2018_FREE()
  USE CMB_SPT3G_TTEEE_2018_clik
  BOK =0
  deallocate(single_lkl)
  deallocate(cl_clik)
  deallocate(CMBParam)
  deallocate(Dataparam)
END SUBROUTINE  SPT3G_ttteee2018_FREE
