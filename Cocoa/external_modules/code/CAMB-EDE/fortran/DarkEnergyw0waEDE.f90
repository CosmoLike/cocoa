module MultiFluidDE
  use DarkEnergyInterface
  use results
  use constants
  use classes
  implicit none

  private
  integer :: max_num_of_fluids = 4
  integer :: max_num_of_params = 10

  type, extends(TCambComponent) :: TMultiFluidDE
    logical :: is_cosmological_constant = .false.
    integer :: num_perturb_equations = 0
    integer :: num_of_fluids = 1
    real(dl), dimension(:,:), allocatable :: de_params
    integer, dimension(:), allocatable :: models
  contains
  procedure :: Init => TMultiFluidDE_Init
  procedure :: w_de => TMultiFluidDE_w_de
  procedure :: grho_de => TMultiFluidDE_grho_de
  procedure :: BackgroundDensityAndPressure => TMultiFluidDE_density
  procedure :: Effective_w_wa => TMultiFluidDE_Effective_w_wa
  procedure :: PerturbedStressEnergy => TMultiFluidDE_PerturbedStressEnergy
  procedure :: PerturbationEvolve => TMultiFluidDE_PerturbationEvolve
  procedure, nopass :: PythonClass => TMultiFluidDE_PythonClass
  procedure, nopass :: SelfPointer => TMultiFluidDE_SelfPointer
  procedure :: ReadParams => TMultiFluidDE_ReadParams
  end type TMultiFluidDE

  subroutine TMultiFluidDE_w_de(this, a, w_de)
    class(TMultiFluidDE), intent(in) :: this
    real(dl), dimension(:), allocatable, intent(out) :: w_de
    real(dl), intent(in) :: a
    integer :: i

    do i = 1, this%num_of_fluids
      if (i == 1) then
        ! Implement DE parametrizations
        if (this%models(i) == 1) then
          ! Implement model 1 for fluid 1
          w_de(i) = -1._dl
        else
          ! Implement error message
      else if (i == 2) then
        ! Implement DE parametrizations
        if (this%models(i) == 1) then
          ! Implement model 1 for fluid 2
          w_de(i) = -1._dl
        else
          ! Implement error message
      else 
        ! Error message: not implemented
      end if
    end do
  end function TMultiFluidDE_w_de

  subroutine TMultiFluidDE_grho_de(this, a, grho_de)
    class(TMultiFluidDE), intent(in) :: this
    real(dl), dimension(:), allocatable, intent(out) :: grho_de
    real(dl), intent(in) :: a
    integer :: i

    do i = 1, this%num_of_fluids
      if (i == 1) then
        ! Implement DE parametrizations
        if (this%models(i) == 1) then
          ! Implement model 1 for fluid 1
          grho_de(i) = 0._dl
        else
          ! Implement error message
      else if (i == 2) then
        ! Implement DE parametrizations
        if (this%models(i) == 1) then
          ! Implement model 1 for fluid 2
          grho_de(i) = 0._dl
        else
          ! Implement error message
      else 
        ! Error message: not implemented
      end if
    end do
  end function TMultiFluidDE_grho_de

  subroutine TMultiFluidDE_density(this, grhov, a, grhov_t, w)
    !Get grhov_t = 8*pi*G*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), dimension(:), allocatable, intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w

    ! grhov = rho_cr * (1 - Omega_m - Omega_r - Omega_nu)

    grhov_t = this%grho_de(a)
    ! VM thinks grhov must be an array (TODO)
    do i = 1, this%num_of_fluids
      ! Ensure a valid result
      if (a > 1e-10) then
          grhov_t(i) = grhov * grhov_t(i) / (a * a)
      else
          grhov_t(i) = 0._dl
      end if
      if (present(w)) then
        w = this%w_de(a)
      end if
    end do
  end subroutine TMultiFluidDE_density

  subroutine TMultiFluidDE_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
        class(TMultiFluidDE), intent(inout) :: this
        real(dl), intent(out) :: dgrhoe, dgqe
        real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
        real(dl), intent(in) :: ay(*)
        real(dl), intent(inout) :: ayprime(*)
        integer, intent(in) :: w_ix

        ! TODO clean variables, change dgrhoe and dgqe to be arrays

        do i=1, this%num_of_fluids
          dgrhoe(i) = ay(w_ix + 2*(i-1)) * grhov_t(i)
          dgqe(i) = ay(w_ix + 1 + 2*(i-1)) * grhov_t(i) * (1 + w(i))
        end do
  end subroutine TMultiFluidDE_PerturbedStressEnergy

  !subroutine w_derivative_lna(this, a)
    ! TODO check the term involving 1/3 * (d(ln(1+w))/lna) in the \delta equation
  !end subroutine

  subroutine TMultiFluidDE_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
        class(TMultiFluidDE), intent(in) :: this
        real(dl), intent(inout) :: ayprime(:)
        real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
        integer, intent(in) :: w_ix
        real(dl) :: Hv3_over_k, loga

        ! Assuming cs_2 = 1 for all fluids
        do i=1, this%num_of_fluids
          Hv3_over_k =  3 * adotoa * y(w_ix + 1 + 2*(i-1)) / k
          !density perturbation
          ayprime(w_ix + 2*(i-1)) = -3 * adotoa * (1._dl - w(i) *  (y(w_ix + 2*(i-1)) + (1 + w(i)) * Hv3_over_k) &
              -  (1 + w(i)) * k * y(w_ix + 1 + 2*(i-1)) - (1 + w(i)) * k * z
          ayprime(w_ix + 2*(i-1)) = ayprime(w_ix + 2*(i-1)) + Hv3_over_k*this%wa*adotoa*a
          ! TODO check the term involving 1/3 * (d(ln(1+w))/lna)
          !velocity
          if (abs(w(i)+1) > 1e-6) then
              ayprime(w_ix + 1 + 2*(i-1)) = -adotoa * (1 - 3 * 1._dl) * y(w_ix + 1 + 2*(i-1)) + &
                  k * 1._dl * y(w_ix + 2*(i-1)) / (1 + w(i))
          else
              ayprime(w_ix + 1 + 2*(i-1)) = 0
          end if
        end do
  end subroutine TMultiFluidDE_PerturbationEvolve

  subroutine TMultiFluidDE_Init(this, State)
        use classes
        class(TMultiFluidDE), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State

        this%num_perturb_equations = 2
        !elseif (this%wa/=0 .and. &
        !    ((1+this%w_lam < -1.e-6_dl) .or. 1+this%w_lam + this%wa < -1.e-6_dl)) then
        !    error stop 'Fluid dark energy model does not allow w crossing -1'
        !end if
  end subroutine TMultiFluidDE_Init

  subroutine TMultiFluidDE_ReadParams(this, Ini)
        use IniObjects
        class(TMultiFluidDE) :: this
        class(TIniFile), intent(in) :: Ini

        !call this%TDarkEnergyModel%ReadParams(Ini)
        !if (Ini%HasKey('AxionEffectiveFluid_a_c')) then
        !    error stop 'AxionEffectiveFluid inputs changed to AxionEffectiveFluid_fde_zc and AxionEffectiveFluid_zc'
        !end if
        !this%w_n  = Ini%Read_Double('AxionEffectiveFluid_w_n')
        !this%fde_zc  = Ini%Read_Double('AxionEffectiveFluid_fde_zc')
        !this%zc  = Ini%Read_Double('AxionEffectiveFluid_zc')
        !call Ini%Read('AxionEffectiveFluid_theta_i', this%theta_i)
  end subroutine TMultiFluidDE_ReadParams

  function TMultiFluidDE_PythonClass()
        character(LEN=:), allocatable :: TMultiFluidDE_PythonClass
        TMultiFluidDE_PythonClass = 'MultiFluidDE'
  end function TMultiFluidDE_PythonClass

  subroutine TMultiFluidDE_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TMultiFluidDE), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType
  end subroutine TMultiFluidDE_SelfPointer
end module MultiFluidDE