    module DarkEnergyFluid
    use DarkEnergyInterface
    use results
    use constants
    use classes
    implicit none

    type, extends(TDarkEnergyEqnOfState) :: TDarkEnergyFluid
        !comoving sound speed is always exactly 1 for quintessence
        !(otherwise assumed constant, though this is almost certainly unrealistic)
    contains
    procedure :: ReadParams => TDarkEnergyFluid_ReadParams
    procedure, nopass :: PythonClass => TDarkEnergyFluid_PythonClass
    procedure, nopass :: SelfPointer => TDarkEnergyFluid_SelfPointer
    procedure :: Init => TDarkEnergyFluid_Init
    procedure :: PerturbedStressEnergy => TDarkEnergyFluid_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TDarkEnergyFluid_PerturbationEvolve
    end type TDarkEnergyFluid

    !Example implementation of fluid model using specific analytic form
    !(approximate effective axion fluid model from arXiv:1806.10608, with c_s^2=1 if n=infinity (w_n=1))
    !This is an example, it's not supposed to be a rigorous model!  (not very well tested)
    type, extends(TDarkEnergyModel) :: TAxionEffectiveFluid
        real(dl) :: w_n = 1._dl !Effective equation of state when oscillating
        real(dl) :: fde_zc = 0._dl ! energy density fraction at a_c (not the same as peak dark energy fraction)
        real(dl) :: zc  !transition redshift (scale factor a_c)
        real(dl) :: theta_i = const_pi/2 !Initial value
        !om is Omega of the early DE component today (assumed to be negligible compared to omega_lambda)
        !omL is the lambda component of the total dark energy omega
        real(dl), private :: a_c, pow, om, omL, acpow, freq, n !cached internally
    contains
    procedure :: ReadParams =>  TAxionEffectiveFluid_ReadParams
    procedure, nopass :: PythonClass => TAxionEffectiveFluid_PythonClass
    procedure, nopass :: SelfPointer => TAxionEffectiveFluid_SelfPointer
    procedure :: Init => TAxionEffectiveFluid_Init
    procedure :: w_de => TAxionEffectiveFluid_w_de
    procedure :: grho_de => TAxionEffectiveFluid_grho_de
    procedure :: PerturbedStressEnergy => TAxionEffectiveFluid_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TAxionEffectiveFluid_PerturbationEvolve
    end type TAxionEffectiveFluid

    type, extends(TDarkEnergyModel) :: TWEDEFluid
        real(dl) :: w_0 = -0.9_dl ! Equation of state for the w constant component
        real(dl) :: wa = 0.0_dl ! w = w0 + wa * (1-a)
        real(dl) :: w_n = 1._dl ! Effective equation of state when oscillating
        real(dl) :: zc  !transition redshift (scale factor a_c)
        real(dl) :: fde_zc = 0._dl ! energy density fraction at a_c (not the same as peak dark energy fraction)
        real(dl) :: theta_i = const_pi/2 !Initial value
        real(dl) :: cs2_de = 1._dl
        real(dl) :: grhode_today
        !om is Omega of the early DE component today (assumed to be negligible compared to omega_w)
        !omW is the w_const component of the total dark energy omega
        real(dl), private :: a_c, pow, om, omW, acpow, freq, n !cached internally
    contains
    procedure :: ReadParams =>  TWEDEFluid_ReadParams
    procedure, nopass :: PythonClass => TWEDEFluid_PythonClass
    procedure, nopass :: SelfPointer => TWEDEFluid_SelfPointer
    procedure :: Init => TWEDEFluid_Init
    procedure :: w_de => TWEDEFluid_w_de
    procedure :: grho_de => TWEDEFluid_grho_de
    procedure :: PerturbedStressEnergy => TWEDEFluid_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TWEDEFluid_PerturbationEvolve 
    procedure :: grho_wwa
    procedure :: omega_axion
    procedure :: w_axion     
    end type TWEDEFluid

    contains

    !----------------- Skeleton subroutines and functions for a fluid model ---------------------

    subroutine TDarkEnergyFluid_ReadParams(this, Ini)
        use IniObjects
        class(TDarkEnergyFluid) :: this
        class(TIniFile), intent(in) :: Ini

        call this%TDarkEnergyEqnOfState%ReadParams(Ini)
        this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)
    end subroutine TDarkEnergyFluid_ReadParams

    function TDarkEnergyFluid_PythonClass()
        character(LEN=:), allocatable :: TDarkEnergyFluid_PythonClass
        TDarkEnergyFluid_PythonClass = 'DarkEnergyFluid'
    end function TDarkEnergyFluid_PythonClass

    subroutine TDarkEnergyFluid_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TDarkEnergyFluid), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TDarkEnergyFluid_SelfPointer

    subroutine TDarkEnergyFluid_Init(this, State)
        use classes
        class(TDarkEnergyFluid), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State

        call this%TDarkEnergyEqnOfState%Init(State)

        if (this%is_cosmological_constant) then
            this%num_perturb_equations = 0
        else
            if (this%use_tabulated_w) then
                if (any(this%equation_of_state%F<-1)) &
                    error stop 'Fluid dark energy model does not allow w crossing -1'
            elseif (this%wa/=0 .and. &
                ((1+this%w_lam < -1.e-6_dl) .or. 1+this%w_lam + this%wa < -1.e-6_dl)) then
                error stop 'Fluid dark energy model does not allow w crossing -1'
            end if
            this%num_perturb_equations = 2
        end if
    end subroutine TDarkEnergyFluid_Init

    subroutine TDarkEnergyFluid_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
        class(TDarkEnergyFluid), intent(inout) :: this
        real(dl), intent(out) :: dgrhoe, dgqe
        real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
        real(dl), intent(in) :: ay(*)
        real(dl), intent(inout) :: ayprime(*)
        integer, intent(in) :: w_ix

        if (this%no_perturbations) then
            dgrhoe=0
            dgqe=0
        else
            dgrhoe = ay(w_ix) * grhov_t
            dgqe = ay(w_ix + 1) * grhov_t * (1 + w)
        end if
    end subroutine TDarkEnergyFluid_PerturbedStressEnergy

    subroutine TDarkEnergyFluid_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
        class(TDarkEnergyFluid), intent(in) :: this
        real(dl), intent(inout) :: ayprime(:)
        real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
        integer, intent(in) :: w_ix
        real(dl) :: Hv3_over_k, loga

        Hv3_over_k =  3*adotoa* y(w_ix + 1) / k
        !density perturbation
        ayprime(w_ix) = -3 * adotoa * (this%cs2_lam - w) *  (y(w_ix) + (1 + w) * Hv3_over_k) &
            -  (1 + w) * k * y(w_ix + 1) - (1 + w) * k * z
        if (this%use_tabulated_w) then
            !account for derivatives of w
            loga = log(a)
            if (loga > this%equation_of_state%Xmin_interp .and. loga < this%equation_of_state%Xmax_interp) then
                ayprime(w_ix) = ayprime(w_ix) - adotoa*this%equation_of_state%Derivative(loga)* Hv3_over_k
            end if
        elseif (this%wa/=0) then
            ayprime(w_ix) = ayprime(w_ix) + Hv3_over_k*this%wa*adotoa*a
        end if
        !velocity
        if (abs(w+1) > 1e-6) then
            ayprime(w_ix + 1) = -adotoa * (1 - 3 * this%cs2_lam) * y(w_ix + 1) + &
                k * this%cs2_lam * y(w_ix) / (1 + w)
        else
            ayprime(w_ix + 1) = 0
        end if
    end subroutine TDarkEnergyFluid_PerturbationEvolve

    !------------------------- Axion Effective Fluid Implementation -------------------------------

    subroutine TAxionEffectiveFluid_ReadParams(this, Ini)
        use IniObjects
        class(TAxionEffectiveFluid) :: this
        class(TIniFile), intent(in) :: Ini

        call this%TDarkEnergyModel%ReadParams(Ini)
        if (Ini%HasKey('AxionEffectiveFluid_a_c')) then
            error stop 'AxionEffectiveFluid inputs changed to AxionEffectiveFluid_fde_zc and AxionEffectiveFluid_zc'
        end if
        this%w_n  = Ini%Read_Double('AxionEffectiveFluid_w_n')
        this%fde_zc  = Ini%Read_Double('AxionEffectiveFluid_fde_zc')
        this%zc  = Ini%Read_Double('AxionEffectiveFluid_zc')
        call Ini%Read('AxionEffectiveFluid_theta_i', this%theta_i)
    end subroutine TAxionEffectiveFluid_ReadParams

    function TAxionEffectiveFluid_PythonClass()
        character(LEN=:), allocatable :: TAxionEffectiveFluid_PythonClass
        TAxionEffectiveFluid_PythonClass = 'AxionEffectiveFluid'
    end function TAxionEffectiveFluid_PythonClass

    subroutine TAxionEffectiveFluid_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TAxionEffectiveFluid), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TAxionEffectiveFluid_SelfPointer

    subroutine TAxionEffectiveFluid_Init(this, State)
        use classes
        class(TAxionEffectiveFluid), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State
        real(dl) :: grho_rad, F, p, mu, xc, n

        select type(State)
        class is (CAMBdata)
            this%is_cosmological_constant = this%fde_zc==0
            this%pow = 3*(1+this%w_n)
            this%a_c = 1/(1+this%zc)
            this%acpow = this%a_c**this%pow
            ! Omega (= rho/rhocr0)  in early de at z=0
            this%om = 2*this%fde_zc/(1-this%fde_zc)*&
                (State%grho_no_de(this%a_c)/this%a_c**4/State%grhocrit + State%Omega_de)/(1 + 1/this%acpow)
            this%omL = State%Omega_de - this%om !Omega_de is total dark energy density today
            this%num_perturb_equations = 2
            if (this%w_n < 0.9999) then
                ! n <> infinity
                !get (very) approximate result for sound speed parameter; arXiv:1806.10608  Eq 30 (but mu may not exactly agree with what they used)
                n = nint((1+this%w_n)/(1-this%w_n))
                !Assume radiation domination, standard neutrino model; H0 factors cancel
                grho_rad = (kappa/c**2*4*sigma_boltz/c**3*State%CP%tcmb**4*Mpc**2*(1+3.046*7._dl/8*(4._dl/11)**(4._dl/3)))
                xc = this%a_c**2/2/sqrt(grho_rad/3)
                F=7./8
                p=1./2
                mu = 1/xc*(1-cos(this%theta_i))**((1-n)/2.)*sqrt((1-F)*(6*p+2)*this%theta_i/n/sin(this%theta_i))
                this%freq =  mu*(1-cos(this%theta_i))**((n-1)/2.)* &
                    sqrt(const_pi)*Gamma((n+1)/(2.*n))/Gamma(1+0.5/n)*2.**(-(n**2+1)/(2.*n))*3.**((1./n-1)/2)*this%a_c**(-6./(n+1)+3) &
                    *( this%a_c**(6*n/(n+1.))+1)**(0.5*(1./n-1))
                this%n = n
            end if
        end select
    end subroutine TAxionEffectiveFluid_Init

    function TAxionEffectiveFluid_w_de(this, a)
        class(TAxionEffectiveFluid) :: this
        real(dl) :: TAxionEffectiveFluid_w_de
        real(dl), intent(IN) :: a
        real(dl) :: rho, apow, acpow

        apow = a**this%pow
        acpow = this%acpow
        rho = this%omL+ this%om*(1+acpow)/(apow+acpow)
        TAxionEffectiveFluid_w_de = this%om*(1+acpow)/(apow+acpow)**2*(1+this%w_n)*apow/rho - 1
    end function TAxionEffectiveFluid_w_de

    function TAxionEffectiveFluid_grho_de(this, a)  ! Relative density (8 pi G a^4 rho_de / grhov)
        class(TAxionEffectiveFluid) :: this
        real(dl) :: TAxionEffectiveFluid_grho_de, apow
        real(dl), intent(IN) :: a

        if(a == 0.d0)then
            TAxionEffectiveFluid_grho_de = 0.d0
        else
            apow = a**this%pow
            TAxionEffectiveFluid_grho_de = (this%omL*(apow+this%acpow) + this%om*(1+this%acpow))*a**4 &
                /((apow+this%acpow)*(this%omL+this%om))
        endif
    end function TAxionEffectiveFluid_grho_de

    subroutine TAxionEffectiveFluid_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
        class(TAxionEffectiveFluid), intent(in) :: this
        real(dl), intent(inout) :: ayprime(:)
        real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
        integer, intent(in) :: w_ix
        real(dl) Hv3_over_k, deriv, apow, acpow, cs2, fac

        if (this%w_n < 0.9999) then
            fac = 2*a**(2-6*this%w_n)*this%freq**2
            cs2 = (fac*(this%n-1) + k**2)/(fac*(this%n+1) + k**2)
        else
            cs2 = 1
        end if
        apow = a**this%pow
        acpow = this%acpow
        Hv3_over_k =  3*adotoa* y(w_ix + 1) / k
        ! dw/dlog a/(1+w)
        deriv  = (acpow**2*(this%om+this%omL)+this%om*acpow-apow**2*this%omL)*this%pow &
            /((apow+acpow)*(this%omL*(apow+acpow)+this%om*(1+acpow)))
        !density perturbation
        ayprime(w_ix) = -3 * adotoa * (cs2 - w) *  (y(w_ix) + Hv3_over_k) &
            -   k * y(w_ix + 1) - (1 + w) * k * z  - adotoa*deriv* Hv3_over_k
        !(1+w)v
        ayprime(w_ix + 1) = -adotoa * (1 - 3 * cs2 - deriv) * y(w_ix + 1) + &
            k * cs2 * y(w_ix)
    end subroutine TAxionEffectiveFluid_PerturbationEvolve

    subroutine TAxionEffectiveFluid_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
        class(TAxionEffectiveFluid), intent(inout) :: this
        real(dl), intent(out) :: dgrhoe, dgqe
        real(dl), intent(in) :: a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
        real(dl), intent(in) :: ay(*)
        real(dl), intent(inout) :: ayprime(*)
        integer, intent(in) :: w_ix

        dgrhoe = ay(w_ix) * grhov_t
        dgqe = ay(w_ix + 1) * grhov_t
    end subroutine TAxionEffectiveFluid_PerturbedStressEnergy

    !-------------------- EDE + wconst implementation --------------------

    subroutine TWEDEFluid_ReadParams(this, Ini)
        use IniObjects
        class(TWEDEFluid) :: this
        class(TIniFile), intent(in) :: Ini

        call this%TDarkEnergyModel%ReadParams(Ini)
        this%w_0 = Ini%Read_Double('wEDE_wconst')
        this%wa = Ini%Read_Double('wEDE_wa')
        this%w_n  = Ini%Read_Double('wEDE_w_n')
        this%fde_zc  = Ini%Read_Double('wEDE_fde_zc')
        this%zc  = Ini%Read_Double('wEDE_zc')
        call Ini%Read('wEDE_theta_i', this%theta_i)
    end subroutine TWEDEFluid_ReadParams

    function TWEDEFluid_PythonClass()
        character(LEN=:), allocatable :: TWEDEFluid_PythonClass
        TWEDEFluid_PythonClass = 'wEDEFluid'
    end function TWEDEFluid_PythonClass

    subroutine TWEDEFluid_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TWEDEFluid), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TWEDEFluid_SelfPointer

    real(dl) function grho_wwa(this, a, grhov)
        ! Returns 8*pi*G*a^4*rho_wconst if using grhov
        ! grhov can be changed if you want other proportionality factors
        class(TWEDEFluid), intent(in) :: this
        real(dl) :: a, grhov
        grho_wwa = grhov * a**(-3*(1+this%w_0 + this%wa)) * a**4
        if (this%wa /= 0) then
            grho_wwa = grho_wwa * exp(-3*this%wa*(1-a))
        end if
    end function grho_wwa

    real(dl) function omega_axion(this, a)
        ! Returns omega_axion = rho_axion/rho_cr0, equation (15) from https://arxiv.org/pdf/1806.10608.pdf
        class(TWEDEFluid), intent(in) :: this
        real(dl), intent(in) :: a
        real(dl) :: z, oma_zc
        z = 1._dl/a - 1._dl
        oma_zc = this%om * (1 + (1+this%zc)**(3*(1+this%w_n)))/2._dl
        omega_axion = (2._dl * oma_zc)/(1 + ( (1+this%zc)/(1+z) )**(3*(1+this%w_n)) )
    end function omega_axion

    function TWEDEFluid_grho_de(this, a)  !relative density (8 pi G a^4 rho_de /grhov)
        class(TWEDEFluid) :: this
        real(dl) :: TWEDEFluid_grho_de, apow
        real(dl), intent(IN) :: a
        real(dl) :: grhow, grhoaxion 

        grhow = this%grho_wwa(a, (this%omW)/(this%omW + this%om))
        grhoaxion = this%omega_axion(a) * a**4 / (this%omW + this%om)

        TWEDEFluid_grho_de = grhow + grhoaxion
    end function TWEDEFluid_grho_de

    real(dl) function w_axion(this, a)
        ! Returns w_axion, equation (16) from https://arxiv.org/pdf/1806.10608.pdf
        class(TWEDEFluid), intent(in) :: this
        real(dl), intent(in) :: a
        real(dl) :: apow
        apow = a**this%pow
        w_axion = (1._dl + this%w_n) / (1._dl + (this%acpow/apow)) - 1._dl
    end function w_axion

    function TWEDEFluid_w_de(this, a)
        ! Here, we need to consider both w_const and axion
        class(TWEDEFluid) :: this
        real(dl) :: TWEDEFluid_w_de
        real(dl), intent(IN) :: a
        real(dl) :: rho, apow, acpow, w_late
        w_late = this%w_0 + this%wa * (1-a)
        rho = this%grho_wwa(a, this%omW) + this%omega_axion(a)
        TWEDEFluid_w_de = (w_late * this%grho_wwa(a, this%omW) + this%w_axion(a)*this%omega_axion(a)) / rho
    end function TWEDEFluid_w_de

    subroutine TWEDEFluid_Init(this, State)
        use classes
        class(TWEDEFluid), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State
        real(dl) :: grho_rad, F, p, mu, xc, n

        select type(State)
        class is (CAMBdata)
            this%grhode_today = State%grhov
            this%is_cosmological_constant = .false.
            this%pow = 3*(1+this%w_n)
            this%a_c = 1/(1+this%zc)
            this%acpow = this%a_c**this%pow
            !Omega in early de at z=0
            this%om = 2*this%fde_zc/(1-this%fde_zc)*&
                (State%grho_no_de(this%a_c)/this%a_c**4/State%grhocrit +&
                 State%Omega_de * this%a_c**(-3*(1+this%w_0)) * exp(-3 * this%wa * (1-this%a_c)))&
                /(1 + 1/this%acpow)
            this%omW = State%Omega_de - this%om ! Omega_de is total dark energy density today
            this%num_perturb_equations = 4
            if (this%w_n < 0.9999) then
                ! n <> infinity
                !get (very) approximate result for sound speed parameter; arXiv:1806.10608  Eq 30 (but mu may not exactly agree with what they used)
                n = nint((1+this%w_n)/(1-this%w_n))
                !Assume radiation domination, standard neutrino model; H0 factors cancel
                grho_rad = (kappa/c**2*4*sigma_boltz/c**3*State%CP%tcmb**4*Mpc**2*(1+3.046*7._dl/8*(4._dl/11)**(4._dl/3)))
                xc = this%a_c**2/2/sqrt(grho_rad/3)
                F=7./8
                p=1./2
                mu = 1/xc*(1-cos(this%theta_i))**((1-n)/2.)*sqrt((1-F)*(6*p+2)*this%theta_i/n/sin(this%theta_i))
                this%freq =  mu*(1-cos(this%theta_i))**((n-1)/2.)* &
                    sqrt(const_pi)*Gamma((n+1)/(2.*n))/Gamma(1+0.5/n)*2.**(-(n**2+1)/(2.*n))*3.**((1./n-1)/2)*this%a_c**(-6./(n+1)+3) &
                    *( this%a_c**(6*n/(n+1.))+1)**(0.5*(1./n-1))
                this%n = n
            end if
        end select
    end subroutine TWEDEFluid_Init

    subroutine TWEDEFluid_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
        class(TWEDEFluid), intent(in) :: this
        real(dl), intent(inout) :: ayprime(:)
        real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
        integer, intent(in) :: w_ix
        real(dl) Hva3_over_k, Hvw3_over_k, deriv, apow, acpow, cs2_axion, fac
        real(dl) :: w_axion, w_late

        ! Indices: w_ix = \delta_a, w_ix + 1 = v_a, w_ix + 2 = \delta_w, w_ix + 3 = v_w 

        w_axion = this%w_axion(a)
        w_late = this%w_0 + this%wa * (1 - a)

        ! -------------- Axion fluid perturbations --------------

        if (this%w_n < 0.9999) then
            fac = 2*a**(2-6*this%w_n)*this%freq**2
            cs2_axion = (fac*(this%n-1) + k**2)/(fac*(this%n+1) + k**2)
        else
            cs2_axion = 1
        end if

        apow = a**this%pow
        acpow = this%acpow
        Hva3_over_k =  3*adotoa* y(w_ix + 1) / k
        ! dw/dlog a/(1+w)
        deriv  = (acpow**2*(this%om+this%omW)+this%om*acpow-apow**2*this%omW)*this%pow &
            /((apow+acpow)*(this%omW*(apow+acpow)+this%om*(1+acpow)))
        ! axion density perturbation
        ! JVR - TESTING, CHANGE AFTER!!
        ayprime(w_ix) = 0
        ayprime(w_ix + 1) = 0
        !ayprime(w_ix) = -3 * adotoa * (cs2_axion - w_axion) *  (y(w_ix) + Hva3_over_k) &
        !    -   k * y(w_ix + 1) - (1 + w_axion) * k * z  - adotoa*deriv * Hva3_over_k
        ! axion (1+w)v
        !ayprime(w_ix + 1) = -adotoa * (1 - 3 * cs2_axion - deriv) * y(w_ix + 1) + &
        !    k * cs2_axion * y(w_ix)

        !------------- w0wa perturbations -----------------
        ! Consider time derivatives!!!!
        ! d(lnw)/d(lna) = a/w * dw/da
        ! if w = w0 + wa(1-a) then dw/da = -wa
        ! so d(lnw)/d(lna) = - (wa * a)/(w0 + wa(1-a))
        ! ayprime(w_ix + 2) = 0
        ! ayprime(w_ix + 3) = 0
        Hvw3_over_k =  3 * adotoa * y(w_ix + 3) / k
        !density perturbation
        ayprime(w_ix + 2) = -3 * adotoa * (this%cs2_de - w_late) * (y(w_ix + 2) + (1 + w_late) * Hvw3_over_k) &
            -  (1 + w_late) * k * y(w_ix + 3) - (1 + w_late) * k * z
        if (this%wa/=0) then
            ayprime(w_ix + 2) = ayprime(w_ix + 2) + Hvw3_over_k*this%wa*adotoa*a
        end if
        !velocity
        if (abs(w_late+1) > 1e-6) then
            ayprime(w_ix + 3) = -adotoa * (1 - 3 * this%cs2_de) * y(w_ix + 3) + &
                k * this%cs2_de * y(w_ix+2) / (1 + w_late)
        else
            ayprime(w_ix + 3) = 0
        end if
    end subroutine TWEDEFluid_PerturbationEvolve

    subroutine TWEDEFluid_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
        class(TWEDEFluid), intent(inout) :: this
        real(dl), intent(out) :: dgrhoe, dgqe
        real(dl), intent(in) :: a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
        real(dl), intent(in) :: ay(*)
        real(dl), intent(inout) :: ayprime(*)
        integer, intent(in) :: w_ix
        real(dl) :: dgrhoe_axion, dgqe_axion, dgrhoe_w, dgqe_w
        real(dl) :: grhow, grhoaxion, faxion, fw, grhov

        ! Again, it's important to talk about indices:
        ! w_ix: \delta_axion
        ! w_ix + 1: v_axion
        ! w_ix + 2: \delta_w
        ! w_ix + 3: v_w

        grhow = this%grhode_today * a ** (1._dl - 3. * this%w_0 - 3. * this%wa) / a**2
        if (this%wa/=0) grhow=grhow*exp(-3. * this%wa * (1._dl - a))
        ! grhow = this%grho_wwa(a, (this%omW)/(this%omW + this%om)) / a**2
        grhoaxion = this%omega_axion(a) * a**2 / (this%omW + this%om)

        !dgrhoe_axion = ay(w_ix) * grhov_t * faxion
        !dgqe_axion = ay(w_ix + 1) * grhov_t * faxion
        dgrhoe_w = ay(w_ix + 2) * grhow
        dgqe_w = ay(w_ix + 3) * grhow * (1 + this%w_0 + this%wa*(1._dl - a))

        dgrhoe = dgrhoe_w
        dgqe = dgqe_w
    end subroutine TWEDEFluid_PerturbedStressEnergy

    end module DarkEnergyFluid