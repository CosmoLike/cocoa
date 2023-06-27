    ! Module for EDE + a component with constant w
    ! Author: João Victor S. Rebouças, December 2021

    module wEDE
        use DarkEnergyInterface
        use ScalarField
        use results
        use constants
        use classes
        implicit none
        private

        real(dl), parameter :: Tpl= sqrt(kappa*hbar/c**5)  ! sqrt(8 pi G hbar/c^5), reduced planck time

        type, extends(TScalarField) :: TEarlyDarkEnergyPlusW
            integer :: which_potential = 1
            real(dl) :: w_const = -0.999 ! Eqn of state for the constant w component
            real(dl) :: n = 5._dl ! Exponent for the EDE potential
            real(dl) :: f = 0.05 ! sqrt(8*pi*G)*f, AxionEDE frequency
            real(dl) :: m = 5d-54 ! AxionEDE mass in reduced Planck mass units
            real(dl) :: theta_i = 3.1_dl ! Initial value of phi/f for the AxionEDE model
            real(dl) :: frac_lambda0 = 1._dl ! fraction of dark energy density that is cosmological constant today
            real(dl) :: V0 = 1.d-115 ! Potential scale parameter for monomial potential
            real(dl) :: initial_phi = 1._dl ! Initial field value for the monomial potential
            logical :: use_zc = .true. ! If you want to adjust EDE parameters to give specified zc/fde 
            real(dl) :: zc = 2000 ! EDE peak redshift
            real(dl) :: fde_zc = 0.05 ! EDE peak energy contribution
            integer :: npoints = 200 ! Baseline number of log a steps; will be increased if needed when there are oscillations
            real(dl) :: zcfdeprec = 1d-3 ! Precision for zc/fde adjustment
            integer :: min_steps_per_osc = 10 ! For adjusting step to sample oscillations
            real(dl), dimension(:), allocatable :: fde, ddfde
        contains
            procedure :: Vofphi => TEarlyDarkEnergyPlusW_VofPhi
            procedure :: Init => TEarlyDarkEnergyPlusW_Init
            procedure :: ReadParams =>  TEarlyDarkEnergyPlusW_ReadParams
            procedure :: BackgroundDensityAndPressure => TEarlyDarkEnergyPlusW_BackgroundDensityAndPressure
            procedure, nopass :: PythonClass => TEarlyDarkEnergyPlusW_PythonClass
            procedure, nopass :: SelfPointer => TEarlyDarkEnergyPlusW_SelfPointer
            procedure, private :: fdeAta
            procedure, private :: fde_peak
            procedure, private :: check_error
            procedure :: calc_zc_fde
            procedure :: grho_w => grho_w
        end type TEarlyDarkEnergyPlusW

        procedure(TClassDverk) :: dverk

        public TEarlyDarkEnergyPlusW
    
    contains

    
    real(dl) function grho_w(this, a)
        ! JVR - calculates the second component energy density as grho_w = 8*pi*G*a^4*rho_w
        class(TEarlyDarkEnergyPlusW), intent(in) :: this
        grho_w = this%state%grhov * a**(-3*(1+this%w_const)) * a**4
    end function grho_w

    ! JVR - OK, but can it be done this way?
    subroutine TEarlyDarkEnergyPlusW_BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
        !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
        class(TEarlyDarkEnergyPlusW), intent(in) :: this
        real(dl), intent(in) :: grhov, a
        real(dl), intent(out) :: grhov_t
        real(dl), optional, intent(out) :: w
        real(dl) V, a2, grhov_lambda, phi, phidot

        ! JVR - can I simply return the sum of both energies from wconst and EDE??

        if (a >= this%astart) then
            a2 = a**2
            call this%ValsAta(a, phi, phidot)
            V = this%Vofphi(phi,0)
            grhov_t = phidot**2/2 + a2*V + this%State%Omega_de * a**(-3*(1+this%w_const)) * a2 + this%grho_w(a) / a2
            if (present(w)) then
                w = (phidot**2/2 - a2*V + this%w_const * this%grho_w(a))/grhov_t
            end if
        else
            grhov_t = this%State%Omega_de * a**(-3*(1+this%w_const)) * a2
            if (present(w)) then
                w = this%w_const
            end if
        end if
    end subroutine TEarlyDarkEnergyPlusW_BackgroundDensityAndPressure

    ! JVR - OK
    subroutine EvolveBackgroundLog(this,num,loga,y,yprime)
        ! Evolve the background equation in terms of loga.
        ! Variables are phi=y(1), a^2 phi' = y(2)
        ! Assume otherwise standard background components
        class(TScalarField) :: this
        integer num
        real(dl) y(num),yprime(num)
        real(dl) loga, a

        a = exp(loga)
        call this%EvolveBackground(num, a, y, yprime)
        yprime = yprime*a
    end subroutine EvolveBackgroundLog

    ! JVR - OK but needs testing
    subroutine EvolveBackground(this,num,a,y,yprime)
        ! Evolve the background equation in terms of a. 
        ! Variables are phi=y(1), a^2 phi' = y(2)
        ! Assume otherwise standard background components
        class(TScalarField) :: this
        integer num
        real(dl) y(num),yprime(num)
        real(dl) a, a2, tot
        real(dl) phi, grhode, phidot, adot

        a2 = a**2
        phi = y(1)
        phidot = y(2)/a2

        grhode = a2 * (0.5d0*phidot**2 + a2 * this%Vofphi(phi,0)) + this%grho_w(a)
        tot = this%state%grho_no_de(a) + grhode

        adot=sqrt(tot/3.0d0)
        yprime(1) = phidot/adot !d phi /d a
        yprime(2) = -a2**2*this%Vofphi(phi,1)/adot
    end subroutine EvolveBackground

    ! JVR - OK
    real(dl) function TScalarField_phidot_start(this,phi)
        class(TScalarField) :: this
        real(dl) :: phi
        TScalarField_phidot_start = 0
    end function TScalarField_phidot_start

    ! JVR - OK
    subroutine ValsAta(this,a,aphi,aphidot)
        class(TScalarField) :: this
        ! Do interpolation for background phi and phidot at a (precomputed in Init)
        real(dl) a, aphi, aphidot
        real(dl) a0,b0,ho2o6,delta,da
        integer ix

        if (a >= 0.9999999d0) then
            aphi= this%phi_a(this%npoints_linear+this%npoints_log)
            aphidot= this%phidot_a(this%npoints_linear+this%npoints_log)
            return
        elseif (a < this%astart) then
            aphi = this%phi_a(1)
            aphidot = 0
            return
        elseif (a > this%max_a_log) then
            delta= a-this%max_a_log
            ix = this%npoints_log + int(delta/this%da)
        else
            delta= log(a)-this%log_astart
            ix = int(delta/this%dloga)+1
        end if
        da = this%sampled_a(ix+1) - this%sampled_a(ix)
        a0 = (this%sampled_a(ix+1) - a)/da
        b0 = 1 - a0
        ho2o6 = da**2/6._dl
        aphi=b0*this%phi_a(ix+1) + a0*(this%phi_a(ix)-b0*((a0+1)*this%ddphi_a(ix)+(2-a0)*this%ddphi_a(ix+1))*ho2o6)
        aphidot=b0*this%phidot_a(ix+1) + a0*(this%phidot_a(ix)-b0*((a0+1)*this%ddphidot_a(ix)+(2-a0)*this%ddphidot_a(ix+1))*ho2o6)
    end subroutine ValsAta

    ! JVR - NOT IMPLEMENTED
    subroutine TEarlyDarkEnergyPlusW_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
        !Get density perturbation and heat flux
        class(TScalarField), intent(inout) :: this
        real(dl), intent(out) :: dgrhoe, dgqe
        real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
        real(dl), intent(in) :: ay(*)
        real(dl), intent(inout) :: ayprime(*)
        integer, intent(in) :: w_ix
        real(dl) phi, phidot, clxq, vq

        call this%ValsAta(a,phi,phidot)
        clxq = ay(w_ix) ! JVR: clxq is \delta\phi, the field perturbation in synchronous gauge.
        vq = ay(w_ix+1) ! vq is \delta\phi', the field perturbation derivative.
        dgrhoe = phidot*vq +clxq*a**2*this%Vofphi(phi,1)
        dgqe = k*phidot*clxq
        ! Check equation (A6) in https://arxiv.org/pdf/astro-ph/9801234.pdf
    end subroutine TEarlyDarkEnergyPlusW_PerturbedStressEnergy

    ! JVR - NOT IMPLEMENTED
    subroutine TEarlyDarkEnergyPlusW_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
        !Get conformal time derivatives of the density perturbation and velocity
        ! Need to set num_perturb_equations and your quantities begin at w_ix
        class(TQuintessence), intent(in) :: this
        real(dl), intent(inout) :: ayprime(:)
        real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
        integer, intent(in) :: w_ix
        real(dl) clxq, vq, phi, phidot

        call this%ValsAta(a,phi,phidot) ! wasting time calling this again..
        clxq=y(w_ix)
        vq=y(w_ix+1)
        ! Check equation (A3) in https://arxiv.org/pdf/astro-ph/9801234.pdf
        ayprime(w_ix)= vq
        ayprime(w_ix+1) = - 2*adotoa*vq - k*z*phidot - k**2*clxq - a**2*clxq*this%Vofphi(phi,2)
    end subroutine TEarlyDarkEnergyPlusW_PerturbationEvolve


    ! JVR - OK
    function TEarlyDarkEnergy_VofPhi(this, phi, deriv) result(VofPhi)
        !The input variable phi is sqrt(8*Pi*G)*psi
        !Returns (8*Pi*G)^(1-deriv/2)*d^{deriv}V(psi)/d^{deriv}psi evaluated at psi
        !return result is in 1/Mpc^2 units [so times (Mpc/c)^2 to get units in 1/Mpc^2]
        class(TEarlyDarkEnergy) :: this
        real(dl) :: phi, Vofphi
        integer :: deriv
        real(dl) :: theta, costheta, V0
        real(dl), parameter :: units = MPC_in_sec**2 /Tpl**2  !convert to units of 1/Mpc^2

        ! JVR - Recipe for implementing potentials:
        ! 1. Write your potential and its derivatives in natural units (assuming all parameters
        ! are in natural units of Mpl)
        ! 2. Multiply the result by this units factor
        
        ! Assume f = sqrt(kappa)*f_theory = f_theory/M_pl
        ! m = m_theory/M_Pl
        select case (this%which_potential)
            case(1) ! Rock 'n' Roll potential, already including the cosmological constant
                V0 = this%V0
                if (deriv == 0) then
                    Vofphi = units * V0 * phi**(2*this%n)
                else if (deriv == 1) then
                    Vofphi = units * V0 * (2 * this%n) * phi**(2*this%n - 1)
                else if (deriv == 2) then
                    Vofphi = units * V0 * (2 * this%n) * (2 * this%n - 1) * phi**(2*this%n - 2)
                end if
            
            case(2) ! AxionEDE
                theta = phi/this%f
                if (deriv==0) then
                   Vofphi = units * this%m**2 * this%f**2 * (1 - cos(theta))**this%n
                else if (deriv ==1) then
                    Vofphi = units*this%m**2*this%f*this%n*(1 - cos(theta))**(this%n-1)*sin(theta)
                else if (deriv ==2) then
                    costheta = cos(theta)
                    Vofphi = units*this%m**2*this%n*(1 - costheta)**(this%n-1)*(this%n*(1+costheta) -1)
                end if
        end select
    end function TEarlyDarkEnergy_VofPhi

    ! JVR - NOT IMPLEMENTED but I think OK???
    subroutine TEarlyDarkEnergyPlusW_Init(this, State)
        use Powell
        class(TEarlyDarkEnergyPlusW), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State
        real(dl) aend, afrom
        integer, parameter ::  NumEqs = 2
        real(dl) c(24),w(NumEqs,9), y(NumEqs)
        integer ind, i, ix
        real(dl), parameter :: splZero = 0._dl
        real(dl) lastsign, da_osc, last_a, a_c
        real(dl) initial_phi, initial_phidot, a2
        real(dl), dimension(:), allocatable :: sampled_a, phi_a, phidot_a, fde
        integer npoints, tot_points, max_ix
        logical has_peak
        real(dl) fzero, xzero
        integer iflag, iter
        Type(TTimer) :: Timer
        Type(TNEWUOA) :: Minimize
        real(dl) log_params(2), param_min(2), param_max(2)
        logical :: outputbackground = .false.
        character(len=50) :: outputbackgroundfilename = "ede_background.dat"
        real(dl) :: fmatter, frad, grhorad, grhomatter, grhoquint, w_phi, quint_to_w_frac

        ! JVR - The Init subroutine is used to initialize all variables in the class that are
        ! needed to compute the background quantities and the perturbation equations
        ! Thus, we need to evolve numerically the field, since rho_phi(a) = phi'(a)^2 / 2a^2 + V(phi(a))
        ! All other components have been initialized at this point

        call this%TScalarField%Init(State)

        if (this%use_zc) then
            call zc_fde_adjust(this)
        end if

        ! Setting initial field value
        if (this%which_potential == 1) then
            initial_phi = this%initial_phi
        else
            initial_phi = this%theta_i*this%f
        end if

        ! Use log spacing in a up to max_a_log, then linear. Switch where step matches
        this%dloga = (-this%log_astart)/(this%npoints-1)
        this%max_a_log = 1.d0/this%npoints/(exp(this%dloga)-1)
        npoints = (log(this%max_a_log)-this%log_astart)/this%dloga + 1

        if (allocated(this%phi_a)) then
            deallocate(this%phi_a,this%phidot_a)
            deallocate(this%ddphi_a,this%ddphidot_a, this%sampled_a)
        end if
        allocate(phi_a(npoints),phidot_a(npoints), sampled_a(npoints), fde(npoints))

        y(1)=initial_phi
        initial_phidot =  this%astart*this%phidot_start(initial_phi)
        y(2)= initial_phidot*this%astart**2

        phi_a(1)=y(1)
        phidot_a(1)=y(2)/this%astart**2
        sampled_a(1)=this%astart
        da_osc = 1
        last_a = this%astart
        max_ix = 0

        ! ----------------- Background Integration -----------------------

        if (outputbackground .eqv. .true.) then
            open(unit = 50, file = outputbackgroundfilename, form = 'formatted', status = 'replace')
            write(50, *) "a     phi     phidot      fde     w     fmatter     frad     grhorad     grhomatter     grhoede   grhow"
        end if

        ind=1
        afrom=this%log_astart
        do i=1, npoints-1
            aend = this%log_astart + this%dloga*i
            ix = i+1
            sampled_a(ix)=exp(aend)
            a2 = sampled_a(ix)**2

            quint_to_w_frac = (phidot_a(i)**2/(2*a2) + this%Vofphi(phi_a(i), 0))/(this%grho_w(sampled_a(i))/a2**2)
            if (quint_to_w_frac < 1e-4) then
                phi_a(ix) = 0
                phidot_a(ix) = 0
                cycle
            end if 
            call dverk(this,NumEqs,EvolveBackgroundLog,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)
            if (.not. this%check_error(exp(afrom), exp(aend))) return
            call EvolveBackgroundLog(this,NumEqs,aend,y,w(:,1))
            phi_a(ix)=y(1)
            phidot_a(ix)=y(2)/a2
            
            ! probing oscillations
            !if (i==1) then
            !    lastsign = y(2)
            !elseif (y(2)*lastsign < 0) then
            !    !derivative has changed sign. Use to probe any oscillation scale:
            !    da_osc = min(da_osc, exp(aend) - last_a)
            !    last_a = exp(aend)
            !    lastsign= y(2)
            !end if

            !Define fde as ratio of early dark energy density to total
            
            !if (max_ix==0 .and. ix > 2 .and. fde(ix)< fde(ix-1)) then
            !    max_ix = ix-1
            !end if
            !if (sampled_a(ix)*(exp(this%dloga)-1)*this%min_steps_per_osc > da_osc) then
            !    !Step size getting too big to sample oscillations well
            !    exit
            !end if

            ! Outputting background data to file
            
            if (outputbackground .eqv. .true.) then
                fde(ix) = 1/((this%state%grho_no_de(sampled_a(ix)) +  this%frac_lambda0*this%State%grhov*a2**2) &
                /(a2*(0.5d0* phidot_a(ix)**2 + a2*this%Vofphi(y(1),0))) + 1)
                grhorad = this%state%grho_radiation(sampled_a(ix))
                grhomatter = this%state%grho_matter(sampled_a(ix))
                grhoquint = a2*(phidot_a(ix)**2/2 + a2*this%Vofphi(phi_a(ix), 0))
                grhow = a2**2 * (this%state%grhov * a)
                write(50, *) sampled_a(ix), phi_a(ix), phidot_a(ix), fde(ix), w_phi, fmatter, frad, grhorad, grhomatter, grhoquint, this%grho_w(sampled_a(ix))
            end if
        end do
        ! Do remaining steps with linear spacing in a, trying to be small enough
        this%npoints_log = ix
        this%max_a_log = sampled_a(ix)
        this%da = min(this%max_a_log *(exp(this%dloga)-1), &
            da_osc/this%min_steps_per_osc, (1- this%max_a_log)/(this%npoints-this%npoints_log))
        this%npoints_linear = int((1- this%max_a_log)/ this%da)+1
        this%da = (1- this%max_a_log)/this%npoints_linear

        tot_points = this%npoints_log+this%npoints_linear
        allocate(this%phi_a(tot_points),this%phidot_a(tot_points))
        allocate(this%ddphi_a(tot_points),this%ddphidot_a(tot_points))
        allocate(this%sampled_a(tot_points), this%fde(tot_points), this%ddfde(tot_points))
        this%sampled_a(1:ix) = sampled_a(1:ix)
        this%phi_a(1:ix) = phi_a(1:ix)
        this%phidot_a(1:ix) = phidot_a(1:ix)
        this%sampled_a(1:ix) = sampled_a(1:ix)
        this%fde(1:ix) = fde(1:ix)

        ind=1
        afrom = this%max_a_log
        do i=1, this%npoints_linear
            ix = this%npoints_log + i
            aend = this%max_a_log + this%da*i
            a2 =aend**2
            this%sampled_a(ix)=aend

            quint_to_w_frac = (this%phidot_a(ix-1)**2/(2*a2) + this%Vofphi(this%phi_a(ix-1), 0))/(this%grho_w(sampled_a(ix-1))/a2**2)
            if (quint_to_w_frac < 1e-4) then
                this%phi_a(ix) = 0
                this%phidot_a(ix) = 0
                cycle
            end if
            call dverk(this,NumEqs,EvolveBackground,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)
            if (.not. this%check_error(afrom, aend)) return
            call EvolveBackground(this,NumEqs,aend,y,w(:,1))
            this%phi_a(ix)=y(1)
            this%phidot_a(ix)=y(2)/a2

            this%fde(ix) = 1/((this%state%grho_no_de(aend) +  this%frac_lambda0*this%State%grhov*a2**2) &
                /(a2*(0.5d0* this%phidot_a(ix)**2 + a2*this%Vofphi(y(1),0))) + 1)
            if (max_ix==0 .and. this%fde(ix)< this%fde(ix-1)) then
                max_ix = ix-1
            end if

            if (outputbackground .eqv. .true.) then
                grhorad = this%state%grho_radiation(this%sampled_a(ix))
                grhomatter = this%state%grho_matter(this%sampled_a(ix))
                grhoquint = a2*(this%phidot_a(ix)**2/2 + a2*this%Vofphi(this%phi_a(ix), 0))
                write(50, *) this%sampled_a(ix), this%phi_a(ix), this%phidot_a(ix), this%fde(ix), w_phi, fmatter, frad, grhorad, grhomatter, grhoquint, this%grho_w(this%sampled_a(ix))
            end if
        end do

        if (outputbackground .eqv. .true.) then
            close(50)
        end if

        call spline(this%sampled_a,this%phi_a,tot_points,splZero,splZero,this%ddphi_a)
        call spline(this%sampled_a,this%phidot_a,tot_points,splZero,splZero,this%ddphidot_a)
        call spline(this%sampled_a,this%fde,tot_points,splZero,splZero,this%ddfde)
        
        ! JVR - do I need this?
        !has_peak = .false.
        !if (max_ix >0) then
        !    ix = max_ix
        !    has_peak = this%fde_peak(a_c, this%sampled_a(ix), this%sampled_a(ix+1), this%fde(ix), &
        !        this%fde(ix+1), this%ddfde(ix), this%ddfde(ix+1))
        !    if (.not. has_peak) then
        !        has_peak = this%fde_peak(a_c, this%sampled_a(ix-1), this%sampled_a(ix), &
        !            this%fde(ix-1), this%fde(ix), this%ddfde(ix-1), this%ddfde(ix))
        !    end if
        !end if
        !if (has_peak) then
        !    this%zc = 1/a_c-1
        !    this%fde_zc = this%fdeAta(a_c)
        !else
        !    if (this%DebugLevel>0) write(*,*) 'TEarlyDarkEnergy: NO PEAK '
        !    this%zc = -1
        !end if
        !if (this%DebugLevel>0) then
        !    write(*,*) 'TEarlyDarkEnergy zc, fde used', this%zc, this%fde_zc
        !end if
    end subroutine TEarlyDarkEnergyPlusW_Init

    logical function check_error(this, afrom, aend)
        class(TEarlyDarkEnergy) :: this
        real(dl) afrom, aend

        if (global_error_flag/=0) then
            write(*,*) 'TEarlyDarkEnergy error integrating', afrom, aend
            write(*,*) this%n, this%f, this%m, this%theta_i
            stop
            check_error = .false.
            return
        end if
        check_error= .true.
    end function check_error

    logical function fde_peak(this, peak, xlo, xhi, Flo, Fhi, ddFlo, ddFhi)
        class(TEarlyDarkEnergy) :: this
        real(dl), intent(out) :: peak
        real(dl) Delta
        real(dl), intent(in) :: xlo, xhi, ddFlo, ddFhi,Flo, Fhi
        real(dl) a, b, c, fac

        !See if derivative has zero in spline interval xlo .. xhi

        Delta = xhi - xlo

        a = 0.5_dl*(ddFhi-ddFlo)/Delta
        b = (xhi*ddFlo-xlo*ddFhi)/Delta
        c = (Fhi-Flo)/Delta+ Delta/6._dl*((1-3*xhi**2/Delta**2)*ddFlo+(3*xlo**2/Delta**2-1)*ddFhi)
        fac = b**2-4*a*c
        if (fac>=0) then
            fac = sqrt(fac)
            peak = (-b + fac)/2/a
            if (peak >= xlo .and. peak <= xhi) then
                fde_peak = .true.
                return
            else
                peak = (-b - fac)/2/a
                if (peak >= xlo .and. peak <= xhi) then
                    fde_peak = .true.
                    return
                end if
            end if
        end if
        fde_peak = .false.
    end function fde_peak

    function match_zc(this, logm)
        class(TEarlyDarkEnergy), intent(inout) :: this
        real(dl), intent(in) :: logm
        real(dl) match_zc, zc, fde_zc

        this%m = exp(logm)
        call this%calc_zc_fde(zc, fde_zc)
        match_zc = zc - this%zc

        end function match_zc

        function match_fde(this, logf)
        class(TEarlyDarkEnergy), intent(inout) :: this
        real(dl), intent(in) :: logf
        real(dl) match_fde, zc, fde_zc

        this%f = exp(logf)
        call this%calc_zc_fde(zc, fde_zc)
        match_fde = fde_zc - this%fde_zc
    end function match_fde

    function match_fde_zc(this, x)
        class(TEarlyDarkEnergy) :: this
        real(dl), intent(in) :: x(:)
        real(dl) match_fde_zc, zc, fde_zc

        if (this%which_potential == 2) then
            this%f = exp(x(1))
            this%m = exp(x(2))
        else
            this%V0 = exp(x(1))
            this%initial_phi = exp(x(2))
        end if
        
        call this%calc_zc_fde(zc, fde_zc)

        match_fde_zc = (log(this%fde_zc)-log(fde_zc))**2 + (log(zc)-log(this%zc))**2
        if (this%DebugLevel>1) then
            write(*,*) 'search f, m, zc, fde_zc, chi2', this%f, this%m, zc, fde_zc, match_fde_zc
        end if
    end function match_fde_zc

    subroutine calc_zc_fde(this, z_c, fde_zc)
        class(TEarlyDarkEnergy), intent(inout) :: this
        real(dl), intent(out) :: z_c, fde_zc
        real(dl) aend, afrom
        integer, parameter ::  NumEqs=2
        real(dl) c(24),w(NumEqs,9), y(NumEqs)
        integer ind, i, ix
        real(dl), parameter :: splZero = 0._dl
        real(dl) a_c
        real(dl) initial_phi, initial_phidot, a2
        real(dl), dimension(:), allocatable :: sampled_a, fde, ddfde
        integer npoints, max_ix
        logical has_peak
        real(dl) a0, b0, da

        ! Get z_c and f_de(z_c) where z_c is the redshift of (first) peak of f_de (de energy fraction)
        ! Do this by forward propagating until peak, then get peak values by cubic interpolation

        if (this%which_potential == 1) then
            initial_phi = this%initial_phi
        else
            initial_phi = this%theta_i*this%f
        end if

        this%log_astart = log(this%astart)
        this%dloga = (-this%log_astart)/(this%npoints-1)

        npoints = this%npoints
        allocate(sampled_a(npoints), fde(npoints), ddfde(npoints))

        y(1) = initial_phi
        initial_phidot = this%astart*this%phidot_start(initial_phi)
        y(2) = initial_phidot*this%astart**2
        sampled_a(1) = this%astart
        max_ix = 0
        ind = 1
        afrom = this%log_astart
        do i = 1, npoints-1
            aend = this%log_astart + this%dloga*i
            ix = i+1
            sampled_a(ix)=exp(aend)
            a2 = sampled_a(ix)**2
            call dverk(this,NumEqs,EvolveBackgroundLog,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)
            if (.not. this%check_error(exp(afrom), exp(aend))) return
            call EvolveBackgroundLog(this,NumEqs,aend,y,w(:,1))
            fde(ix) = 1/((this%state%grho_no_de(sampled_a(ix)) +  this%frac_lambda0*this%State%grhov*a2**2) &
                /((0.5d0*y(2)**2/a2 + a2**2*this%Vofphi(y(1),0))) + 1)
            if (max_ix == 0 .and. ix > 2 .and. fde(ix) < fde(ix-1)) then
                max_ix = ix-1
            end if
            if (max_ix/=0 .and. ix > max_ix+4) exit
        end do

        call spline(sampled_a,fde,ix,splZero,splZero,ddfde)
        has_peak = .false.
        if (max_ix >0) then
            has_peak = this%fde_peak(a_c, sampled_a(max_ix), sampled_a(max_ix+1), fde(max_ix), &
                fde(max_ix+1), ddfde(max_ix), ddfde(max_ix+1))
            if (.not. has_peak) then
                has_peak = this%fde_peak(a_c, sampled_a(max_ix-1), sampled_a(max_ix), &
                    fde(max_ix-1), fde(max_ix), ddfde(max_ix-1), ddfde(max_ix))
            end if
        end if
        if (has_peak) then
            z_c = 1/a_c-1
            ix = int((log(a_c)-this%log_astart)/this%dloga)+1
            da = sampled_a(ix+1) - sampled_a(ix)
            a0 = (sampled_a(ix+1) - a_c)/da
            b0 = 1 - a0
            fde_zc=b0*fde(ix+1) + a0*(fde(ix)-b0*((a0+1)*ddfde(ix)+(2-a0)*ddfde(ix+1))*da**2/6._dl)
        else
            write(*,*) 'calc_zc_fde: NO PEAK'
            z_c = -1
            fde_zc = 0
        end if
    end subroutine calc_zc_fde

    function fdeAta(this,a)
        class(TEarlyDarkEnergy) :: this
        real(dl), intent(in) :: a
        real(dl) fdeAta, aphi, aphidot, a2

        call this%ValsAta(a, aphi, aphidot)
        a2 = a**2
        fdeAta = 1/((this%state%grho_no_de(a) +  this%frac_lambda0*this%State%grhov*a2**2) &
            /(a2*(0.5d0* aphidot**2 + a2*this%Vofphi(aphi,0))) + 1)
    end function fdeAta

    subroutine zc_fde_adjust(this)
        ! Find underlying parameters (V0, initial_phi) for Monomial or (m,f) for AxionEDE
        ! to give specified zc and fde_zc (peak early dark energy fraction)
        ! Input parameters are used as starting values for search, which is done by brute force
        ! (so should generalize easily, but not optimized for this specific potential)
        use Powell
        Type(TTimer) :: Timer
        Type(TNEWUOA) :: Minimize
        real(dl) fzero, xzero
        integer iflag, iter
        type(TEarlyDarkEnergy), intent(inout) :: this
        real(dl) :: log_params(2)

        ! 1 - log_params(*) will be the initial guess for the algorithm
        if (this%which_potential == 1) then
            log_params(1) = log(this%V0)
            log_params(2) = log(this%initial_phi)
        else
            log_params(1) = log(this%f)
            log_params(2) = log(this%m)
        end if

        if (.false.) then
            ! Can just iterate linear optimizations when nearly orthogonal
            call Timer%Start()
            do iter = 1, 2
                call brentq(this,match_fde,log(0.01_dl),log(10._dl), 1d-3,xzero,fzero,iflag)
                if (iflag/=0) print *, 'BRENTQ FAILED f'
                this%f = exp(xzero)
                print *, 'match to m, f =', this%m, this%f, fzero
                call brentq(this,match_zc,log(1d-55),log(1d-52), 1d-3,xzero,fzero,iflag)
                if (iflag/=0) print *, 'BRENTQ FAILED m'
                this%m = exp(xzero)
                print *, 'match to m, f =', this%m, this%f, fzero
                call this%calc_zc_fde(fzero, xzero)
                print *, 'matched outputs', fzero, xzero
            end do
            call Timer%WriteTime('Timing for fitting')
        end if
        
        ! 2 - Could be useful to time
        if (this%DebugLevel>0) call Timer%Start()
        !Minimize in log f, log m
        ! param_min(1) = log(0.001_dl)
        ! param_min(2) = log(1d-58)
        ! param_max(1) = log(1e5_dl)
        ! param_max(2) = log(1d-50)
        ! if (Minimize%BOBYQA(this, match_fde_zc, 2, 5, log_params,param_min, &
        !           param_max, 0.8_dl,1e-4_dl,this%DebugLevel,2000)) then

        ! 3 - NEWUOA is an Unconstrained Optimization Algorithm.
        ! Check its implementation in PowellMinimize.f90
        ! log_params will be updated with values that give desired zc/fde
        ! Arguments: NEWUOA(this - class instance where values are stored;
        !                   match_fde_zc - the function to be optimized;
        !                   2 - the number of variables
        !                   5 - interpolation conditions
        !                   log_params - the initial values)
        if (Minimize%NEWUOA(this, match_fde_zc, 2, 5, log_params,&
            0.8_dl,this%zcfdeprec,this%DebugLevel,500)) then

            if (Minimize%Last_bestfit > 1e-3) then
                global_error_flag = error_darkenergy
                global_error_message = 'TEarlyQuintessence ERROR converging solution for fde, zc'
                write(*,*) 'last-bestfit = ', Minimize%Last_bestfit
                return
            end if

            if (this%which_potential == 1) then
                this%V0 = exp(log_params(1))
                this%initial_phi = exp(log_params(2))
            else
                this%f = exp(log_params(1))
                this%m = exp(log_params(2))
            end if
            
            if (this%DebugLevel>0) then
                call this%calc_zc_fde(fzero, xzero)
                write(*,*) 'matched outputs Bobyqa zc, fde = ', fzero, xzero
            end if
        else
            global_error_flag = error_darkenergy
            global_error_message = 'TEarlyQuintessence ERROR finding solution for fde, zc'
            return
        end if
        
        if (this%DebugLevel>0) call Timer%WriteTime('Timing for parameter fitting')
    end subroutine zc_fde_adjust

    subroutine TEarlyDarkEnergyPlusW_ReadParams(this, Ini)
        use IniObjects
        class(TEarlyDarkEnergy) :: this
        class(TIniFile), intent(in) :: Ini

        call this%TDarkEnergyModel%ReadParams(Ini)

        ! Should CAMB output phi(a)? In which file?
        !this%outputbackground = Ini%Read_Logical('output_background')
        !this%outputbackgroundfilename = trim(Ini%Read_String('output_root')) // trim('_bg.dat')

        ! Reading the w_const eqn of state
        this%w_const = Ini%Read_Double('w_const', -0.999_dl)

        ! Reading EDE potential type and parameters
        this%which_potential = Ini%Read_Int('which_potential', 1)
        this%m = Ini%Read_Double('mass')
        this%f = Ini%Read_Double('freq')
        this%n = Ini%Read_Double('powerlaw', 3._dl)
        this%V0 = Ini%Read_Double('V0', 1d-117)
        this%initial_phi = Ini%Read_Double('initial_phi', 1._dl)
    end subroutine TEarlyDarkEnergyPlusW_ReadParams

    function TEarlyDarkEnergyPlusW_PythonClass()
        character(LEN=:), allocatable :: TEarlyDarkEnergy_PythonClass
        TEarlyDarkEnergy_PythonClass = 'wEDE'
    end function TEarlyDarkEnergyPlusW_PythonClass

    subroutine TEarlyDarkEnergyPlusW_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TEarlyDarkEnergy), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TEarlyDarkEnergyPlusW_SelfPointer

    end module wEDE