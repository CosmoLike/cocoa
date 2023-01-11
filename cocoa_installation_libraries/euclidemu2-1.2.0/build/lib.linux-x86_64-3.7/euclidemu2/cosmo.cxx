/* cosmo.cxx
*  =========
*  This file is part of EuclidEmulator2
*  Copyright (c) 2020 Mischa Knabenhans
*
*  EuclidEmulator2 is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  EuclidEmulator2 is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include "cosmo.h"
#include "units_and_constants.h"

#define EPSCOSMO 1e-6
#define LIMIT 1000

#ifndef PRINT_FLAG
#define PRINT_FLAG 1
#endif

//using namespace planck_units;
using namespace SI_units;

/* CONSTRUCTOR */
Cosmology::Cosmology(double Omega_b, double Omega_m, double Sum_m_nu, double n_s, double h, double w_0, double w_a, double A_s):
	Neff(3.046),
	nTable(101),
	nSteps(101)
	{
	// GSL interpolation object for redshift-to-nStep conversion:
    // The spline is depends on the cosmology but takes redshift as an
    // argument. For this very reason, the spline can computed once when
    // a Cosmology class object is instantiated. From then on, the spline
    // can be evaluated as often as necessary until the Cosmology object
    // is destroyed.
    gsl_wsp = gsl_integration_workspace_alloc(LIMIT);
    acc = gsl_interp_accel_alloc();

	this->cosmo[0] = Omega_b;
    this->cosmo[1] = Omega_m;
    this->cosmo[2] = Sum_m_nu;
    this->cosmo[3] = n_s;
    this->cosmo[4] = h;
    this->cosmo[5] = w_0;
    this->cosmo[6] = w_a;
	this->cosmo[7] = A_s;
	// Compute the present day temperatures of the photon (gamma) and the
	// neutrino (nu) fluids:
    this->T_gamma_0 = T_gamma(1.0);
    this->T_nu_0 = T_nu(1.0);
	// Compute the actual critical density of the Universe in SI units for
    // the given cosmology:
	this->rho_crit = rho_crit_over_h2 * pow(this->cosmo[4],2);
	// Now we can compute the present day values of the density parameters
	// of the neutrino, the photon and the dark energy (DE) fluid:
	this->Omega_nu_0 = Omega_nu(1.0);
	this->Omega_gamma_0 = Omega_gamma(1.0);
	this->Omega_DE_0 = 1 - (this->cosmo[1] + this->Omega_gamma_0 + this->Omega_nu_0);

	//printf("Cosmological parameters assigned successfully\n");
	// Prepare for spline interpolation of z --> nStep mapping:
	t0  = a2t(1.0); // proper time at z = 0 (or equivalently a = 1)
    t10 = a2t(1.0/(10+1)); // proper time at z = 10 (or equivalently a = 0.090909...)
    Delta_t = (t0-t10)/(nSteps-1);

    z2nStep_spline = gsl_spline_alloc(gsl_interp_cspline, nTable);
	// In the next line the spline is actually being computed:
	compute_z2nStep_spline();

	/* Member functions */
	check_parameter_ranges();
	isoprob_tf();
	//convert_to_pkdgrav();
	if(PRINT_FLAG){
	  print_cosmo();
	  print_cosmo_tf();
	}
}

/* DESTRUCTOR */
Cosmology::~Cosmology(){
		gsl_spline_free(z2nStep_spline);
		gsl_integration_workspace_free(gsl_wsp);
		gsl_interp_accel_free(acc);
}

/* CHECK PARAMETER RANGES */
void Cosmology::check_parameter_ranges(){
    for(int i=0 ; i<8; i++) {
        if(cosmo[i] < minima[i] || cosmo[i] > maxima[i] ) {
			std::cout << "Parameter " << i << " is outside allowed range:" << std::endl;
			std::cout << "cosmo[" << i << "] = " << cosmo[i] << std::endl;
			exit(1);
		}
		if(i==5) i++;
    }
	// Check w_a separately because values >0.5 are not
	// allowed although they are inside the range
	if(cosmo[6] < -0.7 || cosmo[6] > 0.5 ) {
        std::cout << "Parameter w_a is outside allowed range:" << std::endl;
        std::cout << " w_a = " << cosmo[6] << std::endl;
        exit(1);
    }
}

/* ISOPROBALISTIC TRANSFORMATION TO UNIT HYPERCUBE */
void Cosmology::isoprob_tf(){
	for (int i=0; i<8; i++){
		cosmo_tf[i] = 2*(cosmo[i] - minima[i])/(maxima[i] - minima[i]) - 1.0;
	}
}

/* COMPUTE a(t) TABLE FOR GIVEN COSMOLOGY */
double Cosmology::Omega_matter(double a){
	/* This function computes the matter density parameter *\
	\* in the Universe at scale factor a.                  */
	return this->cosmo[1] / (a*a*a);
}

double Cosmology::Omega_gamma(double a){
	/* This function computes the photon density parameter in the Universe at scale factor a */

	// Present day photon density corresponding to T_gamma (also in SI units):
	double rho_gamma_0 =  M_PI * M_PI / 15.0
						* pow(kB,4) / (pow(hbar,3)*pow(c,5))
						* pow(Tgamma,4);

	// Scale the photon density and return result
	this->Omega_gamma_0 = rho_gamma_0 / this->rho_crit;
    return this->Omega_gamma_0 / (a*a*a*a);
}

double Cosmology::T_gamma(double a){
	return Tgamma/a;
}

double Cosmology::T_nu(double a){
	return pow(this->Neff/3.0,0.25)*pow(4.0/11.0,1./3.)*Tgamma/a;
}

double Cosmology::rho_nu_i_integrand(double p, void * params){
	/* This function provides the integrand for the function Cosmology::Omega_nu */
	rho_nu_parameters *rho_nu_pars = reinterpret_cast<rho_nu_parameters *>(params);

	double T_nu = rho_nu_pars->csm_instance->T_nu_0;
	double mnui = rho_nu_pars->mnu_i;
	double p2 = p*p;
	double y = p2 / ( exp(c/kB/T_nu * rho_nu_pars->a * p) + 1 );

	return y * sqrt( mnui*mnui*c*c*c*c + p2*c*c);
}

double Cosmology::Omega_nu(double a){
	/* This function computes the neutrino density parameter in the Universe at scale factor a */
	/* REMARK: The assumption of a degenerate neutrino hierarchy is hardcoded. */
	//printf("Entering Omega_nu...\n");
	rho_nu_parameters rho_nu_pars;
	gsl_function F;
	double rho_nu_i, error;

	//the letter "i" in "rho_nu_i" emphasizes that this variable
	//contains the density of only ONE neutrino species

	// Compute the neutrino density of ONE neutrino species of mass m_nu
	// (according to eq. 2.69 in my thesis). Remember that this integral is
	// an integral over the momentum p, so we start by defining an uper bound
	// "pmax" for the integration:
	double pmax = 0.004/a * eV/c;
	double prefactor = 1.0/(pow(M_PI,2) * pow(hbar,3) * pow(c,2));

	// Remember that we always assume degenerate neutrino hierarchy. This is
	// why the mass of a single neutrino species is a third of Sum m_nu (which
	// is stored in this->cosmo[2].
	rho_nu_pars.mnu_i = this->cosmo[2]/3.0 * eV/pow(c,2);
	rho_nu_pars.a = a;
	rho_nu_pars.csm_instance = this;

	//printf("Setting up GSL integrator...\n");
	F.function = &rho_nu_i_integrand;
	F.params = &rho_nu_pars;
	gsl_integration_qag(&F, 0.0, pmax, 0.0,
						EPSCOSMO, LIMIT, GSL_INTEG_GAUSS61,
						gsl_wsp, &rho_nu_i, &error);

	//printf("Omega_nu integration relative error = %.15e\n", error/rho_nu_i);
	rho_nu_i *= prefactor;

	//printf("GSL integration complete.\n");
    // NOTICE: The prefactor 3 comes from the fact that we ALWAYS consider three
    // equal-mass neutrinos while rho_nu_i computes the density of ONE neutrino species
	return 3*rho_nu_i / this->rho_crit;
}

double Cosmology::Omega_DE(double a){
	/* This function computes the dark energy (DE) density parameter of   *\
	|* the Universe. This computation is based on the flatness condition: *|
	|*																	  *|
	|*			Omega_matter + Omega_gamma + Omega_nu + Omega_DE = 1      *|
	|*																	  *|
	\* The assumption of a flat Universe is thus hardcoded.               */
	double w_0 = this->cosmo[5];
	double w_a = this->cosmo[6];

	// The following is equation (2.67) in my thesis
	return this->Omega_DE_0 * pow(a, -3.0 *(1 + w_0 + w_a)) * exp(-3*(1-a)*w_a);
}

double Cosmology::a2Hubble(double a){
	/* This function computes the Hubble parameter at scale factor a */
	H0 = 100*this->cosmo[4] * kilometer/second/Mpc;
	return H0 * sqrt( Cosmology::Omega_matter(a)
					+ Cosmology::Omega_gamma(a)
					+ Cosmology::Omega_nu(a)
					+ Cosmology::Omega_DE(a)
					);
}

double Cosmology::a2t_integrand(double lna, void *params){
	/* This function provides the integrand for the GSL *\
	\* integration in Cosmology::a2t below.             */
	a2t_parameters *a2t_pars = reinterpret_cast<a2t_parameters *>(params);
	return 1./(a2t_pars->csm_instance->a2Hubble(exp(lna)));
}

double Cosmology::a2t(double a){
	/* This function converts a scale factor a to a proper time. *\
	|* The integration is solved by GSL and the integrand is     *|
	\* by Cosmology::a2t_integrand.                              */
	a2t_parameters a2t_params;
	gsl_function F;
	double result, error;
	a2t_params.csm_instance = this;
	F.function = &a2t_integrand;
	F.params = &a2t_params;
	gsl_integration_qag(&F, -15, log(a), 0.0,
						EPSCOSMO, LIMIT, GSL_INTEG_GAUSS61,
						this->gsl_wsp, &result, &error);

	//printf("a2t integration relative error = %.15e\n", error/result);
	return result;
}

void Cosmology::compute_z2nStep_spline(){
	/* This function creates and interpolates the array of (z,nStep)-tuples *\
	|* computed in "Cosmology::compute_time_lookup_table". It returns a     *|
	\* gsl_spline function than can be readily evaluated.                   */

	// Step 1: Setup the array
	double t_current, z;
	double avec[101]={0};
	double frac_nStep[101]={0};

	double z10 = 10.0;
	for(int idx=0; idx<this->nTable; idx++){
		// Loop through redshifts: z \in {10.0, 9.9, ..., 0.1, 0.0}
		z = z10 - idx*0.1;
		// Convert z to a (remember that GSL interpolator expects
		// x-values to be in ascending order)
		avec[idx] = 1.0/(z+1.0);
		// Convert a to t
		t_current = Cosmology::a2t(avec[idx]);
		// Convert t to nStep (fractional)
		frac_nStep[idx] = (t_current - this->t10)/this->Delta_t;
	}

	// Some sanity checks:
	assert(abs(frac_nStep[0]) < EPSCOSMO);
	assert(abs(frac_nStep[this->nTable-1] - (this->nSteps-1)) < EPSCOSMO);

	// Step 2: Interpolate the array
	gsl_spline_init(this->z2nStep_spline, avec, frac_nStep, this->nTable);
}

/* INTERPOLATE OUTPUT STEP */
double Cosmology::compute_step_number(double z){
	/* This function evaluates the spline mapping *\
	\* a redshift to a (fractional) output step.  */
	if(fabs(z) < EPSCOSMO){
		return 100.0;
	}
	else{
		// Now simply evaluate the precomputed spline
		return gsl_spline_eval(this->z2nStep_spline, 1./(z+1.), this->acc);
	}
}

/* PRINT FUNCTIONS */
void Cosmology::print_cosmo(){
	std::cout << "The current cosmology is defined as:" << std::endl;
	std::cout << "\tOmega_b = " << this->cosmo[0] << std::endl;
	std::cout << "\tOmega_m = " << this->cosmo[1] << std::endl;
	std::cout << "\tSum m_nu = " << this->cosmo[2] << " eV" << std::endl;
	std::cout << "\tn_s = " << this->cosmo[3] << std::endl;
	std::cout << "\th = " << this->cosmo[4] << std::endl;
	std::cout << "\tw_0 = " << this->cosmo[5] << std::endl;
	std::cout << "\tw_a = " << this->cosmo[6] << std::endl;
	std::cout << "\tA_s = " << this->cosmo[7] << std::endl;
	std::cout << std::endl;
	std::cout << "Fixed parameters (hardcoded):" << std::endl;
	std::cout << "\tT_gamma = " << this->T_gamma_0 << " K" << std::endl;
    std::cout << std::endl;
	std::cout << "Derived parameters:" << std::endl;
    std::cout << "\tT_nu = " << this->T_nu_0 << " K" << std::endl;
    std::cout << "\tOmega_nu = " << std::setprecision(12) << this->Omega_nu_0 << std::endl;
	std::cout << "\tOmega_gamma = " << this->Omega_gamma_0 << std::endl;
	std::cout << "\tOmega_DE = " << this->Omega_DE_0 << std::endl;
	std::cout << "\trho_crit = " << this->rho_crit/(Msol/pow(Mpc,3)) << " h^2 M_sol/(Mpc^3)" << std::endl;
	std::cout << std::endl;
}

void Cosmology::print_cosmo_tf(){
	std::cout << "The current cosmology is mapped to:" << std::endl;
    std::cout << std::endl;
    std::cout << "\ttf(Omega_b) = " << this->cosmo_tf[0] << std::endl;
    std::cout << "\ttf(Omega_m) = " << this->cosmo_tf[1] << std::endl;
    std::cout << "\ttf(Sum m_nu) = " << this->cosmo_tf[2] << std::endl;
    std::cout << "\ttf(n_s) = " << this->cosmo_tf[3] << std::endl;
    std::cout << "\ttf(h) = " << this->cosmo_tf[4] << std::endl;
    std::cout << "\ttf(w_0) = " << this->cosmo_tf[5] << std::endl;
    std::cout << "\ttf(w_a) = " << this->cosmo_tf[6] << std::endl;
    std::cout << "\ttf(A_s) = " << this->cosmo_tf[7] << std::endl;
	std::cout << std::endl;
}

/* READ COSMOLOGY FROM FILE  */
void read_from_file(char *filename){}
