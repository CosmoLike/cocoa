/* cosmo.h
*  =======
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

#ifndef COSMOLOGY_H
#define COSMOLOGY_H
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

class Cosmology{
	public:
		double cosmo[8], cosmo_tf[8];
		double Omega_gamma_0, Omega_nu_0, Omega_DE_0, rho_crit, T_gamma_0, T_nu_0;

		/* Ranges for cosmological parameters */
        const double minima[8] = {0.04, 0.24, 0.00, 0.92, 0.61, -1.3, -0.7, 1.7e-9};
        const double maxima[8] = {0.06, 0.40, 0.15, 1.00, 0.73, -0.7, 0.7, 2.5e-9};

		/* Member functions */
		Cosmology(double Omega_b, double Omega_m, double Sum_m_nu, double n_s, double h, double w_0, double w_a, double A_s);
		~Cosmology();
		void read_from_file(char *filename);
		void print_cosmo();
		void print_cosmo_tf();
	    double compute_step_number(double z);
		static double rho_nu_i_integrand(double p, void * params);
		static double a2t_integrand(double a, void * params);

	private:
		/* Private members */
		const int nSteps, nTable;
		double t0, t10, Delta_t, Neff, H0;

		typedef struct {
            double mnu_i;
            double a;
            Cosmology * csm_instance;
		} rho_nu_parameters;
		typedef struct {
            Cosmology * csm_instance;
		} a2t_parameters;

		gsl_integration_workspace *gsl_wsp;
		gsl_interp_accel *acc;
    	gsl_spline *z2nStep_spline;

		/* Private member functions*/
		void isoprob_tf();
		void check_parameter_ranges();
		void compute_z2nStep_spline();

		double Hubble(double a);
		double Omega_matter(double a);
		double Omega_gamma(double a);
		double Omega_nu(double a);
		double Omega_DE(double a);
		double T_gamma(double a);
		double T_nu(double a);
		double a2t(double a);
		double a2Hubble(double a);
};
#endif
