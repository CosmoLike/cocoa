#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "utils.h"

double gamma_lanczos_real_cfft(double z) {
/* Lanczos coefficients for g = 7 */
	static double p[] = {
		0.99999999999980993227684700473478,
		676.520368121885098567009190444019,
		-1259.13921672240287047156078755283,
		771.3234287776530788486528258894,
		-176.61502916214059906584551354,
		12.507343278686904814458936853,
		-0.13857109526572011689554707,
		9.984369578019570859563e-6,
		1.50563273514931155834e-7};

	if(z < 0.5) {return M_PI / (sin(M_PI*z)*gamma_lanczos_real_cfft(1. - z));}
	z -= 1;
	double x = p[0];
	for(int n = 1; n < 9; n++){ x += p[n] / (z + (double)(n));}

	double t = z + 7.5;
	return sqrt(2*M_PI) * pow(t, z+0.5) * exp(-t) * x;
}

double lngamma_lanczos_real_cfft(double z) {
/* Lanczos coefficients for g = 7 */
	static double p[] = {
		0.99999999999980993227684700473478,
		676.520368121885098567009190444019,
		-1259.13921672240287047156078755283,
		771.3234287776530788486528258894,
		-176.61502916214059906584551354,
		12.507343278686904814458936853,
		-0.13857109526572011689554707,
		9.984369578019570859563e-6,
		1.50563273514931155834e-7};

	if(z < 0.5) {return log(M_PI) - log(sin(M_PI*z)) - lngamma_lanczos_real_cfft(1. - z);}
	z -= 1;
	double x = p[0];
	for(int n = 1; n < 9; n++){ x += p[n] / (z + (double)(n));}

	double t = z + 7.5;
	return log(2*M_PI) /2.  + (z+0.5)*log(t) -t + log(x);
}


void extrap_log_linear_cfft(double *fk, int N_origin, int N_extra, double *large_fk) {
	double dln_left, dln_right;
	int i;

	dln_left = log(fk[1]/fk[0]);
	// printf("fk[0],fk[1]: %.15e,%.15e,%.15e,%.15e,%.15e\n", fk[0],fk[1],fk[2],fk[3],fk[4]);
	if(fk[0]<=0.) {
		for(i=0; i<N_extra; i++) {
			large_fk[i] = 0.;
		}
	}
	else{
		for(i=0; i<N_extra; i++) {
			large_fk[i] = exp(log(fk[0]) + (i - N_extra) * dln_left);
		}
	}

	for(i=N_extra; i< N_extra+N_origin; i++) {
		large_fk[i] = fk[i - N_extra];
	}

	dln_right = log(fk[N_origin-1]/fk[N_origin-2]);
	if(fk[N_origin-1]<=0.) {
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i] = 0.;
		}
	}
	else {
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i] = exp(log(fk[N_origin-1]) + (i - N_extra - N_origin +1) * dln_right);
		}
	}
}


// void resample_fourier_gauss(double *k, double *fk, config *config, double *k_sample, double *fk_sample) {
// 	long i;
// 	double dlnk = log(k[sizeof(k)-1]/k[0]) / (config->Nk_sample-1.);
// 	for(i=0; i<config->Nk_sample; i++) {
// 		k_sample[i] = k[0] * exp(i*dlnk);
// 		fk_sample[i] = 
// 	}
// }