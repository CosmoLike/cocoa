#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include "utils.h"

double complex gamma_lanczos_cfft(double complex z) {
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

	if(creal(z) < 0.5) {return M_PI / (csin(M_PI*z)*gamma_lanczos_cfft(1. - z));}
	z -= 1;
	double complex x = p[0];
	for(int n = 1; n < 9; n++){ x += p[n] / (z + (double)(n));}

	double complex t = z + 7.5;
	return sqrt(2*M_PI) * cpow(t, z+0.5) * cexp(-t) * x;
}

double complex lngamma_lanczos_cfft(double complex z) {
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

	if(creal(z) < 0.5) {return clog(M_PI) - clog(csin(M_PI*z)) - lngamma_lanczos_cfft(1. - z);}
	z -= 1;
	double complex x = p[0];
	for(int n = 1; n < 9; n++){ x += p[n] / (z + (double)(n));}

	double complex t = z + 7.5;
	return log(2*M_PI) /2.  + (z+0.5)*clog(t) -t + clog(x);
}

void g_l_cfft(double l, double nu, double *eta, double complex *gl, long N) {
/* z = nu + I*eta
Calculate g_l = exp( zln2 + lngamma( (l+nu)/2 + I*eta/2 ) - lngamma( (3+l-nu)/2 - I*eta/2 ) ) */
	long i;
	double complex z;
	for(i=0; i<N; i++) {
		z = nu+I*eta[i];
		// gl[i] = cexp(z*log(2.) + clog(gamma_lanczos((l+z)/2.) ) - clog(gamma_lanczos((3.+l-z)/2.)));
		gl[i] = cexp(z*log(2.) + lngamma_lanczos_cfft((l+z)/2.) - lngamma_lanczos_cfft((3.+l-z)/2.) );		
		// if(isnan(gl[i])) {printf("nan at l,nu,eta, = %lf %lg %lg %lg %lg\n", l,nu, eta[i], gamma_lanczos((l+z)/2.),gamma_lanczos((3.+l-z)/2.));exit(0);}
	}
}

void g_l_1_cfft(double l, double nu, double *eta, double complex *gl1, long N) {
/* z = nu + I*eta
Calculate g_l_1 = exp(zln2 + lngamma( (l+nu-1)/2 + I*eta/2 ) - lngamma( (4+l-nu)/2 - I*eta/2 ) ) */
	long i;
	double complex z;
	for(i=0; i<N; i++) {
		z = nu+I*eta[i];
		gl1[i] = -(z-1.)* cexp((z-1.)*log(2.) + lngamma_lanczos_cfft((l+z-1.)/2.) - lngamma_lanczos_cfft((4.+l-z)/2.));
	}
}

void g_l_2_cfft(double l, double nu, double *eta, double complex *gl2, long N) {
/* z = nu + I*eta
Calculate g_l_1 = exp(zln2 + lngamma( (l+nu-2)/2 + I*eta/2 ) - lngamma( (5+l-nu)/2 - I*eta/2 ) ) */
	long i;
	double complex z;
	for(i=0; i<N; i++) {
		z = nu+I*eta[i];
		gl2[i] = (z-1.)* (z-2.)* cexp((z-2.)*log(2.) + lngamma_lanczos_cfft((l+z-2.)/2.) - lngamma_lanczos_cfft((5.+l-z)/2.));
	}
}

void c_window_cfft(double complex *out, double c_window_width, long halfN) {
	// 'out' is (halfN+1) complex array
	long Ncut;
	Ncut = (long)(halfN * c_window_width);
	long i;
	double W;
	for(i=0; i<=Ncut; i++) { // window for right-side
		W = (double)(i)/Ncut - 1./(2.*M_PI) * sin(2.*i*M_PI/Ncut);
		out[halfN-i] *= W;
	}
}

