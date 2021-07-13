#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include "utils_complex_cfastpt.h"
#include "utils_cfastpt.h"

#include "../log.c/src/log.h"

double complex gamma_lanczos(double complex z) {
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

	if(creal(z) < 0.5) {return M_PI / (csin(M_PI*z)*gamma_lanczos(1. - z));}
	z -= 1;
	double complex x = p[0];
	for(int n = 1; n < 9; n++){ x += p[n] / (z + (double)(n));}

	double complex t = z + 7.5;
	return sqrt(2*M_PI) * cpow(t, z+0.5) * cexp(-t) * x;
}

double complex lngamma_lanczos(double complex z) {
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

	if(creal(z) < 0.5) {return clog(M_PI) - clog(csin(M_PI*z)) - lngamma_lanczos(1. - z);}
	z -= 1;
	double complex x = p[0];
	for(int n = 1; n < 9; n++){ x += p[n] / (z + (double)(n));}

	double complex t = z + 7.5;
	return log(2*M_PI) /2.  + (z+0.5)*clog(t) -t + clog(x);
}

void f_z(double z_real, double *z_imag, double complex *fz, long N) {
/* z = nu + I*eta
Calculate g_l = exp( zln2 + lngamma( (l+nu)/2 + I*eta/2 ) - lngamma( (3+l-nu)/2 - I*eta/2 ) ) */
	long i;
	double complex z;
	g_m_vals(0.5, z_real-0.5, z_imag, fz, N);
	for(i=0; i<N; i++) {
		z = z_real+I*z_imag[i];
		// gl[i] = cexp(z*log(2.) + clog(gamma_lanczos((l+z)/2.) ) - clog(gamma_lanczos((3.+l-z)/2.)));
		fz[i] *= sqrt(M_PI)/2. * cpow(2.,z);		
	}
}

void g_m_vals(double mu, double q_real, double *q_imag, double complex *gm, long N){
	long i;
	for(i=0; i<N; i++) {
		gm[i] = cexp(lngamma_lanczos((mu+1.+q_real+I*q_imag[i])/2.) - lngamma_lanczos((mu+1.-q_real-I*q_imag[i])/2.) );		
	}
}

void gamma_ratios(double l, double nu, double *eta, double complex *gl, long N) {
/* z = nu + I*eta
Calculate g_l = exp( zln2 + lngamma( (l+nu)/2 + I*eta/2 ) - lngamma( (3+l-nu)/2 - I*eta/2 ) ) */
	long i;
	double complex z;
	for(i=0; i<N; i++) {
		z = nu+I*eta[i];
		gl[i] = cexp(lngamma_lanczos((l+z)/2.) - lngamma_lanczos((3.+l-z)/2.) );		
	}
}

void g_l(double l, double nu, double *eta, double complex *gl, long N) {
/* z = nu + I*eta
Calculate g_l = exp( zln2 + lngamma( (l+nu)/2 + I*eta/2 ) - lngamma( (3+l-nu)/2 - I*eta/2 ) ) */
	long i;
	double complex z;
	for(i=0; i<N; i++) {
		z = nu+I*eta[i];
		gl[i] = cexp(z*log(2.) + lngamma_lanczos((l+z)/2.) - lngamma_lanczos((3.+l-z)/2.) );		
	}
}

void g_l_1(double l, double nu, double *eta, double complex *gl1, long N) {
/* z = nu + I*eta
Calculate g_l_1 = exp(zln2 + lngamma( (l+nu-1)/2 + I*eta/2 ) - lngamma( (4+l-nu)/2 - I*eta/2 ) ) */
	long i;
	double complex z;
	for(i=0; i<N; i++) {
		z = nu+I*eta[i];
		gl1[i] = -(z-1.)* cexp((z-1.)*log(2.) + lngamma_lanczos((l+z-1.)/2.) - lngamma_lanczos((4.+l-z)/2.));
	}
}

void g_l_2(double l, double nu, double *eta, double complex *gl2, long N) {
/* z = nu + I*eta
Calculate g_l_1 = exp(zln2 + lngamma( (l+nu-2)/2 + I*eta/2 ) - lngamma( (5+l-nu)/2 - I*eta/2 ) ) */
	long i;
	double complex z;
	for(i=0; i<N; i++) {
		z = nu+I*eta[i];
		gl2[i] = (z-1.)* (z-2.)* cexp((z-2.)*log(2.) + lngamma_lanczos((l+z-2.)/2.) - lngamma_lanczos((5.+l-z)/2.));
	}
}

void c_window(double complex *out, double c_window_width, long halfN) {
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


void fftconvolve(fftw_complex *in1, fftw_complex *in2, long N, fftw_complex *out) {
	long i;
	fftw_complex *a, *b;
	fftw_complex *a1, *b1;
	fftw_complex *c;
	fftw_plan pa, pb, pc;

	long Ntotal;
	if(N%2==1) {
		Ntotal = 2*N - 1;
	}else {
		log_fatal("This fftconvolve doesn't support even size input arrays"); 
		exit(1);
		Ntotal = 2*N;
	}

	a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ntotal );
	b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ntotal );
	a1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ntotal );
	b1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ntotal );

	c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ntotal );

	for(i=0;i<N;i++){
		a[i] = in1[i];
		b[i] = in2[i];
	}
	for( ;i<Ntotal;i++){
		a[i] = 0.;
		b[i] = 0.;
	}

	pa = fftw_plan_dft_1d(Ntotal, a, a1, FFTW_FORWARD, FFTW_ESTIMATE);
	pb = fftw_plan_dft_1d(Ntotal, b, b1, FFTW_FORWARD, FFTW_ESTIMATE);
	#pragma omp parallel
	#pragma omp single 
	{
		#pragma omp task
		fftw_execute(pa);
		#pragma omp task
		fftw_execute(pb);
		#pragma omp taskwait
	}

	for(i=0;i<Ntotal;i++){
		a1[i] *= b1[i];
	}
	pc = fftw_plan_dft_1d(Ntotal, a1, c, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(pc);

	for(i=0;i<Ntotal; i++){
		out[i] = c[i]/(double complex)Ntotal;
	}

	fftw_destroy_plan(pa);
	fftw_destroy_plan(pb);
	fftw_destroy_plan(pc);
	fftw_free(a);
	fftw_free(b);
	fftw_free(a1);
	fftw_free(b1);
	fftw_free(c);
}

void fftconvolve_optimize(fftw_complex *in1, fftw_complex *in2, long N, fftw_complex *out, fftw_complex *a, 
fftw_complex *b, fftw_complex *a1, fftw_complex *b1, fftw_complex *c, fftw_plan pa, fftw_plan pb, fftw_plan pc) 
{
	long i;

	long Ntotal;
	if(N%2==1) 
	{
		Ntotal = 2*N - 1;
	}
	else 
	{
		log_fatal("This fftconvolve doesn't support even size input arrays"); 
		exit(1);
		Ntotal = 2*N;
	}

	for(i=0;i<N;i++){
		a[i] = in1[i];
		b[i] = in2[i];
	}
	for( ;i<Ntotal;i++){
		a[i] = 0.;
		b[i] = 0.;
	}

//	#pragma omp parallel
//	#pragma omp single 
	{
//		#pragma omp task
		fftw_execute(pa);
//		#pragma omp task
		fftw_execute(pb);
//		#pragma omp taskwait
	}

	for(i=0;i<Ntotal;i++)
	{
		a1[i] *= b1[i];
	}
	
	fftw_execute(pc);

	for(i=0;i<Ntotal; i++)
	{
		out[i] = c[i]/(double complex)Ntotal;
	}
}

void fftconvolve_real(double *in1, double *in2, long N1, long N2, double *out) 
{
	long i;
	double *a, *b;
	fftw_complex *a1, *b1;
	double *c;
	fftw_plan pa, pb, pc;

	long Ntotal = N1+N2-1;
	long Ncomplex;
	if(Ntotal%2==1){Ncomplex = (Ntotal+1)/2;}
	else{Ncomplex = Ntotal/2+1;}

	a = (double*) fftw_malloc(sizeof(double) * Ntotal );
	b = (double*) fftw_malloc(sizeof(double) * Ntotal );
	a1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ncomplex );
	b1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ncomplex );

	c = (double*) fftw_malloc(sizeof(double) * Ntotal );

	for(i=0;i<N1;i++)
	{
		a[i] = in1[i];
	}
	for( ;i<Ntotal;i++)
	{
		a[i] = 0.;
	}
	for(i=0;i<N2;i++)
	{
		b[i] = in2[i];
	}
	for( ;i<Ntotal;i++)
	{
		b[i] = 0.;
	}

	pa = fftw_plan_dft_r2c_1d(Ntotal, a, a1, FFTW_ESTIMATE);
	pb = fftw_plan_dft_r2c_1d(Ntotal, b, b1, FFTW_ESTIMATE);
	
	#pragma omp parallel
	#pragma omp single 
	{
		#pragma omp task
		fftw_execute(pa);
		#pragma omp task
		fftw_execute(pb);
		#pragma omp taskwait
	}
	
	for(i=0;i<Ncomplex;i++){
		a1[i] *= b1[i];
	}
	pc = fftw_plan_dft_c2r_1d(Ntotal, a1, c, FFTW_ESTIMATE);
	fftw_execute(pc);

	for(i=0;i<Ntotal; i++){
		out[i] = c[i]/(double)Ntotal;
	}

	fftw_destroy_plan(pa);
	fftw_destroy_plan(pb);
	fftw_destroy_plan(pc);
	fftw_free(a);
	fftw_free(b);
	fftw_free(a1);
	fftw_free(b1);
	fftw_free(c);
}