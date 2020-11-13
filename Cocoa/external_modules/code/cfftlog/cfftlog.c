#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include "cfftlog.h"
#include "utils.h"
#include "utils_complex.h"


void cfftlog(double *x, double *fx, long N, config *config, int ell, double *y, double *Fy) {	
	long N_original = N;
	long N_pad = config->N_pad;
	N += 2*N_pad;

	if(N % 2) {printf("Please use even number of x !\n"); exit(0);}
	long halfN = N/2;

	double x0, y0;
	x0 = x[0];

	double dlnx;
	dlnx = log(x[1]/x0);

	// Only calculate the m>=0 part
	double eta_m[halfN+1];
	long i;
	for(i=0; i<=halfN; i++) {eta_m[i] = 2*M_PI / dlnx / N * i;}

	double complex gl[halfN+1];
	
	switch(config->derivative) {
		case 0: g_l_cfft((double)ell, config->nu, eta_m, gl, halfN+1); break;
		case 1: g_l_1_cfft((double)ell, config->nu, eta_m, gl, halfN+1); break;
		case 2: g_l_2_cfft((double)ell, config->nu, eta_m, gl, halfN+1); break;
		default: printf("Integral Not Supported! Please choose config->derivative from [0,1,2].\n");
	}
	// printf("g2[0]: %.15e+I*(%.15e)\n", creal(g2[0]),cimag(g2[0]));

	// calculate y arrays
	for(i=0; i<N_original; i++) {y[i] = (ell+1.) / x[N_original-1-i];}
	y0 = y[0];

	// biased input func
	double *fb;
	fb = malloc(N* sizeof(double));
	for(i=0; i<N_pad; i++) {
		fb[i] = 0.;
		fb[N-1-i] = 0.;
	}
	for(i=N_pad; i<N_pad+N_original; i++) {
		fb[i] = fx[i-N_pad] / pow(x[i-N_pad], config->nu) ;
	}

	fftw_complex *out;
	fftw_plan plan_forward, plan_backward;
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (halfN+1) );
	plan_forward = fftw_plan_dft_r2c_1d(N, fb, out, FFTW_ESTIMATE);

	fftw_execute(plan_forward);

	c_window_cfft(out, config->c_window_width, halfN);
	// printf("out[1]:%.15e+i*(%.15e)\n", creal(out[1]), cimag(out[1]));

	for(i=0; i<=halfN; i++) {
		out[i] *= cpow(x0*y0/exp(2*N_pad*dlnx), -I*eta_m[i]) * gl[i] ;
		out[i] = conj(out[i]);
	}

	// printf("out[1]:%.15e+i*(%.15e)\n", creal(out[1]), cimag(out[1]));
	double *out_ifft;
	out_ifft = malloc(sizeof(double) * N );
	plan_backward = fftw_plan_dft_c2r_1d(N, out, out_ifft, FFTW_ESTIMATE);

	fftw_execute(plan_backward);

	for(i=0; i<N_original; i++) {
		Fy[i] = out_ifft[i-N_pad] * sqrt(M_PI) / (4.*N * pow(y[i], config->nu));
	}

	fftw_destroy_plan(plan_forward);
	fftw_destroy_plan(plan_backward);
	fftw_free(out);
	free(out_ifft);
	free(fb);
}

void cfftlog_ells(double *x, double *fx, long N, config *config, int* ell, long Nell, double **y, double **Fy) {	
	long N_original = N;
	long N_pad = config->N_pad;
	long N_extrap_low = config->N_extrap_low;
	long N_extrap_high = config->N_extrap_high;
	N += (2*N_pad + N_extrap_low+N_extrap_high);

	if(N % 2) {
		printf("Please use even number of x !\n"); 
		exit(0);
	}
	long halfN = N/2;

	double x0;
	x0 = x[0];

	double dlnx;
	dlnx = log(x[1]/x0);

	// Only calculate the m>=0 part
	double eta_m[halfN+1];
	for(int i=0; i<=halfN; i++) {
		eta_m[i] = 2*M_PI / dlnx / N * i;
	}

	// biased input func
	double *fb = malloc(N* sizeof(double));
	for(int i=0; i<N_pad; i++) {
		fb[i] = 0.;
		fb[N-1-i] = 0.;
	}
	double xi;
	int sign;
	if(N_extrap_low) {
		if(fx[0]==0) {
			printf("Can't log-extrapolate zero on the low side!\n");
			exit(1);
		}
		else if(fx[0]>0) {sign = 1;}
		else {sign=-1;}
		if(fx[1]/fx[0]<=0) {printf("Log-extrapolation on the low side fails due to sign change!\n"); exit(1);}
		double dlnf_low = log(fx[1]/fx[0]);
		for(int i=N_pad; i<N_pad+N_extrap_low; i++) {
			xi = exp(log(x0) + (i-N_pad - N_extrap_low)*dlnx);
			fb[i] = sign * exp(log(fx[0]*sign) + (i- N_pad - N_extrap_low)*dlnf_low) / pow(xi, config->nu);
		}
	}
	for(int i=N_pad+N_extrap_low; i<N_pad+N_extrap_low+N_original; i++) {
		fb[i] = fx[i-N_pad-N_extrap_low] / pow(x[i-N_pad-N_extrap_low], config->nu) ;
	}
	if(N_extrap_high) {
		if(fx[N_original-1]==0) {
			printf("Can't log-extrapolate zero on the high side!\n");
			exit(1);
		}
		else if(fx[N_original-1]>0) {sign = 1;}
		else {sign=-1;}
		if(fx[N_original-1]/fx[N_original-2]<=0) {printf("Log-extrapolation on the high side fails due to sign change!\n"); exit(1);}
		double dlnf_high = log(fx[N_original-1]/fx[N_original-2]);
		for(int i=N-N_pad-N_extrap_high; i<N-N_pad; i++) {
			xi = exp(log(x[N_original-1]) + (i-N_pad - N_extrap_low- N_original)*dlnx);
			fb[i] = sign * exp(log(fx[N_original-1]*sign) + (i- N_pad - N_extrap_low- N_original)*dlnf_high) / pow(xi, config->nu);
		}
	}

	fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (halfN+1) );
	fftw_plan plan_forward;
	
	plan_forward = fftw_plan_dft_r2c_1d(N, fb, out, FFTW_ESTIMATE);
	fftw_execute(plan_forward);
	c_window_cfft(out, config->c_window_width, halfN);
	

	double **out_ifft = malloc(sizeof(double*) * Nell);
	fftw_complex **out_vary = malloc(sizeof(fftw_complex*) * Nell);
	fftw_plan* plan_backward = malloc(sizeof(fftw_plan*) * Nell); 
	for(int j=0; j<Nell; j++) {
		out_ifft[j] = malloc(sizeof(double) * N);
		out_vary[j] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (halfN+1) );
		plan_backward[j] = fftw_plan_dft_c2r_1d(N, out_vary[j], out_ifft[j], FFTW_ESTIMATE);
	}

	#pragma omp parallel for
	for (int j=0; j<Nell; j++) {
		double complex gl[halfN+1];

		switch(config->derivative) {
			case 0: g_l_cfft((double)ell[j], config->nu, eta_m, gl, halfN+1); break;
			case 1: g_l_1_cfft((double)ell[j], config->nu, eta_m, gl, halfN+1); break;
			case 2: g_l_2_cfft((double)ell[j], config->nu, eta_m, gl, halfN+1); break;
			default: printf("Integral Not Supported! Please choose config->derivative from [0,1,2].\n");
		}

		// calculate y arrays
		for(int i=0; i<N_original; i++) {
			y[j][i] = (ell[j]+1.) / x[N_original-1-i];
		}
		const double y0 = y[j][0];

		for(int i=0; i<=halfN; i++) {
			out_vary[j][i] = conj(out[i] * cpow(x0*y0/exp((N-N_original)*dlnx), -I*eta_m[i]) * gl[i]) ;
		}

		fftw_execute(plan_backward[j]);

		for(int i=0; i<N_original; i++) {
			Fy[j][i] = out_ifft[j][i+N_pad+N_extrap_high] * sqrt(M_PI) / (4.*N * pow(y[j][i], config->nu));
		}
	}

	for (int j=0; j<Nell; j++) {
		fftw_destroy_plan(plan_backward[j]);
		fftw_free(out_vary[j]);
		free(out_ifft[j]);
	}
	free(plan_backward);
	free(out_vary);
	free(out_ifft);
	fftw_destroy_plan(plan_forward);
	fftw_free(out);
	free(fb);
}

void cfftlog_ells_increment(double *x, double *fx, long N, config *config, int* ell, long Nell, double **y, double **Fy) {

	long N_original = N;
	long N_pad = config->N_pad;
	long N_extrap_low = config->N_extrap_low;
	long N_extrap_high = config->N_extrap_high;
	N += (2*N_pad + N_extrap_low+N_extrap_high);

	if(N % 2) {
		printf("Please use even number of x !\n"); exit(0);
	}
	long halfN = N/2;

	double x0;
	x0 = x[0];

	double dlnx;
	dlnx = log(x[1]/x0);

	// Only calculate the m>=0 part
	double eta_m[halfN+1];
	for(int i=0; i<=halfN; i++) {
		eta_m[i] = 2*M_PI / dlnx / N * i;
	}

	// biased input func
	double *fb;
	fb = malloc(N* sizeof(double));
	for(int i=0; i<N_pad; i++) {
		fb[i] = 0.;
		fb[N-1-i] = 0.;
	}
	double xi;
	int sign;
	if(N_extrap_low) {
		if(fx[0]==0) {
			printf("Can't log-extrapolate zero on the low side!\n");
			exit(1);
		}
		else if(fx[0]>0) {sign = 1;}
		else {sign=-1;}
		if(fx[1]/fx[0]<=0) {printf("Log-extrapolation on the low side fails due to sign change!\n"); exit(1);}
		double dlnf_low = log(fx[1]/fx[0]);
		for(int i=N_pad; i<N_pad+N_extrap_low; i++) {
			xi = exp(log(x0) + (i-N_pad - N_extrap_low)*dlnx);
			fb[i] = sign * exp(log(fx[0]*sign) + (i- N_pad - N_extrap_low)*dlnf_low) / pow(xi, config->nu);
		}
	}
	for(int i=N_pad+N_extrap_low; i<N_pad+N_extrap_low+N_original; i++) {
		fb[i] = fx[i-N_pad-N_extrap_low] / pow(x[i-N_pad-N_extrap_low], config->nu) ;
	}
	if(N_extrap_high) {
		if(fx[N_original-1]==0) {
			printf("Can't log-extrapolate zero on the high side!\n");
			exit(1);
		}
		else if(fx[N_original-1]>0) {sign = 1;}
		else {sign=-1;}
		if(fx[N_original-1]/fx[N_original-2]<=0) {printf("Log-extrapolation on the high side fails due to sign change!\n"); exit(1);}
		double dlnf_high = log(fx[N_original-1]/fx[N_original-2]);
		for(int i=N-N_pad-N_extrap_high; i<N-N_pad; i++) {
			xi = exp(log(x[N_original-1]) + (i-N_pad - N_extrap_low- N_original)*dlnx);
			fb[i] = sign * exp(log(fx[N_original-1]*sign) + (i- N_pad - N_extrap_low- N_original)*dlnf_high) / pow(xi, config->nu);
		}
	}

	fftw_complex *out;
	fftw_plan plan_forward;
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (halfN+1) );
	plan_forward = fftw_plan_dft_r2c_1d(N, fb, out, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	c_window_cfft(out, config->c_window_width, halfN);

	double **out_ifft = malloc(sizeof(double*) * Nell);
	fftw_complex **out_vary = malloc(sizeof(fftw_complex*) * Nell);
	fftw_plan* plan_backward = malloc(sizeof(fftw_plan*) * Nell); 
	for(int j=0; j<Nell; j++) {
		out_ifft[j] = malloc(sizeof(double) * N);
		out_vary[j] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (halfN+1) );
		plan_backward[j] = fftw_plan_dft_c2r_1d(N, out_vary[j], out_ifft[j], FFTW_ESTIMATE);
	}

	#pragma omp parallel for
	for(int j=0; j<Nell; j++){
		double complex gl[halfN+1];
		switch(config->derivative) {
			case 0: g_l_cfft((double)ell[j], config->nu, eta_m, gl, halfN+1); break;
			case 1: g_l_1_cfft((double)ell[j], config->nu, eta_m, gl, halfN+1); break;
			case 2: g_l_2_cfft((double)ell[j], config->nu, eta_m, gl, halfN+1); break;
			default: printf("Integral Not Supported! Please choose config->derivative from [0,1,2].\n");
		}

		// calculate y arrays
		for(int i=0; i<N_original; i++) {
			y[j][i] = (ell[j]+1.) / x[N_original-1-i];
		}
		const double y0 = y[j][0];

		for(int i=0; i<=halfN; i++) {
			out_vary[j][i] = conj(out[i] * cpow(x0*y0/exp((N-N_original)*dlnx), -I*eta_m[i]) * gl[i]) ;
		}

		fftw_execute(plan_backward[j]);

		for(int i=0; i<N_original; i++) {
			Fy[j][i] += out_ifft[j][i+N_pad+N_extrap_high] * sqrt(M_PI) / (4.*N * pow(y[j][i], config->nu));
		}
	}

	for (int j=0; j<Nell; j++) {
		fftw_destroy_plan(plan_backward[j]);
		fftw_free(out_vary[j]);
		free(out_ifft[j]);
	}
	free(plan_backward);
	free(out_vary);
	free(out_ifft);
	fftw_destroy_plan(plan_forward);
	fftw_free(out);
	free(fb);
}
