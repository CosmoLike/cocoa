#include <complex.h>
#include <fftw3.h>

#ifndef __CFFTLOG_UTILS_COMPLEX_H
#define __CFFTLOG_UTILS_COMPLEX_H
#ifdef __cplusplus
extern "C" {
#endif

void g_l_cfft(double l, double nu, double* eta, double complex* gl, long N);

void g_l_1_cfft(double l, double nu, double* eta, double complex* gl1, long N);

void g_l_2_cfft(double l, double nu, double* eta, double complex* gl2, long N);

void c_window_cfft(double complex *out, double c_window_width, long halfN);

double complex gamma_lanczos_cfft(double complex z);

double complex lngamma_lanczos_cfft(double complex z);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD