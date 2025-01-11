#include <fftw3.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include "structs.h"

#ifndef __COSMOLIKE_BASICS_H
#define __COSMOLIKE_BASICS_H
#ifdef __cplusplus
extern "C" {
#endif

#define NR_END 1
#define FREE_ARG char *

gsl_integration_glfixed_table* malloc_gslint_glfixed(const int n);
void* malloc1d(const int nx);
void** malloc2d(const int nx, const int ny);
void* calloc1d(const int nx);

int fdiff(const double a, const double b);

double fmin(const double a, const double b);

double fmax(const double a, const double b);

bin_avg set_bin_average(
    const int i_theta, 
    const int j_L
  );

double int_gsl_integrate_high_precision(double (*func)(double, void *),
  void *arg, double a, double b, double *error, int niter);

double int_gsl_integrate_medium_precision(double (*func)(double, void *),
  void *arg, double a, double b, double *error, int niter);

double int_gsl_integrate_low_precision(double (*func)(double, void *),
  void *arg, double a, double b, double *error, int niter);

double interpol2d(double **f, int nx, double ax, double bx, double dx, double x,
  int ny, double ay, double by, double dy, double y, double lower, double upper);

double interpol2d_fitslope(double **f, int nx, double ax, double bx, double dx,
  double x, int ny, double ay, double by, double dy, double y, double lower);

double interpol(const double* const f, const int n, const double a, 
const double b, const double dx, const double x, const double lower, 
const double upper);

double interpol_fitslope(double *f, int n, double a, double b, double dx,
  double x, double lower);

int line_count(char *filename);

void error(char *s);

void hankel_kernel_FT(double x, fftw_complex *res, double *arg,
int argc __attribute__((unused)));

void cdgamma(fftw_complex x, fftw_complex *res);

void hankel_kernel_FT_3D(double x, fftw_complex *res, double *arg, int argc);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD
