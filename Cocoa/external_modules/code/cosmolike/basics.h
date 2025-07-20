#include <fftw3.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
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

typedef struct 
{
  double xmin;
  double xmax;
  double Pmin;
  double Pmax;
  double dPmin;
  double dPmax;
} bin_avg;

gsl_interp* malloc_gsl_interp(const int n);

gsl_spline* malloc_gsl_spline(const int n);

gsl_integration_glfixed_table* malloc_gslint_glfixed(const int n);

void* malloc1d(const int nx);

void** malloc2d(const int nx, const int ny);

void*** malloc3d(const int nx, const int ny, const int nz);

void* calloc1d(const int nx);

int fdiff(const double a, const double b);

double fmin(const double a, const double b);

double fmax(const double a, const double b);

bin_avg set_bin_average(const int i_theta, const int j_L);

double interpol1d(
  const double* const f, 
  const int n, 
  const double a, 
  const double b, 
  const double dx, 
  const double x); 

double interpol2d(double **f, int nx, double ax, double bx, double dx, double x,
  int ny, double ay, double by, double dy, double y);

int line_count(char *filename);

void hankel_kernel_FT(double x, fftw_complex *res, double *arg,
int argc __attribute__((unused)));

void cdgamma(fftw_complex x, fftw_complex *res);

void hankel_kernel_FT_3D(double x, fftw_complex *res, double *arg, int argc);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD
