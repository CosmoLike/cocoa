#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#include "structs.h"
#include "basics.h"

#include "log.c/src/log.h"

double fmin(const double a, const double b)
{
  return a < b ? a : b;
}

double fmax(const double a, const double b)
{
  return a > b ? a : b;
}

bin_avg set_bin_average(int i_theta, int j_L)
{
  static double** Pmin  = 0;
  static double** Pmax  = 0;
  static double** dPmin = 0;
  static double** dPmax = 0;
  static double* xmin = 0;
  static double* xmax = 0;
  
  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  if (!(i_theta < like.Ntheta))
   {
     log_fatal("bad i_theta index");
     exit(1);
   }
   if (!(j_L < limits.LMAX))
   {
     log_fatal("bad j_L index");
     exit(1);
   }

  if (Pmin == 0)
  {
    Pmin = (double**) malloc(like.Ntheta*sizeof(double*));
    Pmax = (double**) malloc(like.Ntheta*sizeof(double*));
    dPmin = (double**) malloc(like.Ntheta*sizeof(double*));
    dPmax = (double**) malloc(like.Ntheta*sizeof(double*));
    for(int i=0; i<like.Ntheta ; i++)
    {
      Pmin[i] = (double*) calloc(limits.LMAX, sizeof(double));
      Pmax[i] = (double*) calloc(limits.LMAX, sizeof(double));
      dPmin[i] = (double*) calloc(limits.LMAX, sizeof(double));
      dPmax[i] = (double*) calloc(limits.LMAX, sizeof(double));
    }
    xmin = (double*) calloc(like.Ntheta, sizeof(double));
    xmax = (double*) calloc(like.Ntheta, sizeof(double));
    const double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    for(int i=0; i<like.Ntheta ; i++)
    {
      xmin[i] = cos(exp(log(like.vtmin) + (i + 0.0)*logdt));
      xmax[i] = cos(exp(log(like.vtmin) + (i + 1.0)*logdt));
    }
    #pragma omp parallel for
    for (int i=0; i<like.Ntheta; i++)
    {
      int status = gsl_sf_legendre_Pl_deriv_array(limits.LMAX, xmin[i], Pmin[i], dPmin[i]);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
      status = gsl_sf_legendre_Pl_deriv_array(limits.LMAX, xmax[i], Pmax[i], dPmax[i]);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      } 
    }
  }
  bin_avg r;
  r.xmin = xmin[i_theta];
  r.xmax = xmax[i_theta];
  r.Pmin = Pmin[i_theta][j_L];
  r.Pmax = Pmax[i_theta][j_L];
  r.dPmin = dPmin[i_theta][j_L];
  r.dPmax = dPmax[i_theta][j_L];
  return r;
}

/*
double int_gsl_integrate_high_precision(double (*func)(double, void*),
void* arg, double a, double b, double* error, int niter)
{
  double res, err;
  gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(niter);
  gsl_function F;
  F.function = func;
  F.params = arg;
  gsl_integration_cquad(&F, a, b, 0, precision.high, w, &res, &err, 0);
  if (NULL != error)
    *error = err;
  gsl_integration_cquad_workspace_free(w);
  return res;
}
*/

double int_gsl_integrate_high_precision(double (*func)(double, void*),
void* arg, double a, double b, double* error, int niter)
{
  double res, err;
  
  // Hacked from GSL souce code (cquad.c) - stack (faster) vs heap allocation
  size_t heap[niter];
  gsl_integration_cquad_ival ival[niter];
  gsl_integration_cquad_workspace w;
  w.size = niter;
  w.ivals = ival;
  w.heap = heap;  

  gsl_function F;
  F.function = func;
  F.params = arg;
  gsl_integration_cquad(&F, a, b, 0, precision.high, &w, &res, &err, 0);
  if (NULL != error) 
  {
    *error = err;
  }
  return res;
}

/*
double int_gsl_integrate_medium_precision(double (*func)(double, void*),
void* arg, double a, double b, double* error, int niter)
{
  double res, err;
  gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(niter);
  gsl_function F;
  F.function = func;
  F.params = arg;
  gsl_integration_cquad(&F, a, b, 0, precision.medium, w, &res, &err, 0);
  if (NULL != error) 
  {
    *error = err;
  }
  gsl_integration_cquad_workspace_free(w);
  return res;
}
*/


double int_gsl_integrate_medium_precision(double (*func)(double, void*),
void* arg, double a, double b, double* error, int niter)
{
  double res, err;
  
  // Hacked from GSL souce code (cquad.c) - stack (faster) vs heap allocation
  size_t heap[niter];
  gsl_integration_cquad_ival ival[niter];
  gsl_integration_cquad_workspace w;
  w.size = niter;
  w.ivals = ival;
  w.heap = heap;

  gsl_function F;
  F.function = func;
  F.params = arg;
  gsl_integration_cquad(&F, a, b, 0, precision.medium, &w, &res, &err, 0);
  if (NULL != error) 
  {
    *error = err;
  }
  return res;
}

double int_gsl_integrate_low_precision(double (*func)(double, void*),
void* arg, double a, double b, double* error __attribute__((unused)),
int niter __attribute__((unused))) 
{
  gsl_set_error_handler_off();
  double res, err;
  gsl_function F;
  F.function = func;
  F.params = arg;
  size_t neval;
  gsl_integration_qng(&F, a, b, 0, precision.medium, &res, &err, &neval);
  return res;
}

void error(char *s) {
  printf("error:%s\n ", s);
  exit(1);
}

// ============================================================
// Interpolates f at the value x, where f is a double[n] array,
// representing a function between a and b, stepwidth dx.
// 'lower' and 'upper' are powers of a logarithmic power law
// extrapolation. If no	extrapolation desired, set these to 0
// ============================================================

double interpol(double *f, int n, double a, double b, double dx, double x,
double lower, double upper) {
  double r;
  int i;
  if (x < a) {
    if (lower == 0.) {
      return 0.0;
    }
    return f[0] + lower * (x - a);
  }
  r = (x - a) / dx;
  i = (int)(floor(r));
  if (i + 1 >= n) {
    if (upper == 0.0) {
      if (i + 1 == n) {
        return f[i]; // constant extrapolation
      } else {
        return 0.0;
      }
    } else {
      return f[n - 1] + upper * (x - b); /* linear extrapolation */
    }
  } else {
    return (r - i) * (f[i + 1] - f[i]) + f[i]; /* interpolation */
  }
}

int line_count(char *filename) {
  FILE *n;
  n = fopen(filename, "r");
  if (!n) {
    printf("line_count: %s not found!\nEXIT!\n", filename);
    exit(1);
  }
  int ch = 0, prev = 0, number_of_lines = 0;

  do {
    prev = ch;
    ch = fgetc(n);
    if (ch == '\n')
      number_of_lines++;
  } while (ch != EOF);
  fclose(n);
  // last line might not end with \n, but if previous character does, last line
  // is empty
  if (ch != '\n' && prev != '\n' && number_of_lines != 0)
    number_of_lines++;
  return number_of_lines;
}

double interpol_fitslope(double *f, int n, double a, double b, double dx,
                         double x, double lower) {
  double r;
  int i, fitrange;
  if (x < a) {
    if (lower == 0.) {
      return 0.0;
    }
    return f[0] + lower * (x - a);
  }
  r = (x - a) / dx;
  i = (int)(floor(r));
  if (i + 1 >= n) {
    if (n > 50) {
      fitrange = 5;
    } else {
      fitrange = (int)floor(n / 10);
    }
    double upper = (f[n - 1] - f[n - 1 - fitrange]) / (dx * fitrange);
    return f[n - 1] + upper * (x - b); /* linear extrapolation */

  } else {
    return (r - i) * (f[i + 1] - f[i]) + f[i]; /* interpolation */
  }
}

/* ============================================================ *
 * like interpol, but f beeing a 2d-function			*
 * 'lower' and 'upper' are the powers of a power law extra-	*
 * polation in the second argument				*
 * ============================================================ */
double interpol2d(double **f, int nx, double ax, double bx, double dx, double x,
int ny, double ay, double by, double dy, double y, double lower, double upper) {
  double t, dt, s, ds;
  int i, j;
  if (x < ax) {
    return 0.;
  }
  if (x > bx) {
    return 0.;
  }
  t = (x - ax) / dx;
  i = (int)(floor(t));
  dt = t - i;
  if (y < ay) {
    return ((1. - dt) * f[i][0] + dt * f[i + 1][0]) + (y - ay) * lower;
  } else if (y > by) {
    return ((1. - dt) * f[i][ny - 1] + dt * f[i + 1][ny - 1]) +
           (y - by) * upper;
  }
  s = (y - ay) / dy;
  j = (int)(floor(s));
  ds = s - j;
  if ((i + 1 == nx) && (j + 1 == ny)) {
    return (1. - dt) * (1. - ds) * f[i][j];
  }
  if (i + 1 == nx) {
    return (1. - dt) * (1. - ds) * f[i][j] + (1. - dt) * ds * f[i][j + 1];
  }
  if (j + 1 == ny) {
    return (1. - dt) * (1. - ds) * f[i][j] + dt * (1. - ds) * f[i + 1][j];
  }
  return (1. - dt) * (1. - ds) * f[i][j] + (1. - dt) * ds * f[i][j + 1] +
         dt * (1. - ds) * f[i + 1][j] + dt * ds * f[i + 1][j + 1];
}

double interpol2d_fitslope(double **f, int nx, double ax, double bx, double dx,
double x, int ny, double ay, double by, double dy, double y, double lower) {
  double t, dt, s, ds, upper;
  int i, j, fitrange;
  if (x < ax) {
    return 0.;
  }

  if (x > bx) {
    return 0.;
  }
  t = (x - ax) / dx;
  i = (int)(floor(t));
  dt = t - i;
  if (y < ay) {
    return ((1. - dt) * f[i][0] + dt * f[i + 1][0]) + (y - ay) * lower;
  } else if (y > by) {
    if (ny > 25) {
      fitrange = 5;
    } else {
      fitrange = (int)floor(ny / 5);
    }
    upper = ((1. - dt) * (f[i][ny - 1] - f[i][ny - 1 - fitrange]) +
             dt * (f[i + 1][ny - 1] - f[i + 1][ny - 1 - fitrange])) /
            (dy * fitrange);
    return ((1. - dt) * f[i][ny - 1] + dt * f[i + 1][ny - 1]) +
           (y - by) * upper;
  }
  s = (y - ay) / dy;
  j = (int)(floor(s));
  ds = s - j;
  if ((i + 1 == nx) && (j + 1 == ny)) {
    return (1. - dt) * (1. - ds) * f[i][j];
  }
  if (i + 1 == nx) {
    return (1. - dt) * (1. - ds) * f[i][j] + (1. - dt) * ds * f[i][j + 1];
  }
  if (j + 1 == ny) {
    return (1. - dt) * (1. - ds) * f[i][j] + dt * (1. - ds) * f[i + 1][j];
  }
  return (1. - dt) * (1. - ds) * f[i][j] + (1. - dt) * ds * f[i][j + 1] +
         dt * (1. - ds) * f[i + 1][j] + dt * ds * f[i + 1][j + 1];
}

void hankel_kernel_FT(double x, fftw_complex *res, double *arg,
int argc __attribute__((unused))) {
  fftw_complex a1, a2, g1, g2;

  // arguments for complex gamma
  const double q = arg[0];
  const int mu = (int)(arg[1] + 0.1);
  a1[0] = 0.5 * (1.0 + mu + q);
  a2[0] = 0.5 * (1.0 + mu - q);
  a1[1] = 0.5 * x;
  a2[1] = -a1[1];

  cdgamma(a1, &g1);
  cdgamma(a2, &g2);

  const double xln2 = x * M_LN2;
  const double si = sin(xln2);
  const double co = cos(xln2);
  const double d1 = g1[0] * g2[0] + g1[1] * g2[1]; /* Re */
  const double d2 = g1[1] * g2[0] - g1[0] * g2[1]; /* Im */
  const double mod = g2[0] * g2[0] + g2[1] * g2[1];
  const double pref = exp(M_LN2 * q) / mod;

  (*res)[0] = pref * (co * d1 - si * d2);
  (*res)[1] = pref * (si * d1 + co * d2);
}

void cdgamma(fftw_complex x, fftw_complex *res) {
  double xr, xi, wr, wi, ur, ui, vr, vi, yr, yi, t;

  xr = (double)x[0];
  xi = (double)x[1];

  if (xr < 0) {
    wr = 1 - xr;
    wi = -xi;
  } else {
    wr = xr;
    wi = xi;
  }

  ur = wr + 6.00009857740312429;
  vr = ur * (wr + 4.99999857982434025) - wi * wi;
  vi = wi * (wr + 4.99999857982434025) + ur * wi;
  yr = ur * 13.2280130755055088 + vr * 66.2756400966213521 +
       0.293729529320536228;
  yi = wi * 13.2280130755055088 + vi * 66.2756400966213521;
  ur = vr * (wr + 4.00000003016801681) - vi * wi;
  ui = vi * (wr + 4.00000003016801681) + vr * wi;
  vr = ur * (wr + 2.99999999944915534) - ui * wi;
  vi = ui * (wr + 2.99999999944915534) + ur * wi;
  yr += ur * 91.1395751189899762 + vr * 47.3821439163096063;
  yi += ui * 91.1395751189899762 + vi * 47.3821439163096063;
  ur = vr * (wr + 2.00000000000603851) - vi * wi;
  ui = vi * (wr + 2.00000000000603851) + vr * wi;
  vr = ur * (wr + 0.999999999999975753) - ui * wi;
  vi = ui * (wr + 0.999999999999975753) + ur * wi;
  yr += ur * 10.5400280458730808 + vr;
  yi += ui * 10.5400280458730808 + vi;
  ur = vr * wr - vi * wi;
  ui = vi * wr + vr * wi;
  t = ur * ur + ui * ui;
  vr = yr * ur + yi * ui + t * 0.0327673720261526849;
  vi = yi * ur - yr * ui;
  yr = wr + 7.31790632447016203;
  ur = log(yr * yr + wi * wi) * 0.5 - 1;
  ui = atan2(wi, yr);
  yr = exp(ur * (wr - 0.5) - ui * wi - 3.48064577727581257) / t;
  yi = ui * (wr - 0.5) + ur * wi;
  ur = yr * cos(yi);
  ui = yr * sin(yi);
  yr = ur * vr - ui * vi;
  yi = ui * vr + ur * vi;
  if (xr < 0) {
    wr = xr * 3.14159265358979324;
    wi = exp(xi * 3.14159265358979324);
    vi = 1 / wi;
    ur = (vi + wi) * sin(wr);
    ui = (vi - wi) * cos(wr);
    vr = ur * yr + ui * yi;
    vi = ui * yr - ur * yi;
    ur = 6.2831853071795862 / (vr * vr + vi * vi);
    yr = ur * vr;
    yi = ur * vi;
  }

  (*res)[0] = yr;
  (*res)[1] = yi;
}

void hankel_kernel_FT_3D(double x, fftw_complex *res, double *arg, int argc __attribute__((unused)))
{
      fftw_complex a1, a2, g1, g2;
      double           mu;
      double        mod, xln2, si, co, d1, d2, pref, q;
      q = arg[0];
      mu = arg[1];

      /* arguments for complex gamma */
      a1[0] = 0.5*(1.0+mu+q);
      a2[0] = 0.5*(1.0+mu-q);
      a1[1] = 0.5*x; a2[1]=-a1[1];
      cdgamma(a1,&g1);
      cdgamma(a2,&g2);
      xln2 = x*M_LN2;
      si   = sin(xln2);
      co   = cos(xln2);
      d1   = g1[0]*g2[0]+g1[1]*g2[1]; /* Re */
      d2   = g1[1]*g2[0]-g1[0]*g2[1]; /* Im */
      mod  = g2[0]*g2[0]+g2[1]*g2[1];
      pref = exp(M_LN2*q)/mod;

      (*res)[0] = pref*(co*d1-si*d2);
      (*res)[1] = pref*(si*d1+co*d2);
}
