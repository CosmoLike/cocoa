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

gsl_spline* malloc_gsl_spline(const int n)
{
  gsl_spline* result;

  if (Ntable.photoz_interpolation_type == 0)
    result = gsl_spline_alloc(gsl_interp_cspline, n);
  else if (Ntable.photoz_interpolation_type == 1)
    result = gsl_spline_alloc(gsl_interp_linear, n);
  else
    result = gsl_spline_alloc(gsl_interp_steffen, n);

  if (result == NULL)
  {
    log_fatal("array allocation failed");
    exit(1);
  }
  return result;
}

gsl_integration_glfixed_table* malloc_gslint_glfixed(const int n)
{
  gsl_integration_glfixed_table* w = gsl_integration_glfixed_table_alloc(n);
  if (w == NULL)
  {
    log_fatal("array allocation failed");
    exit(1);
  }
  return w;
}

void*** malloc3d(const int nx, const int ny, const int nz)
{
  double*** tab = (double***) malloc(sizeof(double**)*nx +
                                     sizeof(double*)*nx*ny + 
                                     sizeof(double)*nx*ny*nz);
  if (tab == NULL)
  {
    log_fatal("array allocation failed");
    exit(1);
  }

  #pragma omp parallel for collapse(2)
  for (int i=0; i<nx; i++)
  {
    for (int j=0; j<ny; j++)
    {
      tab[i] = (double**)(tab + nx) + ny*i;
      tab[i][j] = (double*)((double**)(tab+nx) + nx*ny) + ny*nz*i + nz*j;
    }
  }

  return (void***) tab;
}

void** malloc2d(const int nx, const int ny)
{
  double** tab = (double**) malloc(sizeof(double*)*nx + 
                                   sizeof(double)*nx*ny);
  if (tab == NULL)
  {
    log_fatal("array allocation failed");
    exit(1);
  }

  #pragma omp parallel for
  for (int i=0; i<nx; i++)
    tab[i] = (double*)(tab + nx) + ny*i;
  
  return (void**) tab;
}

void* malloc1d(const int nx)
{
  double* vec = (double*) malloc(sizeof(double)*nx);
  if (vec == NULL)
  {
    log_fatal("array allocation failed");
    exit(1);
  }
  return (void*) vec;
}

void* calloc1d(const int nx)
{
  double* vec = (double*) calloc(nx, sizeof(double));
  if (vec == NULL)
  {
    log_fatal("array allocation failed");
    exit(1);
  }
  return (void*) vec;
}

int fdiff(const double a, const double b)
{
  return (fabs(a-b) < 1.0e-13 * fabs(a+b) || fabs(a-b) < 2.0e-38) ? 0 : 1;
}

double fmin(const double a, const double b)
{
  return a < b ? a : b;
}

double fmax(const double a, const double b)
{
  return a > b ? a : b;
}

bin_avg set_bin_average(const int i_theta, const int j_L)
{
  static double*** P  = NULL;
  static double** xminmax = NULL;
  static int ntheta = 0;
  
  if (Ntable.Ntheta == 0)
  {
    log_fatal("Ntable.Ntheta not initialized");
    exit(1);
  }
  if (!(i_theta < Ntable.Ntheta))
  {
    log_fatal("bad i_theta index");
    exit(1);
  }
  if (j_L > limits.LMAX)
  {
    log_fatal("bad j_L index");
    exit(1);
  }

  if (P == NULL || (ntheta != Ntable.Ntheta))
  {
    if (P != NULL) free(P);
    if (xminmax != NULL) free(xminmax);

    // Legendre computes l=0,...,lmax (inclusive)
    P  = (double***) malloc3d(4, Ntable.Ntheta, limits.LMAX+1);
    double** Pmin  = P[0]; double** Pmax  = P[1];
    double** dPmin = P[2]; double** dPmax = P[3];

    xminmax = (double**) malloc2d(2, Ntable.Ntheta);

    const double logdt = 
      (log(Ntable.vtmax)-log(Ntable.vtmin))/ (double) Ntable.Ntheta;
    
    for(int i=0; i<Ntable.Ntheta ; i++)
    {
      xminmax[0][i] = cos(exp(log(Ntable.vtmin) + (i + 0.0)*logdt));
      xminmax[1][i] = cos(exp(log(Ntable.vtmin) + (i + 1.0)*logdt));
    }

    #pragma omp parallel for
    for (int i=0; i<Ntable.Ntheta; i++)
    {
      if (abs(xminmax[0][i]) > 1)
      {
        log_fatal("logical error: Legendre argument xmin = %.3e>1", xminmax[0][i]);
        exit(1);
      }
      if (abs(xminmax[1][i]) > 1)
      {
        log_fatal("logical error: Legendre argument xmax = %.3e>1", xminmax[1][i]);
        exit(1);
      }
      
      int status = 
      gsl_sf_legendre_Pl_deriv_array(limits.LMAX, xminmax[0][i], Pmin[i], dPmin[i]);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }

      status = 
      gsl_sf_legendre_Pl_deriv_array(limits.LMAX, xminmax[1][i], Pmax[i], dPmax[i]);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      } 
    }

    ntheta = Ntable.Ntheta;
  }

  bin_avg r;
  r.xmin = xminmax[0][i_theta];
  r.xmax = xminmax[1][i_theta];
  r.Pmin = P[0][i_theta][j_L];
  r.Pmax = P[1][i_theta][j_L];
  r.dPmin = P[2][i_theta][j_L];
  r.dPmax = P[3][i_theta][j_L];
  return r;
}

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
  gsl_integration_cquad(&F, a, b, 0, precision.medium/2, &w, &res, &err, 0);
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

void error(char *s) 
{
  printf("error:%s\n ", s);
  exit(1);
}

// ============================================================
// Interpolates f at the value x, where f is a double[n] array,
// representing a function between a and b, stepwidth dx.
// 'lower' and 'upper' are powers of a logarithmic power law
// extrapolation. If no	extrapolation desired, set these to 0
// ============================================================

double interpol1d(
  const double* const f, 
  const int n, 
  const double a, 
  const double b, 
  const double dx, 
  const double x) 
{
  double ans;
  if (x < a) 
    ans = f[0] + (x - a);
  else
  {
    const double r = (x - a) / dx;
    const int i = (int) floor(r);
    if (i + 1 >= n) 
      ans = f[n - 1] + (x - b);
    else 
      ans = (r - i) * (f[i + 1] - f[i]) + f[i];
  }
  return ans;
}

double interpol(const double* const f, const int n, const double a, 
const double b, const double dx, const double x, const double lower, 
const double upper) 
{
  if (x < a) 
  {
    if (lower == 0.) 
    {
      return 0.0;
    }
    return f[0] + lower * (x - a);
  }

  const double r = (x - a) / dx;
  const int i = (int) floor(r);
  
  if (i + 1 >= n) 
  {
    if (upper == 0.0) 
    {
      if (i + 1 == n) 
      {
        return f[i]; // constant extrapolation
      } 
      else 
      {
        return 0.0;
      }
    } 
    else 
    {
      return f[n - 1] + upper * (x - b); // linear extrapolation
    }
  } 
  else 
  {
    return (r - i) * (f[i + 1] - f[i]) + f[i]; // interpolation
  }
}


int line_count(char* filename) 
{  
  FILE* ein = fopen(filename, "r");
  if (ein == NULL) 
  {
    log_fatal("File not open (%s)", filename);
    exit(1);
  }
  
  int ch = 0; 
  int prev = 0; 
  int nlines = 0;

  do 
  {
    prev = ch;
    
    ch = fgetc(ein);
    
    if (ch == '\n')
    {
      nlines++;
    }
  } while (ch != EOF);
  
  fclose(ein);
  
  // last line might not end with "\n". 
  // However, if previous character does, then the last line is empty
  if (ch != '\n' && prev != '\n' && nlines != 0) 
  {
    nlines++;
  }
  return nlines;
}

double interpol_fitslope(double *f, int n, double a, double b, double dx,
                         double x, double lower) 
{
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
int ny, double ay, double by, double dy, double y, double lower, double upper) 
{
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

double interpol2d_fitslope(
    double **f, 
    int nx, 
    double ax, 
    double bx, 
    double dx,
    double x, 
    int ny, 
    double ay, 
    double by, 
    double dy, 
    double y, 
    double lower
  ) 
{  
  if (x < ax || x > bx)
    return 0.;
  
  const double t = (x - ax) / dx;
  int i = (int)(floor(t));
  const double dt = (x - ax) / dx - i;
  
  if (y < ay) 
    return ((1. - dt) * f[i][0] + dt * f[i + 1][0]) + (y - ay) * lower;
  else if (y > by) 
  {
    int fitrange;

    if (ny > 25) 
      fitrange = 5;
    else 
      fitrange = (int) floor(ny / 5);

    const double upper = ((1. - dt) * (f[i][ny - 1] - f[i][ny - 1 - fitrange]) + 
      dt * (f[i + 1][ny - 1] - f[i + 1][ny - 1 - fitrange])) /(dy * fitrange);
    
    return ((1. - dt) * f[i][ny - 1] + dt * f[i + 1][ny - 1]) + (y - by) * upper;
  }
  
  const double s = (y - ay) / dy;
  const int j = (int)(floor(s));
  const double ds = s - j;
  
  if ((i + 1 == nx) && (j + 1 == ny)) 
    return (1. - dt) * (1. - ds) * f[i][j];
  
  if (i + 1 == nx) 
    return (1. - dt) * (1. - ds) * f[i][j] + (1. - dt) * ds * f[i][j + 1];

  if (j + 1 == ny)
    return (1. - dt) * (1. - ds) * f[i][j] + dt * (1. - ds) * f[i + 1][j];
  
  return (1. - dt) * (1. - ds) * f[i][j] + (1. - dt) * ds * f[i][j + 1] +
         dt * (1. - ds) * f[i + 1][j] + dt * ds * f[i + 1][j + 1];
}

void hankel_kernel_FT(double x, fftw_complex *res, double *arg,
int argc __attribute__((unused))) 
{
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

void cdgamma(fftw_complex x, fftw_complex *res) 
{
  double xr, xi, wr, wi, ur, ui, vr, vi, yr, yi, t;

  xr = (double) x[0];
  xi = (double) x[1];

  if (xr < 0) 
  {
    wr = 1 - xr;
    wi = -xi;
  } else 
  {
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
