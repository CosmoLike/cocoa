#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h> 
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#include "structs.h"
#include "basics.h"

#include "log.c/src/log.h"

gsl_interp* malloc_gsl_interp(const int n)
{
  gsl_interp* result;

  if (Ntable.photoz_interpolation_type == 0)
    result = gsl_interp_alloc(gsl_interp_cspline, n);
  else if (Ntable.photoz_interpolation_type == 1)
    result = gsl_interp_alloc(gsl_interp_linear, n);
  else
    result = gsl_interp_alloc(gsl_interp_steffen, n);

  if (result == NULL) {
    log_fatal("array allocation failed");
    exit(1);
  }
  return result;
}

gsl_spline* malloc_gsl_spline(const int n)
{
  gsl_spline* result;

  if (Ntable.photoz_interpolation_type == 0)
    result = gsl_spline_alloc(gsl_interp_cspline, n);
  else if (Ntable.photoz_interpolation_type == 1)
    result = gsl_spline_alloc(gsl_interp_linear, n);
  else
    result = gsl_spline_alloc(gsl_interp_steffen, n);

  if (result == NULL) {
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

/*
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
*/

void*** malloc3d(const int nx, const int ny, const int nz)
{ // Got help from ChatGPT to do the align version of my previous malloc func
  void* raw_block = NULL;
  if (posix_memalign(&raw_block, 64, nx * sizeof(double**) + 
                                     nx * ny * sizeof(double*) + 
                                     nx * ny * nz * sizeof(double)) != 0) {
    log_fatal("posix_memalign failed in malloc3d_aligned64");
    exit(EXIT_FAILURE);
  }

  double*** tab = (double***) raw_block;

  #pragma omp parallel for collapse(2)
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      tab[i]    = (double**) ((char*) raw_block + nx*sizeof(double**)) + ny*i;
      tab[i][j] = (double*)  ((char*) raw_block + 
                              nx*sizeof(double**) + 
                              nx*ny*sizeof(double*)) + nz*(ny*i + j);
    }
  }

  return (void***) tab;
}

/*
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
  for (int i=0; i<nx; i++) {
    tab[i] = (double*)(tab + nx) + ny*i;
  }
  
  return (void**) tab;
}
*/

void** malloc2d(const int nx, const int ny)
{ // Got help from ChatGPT to do the align version of my previous malloc func
  const size_t pointer_bytes = nx * sizeof(double*);
  const size_t data_bytes    = nx * ny * sizeof(double);
  const size_t total_bytes   = pointer_bytes + data_bytes;

  void* raw_block = NULL;
  if (posix_memalign(&raw_block, 
                     64, 
                     sizeof(double*)*nx+sizeof(double)*nx*ny) != 0) {
    log_fatal("posix_memalign failed for malloc2d");
    exit(1);
  }

  double** tab = (double**) raw_block;

  #pragma omp parallel for
  for (int i = 0; i < nx; ++i) {
    tab[i] = (double*) ((char*) raw_block + nx * sizeof(double*)) + i * ny;
  }
  return (void**) tab;
}

/*
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
*/

void* malloc1d(const int nx)
{ // Got help from ChatGPT to do the align version of my previous malloc func
  void* vec = NULL;
  if (posix_memalign(&vec, 64, sizeof(double) * nx) != 0) {
    log_fatal("array allocation failed (malloc1d)");
    exit(1); 
  }
  return vec;
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
  static double cache[MAX_SIZE_ARRAYS];

  if (Ntable.Ntheta == 0) {
    log_fatal("Ntable.Ntheta not initialized");
    exit(1);
  }
  if (P == NULL || (ntheta != Ntable.Ntheta) || fdiff(cache[0], Ntable.random))
  {
    if (P != NULL) {
      free(P);
    }
    if (xminmax != NULL) {
      free(xminmax);
    }

    // Legendre computes l=0,...,lmax (inclusive)
    P  = (double***) malloc3d(4, Ntable.Ntheta, Ntable.LMAX+1);
    double** Pmin  = P[0]; double** Pmax  = P[1];
    double** dPmin = P[2]; double** dPmax = P[3];

    xminmax = (double**) malloc2d(2, Ntable.Ntheta);

    const double logdt = (log(Ntable.vtmax)-log(Ntable.vtmin))/ Ntable.Ntheta;
    for(int i=0; i<Ntable.Ntheta ; i++) {
      xminmax[0][i] = cos(exp(log(Ntable.vtmin) + (i + 0.)*logdt));
      xminmax[1][i] = cos(exp(log(Ntable.vtmin) + (i + 1.)*logdt));
    }

    #pragma omp parallel for
    for (int i=0; i<Ntable.Ntheta; i++) {
      if (abs(xminmax[0][i]) > 1) {
        log_fatal("logical error: Legendre argument xmin = %.3e>1", xminmax[0][i]);
        exit(1);
      }
      if (abs(xminmax[1][i]) > 1) {
        log_fatal("logical error: Legendre argument xmax = %.3e>1", xminmax[1][i]);
        exit(1);
      }
      
      int status = 
      gsl_sf_legendre_Pl_deriv_array(Ntable.LMAX, xminmax[0][i], Pmin[i], dPmin[i]);
      if (status) {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
      status = 
      gsl_sf_legendre_Pl_deriv_array(Ntable.LMAX, xminmax[1][i], Pmax[i], dPmax[i]);
      if (status) {
        log_fatal(gsl_strerror(status));
        exit(1);
      } 
    }
    ntheta = Ntable.Ntheta;
    cache[0] = Ntable.random;
  }
  if (!(i_theta < Ntable.Ntheta)) {
    log_fatal("bad i_theta index");
    exit(1);
  }
  if (j_L > Ntable.LMAX) {
    log_fatal("bad j_L index");
    exit(1);
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
  {  
    ans = f[0]; // constant extrapolation
  }
  else
  {
    const double r = (x - a) / dx;
    const int i = (int) floor(r);
    if (i + 1 >= n)
    {
      ans = f[n-1]; // constant extrapolation
    }
    else 
      ans = (r - i) * (f[i + 1] - f[i]) + f[i];
  }
  return ans;
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

double interpol2d(double **f, int nx, double ax, double bx, double dx, double x,
int ny, double ay, double by, double dy, double y) 
{
  double t, dt, s, ds;
  int i, j;
  
  if (x < ax) 
    return 0.;
  if (x > bx) 
    return 0.;

  t = (x - ax) / dx;
  i = (int)(floor(t));
  dt = t - i;
  
  if (y < ay) 
  {
    return ((1. - dt) * f[i][0] + dt * f[i + 1][0]) + (y - ay);
  } 
  else if (y > by) 
  {
    return ((1. - dt) * f[i][ny - 1] + dt * f[i + 1][ny - 1]) + (y - by);
  }
  s = (y - ay) / dy;
  j = (int)(floor(s));
  ds = s - j;
  
  if ((i + 1 == nx) && (j + 1 == ny)) 
  {
    return (1. - dt) * (1. - ds) * f[i][j];
  }
  if (i + 1 == nx) 
  {
    return (1. - dt) * (1. - ds) * f[i][j] + (1. - dt) * ds * f[i][j + 1];
  }
  if (j + 1 == ny) 
  {
    return (1. - dt) * (1. - ds) * f[i][j] + dt * (1. - ds) * f[i + 1][j];
  }
  
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
