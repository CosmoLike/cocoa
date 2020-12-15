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
#include "basics.h"

#include "log.c/src/log.h"

con constants = {
    3.14159265358979323846, // pi
    9.86960440108935861883, // pisqr
    6.28318530717958647693, // twopi
    0.69314718,
    2.90888208665721580e-4, // arcmin
    299792.458              // speed of light km/s
};

pre precision = {
    1e-2, // low
    1e-3, // medium
    1e-5, // high
    1e-7  // insane
};

lim limits = {
    //	0.19, //a_min (in order to compute z=4 WFIRST)
    1. / (1. + 10.), // a_min (z=10, needed for CMB lensing)
    6.667e-6,        // k_min_mpc
    1.e3,            // k_max_mpc
    50.,             // k_max_mpc_class
    2.e-2,           // k_min_cH0
    3.e+6,           // k_max_cH0
    0.1,             // P_2_s_min
    1.0e5,           // P_2_s_max
    3.0e-7,          // xi_via_hankel_theta_min
    0.12,            // xi_via_hankel_theta_max
    1.e-5,           // xi_3d_rmin
    1.e+0,           // xi_3d_rmax
    1.e+6,           // M_min
    1.e+17,          // M_max
};

Ntab Ntable = {
    100,  // N_a !!!!!!DO NOT F!@#$ING DECREASE THIS NUMBER UNLESS YOU ARE
          // ELISABETH AND HAVE ASKED TIM BEFORE!!!!
    500,  // N_k_lin
    500,  // N_k_nlin
    200,  // N_ell
    200,  // N_theta
    2048, // N_theta for Hankel
    1000, // N_S2
    1000, // N_DS
    50,   // N_norm
    50,   // N_r_3d
    25,   // N_k_3d
    20,   // N_a_halo
};

void SVD_inversion(gsl_matrix *cov, gsl_matrix *inverseSVD, int Nmatrix) {
  int i, j;
  gsl_matrix *V = gsl_matrix_calloc(Nmatrix, Nmatrix);
  gsl_matrix *U = gsl_matrix_calloc(Nmatrix, Nmatrix);
  gsl_vector *S = gsl_vector_calloc(Nmatrix);
  gsl_vector *work = gsl_vector_calloc(Nmatrix);

  gsl_matrix_memcpy(U, cov);

  gsl_linalg_SV_decomp(U, V, S, work);

  for (i = 0; i < Nmatrix; i++) {
    gsl_vector *b = gsl_vector_calloc(Nmatrix);
    gsl_vector *x = gsl_vector_calloc(Nmatrix);
    gsl_vector_set(b, i, 1.0);
    gsl_linalg_SV_solve(U, V, S, b, x);
    for (j = 0; j < Nmatrix; j++) {
      gsl_matrix_set(inverseSVD, i, j, gsl_vector_get(x, j));
    }
    gsl_vector_free(b);
    gsl_vector_free(x);
  }
  gsl_matrix_free(V);
  gsl_matrix_free(U);
  gsl_vector_free(S);
  gsl_vector_free(work);
}

void invert_matrix_colesky(gsl_matrix *A) {
  gsl_linalg_cholesky_decomp(A); // Adummy will be overwritten  */
  gsl_linalg_cholesky_invert(A);
}

double int_gsl_integrate_high_precision(double (*func)(double, void *),
void *arg, double a, double b, double *error, int niter) {
  double res, err;
  gsl_integration_cquad_workspace *w =
      gsl_integration_cquad_workspace_alloc(niter);
  gsl_function F;
  F.function = func;
  F.params = arg;
  gsl_integration_cquad(&F, a, b, 0, precision.high, w, &res, &err, 0);
  if (NULL != error)
    *error = err;
  gsl_integration_cquad_workspace_free(w);
  return res;
}

double int_gsl_integrate_medium_precision(double (*func)(double, void *),
void *arg, double a, double b, double *error, int niter) {
  double res, err;
  gsl_integration_cquad_workspace *w =
    gsl_integration_cquad_workspace_alloc(niter);
  gsl_function F;
  F.function = func;
  F.params = arg;
  gsl_integration_cquad(&F, a, b, 0, precision.medium, w, &res, &err, 0);
  if (NULL != error) {
    *error = err;
  }
  gsl_integration_cquad_workspace_free(w);
  return res;
}

double int_gsl_integrate_low_precision(double (*func)(double, void *),
void *arg, double a, double b, double *error,
int niter __attribute__((unused))) {
  gsl_set_error_handler_off();
  double res, err;
  gsl_function F;
  F.function = func;
  F.params = arg;
  size_t neval;
  gsl_integration_qng(&F, a, b, 0, precision.medium, &res, &err, &neval);
  return res;
}

// allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]
double **create_double_matrix(long nrl, long nrh, long ncl, long nch) {
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double **m;

  // allocate pointers to rows
  m = (double **)calloc(nrow + NR_END, sizeof(double *));
  if (!m) {
    log_fatal("allocation failure 1 in create_double_matrix()");
    exit(1);
  }
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (double *)calloc(nrow * ncol + NR_END, sizeof(double));
  if (!m[nrl]) {
    log_fatal("allocation failure 2 in create_double_matrix()");
    exit(1);
  }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;

  // return pointer to array of pointers to rows
  return m;
}

// free a double matrix allocated by create_double_matrix()
void free_double_matrix(double **m, long nrl, long nrh __attribute__((unused)),
long ncl, long nch __attribute__((unused))) {
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

// allocate a double vector with subscript range v[nl..nh]
double *create_double_vector(long nl, long nh) {
  double *v;
  v = (double *) calloc(nh - nl + 1 + NR_END, sizeof(double));
  if (!v) {
    log_fatal("allocation failure in double vector()");
    exit(1);
  }
  return v - nl + NR_END;
}

// allocate a int vector with subscript range v[nl..nh]
int *int_vector(long nl, long nh) {
  int *v;
  v = (int *)calloc(nh - nl + 1 + NR_END, sizeof(int));
  if (!v) {
    log_fatal("allocation failure in int vector()");
    exit(1);
  }
  return v - nl + NR_END;
}

// allocate a int vector with subscript range v[nl..nh]
long *long_vector(long nl, long nh) {
  long *v;
  v = (long *)calloc(nh - nl + 1 + NR_END, sizeof(long));
  if (!v) {
    log_fatal("allocation failure in int vector()");
    exit(1);
  }
  return v - nl + NR_END;
}

// free a double vector allocated with vector()
void free_double_vector(double *v, long nl, long nh __attribute__((unused))) {
  free((FREE_ARG)(v + nl - NR_END));
}

void error(char *s) {
  printf("error:%s\n ", s);
  exit(1);
}

/* ============================================================ *
 * Interpolates f at the value x, where f is a double[n] array,	*
 * representing a function between a and b, stepwidth dx.	*
 * 'lower' and 'upper' are powers of a logarithmic power law	*
 * extrapolation. If no	extrapolation desired, set these to 0	*
 * ============================================================ */
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
        return f[i]; /* constant extrapolation */
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
                  int ny, double ay, double by, double dy, double y,
                  double lower, double upper) {
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
    // printf("%d %d\n",i+1,j+1);
    return (1. - dt) * (1. - ds) * f[i][j];
  }
  if (i + 1 == nx) {
    // printf("%d %d\n",i+1,j+1);
    return (1. - dt) * (1. - ds) * f[i][j] + (1. - dt) * ds * f[i][j + 1];
  }
  if (j + 1 == ny) {
    // printf("%d %d\n",i+1,j+1);
    return (1. - dt) * (1. - ds) * f[i][j] + dt * (1. - ds) * f[i + 1][j];
  }
  return (1. - dt) * (1. - ds) * f[i][j] + (1. - dt) * ds * f[i][j + 1] +
         dt * (1. - ds) * f[i + 1][j] + dt * ds * f[i + 1][j + 1];
}

double interpol2d_fitslope(double **f, int nx, double ax, double bx, double dx,
                           double x, int ny, double ay, double by, double dy,
                           double y, double lower) {
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
    // printf("%d %d\n",i+1,j+1);
    return (1. - dt) * (1. - ds) * f[i][j];
  }
  if (i + 1 == nx) {
    // printf("%d %d\n",i+1,j+1);
    return (1. - dt) * (1. - ds) * f[i][j] + (1. - dt) * ds * f[i][j + 1];
  }
  if (j + 1 == ny) {
    // printf("%d %d\n",i+1,j+1);
    return (1. - dt) * (1. - ds) * f[i][j] + dt * (1. - ds) * f[i + 1][j];
  }
  return (1. - dt) * (1. - ds) * f[i][j] + (1. - dt) * ds * f[i][j + 1] +
         dt * (1. - ds) * f[i + 1][j] + dt * ds * f[i + 1][j + 1];
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

/////// Hankel transform routine used in cosmo2D_xx/////////////////
void hankel_kernel_FT(double x, fftw_complex *res, double *arg, int argc __attribute__((unused))) {
  fftw_complex a1, a2, g1, g2;
  int mu;
  double mod, xln2, si, co, d1, d2, pref, q;
  q = arg[0];
  mu = (int)(arg[1] + 0.1);

  /* arguments for complex gamma */
  a1[0] = 0.5 * (1.0 + mu + q);
  a2[0] = 0.5 * (1.0 + mu - q);
  a1[1] = 0.5 * x;
  a2[1] = -a1[1];
  cdgamma(a1, &g1);
  cdgamma(a2, &g2);
  xln2 = x * constants.ln2;
  si = sin(xln2);
  co = cos(xln2);
  d1 = g1[0] * g2[0] + g1[1] * g2[1]; /* Re */
  d2 = g1[1] * g2[0] - g1[0] * g2[1]; /* Im */
  mod = g2[0] * g2[0] + g2[1] * g2[1];
  pref = exp(constants.ln2 * q) / mod;

  (*res)[0] = pref * (co * d1 - si * d2);
  (*res)[1] = pref * (si * d1 + co * d2);
}
