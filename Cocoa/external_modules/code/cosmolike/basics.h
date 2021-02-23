#include <fftw3.h>
#include <gsl/gsl_matrix.h>

#ifndef __COSMOLIKE_BASICS_H
#define __COSMOLIKE_BASICS_H
#ifdef __cplusplus
extern "C" {
#endif

#define NR_END 1
#define FREE_ARG char *
#define EXIT_MISSING_FILE(ein, purpose, filename)                              \
  if (!ein) {                                                                  \
    fprintf(stderr, "Could not find %s file %s\n", purpose, filename);         \
    exit(1);                                                                   \
  }

static double darg __attribute__((unused)), maxarg1 __attribute__((unused)),
    maxarg2 __attribute__((unused));

#define FMAX(a, b)                                                             \
  (maxarg1 = (a), maxarg2 = (b), (maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define FMIN(a, b)                                                             \
  (maxarg1 = (a), maxarg2 = (b), (maxarg1) < (maxarg2) ? (maxarg1) : (maxarg2))

typedef struct {
  double pi;
  double pi_sqr;
  double twopi;
  double ln2;
  double arcmin;
  double lightspeed;
} con;

typedef struct {
  double low;
  double medium;
  double high;
  double insane;
} pre;

typedef struct {
  double a_min;
  double k_min_cH0;
  double k_max_cH0;
  double M_min;
  double M_max;
  int LMIN_tab; // Cosmo2D
  int LMAX; // Cosmo2D
  int LMAX_NOLIMBER; // Cosmo2D
} lim;

typedef struct {
  int N_a;
  int N_k_lin;
  int N_k_nlin;
  int N_ell;
  int N_theta;
  int N_thetaH;
  int N_S2;
  int N_DS;
  int N_ell_TATT; // Cosmo2D
} Ntab;

extern con constants;

extern pre precision;

extern lim limits;

extern Ntab Ntable;

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

double interpol(double *f, int n, double a, double b, double dx, double x,
  double lower, double upper);

double interpol_fitslope(double *f, int n, double a, double b, double dx,
  double x, double lower);

void free_double_vector(double *v, long nl, long nh);

double *create_double_vector(long nl, long nh);

void free_double_matrix(double **m, long nrl, long nrh, long ncl, long nch);

double **create_double_matrix(long nrl, long nrh, long ncl, long nch);

int line_count(char *filename);

void error(char *s);

void hankel_kernel_FT(double x, fftw_complex *res, double *arg,
int argc __attribute__((unused)));

void cdgamma(fftw_complex x, fftw_complex *res);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD