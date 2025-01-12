#include <assert.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basics.h"
#include "cosmo3D.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// integration boundary routines
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double amin_source(int ni) 
{
  if (ni < -1 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  return 1. / (redshift.shear_zdist_zmax_all + 1.);
}

double amax_source(int i __attribute__((unused))) 
{
  return 1. / (1. + fmax(redshift.shear_zdist_zmin_all, 0.001));
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double amax_source_IA(int ni) 
{
  if (ni < -1 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  return 1. / (1. + fmax(redshift.shear_zdist_zmin_all, 0.001));
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double amin_lens(int ni) 
{
  if (ni < -1 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  
  if (ni == -1) 
  {
    return 1. / (redshift.clustering_zdist_zmax_all + 1.);
  }

  return 1. / (1 + redshift.clustering_zdist_zmax[ni] + 
                 2. * fabs(nuisance.bias_zphot_clustering[ni]));
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double amax_lens(int ni) 
{
  if (ni < -1 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  if (gbias.b_mag[ni] != 0) 
  {
    return 1. / (1. + fmax(redshift.shear_zdist_zmin_all, 0.001));
  }
  
  if (ni == -1) 
  {
    return 1. / (1. + fmax(redshift.clustering_zdist_zmin_all, 0.001));
  }
  
  return 1. / (1 + fmax(redshift.clustering_zdist_zmin[ni] -
      2. * fabs(nuisance.bias_zphot_clustering[ni]), 0.001));
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int test_kmax(double l, int ni) // return 1 if true, 0 otherwise
{ // test whether the (l, ni) bin is in the linear clustering regime
  static double chiref[10] = {-1.};
    
  if (chiref[0] < 0)
  {
    for (int i=0; i<redshift.clustering_nbin; i++)
    {
      chiref[i] = chi(1.0/(1. + 0.5 * (redshift.clustering_zdist_zmin[i] + 
                                       redshift.clustering_zdist_zmax[i])));
    }
  }
  
  if (ni < -1 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  
  const double R_min = like.Rmin_bias; // set minimum scale to which
                                       // we trust our bias model, in Mpc/h
  const double kmax = 2.0*M_PI / R_min * cosmology.coverH0;
  
  int res = 0.0;
  if ((l + 0.5) / chiref[ni] < kmax)
  {
    res = 1;
  }
  return res;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int test_zoverlap(int ni, int nj) // test whether source bin nj is behind lens bin ni
{ // Galaxy-Galaxy Lensing bins (redshift overlap tests)
  if (ni < 0 || ni > redshift.clustering_nbin - 1 || 
      nj < 0 || nj > redshift.shear_nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  return 1;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int ZL(int ni) 
{
  static int N[MAX_SIZE_ARRAYS*MAX_SIZE_ARRAYS] = {-42};
  if (N[0] < -1) 
  {
    int n = 0;
    for (int i=0; i<redshift.clustering_nbin; i++) 
    {
      for (int j=0; j<redshift.shear_nbin; j++) 
      {
        if (test_zoverlap(i, j)) 
        {
          N[n] = i;
          n++;
        }
      }
    }
  }
  
  if (ni < 0 || ni > tomo.ggl_Npowerspectra - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", ni, tomo.ggl_Npowerspectra);
    exit(1);
  }
  return N[ni];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int ZS(int nj) 
{
  static int N[MAX_SIZE_ARRAYS*MAX_SIZE_ARRAYS] = {-42};
  if (N[0] < -1) 
  {
    int n = 0;
    for (int i = 0; i < redshift.clustering_nbin; i++) 
    {
      for (int j = 0; j < redshift.shear_nbin; j++) 
      {
        if (test_zoverlap(i, j)) 
        {
          N[n] = j;
          n++;
        }
      }
    }
  }
  
  if (nj < 0 || nj > tomo.ggl_Npowerspectra - 1)
  {
    log_fatal("invalid bin input nj = %d (max %d)", nj, tomo.ggl_Npowerspectra);
    exit(1);
  }
  return N[nj];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int N_ggl(int ni, int nj) 
{ // ni = redshift bin of the lens, nj = redshift bin of the source
  static int N[MAX_SIZE_ARRAYS][MAX_SIZE_ARRAYS] = {{-42}};
  if (N[0][0] < 0) 
  {
    int n = 0;
    for (int i=0; i<redshift.clustering_nbin; i++) 
    {
      for (int j=0; j<redshift.shear_nbin; j++) 
      {
        if (test_zoverlap(i, j)) 
        {
          N[i][j] = n;
          n++;
        } 
        else 
        {
          N[i][j] = -1;
        }
      }
    }
  }
  
  if (ni < 0 || ni > redshift.clustering_nbin - 1 || 
      nj < 0 || nj > redshift.shear_nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
  return N[ni][nj];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int Z1(int ni) 
{ // find z1 of tomography combination (z1, z2) constituting shear tomo bin Nbin
  static int N[MAX_SIZE_ARRAYS*MAX_SIZE_ARRAYS] = {-42};
  if (N[0] < -1) 
  {
    int n = 0;
    for (int i=0; i < redshift.shear_nbin; i++) 
    {
      for (int j=i; j < redshift.shear_nbin; j++) 
      {
        N[n] = i;
        n++;
      }
    }
  }
  
  if (ni < 0 || ni > tomo.shear_Npowerspectra - 1)
  {
    log_fatal("invalid bin input ni = %d (max = %d)", 
      ni, tomo.shear_Npowerspectra);
    exit(1);
  }
  return N[ni];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int Z2(int nj) 
{ // find z2 of tomography combination (z1,z2) constituting shear tomo bin Nbin
  static int N[MAX_SIZE_ARRAYS*MAX_SIZE_ARRAYS] = {-42};
  if (N[0] < -1) 
  {
    int n = 0;
    for (int i=0; i<redshift.shear_nbin; i++) 
    {
      for (int j=i; j<redshift.shear_nbin; j++) 
      {
        N[n] = j;
        n++;
      }
    }
  }
  
  if (nj < 0 || nj > tomo.shear_Npowerspectra - 1)
  {
    log_fatal("invalid bin input nj = %d (max = %d)", 
      nj, tomo.shear_Npowerspectra);
    exit(1);
  }
  return N[nj];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int N_shear(int ni, int nj) 
{ // find shear tomography bin number N_shear of tomography combination (z1, z2)
  static int N[MAX_SIZE_ARRAYS][MAX_SIZE_ARRAYS] = {{-42}};
  if (N[0][0] < -1) 
  {
    int n = 0;
    for (int i=0; i<redshift.shear_nbin; i++) 
    {
      for (int j=i; j<redshift.shear_nbin; j++) 
      {
        N[i][j] = n;
        N[j][i] = n;
        n++;
      }
    }
  }

  const int ntomo = redshift.shear_nbin;
  if (ni < 0 || ni > ntomo - 1 || nj < 0 || nj > ntomo - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d) (max = %d)", ni, nj, ntomo);
    exit(1);
  }
  return N[ni][nj];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int ZCL1(int ni) 
{ // find ZCL1 of tomography combination (zcl1, zcl2) constituting tomo bin Nbin
  static int N[MAX_SIZE_ARRAYS*MAX_SIZE_ARRAYS] = {-42};
  if (N[0] < -1) 
  {
    int n = 0;
    for (int i=0; i<redshift.clustering_nbin; i++) 
    {
      for (int j=i; j<redshift.clustering_nbin; j++) 
      {
        N[n] = i;
        n++;
      }
    }
  }

  if (ni < 0 || ni > tomo.clustering_Npowerspectra - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", ni, 
      tomo.clustering_Npowerspectra);
    exit(1);
  }
  return N[ni];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int ZCL2(int nj) 
{ // find ZCL2 of tomography combination (zcl1, zcl2) constituting tomo bin Nbin
  static int N[MAX_SIZE_ARRAYS*MAX_SIZE_ARRAYS] = {-42};
  if (N[0] < -1) 
  {
    int n = 0;
    for (int i=0; i<redshift.clustering_nbin; i++) 
    {
      for (int j=i; j<redshift.clustering_nbin; j++) 
      {
        N[n] = j;
        n++;
      }
    }
  }

  if (nj < 0 || nj > tomo.clustering_Npowerspectra - 1)
  {
    log_fatal("invalid bin input nj = %d (max %d)", nj, 
      tomo.clustering_Npowerspectra);
    exit(1);
  }
  return N[nj];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int N_CL(int ni, int nj) 
{
  static int N[MAX_SIZE_ARRAYS][MAX_SIZE_ARRAYS] = {{-42}};
  if (N[0][0] < -1) 
  {
    int n = 0;
    for (int i=0; i<redshift.clustering_nbin; i++) 
    {
      for (int j=i; j<redshift.clustering_nbin; j++) 
      {
        N[i][j] = n;
        N[j][i] = n;
        n++;
      }
    }
  }

  const int ntomo = redshift.clustering_nbin;
  if (ni < 0 || ni > ntomo - 1 ||  nj < 0 || nj > ntomo - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d) (max = %d)", ni, nj, ntomo);
    exit(1);
  }
  return N[ni][nj];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Shear routines for redshift distributions
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double zdistr_histo_n(double z, const int ni)
{
  if (redshift.shear_zdist_table == NULL) 
  {
    log_fatal("redshift n(z) not loaded");
    exit(1);
  } 
  
  double res = 0.0;
  if ((z >= redshift.shear_zdist_zmin_all) && 
      (z < redshift.shear_zdist_zmax_all)) 
  {
    // alias
    const int ntomo = redshift.shear_nbin;
    const int nzbins = redshift.shear_nzbins;
    double** tab = redshift.shear_zdist_table;
    double* z_v = redshift.shear_zdist_table[ntomo];
    
    const double dz_histo = (z_v[nzbins - 1] - z_v[0]) / ((double) nzbins - 1.);
    const double zhisto_min = z_v[0];
    const double zhisto_max = z_v[nzbins - 1] + dz_histo;

    const int nj = (int) floor((z - zhisto_min) / dz_histo);
    
    if (ni < 0 || ni > ntomo - 1 || nj < 0 || nj > nzbins - 1)
    {
      log_fatal("invalid bin input (zbin = ni, bin = nj) = (%d, %d)", ni, nj);
      exit(1);
    } 
    res = tab[ni][nj];
  }
  return res;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double zdistr_photoz(double zz, const int nj) 
{
  static double** table = NULL;
  static gsl_spline* photoz_splines[MAX_SIZE_ARRAYS+1];
  
  if (table == NULL) 
  {  
    const int ntomo  = redshift.shear_nbin;
    const int nzbins = redshift.shear_nzbins;

    if (table != NULL)
    {
      free(table);
    }
    table = (double**) malloc2d(ntomo + 2, nzbins);
    
    const double zmin = redshift.shear_zdist_zmin_all;
    const double zmax = redshift.shear_zdist_zmax_all;
    const double dz_histo = (zmax - zmin) / ((double) nzbins);  
    for (int k=0; k<nzbins; k++) 
    { // redshift stored at zv = table[ntomo+1]
      table[ntomo+1][k] = zmin + (k + 0.5) * dz_histo;
    }
    
    double NORM[MAX_SIZE_ARRAYS];    
    double norm = 0; 
    #pragma omp parallel for reduction( + : norm )
    for (int i=0; i<ntomo; i++) 
    {
      NORM[i] = 0.0;
      for (int k=0; k<nzbins; k++) 
      {    
        const double z = table[ntomo+1][k];  
        NORM[i] += zdistr_histo_n(z, i) * dz_histo;
      }
      norm += NORM[i];
    }
        
    #pragma omp parallel for
    for (int k=0; k<nzbins; k++) 
    { 
      table[0][k] = 0; // store normalization in table[0][:]
      for (int i=0; i<ntomo; i++) 
      {
        const double z = table[ntomo+1][k];
        table[i + 1][k] = zdistr_histo_n(z, i)/NORM[i];
        table[0][k] += table[i+1][k] * NORM[i] / norm;
      }
    }

    for (int i=0; i<ntomo+1; i++) 
    {
      photoz_splines[i] = malloc_gsl_spline(nzbins);
    }

    #pragma omp parallel for
    for (int i=0; i<ntomo+1; i++) 
    {
      int status = gsl_spline_init(photoz_splines[i], 
                                   table[ntomo+1], // z_v = table[ntomo+1]
                                   table[i], 
                                   nzbins);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
    }
  }
  
  const int ntomo  = redshift.shear_nbin;
  const int nzbins = redshift.shear_nzbins;

  if (nj < -1 || nj > ntomo - 1) 
  {
    log_fatal("nj = %d bin outside range (max = %d)", nj, redshift.shear_nbin);
    exit(1);
  }
  
  zz = zz - nuisance.bias_zphot_shear[nj];
  
  double res; 
  if (zz <= table[ntomo+1][0] || zz >= table[ntomo+1][nzbins - 1])
  { // z_v = table[ntomo+1]
    res = 0.0;
  }
  else
  {
    int status = gsl_spline_eval_e(photoz_splines[nj+1], zz, NULL, &res);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
  }
  return res;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double int_for_zmean_source(double z, void* params) 
{
  double* ar = (double*) params;
  const int ni = (int) ar[0];
  
  if (ni < 0 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return z * zdistr_photoz(z, ni);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double zmean_source(int ni) 
{ // mean true redshift of source galaxies in tomography bin j
  static double* table = NULL;
  static Ntab numtable;

  if (table == NULL || recompute_table(numtable))
  {
    if (table != NULL)  free(table);
    table = (double*) malloc1d(redshift.shear_nbin);
   
    const size_t szint = 300 + 50 * (Ntable.high_def_integration);
    gsl_integration_glfixed_table* w = malloc_gslint_glfixed(szint);

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
    {
      double init = zdistr_photoz(0., 0);
    }
    #pragma GCC diagnostic pop
    #pragma GCC diagnostic pop

    #pragma omp parallel for
    for (int i=0; i<redshift.shear_nbin; i++) 
    {
      double ar[1] = {(double) i};
      gsl_function F;
      F.params = ar;
      F.function = int_for_zmean_source;
      
      table[i] = gsl_integration_glfixed(&F, 
        redshift.shear_zdist_zmin[i], redshift.shear_zdist_zmax[i], w);
    }

    gsl_integration_glfixed_table_free(w);
  }

  if (ni < 0 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  return table[ni];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Lenses routines for redshift distributions
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double pf_histo_n(double z, const int ni) 
{ // based file with structure z[i] nz[0][i] .. nz[redshift.clustering_nbin-1][i]
  if (redshift.clustering_zdist_table == NULL) 
  {
    log_fatal("redshift n(z) not loaded");
    exit(1);
  } 
  
  double res = 0.0;
  if ((z >= redshift.clustering_zdist_zmin_all) && 
      (z < redshift.clustering_zdist_zmax_all)) 
  {
    
    const int ntomo = redshift.clustering_nbin;                 // alias
    const int nzbins = redshift.clustering_nzbins;          // alias
    double** tab = redshift.clustering_zdist_table;         // alias
    double* z_v = redshift.clustering_zdist_table[ntomo];   // alias
    
    const double dz_histo = (z_v[nzbins - 1] - z_v[0]) / ((double) nzbins - 1.);
    const double zhisto_min = z_v[0];
    const double zhisto_max = z_v[nzbins - 1] + dz_histo;

    const int nj = (int) floor((z - zhisto_min) / dz_histo);
    
    if (ni < 0 || ni > ntomo - 1 || nj < 0 || nj > nzbins - 1)
    {
      log_fatal("invalid bin input (zbin = ni, bin = nj) = (%d, %d)", ni, nj);
      exit(1);
    } 
    res = tab[ni][nj];
  }
  return res;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double pf_photoz(double zz, int nj) 
{
  static double** table = NULL;
  static gsl_spline* photoz_splines[MAX_SIZE_ARRAYS+1];

  if (table == NULL) 
  {  
    const int ntomo  = redshift.clustering_nbin;          // alias
    const int nzbins = redshift.clustering_nzbins;    // alias
    
    if (table != NULL) free(table);
    table = (double**) malloc2d(ntomo + 2, nzbins);

    const double zmin = redshift.clustering_zdist_zmin_all;
    const double zmax = redshift.clustering_zdist_zmax_all;
    const double dz_histo = (zmax - zmin) / ((double) nzbins);  
    for (int k=0; k<nzbins; k++) 
    { // redshift stored at zv = table[ntomo+1]
      table[ntomo+1][k] = zmin + (k + 0.5) * dz_histo;
    }
        
    double NORM[MAX_SIZE_ARRAYS];
    double norm = 0;
    #pragma omp parallel for reduction( + : norm )
    for (int i=0; i<ntomo; i++) 
    {
      NORM[i] = 0.0;
      for (int k=0; k<nzbins; k++) 
      {    
        const double z = table[ntomo+1][k];  
        NORM[i] += pf_histo_n(z, i) * dz_histo;
      }
      norm += NORM[i];
    }

    #pragma omp parallel for
    for (int k=0; k<nzbins; k++) 
    { 
      table[0][k] = 0; // store normalization in table[0][:]
      for (int i=0; i<ntomo; i++) 
      {
        const double z = table[ntomo+1][k];
        table[i + 1][k] = pf_histo_n(z, i)/NORM[i];
        table[0][k] += table[i+1][k] * NORM[i] / norm;
      }
    }

    for (int i=0; i<ntomo+1; i++) 
    {
      photoz_splines[i] = malloc_gsl_spline(nzbins);
    }

    #pragma omp parallel for
    for (int i=0; i<ntomo+1; i++) 
    {
      int status = gsl_spline_init(photoz_splines[i], 
                                   table[ntomo+1], // z_v = table[ntomo+1]
                                   table[i], 
                                   nzbins);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
    }
  }
  
  const int ntomo  = redshift.clustering_nbin;
  const int nzbins = redshift.clustering_nzbins;

  if (nj < 0 || nj > ntomo - 1) 
  {
    log_fatal("nj = %d bin outside range (max = %d)", nj, ntomo);
    exit(1);
  }
  
  zz = zz - nuisance.bias_zphot_clustering[nj];
  
  double res; 
  if (zz <= table[ntomo+1][0] || zz >= table[ntomo+1][nzbins - 1])
  { // z_v = table[ntomo+1]
    res = 0.0;
  }
  else
  {
    int status = gsl_spline_eval_e(photoz_splines[nj+1], zz, NULL, &res);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
  }
  return res;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double int_for_zmean(double z, void* params) 
{
  double* ar = (double*) params;
  const int ni = (int) ar[0];
  
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return z * pf_photoz(z, ni);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double norm_for_zmean(double z, void* params) 
{
  double* ar = (double*) params;
  const int ni = (int) ar[0];
  
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return pf_photoz(z, ni);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double zmean(const int ni)
{ // mean true redshift of galaxies in tomography bin j
  static double* table = NULL;
  static Ntab numtable;

  if (table == NULL || recompute_table(numtable))
  {
    if (table != NULL) free(table);
    table = (double*) malloc1d(redshift.clustering_nbin+1);

    const size_t szint = 300 + 50 * (Ntable.high_def_integration);
    gsl_integration_glfixed_table* w = malloc_gslint_glfixed(szint);

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // COCOA: init static variables.
      double init = pf_photoz(0., 0);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for
    for (int i=0; i<redshift.clustering_nbin; i++) 
    {
      double ar[1] = {(double) i};
      gsl_function F;
      F.params = ar;
      
      F.function = int_for_zmean;
      const double num = gsl_integration_glfixed(&F, 
        redshift.clustering_zdist_zmin[i], redshift.clustering_zdist_zmax[i], w);
      
      F.function = norm_for_zmean;
      const double den = gsl_integration_glfixed(&F, 
        redshift.clustering_zdist_zmin[i], redshift.clustering_zdist_zmax[i], w);
      
      table[i] = num/den;
    }

    gsl_integration_glfixed_table_free(w);
    update_table(&numtable);
  }

  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  
  return table[ni];
}

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// Bin-averaged lens efficiencies
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

double int_for_g_tomo(double aprime, void* params) 
{
  if (!(aprime>0) || !(aprime<1))
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double *ar = (double*)params;
  
  const int ni = (int) ar[0];
  
  if (ni < 0 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  
  const double chi1 = ar[1];
  const double chi_prime = chi(aprime);
  return zdistr_photoz(1. / aprime - 1., ni) * f_K(chi_prime - chi1) /
    f_K(chi_prime) / (aprime * aprime);
}

double g_tomo(double ainput, int ni) // for tomography bin ni 
{  
  static nuisancepara N;
  static cosmopara C;
  static Ntab numtable;
  static double** table = NULL;

  const double amin = 1.0/(redshift.shear_zdist_zmax_all + 1.0);
  const double amax = 0.999999;
  const double da = (amax - amin)/((double) Ntable.N_a - 1.0);

  if (table == NULL || recompute_table(numtable)) 
  {
    if (table != NULL) free(table);
    table = (double**) malloc2d(redshift.shear_nbin, Ntable.N_a);
  }

  if (recompute_zphot_shear(N) || 
      recompute_cosmo3D(C) || 
      recompute_table(numtable)) 
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // COCOA: init static variables - allows the OpenMP on the next loop
      double ar[2];
      ar[0] = (double) 0;
      ar[1] = chi(amin);
      double trash = int_for_g_tomo(amin, (void*) ar);
    }
    #pragma GCC diagnostic pop
    
    const size_t szint = 250 + 50 * (Ntable.high_def_integration);
    gsl_integration_glfixed_table* w = malloc_gslint_glfixed(szint);

    #pragma omp parallel for collapse(2)
    for (int j=0; j<redshift.shear_nbin; j++) 
    {
      for (int i=0; i<Ntable.N_a; i++) 
      {
        const double a = amin + i*da;
        
        double ar[2];
        ar[0] = (double) j;          
        ar[1] = chi(a);

        gsl_function F;
        F.params = ar;
        F.function = int_for_g_tomo;
        
        table[j][i] = gsl_integration_glfixed(&F, amin, a, w);
      }
    }

    gsl_integration_glfixed_table_free(w);
    update_cosmopara(&C);
    update_nuisance(&N);
    update_table(&numtable);
  }

  if (ni < 0 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 

  double res;
  if (ainput <= amin || ainput > 1.0 - da) 
  {
    res = 0.0;
  }
  else
  {
    res = interpol(table[ni], Ntable.N_a, amin, amax, da, ainput, 1.0, 1.0); 
  }
  return res;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double int_for_g2_tomo(double aprime, void* params) 
{ // \int n(z') W(z,z')^2 routines for source clustering
  if (!(aprime>0) || !(aprime<1))
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double *ar = (double*) params;
  
  const int ni = (int) ar[0];
  
  if (ni < 0 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  
  const double chi1 = ar[1];
  const double chi_prime = chi(aprime);
  
  return zdistr_photoz(1./aprime-1., ni) * f_K(chi_prime - chi1)/
    f_K(chi_prime - chi1)/(f_K(chi_prime)*f_K(chi_prime))/(aprime*aprime);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double g2_tomo(double a, int ni) 
{ // for tomography bin ni
  static nuisancepara N;
  static cosmopara C;
  static Ntab numtable;
  static double** table = NULL;

  const double amin = 1.0/(redshift.shear_zdist_zmax_all + 1.0);
  const double amax = 0.999999;
  const double da = (amax - amin)/((double) Ntable.N_a - 1.0);

  if (table == 0 || recompute_table(numtable)) 
  {
    if (table != NULL) free(table);
    table = (double**) malloc2d(redshift.shear_nbin, Ntable.N_a);
  }

  if (recompute_zphot_shear(N) || 
      recompute_cosmo3D(C) || 
      recompute_table(numtable)) 
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // init static variables
      double ar[2];
      ar[0] = (double) 0;
      ar[1] = chi(amin);
      double trash = int_for_g2_tomo(amin, (void*) ar);
    }
    #pragma GCC diagnostic pop

    const size_t szint = 250 + 50 * (Ntable.high_def_integration);
    gsl_integration_glfixed_table* w = malloc_gslint_glfixed(szint);

    #pragma omp parallel for collapse(2)
    for (int j=0; j<redshift.shear_nbin; j++) 
    {
      for (int i=0; i<Ntable.N_a; i++) 
      {
        const double a = amin + i*da;
        
        double ar[2];
        ar[0] = (double) j;          
        ar[1] = chi(a);

        gsl_function F;
        F.params = ar;
        F.function = int_for_g2_tomo;

        table[j][i] = gsl_integration_glfixed(&F, amin, a, w);
      } 
    }

    gsl_integration_glfixed_table_free(w);   
    update_nuisance(&N);
    update_cosmopara(&C);
    update_table(&numtable);
  }

  if (ni < 0 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 

  double res;
  if (a <= amin || a > 1.0 - da) 
  {
    res = 0.0;
  }
  else 
  {
    res = interpol(table[ni], Ntable.N_a, amin, amax, da, a, 1.0, 1.0);
  }  
  return res;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double int_for_g_lens(double aprime, void* params) 
{
  double *ar = (double*) params;

  const int ni = (int) ar[0];
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  const double chi1 = ar[1];
  const double chi_prime = chi(aprime);
  
  return pf_photoz(1. / aprime - 1., ni) * f_K(chi_prime - chi1) /
        f_K(chi_prime) / (aprime * aprime);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double g_lens(double a, int ni) 
{ // for *lens* tomography bin ni
  static nuisancepara N;
  static cosmopara C;
  static Ntab numtable;
  static double** table = NULL;

  const double amin = 1.0/(redshift.clustering_zdist_zmax_all + 1.0);
  const double amax = 0.999999;
  const double da = (amax - amin)/((double) Ntable.N_a - 1.0);
  const double amin_shear = 1. / (redshift.shear_zdist_zmax_all + 1.);

  if (table == NULL || recompute_table(numtable)) 
  {
    if (table != NULL)  free(table);
    table = (double**) malloc2d(redshift.clustering_nbin , Ntable.N_a);
  }

  if (recompute_zphot_clustering(N) || 
      recompute_cosmo3D(C) || 
      recompute_table(numtable)) 
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // COCOA: init static variables - allows the OpenMP on the next loop
      double ar[2];
      ar[0] = (double) 0;
      ar[1] = chi(amin);
      double trash = int_for_g_lens(amin_shear, (void*) ar);
    }
    #pragma GCC diagnostic pop 

    const size_t szint = 275 + 50 * (Ntable.high_def_integration);
    gsl_integration_glfixed_table* w = malloc_gslint_glfixed(szint);

    #pragma omp parallel for collapse(2)
    for (int j=0; j<redshift.clustering_nbin; j++) 
    {
      for (int i=0; i<Ntable.N_a; i++) 
      {
        const double a =  amin + i*da;
        double ar[2];
        ar[0] = (double) j;
        ar[1] = chi(a);

        gsl_function F;
        F.params = ar;
        F.function = int_for_g_lens;

        table[j][i] = gsl_integration_glfixed(&F, amin_shear, a, w);
      }
    }

    gsl_integration_glfixed_table_free(w);
    update_table(&numtable);
    update_nuisance(&N);
    update_cosmopara(&C);
  }

  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  double res;
  if (a < amin || a > 1.0 - da) 
  {
    res = 0.0;
  }
  else 
  {
    res = interpol(table[ni], Ntable.N_a, amin, amax, da, a, 1.0, 1.0); 
  }
  return res;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double g_cmb(double a) 
{
  static cosmopara C;
  static double chi_cmb = 0.;
  static double fchi_cmb = 0.;
  if (recompute_cosmo3D(C)) 
  {
    chi_cmb = chi(1./1091.);
    fchi_cmb = f_K(chi_cmb);
    update_cosmopara(&C);
  }
  return f_K(chi_cmb - chi(a))/fchi_cmb;
}





















































/*


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// TODO: CLUSTER ROUTINES BELOW ARE IN ALPHA STAGE
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// cluster routines (for redshift distributions, including photo-zs)
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double pf_cluster_histo_n(double z,  void* params) 
{
  static int njmax;
  static double** table;
  static double zhisto_max, zhisto_min, dz;

  if (table == 0)
  {
    const int zbins = line_count(redshift.clusters_REDSHIFT_FILE);
    njmax = zbins;
    table = (double**) malloc(sizeof(double*)*redshift.cluster_nbin);
    for (int i=0; i<redshift.cluster_nbin; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*zbins);
    }
    double* z_v = (double*) malloc(sizeof(double)*zbins);
    
    FILE* ein = fopen(redshift.clusters_REDSHIFT_FILE, "r");
    if (ein == NULL)
    {
      log_fatal("file not opened");
      exit(1);
    }
    int p = 0;
    for (int i=0; i<zbins; i++)
    {
      fscanf(ein, "%le", &z_v[i]);
      p++;
      if (i > 0 && z_v[i] < z_v[i-1]) 
      {
        break;
      }
      for (int k=0; k<redshift.cluster_nbin; k++)
      {
        fscanf(ein," %le", &table[k][i]);
      }
    }
    fclose(ein);

    dz = (z_v[p - 1] - z_v[0]) / ((double) p - 1.0);
    zhisto_max = z_v[p - 1] + dz;
    zhisto_min = z_v[0];
    
    for (int k=0; k<redshift.cluster_nbin; k++)
    { // now, set tomography bin boundaries
      double max = table[k][0];
      for (int i=1; i<zbins; i++)
      {
        if (table[k][i]> max)
        {
          max = table[k][i];
        }
      }
      {
        int i = 0;
        while (table[k][i] <1.e-8*max && i < zbins-2)
        {
          i++;
        }
        tomo.cluster_zmin[k] = z_v[i];
      }
      {
        int i = zbins-1;
        while (table[k][i] <1.e-8*max && i > 0)
        {
          i--;
        }
        tomo.cluster_zmax[k] = z_v[i];
      }
      log_info("tomo.cluster_zmin[%d] = %.3f,tomo.cluster_zmax[%d] = %.3f",
        k, tomo.cluster_zmin[k], k, tomo.cluster_zmax[k]);
    }
    free(z_v);
    
    if (zhisto_max < tomo.cluster_zmax[redshift.cluster_nbin - 1] || zhisto_min > tomo.cluster_zmin[0])
    {
      log_fatal("%e %e %e %e", zhisto_min,tomo.cluster_zmin[0], 
        zhisto_max,tomo.cluster_zmax[redshift.cluster_nbin-1]);
      log_fatal("pf_cluster_histo_n: redshift file = %s is incompatible with bin choice", 
        redshift.clusters_REDSHIFT_FILE);
      exit(1);
    }
  }

  double res = 0.0;
  if ((z >= zhisto_min) && (z < zhisto_max)) 
  {
    double *ar = (double*) params;
    const int ni = (int) ar[0];
    const int nj = (int) floor((z - zhisto_min)/dz);

    if (ni < 0 || ni > redshift.cluster_nbin - 1 || nj < 0 || nj > njmax - 1)
    {
      log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
      exit(1);
    } 
    res = table[ni][nj];
  }
  return res;
}

double pz_cluster(const double zz, const int nz)
{
  static double** table = 0;
  static double* z_v = 0;
  static int zbins = -1;
  static gsl_spline* photoz_splines[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = 0;

  if (nz < -1 || nz > redshift.cluster_nbin - 1)
  {
    log_fatal("invalid bin input nz = %d (max %d)", nz, redshift.cluster_nbin);
    exit(1);
  }

  if (redshift.clusters_photoz == 0)
  {
    if ((zz >= tomo.cluster_zmin[nz]) & (zz <= tomo.cluster_zmax[nz])) 
    {
      return 1;
    }
    else
    { 
      return 0;
    }
  }

  if (table == 0)
  {
    const size_t nsize_integration = 2250 + 500 * (Ntable.high_def_integration);
    w = gsl_integration_glfixed_table_alloc(nsize_integration);

    zbins = line_count(redshift.clusters_REDSHIFT_FILE);
 
    table = (double**) malloc(sizeof(double*)*(redshift.cluster_nbin+1));
    for (int i=0; i<(redshift.cluster_nbin+1); i++)
    {
      table[i] = (double*) malloc(sizeof(double)*zbins);
    }
    z_v = (double*) malloc(sizeof(double)*zbins);

    for (int i=0; i<redshift.clustering_nbin+1; i++) 
    {
      if (Ntable.photoz_interpolation_type == 0)
      {
        photoz_splines[i] = gsl_spline_alloc(gsl_interp_cspline, zbins);
      }
      else if (Ntable.photoz_interpolation_type == 1)
      {
        photoz_splines[i] = gsl_spline_alloc(gsl_interp_linear, zbins);
      }
      else
      {
        photoz_splines[i] = gsl_spline_alloc(gsl_interp_steffen, zbins);
      }
      if (photoz_splines[i] == NULL)
      {
        log_fatal("fail allocation");
        exit(1);
      }
    }

    FILE* ein = fopen(redshift.clusters_REDSHIFT_FILE,"r");
    int p = 0;
    for (int i=0; i<zbins; i++)
    {
      fscanf(ein, "%le", &z_v[i]);
      p++;
      if (i > 0 && z_v[i] < z_v[i-1]) 
      {
        break;
      }
      for (int k=0; k<redshift.cluster_nbin; k++)
      {
        double space;
        fscanf(ein,"%le",&space);
      }
    }
    fclose(ein);

    {
      const double zhisto_min = fmax(z_v[0], 1.e-5);
      const double zhisto_max = z_v[p - 1] + (z_v[p - 1] - z_v[0])/((double) p - 1.0);  
      const double da = (zhisto_max - zhisto_min)/((double) zbins);
      for (int i=0; i<zbins; i++)
      { 
        z_v[i] = zhisto_min + (i + 0.5)*da;
      }
    }
    
    { // init the function pf_cluster_histo_n
      double ar[1] = {(double) 0};
      pf_cluster_histo_n(0., (void*) ar);
    }

    double NORM[MAX_SIZE_ARRAYS]; 
    
    #pragma omp parallel for
    for (int i=0; i<redshift.cluster_nbin; i++)
    {
      double ar[1] = {(double) i};

      gsl_function F;
      F.params = (void*) ar;
      F.function = pf_cluster_histo_n;
      const double norm = gsl_integration_glfixed(&F, 1E-5, tomo.cluster_zmax[i] + 1.0, w) / 
        (tomo.cluster_zmax[i] - tomo.cluster_zmin[i]);
      
      if (norm == 0) 
      {
        log_fatal("pz_cluster: norm(nz = %d) = 0", i);
        exit(1);
      }

      for (int k=0; k<zbins; k++)
      { 
        table[i + 1][k] = pf_cluster_histo_n(z_v[k], (void*) ar)/norm;
      }
      NORM[i] = norm;
    }
    
    double norm = 0;
    for (int i=0; i<redshift.cluster_nbin; i++)
    { // calculate normalized overall redshift distribution (without bins), store in table[0][:]
      norm += NORM[i];
    }

    for (int k=0; k<zbins; k++)
    {
      table[0][k] = 0; 
      for (int i=0; i<redshift.cluster_nbin; i++)
      {
        table[0][k] += table[i+1][k]*NORM[i]/norm;
      }
    }

    #pragma omp parallel for
    for (int i=-1; i<redshift.cluster_nbin; i++) 
    {
      int status = gsl_spline_init(photoz_splines[i + 1], z_v, table[i + 1], zbins);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
    }
  }

  if (nz > redshift.cluster_nbin - 1 || nz < -1)
  {
    log_fatal("pz_cluster(z, %d) outside redshift.cluster_nbin range", nz);
    exit(1);
  }

  double res;
  if (zz <= z_v[0] || zz >= z_v[zbins - 1]) 
  {
    res = 0.0;
  }
  else
  {
    double result = 0.0;
    int status = gsl_spline_eval_e(photoz_splines[nz + 1], zz, NULL, &result);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    res = result;
  }
  return res;
}

double dV_cluster(double z, void* params)
{
  double* ar = (double*) params;
  const int nz = ar[0];
  const double a = 1.0/(1.0 + z);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  return (fK*fK/hoverh0)*pz_cluster(z, nz);
}

double norm_z_cluster(const int nz)
{
  static cosmopara C;
  static double* table;
  static gsl_integration_glfixed_table* w = 0;

  const int N_z = redshift.clustering_nbin;
  const double zmin = tomo.cluster_zmin[nz];
  const double zmax = tomo.cluster_zmax[nz];
  
  if (table == 0) 
  {
    table = (double*) malloc(sizeof(double)*N_z);
    const size_t nsize_integration = 2500 + 50 * (Ntable.high_def_integration);
    w = gsl_integration_glfixed_table_alloc(nsize_integration);
  }
  if (recompute_cosmo3D(C))
  {
    { // init static vars only 
      const int i = 0;
      double params[1] = {0.0};
      dV_cluster(tomo.cluster_zmin[i], (void*) params);
    }
    #pragma omp parallel for
    for (int i=0; i<N_z; i++)
    {
      double ar[1] = {(double) i};

      gsl_function F;
      F.params = (void*) ar;
      F.function = dV_cluster;
      table[i] = gsl_integration_glfixed(&F, zmin, zmax, w);
    }
    update_cosmopara(&C);
  }
  return table[nz];
}

double zdistr_cluster(const int nz, const double z)
{ //simplfied selection function, disregards evolution of N-M relation+mass function within z bin
  double params[1] = {nz};
  return dV_cluster(z, (void*) params)/norm_z_cluster(nz);
}


double int_for_g_lens_cl(double aprime, void* params)
{
  double* ar = (double*) params;
  const int ni = (int) ar[0];
  const double a = ar[1];
  const int nl = (int) ar[2];
  
  if (ni < -1 || ni > redshift.cluster_nbin - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", ni, redshift.cluster_nbin);
    exit(1);
  }
  if (nl < 0 || nl > Cluster.N200_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", nl, Cluster.N200_Nbin);
    exit(1);
  } 
  
  const double zprime = 1.0/aprime - 1.0;
  const double chi1 = chi(a);
  const double chiprime = chi(aprime);
  return zdistr_cluster(zprime, ni)*f_K(chiprime - chi1)/f_K(chiprime)/(aprime*aprime);
}

double g_lens_cluster(const double a, const int nz, const int nl)
{ 
  static cosmopara C;
  static double** table = 0;
  static gsl_integration_glfixed_table* w = 0;

  const int N_z = redshift.cluster_nbin;
  const int N_l = Cluster.N200_Nbin;
  if (nl < 0 || nl > N_l - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", nl, N_l);
    exit(1);
  } 
  
  const double amin = 1.0/(tomo.cluster_zmax[redshift.cluster_nbin - 1] + 1.);
  const double amax = 0.999999;
  const double da = (amax - amin)/((double) Ntable.N_a - 1.0);
  
  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*(redshift.cluster_nbin + 1));
    for (int i=0; i<redshift.cluster_nbin+1; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*Ntable.N_a);
    }

    const size_t nsize_integration = 300 + 50 * (Ntable.high_def_integration);
    w = gsl_integration_glfixed_table_alloc(nsize_integration);
  }

  if (recompute_cosmo3D(C)) // there is no nuisance bias/sigma parameters yet
  {
    { // init static vars only, if j=-1, no tomography is being done
      double ar[3];
      ar[0] = (double) -1; // j = -1 no tomography is being done
      ar[1] = amin;
      ar[2] = (double) 0;
      int_for_g_lens_cl(amin, (void*) ar);
      if (N_z > 0)
      {
        ar[0] = (double) 0;
        int_for_g_lens_cl(amin, (void*) ar);
      }
    }

    #pragma omp parallel for collapse(2)
    for (int j=-1; j<N_z; j++) 
    { 
      for (int i=0; i<Ntable.N_a; i++) 
      {
        const double aa = amin + i*da;
        double ar[3];
        ar[0] = (double) j; 
        ar[1] = aa;
        ar[2] = nl;

        gsl_function F;
        F.params = ar;
        F.function = int_for_g_lens_cl;

        table[j + 1][i] = 
          gsl_integration_glfixed(&F, 1.0/(redshift.shear_zdist_zmax_all + 1.0), aa, w);
      }      
    } 

    update_cosmopara(&C);
  }

  if (nz < -1 || nz > N_z - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", nz, N_z);
    exit(1);
  }
  if (a < amin || a > amax)
  {
    log_fatal("a = %e outside look-up table range [%e,%e]", a, amin, amax);
    exit(1);
  }
  return interpol(table[nz + 1], Ntable.N_a, amin, amax, da, a, 1.0, 1.0); 
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Cluster-Galaxy Cross Clustering
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int ZCGCL1(const int ni)
{
  if (ni < 0 || ni > tomo.cg_clustering_Npowerspectra - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", ni, tomo.cg_clustering_Npowerspectra);
    exit(1);
  }
  // we assume cluster bin = galaxy bin (no cross)
  return tomo.external_selection_cg_clustering[ni];
}

int ZCGCL2(const int nj)
{
  if (nj < 0 || nj > tomo.cg_clustering_Npowerspectra - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", nj, tomo.cg_clustering_Npowerspectra);
    exit(1);
  }
  // we assume cluster bin = galaxy bin (no cross)
  return tomo.external_selection_cg_clustering[nj];
}

int N_CGCL(const int ni, const int nj)
{ // ni = Cluster Nbin, nj = Galaxy Nbin
  static int N[MAX_SIZE_ARRAYS][MAX_SIZE_ARRAYS] = {{-42}};
  if (N[0][0] < 0)
  {
    int n = 0;
    for (int i=0; i<redshift.cluster_nbin; i++)
    {
      for (int j=0; j<redshift.clustering_nbin; j++)
      {
        if (i == j) // we are not considering cross spectrum
        {
          for (int k=0; k<tomo.cg_clustering_Npowerspectra; k++)
          {
            if (i == tomo.external_selection_cg_clustering[k])
            {
              N[i][j] = n;
              n++;
            }
            else
            {
              N[i][j] = -1;
            }
          }
        }
        else
        {
          N[i][j] = -1;
        }
      }
    }
  }
  if (ni < 0 || ni > redshift.cluster_nbin - 1 || nj < 0 || nj > redshift.clustering_nbin - 1)
  {
    log_fatal("invalid bin input (ni (cluster nbin), nj (galaxy nbin)) = (%d, %d)", ni, nj);
    exit(1);
  }
  return N[ni][nj];
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Cluster-Galaxy Lensing bins (redshift overlap tests)
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int test_zoverlap_c(int zc, int zs) // test whether source bin zs is behind lens bin zl
{
  if (redshift.shear_photoz < 4 && tomo.cluster_zmax[zc] <= tomo.shear_zmin[zs]) 
  {
    return 1;
  }
  if (redshift.shear_photoz == 4 && tomo.cluster_zmax[zc] < zmean_source(zs)) 
  {
    return 1;
  }
  return 0;
}

int ZCL(int ni) 
{
  static int N[MAX_SIZE_ARRAYS*MAX_SIZE_ARRAYS] = {-42};
  if (N[0] < -1) 
  {
    int n = 0;
    for (int i=0; i<redshift.cluster_nbin; i++) 
    {
      for (int j=0; j<redshift.shear_nbin; j++) 
      {
        if (test_zoverlap_c(i, j)) 
        {
          N[n] = i;
          n++;
        }
      }
    }
  }
  if (ni < 0 || ni > tomo.cgl_Npowerspectra - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", ni, tomo.cgl_Npowerspectra);
    exit(1);
  }
  return N[ni];
}

int ZCS(int nj) 
{
  static int N[MAX_SIZE_ARRAYS*MAX_SIZE_ARRAYS] = {-42};
  if (N[0] < -1) 
  {
    int n = 0;
    for (int i=0; i<redshift.cluster_nbin; i++) 
    {
      for (int j=0; j<redshift.shear_nbin; j++) 
      {
        if (test_zoverlap_c(i, j)) 
        {
          N[n] = j;
          n++;
        }
      }
    }
  }
  if (nj < 0 || nj > tomo.cgl_Npowerspectra - 1)
  {
    log_fatal("invalid bin input nj = %d (max %d)", nj, tomo.cgl_Npowerspectra);
    exit(1);
  }
  return N[nj];
}

int N_cgl(int ni, int nj) 
{
  static int N[MAX_SIZE_ARRAYS][MAX_SIZE_ARRAYS] = {{-42}};
  if (N[0][0] < 0) 
  {
    int n = 0;
    for (int i=0; i<redshift.cluster_nbin; i++) 
    {
      for (int j=0; j<redshift.shear_nbin; j++) 
      {
        if (test_zoverlap_c(i, j)) 
        {
          N[i][j] = n;
          n++;
        } 
        else 
        {
          N[i][j] = -1;
        }
      }
    }
  }
  if (ni < 0 || ni > redshift.cluster_nbin - 1 || nj < 0 || nj > redshift.shear_nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
  return N[ni][nj];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Cluster Clustering (redshift overlap tests)
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int ZCCL1(const int ni)
{ // find ZCCL1 of tomography combination (zccl1, zccl2) constituting c-clustering tomo bin Nbin
  static int N[MAX_SIZE_ARRAYS*MAX_SIZE_ARRAYS] = {-42};
  if (N[0] < -1) 
  {
    int n = 0;
    for (int i=0; i<redshift.cluster_nbin; i++) 
    {
      for (int j=i; j<redshift.cluster_nbin; j++) 
      {
        N[n] = i;
        n++;
      }
    }
  }
  if (ni < 0 || ni > tomo.cc_clustering_Npowerspectra - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", ni, tomo.cc_clustering_Npowerspectra);
    exit(1);
  }
  return N[ni];
}

int ZCCL2(const int nj)
{ // find ZCCL2 of tomography combination (zcl1, zcl2) constituting c-clustering tomo bin Nbin
  static int N[MAX_SIZE_ARRAYS*MAX_SIZE_ARRAYS] = {-42};
  if (N[0] < -1) 
  {
    int n = 0;
    for (int i=0; i<redshift.cluster_nbin; i++) 
    {
      for (int j=i; j<redshift.cluster_nbin; j++) 
      {
        N[n] = j;
        n++;
      }
    }
  }
  if (nj < 0 || nj > tomo.cc_clustering_Npowerspectra - 1)
  {
    log_fatal("invalid bin input nj = %d (max %d)", nj, tomo.cc_clustering_Npowerspectra);
    exit(1);
  }
  return N[nj];
}

int N_CCL(const int ni, const int nj)
{
  static int N[MAX_SIZE_ARRAYS][MAX_SIZE_ARRAYS] = {{-42}};
  if (N[0][0] < -1) 
  {
    int n = 0;
    for (int i=0; i<redshift.cluster_nbin; i++) 
    {
      for (int j=i; j<redshift.cluster_nbin; j++) 
      {
        N[i][j] = n;
        N[j][i] = n;
        n++;
      }
    }
  }
  if (ni < 0 || ni > redshift.cluster_nbin - 1 || nj < 0 || nj > redshift.cluster_nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
  return N[ni][nj];
}

*/