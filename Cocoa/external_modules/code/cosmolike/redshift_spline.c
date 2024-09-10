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

// In lens = source LSST systematic paper, Supranta made a case where photo-z = 4 
// with one important change: the test_overlap would be given by 
// tomo.clustering_zmax[ni] <= tomo.shear_zmin[nj]. We will call this IA = 7

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// integration boundary routines
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double amin_source(int ni) 
{
  if(ni < -1 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  return 1. / (redshift.shear_zdistrpar_zmax + 1.);
}

double amax_source(int i __attribute__((unused))) 
{
  return 1. / (1. + fmax(redshift.shear_zdistrpar_zmin, 0.001));
}

double amax_source_IA(int ni) 
{
  if(ni < -1 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  return 1. / (1. + fmax(redshift.shear_zdistrpar_zmin, 0.001));
}

double amin_lens(int ni) 
{
  if(ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  
  if (ni == -1) 
  {
    return 1. / (redshift.clustering_zdistrpar_zmax + 1.);
  }

  return 1. / (1 + tomo.clustering_zmax[ni] + 
                 2. * fabs(nuisance.bias_zphot_clustering[ni]));
}

double amax_lens(int ni) 
{
  if(ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  if (gbias.b_mag[ni] != 0) 
  {
    return 1. / (1. + fmax(redshift.shear_zdistrpar_zmin, 0.001));
  }
  
  if (ni == -1) 
  {
    return 1. / (1. + fmax(redshift.clustering_zdistrpar_zmin, 0.001));
  }
  
  return 1. / (1 + fmax(tomo.clustering_zmin[ni] -
      2. * fabs(nuisance.bias_zphot_clustering[ni]), 0.001));
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// redshift overlap tests, allowed tomography combinations
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int test_kmax(double l, int ni) // return 1 if true, 0 otherwise
{ // test whether the (l, ni) bin is in the linear clustering regime
  static double chiref[10] = {-1.};
    
  if (chiref[0] < 0)
  {
    for (int i=0; i<tomo.clustering_Nbin; i++)
    {
      chiref[i] = chi(1.0/(1. + 0.5 * (tomo.clustering_zmin[i] + 
                                       tomo.clustering_zmax[i])));
    }
  }
  
  if(ni < -1 || ni > tomo.clustering_Nbin - 1)
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
// Galaxy-Galaxy Lensing bins (redshift overlap tests)
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int test_zoverlap(int ni, int nj) // test whether source bin nj is behind lens bin ni
{ 
  if(ni < 0 || ni > tomo.clustering_Nbin - 1 || nj < 0 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  return 1;

  /*if (like.use_ggl_efficiency_zoverlap == 1)
  {
    if (ggl_efficiency(ni, nj) > survey.ggl_overlap_cut) 
    {
      return 1;
    }
    else
    {
      return 0;
    }
  }
  else
  {
    const int lphotoz = redshift.clustering_photoz;
    const int sphotoz = redshift.shear_photoz;

    if (sphotoz < 4 || sphotoz == 5 || sphotoz == 6) 
    {
      if (tomo.clustering_zmax[ni] <= tomo.shear_zmin[nj])
      {
        return 1;
      }
      else
      {
        return 0;
      }
    }
    
    if (sphotoz == 7) 
    {
      if (tomo.shear_zmax[nj] >= tomo.clustering_zmin[ni])
      {
        return 1;
      }
      else
      {
        return 0;
      }
    }

    const double zmeanlens = zmean(ni);
    const double zmeansource = zmean_source(nj);

    if (sphotoz == 4 && lphotoz != 4 && tomo.clustering_zmax[ni] < zmeansource) 
    {
      return 1;
    }
    else
    {
      return 0;
    }
    
    if (sphotoz == 4 && lphotoz == 4 && zmeanlens + 0.1 < zmeansource) 
    {
      return 1;
    }
    else
    {
      return 0;
    }

    return 0;
  }*/
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
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    {
      for (int j=0; j<tomo.shear_Nbin; j++) 
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
    for (int i = 0; i < tomo.clustering_Nbin; i++) 
    {
      for (int j = 0; j < tomo.shear_Nbin; j++) 
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
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    {
      for (int j=0; j<tomo.shear_Nbin; j++) 
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
  
  if (ni < 0 || 
      ni > tomo.clustering_Nbin - 1 || 
      nj < 0 || 
      nj > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
  return N[ni][nj];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Cosmic Shear bins
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int Z1(int ni) 
{ // find z1 of tomography combination (z1, z2) constituting shear tomo bin Nbin
  static int N[MAX_SIZE_ARRAYS*MAX_SIZE_ARRAYS] = {-42};
  if (N[0] < -1) 
  {
    int n = 0;
    for (int i=0; i < tomo.shear_Nbin; i++) 
    {
      for (int j=i; j < tomo.shear_Nbin; j++) 
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
    for (int i=0; i<tomo.shear_Nbin; i++) 
    {
      for (int j=i; j<tomo.shear_Nbin; j++) 
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
    for (int i=0; i<tomo.shear_Nbin; i++) 
    {
      for (int j=i; j<tomo.shear_Nbin; j++) 
      {
        N[i][j] = n;
        N[j][i] = n;
        n++;
      }
    }
  }
  if (ni < 0 || 
      ni > tomo.shear_Nbin - 1 || 
      nj < 0 || 
      nj > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d) (max = %d)", 
      ni, nj, tomo.shear_Nbin);
    exit(1);
  }
  return N[ni][nj];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Galaxy Clustering bins
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int ZCL1(int ni) 
{ // find ZCL1 of tomography combination (zcl1, zcl2) constituting tomo bin Nbin
  static int N[MAX_SIZE_ARRAYS*MAX_SIZE_ARRAYS] = {-42};
  if (N[0] < -1) 
  {
    int n = 0;
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    {
      for (int j=i; j<tomo.clustering_Nbin; j++) 
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
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    {
      for (int j=i; j<tomo.clustering_Nbin; j++) 
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
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    {
      for (int j=i; j<tomo.clustering_Nbin; j++) 
      {
        N[i][j] = n;
        N[j][i] = n;
        n++;
      }
    }
  }
  if (ni < 0 || 
      ni > tomo.clustering_Nbin - 1 || 
      nj < 0 || 
      nj > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
  return N[ni][nj];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Shear routines for redshift distributions, including photo-zs * (optional)
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double bias_zphot_shear(double z, int ni) 
{
  if (ni < 0 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  return (1.0 + z) * nuisance.bias_zphot_shear[ni];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double sigma_zphot_shear(double z, int ni) 
{
  if (ni < 0 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  return (1.0 + z) * nuisance.sigma_zphot_shear[ni];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double zdistr_histo_n(double z, void* params)
{ // return nz(z, j) based on redshift file with structure 
  // z[i] nz[0][i] .. nz[tomo.shear_Nbin-1][i]
  static double** tab = NULL;
  static int zbins;
  static double zhisto_max; 
  static double zhisto_min; 
  static double dz_histo;
  
  if (tab == NULL) 
  {
    zbins = line_count(redshift.shear_REDSHIFT_FILE);
    
    tab = (double**) malloc(sizeof(double*)*tomo.shear_Nbin + 
                            sizeof(double)*tomo.shear_Nbin*zbins);
    if (tab == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }
    for (int i=0; i<tomo.shear_Nbin; i++)
    {
      tab[i] = ((double*)(tab + tomo.shear_Nbin) + zbins*i);
    }
    
    double* z_v = (double*) malloc(sizeof(double)*zbins);
    if (z_v == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }

    FILE* ein = fopen(redshift.shear_REDSHIFT_FILE, "r");
    if (ein == NULL)
    {
      log_fatal("file not opened");
      exit(1);
    }

    for (int i=0; i<zbins; i++) 
    {
      int count = fscanf(ein, "%le", &z_v[i]);
      if(count != 1)
      {
        log_fatal("fscanf failed to read the file");
        exit(1);
      }

      if (i > 0) 
      {
        if (z_v[i] < z_v[i - 1])
        {
          log_fatal("bad n(z) file (dz not monotonic)");
          exit(1);
        }
      }

      for (int k=0; k<tomo.shear_Nbin; k++) 
      {
        count = fscanf(ein, "%le", &tab[k][i]);
        if(count != 1)
        {
          log_fatal("fscanf failed to read the file");
          exit(1);
        }
      }
    }
    
    dz_histo = (z_v[zbins - 1] - z_v[0]) / ((double) zbins - 1.);
    zhisto_min = z_v[0];
    zhisto_max = z_v[zbins - 1] + dz_histo;

    // VM: added lines below (JULY 2003)
    redshift.shear_zdistrpar_zmin = fmax(zhisto_min, 1.e-5);
    redshift.shear_zdistrpar_zmax = zhisto_max;

    // now, set tomography bin boundaries
    for (int k=0; k<tomo.shear_Nbin; k++) 
    {
      double max = tab[k][0];
      for (int i=1; i<zbins; i++) 
      {
        if (tab[k][i] > max) 
        {
          max = tab[k][i];
        }
      }
      {
        int i = 0;
        while (tab[k][i] < 1.e-8 * max && i < zbins - 2) 
        {
          i++;
        }
        tomo.shear_zmin[k] = z_v[i];
      }
      {
        int i = zbins - 1;
        while (tab[k][i] < 1.e-8 * max && i > 0) 
        {
          i--;
        }
        tomo.shear_zmax[k] = z_v[i];
      }

      //log_debug("tomo.shear_zmin[%d] = %.3f,tomo.shear_zmax[%d] = %.3f", k,
      //  tomo.shear_zmin[k], k, tomo.shear_zmax[k]);
    }
    
    fclose(ein);
    free(z_v);
    
    if (zhisto_max < tomo.shear_zmax[tomo.shear_Nbin - 1] || 
        zhisto_min > tomo.shear_zmin[0]) 
    {
      log_fatal("zhisto_min = %e,zhisto_max = %e", zhisto_min, zhisto_max);
      
      log_fatal("tomo.shear_zmin[0] = %e, tomo.shear_zmax[N-1] = %e", 
        tomo.shear_zmin[0], tomo.shear_zmax[tomo.shear_Nbin - 1]);
      
      log_fatal("%s parameters incompatible with tomo.shear bin choice", 
        redshift.shear_REDSHIFT_FILE);
      
      exit(1);
    }
  }

  double res = 0.0;
  if ((z >= zhisto_min) && (z < zhisto_max)) 
  {
    double* ar = (double*) params; 
    const int ni = (int) ar[0];
    const int nj = (int) floor((z - zhisto_min) / dz_histo);
    
    if(ni < 0 || 
       ni > tomo.shear_Nbin - 1 || 
       nj < 0 || 
       nj > zbins - 1)
    {
      log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
      exit(1);
    } 
    res = tab[ni][nj];
  }
  return res;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double n_of_z(double z, int nz) 
{
  double sigma = 0.15;
  double x = 1.0 + 0.5 * nz - z;
  return 1./sqrt(2.*M_PI*sigma*sigma) * exp(-x * x / (2. * sigma * sigma));
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double zdistr_photoz(double zz, int nj) 
{
  static double** table = NULL;
  static double* z_v = 0;
  static int zbins = -1;
  static gsl_spline* photoz_splines[MAX_SIZE_ARRAYS+1];
  static int zbins_file = -1;
  static double* z_v_file = 0;

  if (redshift.shear_photoz == -1) 
  {
    return n_of_z(zz, nj);
  } 
  
  if (table == NULL) 
  {
    zbins_file = line_count(redshift.shear_REDSHIFT_FILE);
    
    z_v_file = (double*) malloc(sizeof(double)*zbins_file);
    if (z_v_file == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }

    FILE* ein = fopen(redshift.shear_REDSHIFT_FILE, "r");
    if (ein == NULL)
    {
      log_fatal("file not opened");
      exit(1);
    }

    for (int i=0; i<zbins_file; i++) 
    {
      int count = fscanf(ein, "%le", &z_v_file[i]);
      if(count != 1)
      {
        log_fatal("fscanf failed to read the file");
        exit(1);
      }

      if (i > 0) 
      {
        if (z_v_file[i] < z_v_file[i - 1])
        {
          log_fatal("bad n(z) file (dz not monotonic)");
          exit(1);
        }
      }

      for (int k=0; k<tomo.shear_Nbin; k++) 
      {
        double space;
        count = fscanf(ein, "%le", &space);
        if(count != 1)
        {
          log_fatal("fscanf failed to read the file");
          exit(1);
        }
      }
    }
    fclose(ein);

    const double zmin_file = z_v_file[0];
    const double zmax_file = z_v_file[zbins_file - 1];
    
    const double dz_file = (zmax_file - zmin_file) / ((double) zbins_file-1.0);
    redshift.shear_zdistrpar_zmin = fmax(zmin_file, 1.e-5);
    redshift.shear_zdistrpar_zmax = zmax_file + dz_file;
 
    zbins = zbins_file * Ntable.acc_boost_photoz_sampling;

    table = (double**) malloc(sizeof(double*)*(tomo.shear_Nbin+1) + 
                              sizeof(double)*(tomo.shear_Nbin+1)*zbins);
    if (table == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }
    for (int i=0; i<(tomo.shear_Nbin+1); i++)
    {
      table[i]  = ((double*)(table + (tomo.shear_Nbin+1)) + zbins*i);
      for (int j=0; j<zbins; j++)
      {
        table[i][j] = 0.0;
      }
    }
  
    z_v = (double*) malloc(sizeof(double)*zbins);
    if (z_v == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }

    for (int i=0; i<tomo.shear_Nbin+1; i++) 
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
        log_fatal("array allocation failed");
        exit(1);
      }
    }

    const double zhisto_min = redshift.shear_zdistrpar_zmin;
    const double zhisto_max = redshift.shear_zdistrpar_zmax;
    const double dz_histo = (zhisto_max - zhisto_min) / ((double) zbins);  
    
    for (int i=0; i<zbins; i++) 
    {
      z_v[i] = zhisto_min + (i + 0.5) * dz_histo;
    }

    double NORM[MAX_SIZE_ARRAYS];
    {
      const int i = 0;
      double ar[4]; 
      ar[0] = (double) i;
      zdistr_histo_n(0.0, (void*) ar);
    }
    #pragma omp parallel for
    for (int i=0; i<tomo.shear_Nbin; i++) 
    {
      NORM[i] = 0.0;
      
      double ar[1] = {(double) i};

      for (int k=0; k<zbins; k++) 
      {
        table[i + 1][k] = zdistr_histo_n(z_v[k], (void*) ar); 
        NORM[i] += table[i + 1][k] * dz_histo;
      }
      for (int k=0; k<zbins; k++) 
      {
        table[i + 1][k] /= NORM[i];
      }
    }
        
    double norm = 0;
    for (int i=0; i<tomo.shear_Nbin; i++) 
    { // calculate normalized overall redshift distribution
      norm += NORM[i];
    }
    #pragma omp parallel for
    for (int k=0; k<zbins; k++) 
    { // store normalization in table[0][:]
      table[0][k] = 0;
      for (int i=0; i<tomo.shear_Nbin; i++) 
      {
        table[0][k] += table[i+1][k] * NORM[i] / norm;
      }
    }

    #pragma omp parallel for
    for (int i=-1; i<tomo.shear_Nbin; i++) 
    {
      int status = gsl_spline_init(photoz_splines[i+1], z_v, table[i+1], zbins);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
    }
  }
  
  if (nj < -1 || nj > tomo.shear_Nbin - 1) 
  {
    log_fatal("nj = %d bin outside range (max = %d)", nj, tomo.shear_Nbin);
    exit(1);
  }
  
  zz = zz - nuisance.bias_zphot_shear[nj];
  //hello
  double res; 
  if (zz <= z_v[0] || zz >= z_v[zbins - 1])
  {
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
// Clustering routines for redshift distributions
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double bias_zphot_clustering(double z, int ni) 
{
  if (ni < 0 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return (1. + z) * nuisance.bias_zphot_clustering[ni];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double sigma_zphot_clustering(double z, int ni) 
{
  if (ni < 0 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return (1. + z) * nuisance.sigma_zphot_clustering[ni];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double pf_histo_n(double z, void* params) 
{ // based file with structure z[i] nz[0][i] .. nz[tomo.clustering_Nbin-1][i]
  static double** table = NULL;
  static int zbins = 0; 
  static double zhisto_min = 0;
  static double zhisto_max = 0; 
  static double dz_histo = 0;

  if (table == NULL) 
  {
    zbins = line_count(redshift.clustering_REDSHIFT_FILE);
    
    table = (double**) malloc(sizeof(double*)*tomo.clustering_Nbin + 
                              sizeof(double)*tomo.clustering_Nbin*zbins);
    if (table == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }
    for (int i=0; i<tomo.clustering_Nbin; i++)
    {
      table[i]  = ((double*)(table + tomo.clustering_Nbin) + zbins*i);
      for (int j=0; j<zbins; j++)
      {
        table[i][j] = 0.0;
      }
    }

    double* z_v = (double*) malloc(sizeof(double)*zbins);
    if (z_v == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }

    FILE* ein = fopen(redshift.clustering_REDSHIFT_FILE, "r");
    if (ein == NULL)
    {
      log_fatal("file not opened");
      exit(1);
    }

    for (int i=0; i<zbins; i++) 
    {
      int count = fscanf(ein, "%le", &z_v[i]);
      if(count != 1)
      {
        log_fatal("fscanf failed to read the file");
        exit(1);
      }
      if (i > 0) 
      {
        if (z_v[i] < z_v[i - 1])
        {
          log_fatal("bad n(z) file (redshift not monotonic)");
          exit(1);
        }
      }
      for (int k=0; k < tomo.clustering_Nbin; k++) 
      {
        count = fscanf(ein, "%le", &table[k][i]);
        if(count != 1)
        {
          log_fatal("fscanf failed to read the file");
          exit(1);
        }
      }
    }
    fclose(ein);

    dz_histo = (z_v[zbins - 1] - z_v[0]) / ((double) zbins - 1.0);
    zhisto_min = z_v[0];
    zhisto_max = z_v[zbins - 1] + dz_histo;
    
    // VM: added lines below (JULY 2003).
    redshift.clustering_zdistrpar_zmin = fmax(zhisto_min, 1.e-5);
    redshift.clustering_zdistrpar_zmax = zhisto_max;

    // now, set tomography bin boundaries
    for (int k=0; k<tomo.clustering_Nbin; k++) 
    {
      double max = table[k][0];
      for (int i = 1; i < zbins; i++) 
      {
        if (table[k][i] > max) 
        {
          max = table[k][i];
        }
      }
      {
        int i = 0;
        while (table[k][i] < 1.e-8 * max && i < zbins - 2) 
        {
          i++;
        }
        tomo.clustering_zmin[k] = z_v[i];
      }
      {
        int i = zbins - 1;
        while (table[k][i] < 1.e-8 * max && i > 0) 
        {
          i--;
        }
        tomo.clustering_zmax[k] = z_v[i]; 
      }

      //log_debug("tomo.clustering_zmin[%d] = %.3f,tomo.clustering_zmax[%d] = %.3f",
      //    k, tomo.clustering_zmin[k], k, tomo.clustering_zmax[k]);
    }

    free(z_v);
    
    if (zhisto_max < tomo.clustering_zmax[tomo.clustering_Nbin - 1] ||
        zhisto_min > tomo.clustering_zmin[0]) 
    {
      log_fatal("%e %e   %e %e", 
                zhisto_min, 
                tomo.clustering_zmin[0], 
                zhisto_max, 
                tomo.clustering_zmax[tomo.clustering_Nbin - 1]);
      log_fatal("%s parameters incompatible with tomo.clustering bin choice", 
                redshift.clustering_REDSHIFT_FILE);
      exit(1);
    }
  }

  double res = 0.0;
  if ((z >= zhisto_min) && (z < zhisto_max)) 
  {
    double* ar = (double*) params;
    const int ni = (int) ar[0];
    const int nj = (int) floor((z - zhisto_min) / dz_histo);
    
    if(ni < 0 || ni > tomo.clustering_Nbin - 1 || nj < 0 || nj > zbins - 1)
    {
      log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
      exit(1);
    } 
    res = table[ni][nj];
  }
  return res;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double pf_photoz(double zz, int nj) 
{ // works only with binned distributions; nj =-1 -> no tomo; nj>= 0 -> tomo bin nj
  static int zbins = -1;
  static gsl_spline* photoz_splines[MAX_SIZE_ARRAYS+1];
  static double* z_v = NULL;
  static double* z_v_file = NULL;
  static double** table = NULL;
  static gsl_integration_glfixed_table* gslw = NULL;
  static Ntab numtable;

  if (redshift.clustering_photoz == -1) 
  {
    return n_of_z(zz, nj);
  }

  if (table == NULL || recompute_table(numtable)) 
  {
    if (table == NULL)
    {
      const int zbins_file = line_count(redshift.clustering_REDSHIFT_FILE);
      
      z_v_file = (double*) malloc(sizeof(double)*zbins_file);
      if (z_v_file == NULL)
      {
        log_fatal("array allocation failed");
        exit(1);
      }

      FILE* ein = fopen(redshift.clustering_REDSHIFT_FILE, "r");
      if (ein == NULL)
      {
        log_fatal("file not opened");
        exit(1);
      }

      for (int i=0; i<zbins_file; i++) 
      {
        int count = fscanf(ein, "%le", &z_v_file[i]);
        if(count != 1)
        {
          log_fatal("fscanf failed to read the file");
          exit(1);
        }
        
        if (i > 0) 
        {
          if (z_v_file[i] < z_v_file[i - 1])
          {
            log_fatal("bad n(z) file (redshift not monotonic)");
            exit(1);
          }
        }

        for (int k=0; k<tomo.clustering_Nbin; k++) 
        {
          double space;
          int count = fscanf(ein, "%le", &space);
          if(count != 1)
          {
            log_fatal("fscanf failed to read the file");
            exit(1);
          }
        }
        
        const double zmin_file = z_v_file[0];
        const double zmax_file = z_v_file[zbins_file - 1];
        const double dz_file = (zmax_file - zmin_file)/((double) zbins_file-1.0);
        redshift.clustering_zdistrpar_zmin = fmax(zmin_file, 1.e-5);
        redshift.clustering_zdistrpar_zmax = zmax_file + dz_file;
      }

      fclose(ein);

      zbins = zbins_file * Ntable.acc_boost_photoz_sampling; //VM: upsampling if needed
      
      table = (double**) malloc(sizeof(double*)*(tomo.clustering_Nbin+1) + 
                                sizeof(double)*(tomo.clustering_Nbin+1)*zbins);
      if (table == NULL)
      {
        log_fatal("array allocation failed");
        exit(1);
      }
      for (int i=0; i<(tomo.clustering_Nbin+1); i++)
      {
        table[i]  = ((double*)(table + (tomo.clustering_Nbin+1)) + zbins*i);
        for (int j=0; j<zbins; j++)
        {
          table[i][j] = 0.0;
        }
      }

      z_v = (double*) malloc(sizeof(double)*zbins);
      if (z_v == NULL)
      {
        log_fatal("array allocation failed");
        exit(1);
      }

      for (int i=0; i<tomo.clustering_Nbin+1; i++) 
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

      for (int i=0; i<zbins; i++) 
      {
        const double zhisto_min = redshift.clustering_zdistrpar_zmin;
        const double zhisto_max = redshift.clustering_zdistrpar_zmax;
        const double dz_histo = (zhisto_max - zhisto_min) / ((double) zbins);
        z_v[i] = zhisto_min + (i + 0.5) * dz_histo;
      }
    }

    const size_t nsize_integration = 300 + 50 * (Ntable.high_def_integration);
    if (gslw != NULL)
    {
      gsl_integration_glfixed_table_free(gslw);    
    }
    gslw = gsl_integration_glfixed_table_alloc(nsize_integration);

    double NORM[MAX_SIZE_ARRAYS]; 
    { // initialization that allows openmp threading of the next loop
      const int i = 0;
      double ar[1];        
      ar[0] = (double) i;
      pf_histo_n(0., (void*) ar);
    } 
    #pragma omp parallel for
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    {
      double ar[1] = {(double) i};
      gsl_function F;
      F.params = (void*) ar;
      F.function = pf_histo_n;

      NORM[i] = gsl_integration_glfixed(&F, tomo.clustering_zmin[i], 
                                        tomo.clustering_zmax[i], gslw);
      if (NORM[i] == 0) 
      {
        log_fatal("pf_photoz: norm(nz = %d) = 0", i);
        exit(1);
      }

      for (int k=0; k<zbins; k++) 
      {
        table[i+1][k] = pf_histo_n(z_v[k], (void*) ar)/NORM[i];
      }
    }

    double norm = 0;
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    { // calculate normalized overall redshift distribution
      norm += NORM[i];
    }
    for (int k=0; k<zbins; k++) 
    { // store overall normalization in table[0][:]
      table[0][k] = 0;
      for (int i=0; i<tomo.clustering_Nbin; i++) 
      {
        table[0][k] += table[i + 1][k] * NORM[i] / norm;
      }
    }

    #pragma omp parallel for
    for (int i=-1; i<tomo.clustering_Nbin; i++) 
    {
      int status = gsl_spline_init(photoz_splines[i+1], z_v, table[i+1], zbins);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
    }

    update_table(&numtable);
  }
  
  if (nj > tomo.clustering_Nbin - 1 || nj < -1)  
  {
    log_fatal("pf_photoz(z, %d) outside tomo.clustering_Nbin range", nj);
    exit(1);
  }
  
  zz = zz - nuisance.bias_zphot_clustering[nj];
  
  double res;
  if (zz <= z_v[0] || zz >= z_v[zbins - 1]) 
  {
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
// routines calculating the number of source and lens galaxies per bin
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double int_nsource(double z, void* param __attribute__((unused))) 
{
  return zdistr_photoz(z, -1);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double nsource(const int ni) // ni =-1 -> no tomography; ni>= 0 -> tomography bin ni
{ // returns n_gal for shear tomography bin j, works only with binned
  static double* table = 0;

  if (table == 0) 
  {
    table = (double*) malloc(sizeof(double)*(tomo.shear_Nbin+1));
    
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = zdistr_photoz(0.0, 0.0);
      
    }
    #pragma GCC diagnostic pop

    for (int i=0; i<tomo.shear_Nbin; i++) 
    {
      table[i + 1] = tomo.n_source[i];
      
      if (tomo.n_source[i] < 0.01 || tomo.n_source[i] > 100.) 
      {
        log_fatal("nsource(%d) = %e\n outside expected range", i, tomo.n_source[i]);
        exit(0);
      }
    }
    
    table[0] = survey.n_gal;
  }
  
  if (ni < -1 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  
  return table[ni + 1];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double int_nlens(double z, void *param __attribute__((unused))) 
{
  return pf_photoz(z, -1);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double nlens(int ni) // ni =-1 -> no tomography; ni>= 0 -> tomography bin ni
{ // returns n_gal for clustering tomography bin ni
  static double* table = 0;
  
  if (table == 0) 
  {
    table = (double*) malloc(sizeof(double)*(tomo.clustering_Nbin+1));

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = pf_photoz(0., 0);
    }
    #pragma GCC diagnostic pop
         
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    {
      table[i + 1] = tomo.n_lens[i];
    }
    table[0] = survey.n_lens;
  }
  
  if (ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  
  return table[ni + 1];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double int_for_zmean(double z, void* params) 
{
  double* ar = (double*) params;
  const int ni = (int) ar[0];
  if (ni < -1 || ni > tomo.clustering_Nbin - 1)
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
  
  if (ni < -1 || ni > tomo.clustering_Nbin - 1)
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
  static gsl_integration_glfixed_table* gslw = 0;
  static Ntab numtable;

  if (table == NULL || recompute_table(numtable))
  {
    if (table == NULL) 
    {
      table = (double*) malloc(sizeof(double)*(tomo.clustering_Nbin+1));
      if (table == NULL)
      {
        log_fatal("array allocation failed");
        exit(1);
      }
    }

    const size_t nsize_integration = 300 + 50 * (Ntable.high_def_integration);
    if (gslw != NULL)
    {
      gsl_integration_glfixed_table_free(gslw);    
    }
    gslw = gsl_integration_glfixed_table_alloc(nsize_integration);

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
    {
      double init = pf_photoz(0., 0);
    }
    #pragma GCC diagnostic pop
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    {
      double ar[1] = {(double) i};
      gsl_function F;
      F.params = ar;
      
      F.function = int_for_zmean;
      const double num = gsl_integration_glfixed(&F, tomo.clustering_zmin[i], 
                                                 tomo.clustering_zmax[i], gslw);
      F.function = norm_for_zmean;
      const double den = gsl_integration_glfixed(&F, tomo.clustering_zmin[i], 
                                                 tomo.clustering_zmax[i], gslw);
      table[i] = num/den;
    }

    update_table(&numtable);
  }

  if (ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  
  return table[ni];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double int_for_zmean_source(double z, void* params) 
{
  double* ar = (double*) params;
  const int ni = (int) ar[0];
  
  if (ni < -1 || ni > tomo.shear_Nbin - 1)
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
  static gsl_integration_glfixed_table* gslw = 0;
  static double* table = NULL;
  static Ntab numtable;

  if (table == NULL || recompute_table(numtable))
  {
    if (table == NULL) 
    {
      table = (double*) malloc(sizeof(double)*tomo.shear_Nbin);
      if (table == NULL)
      {
        log_fatal("array allocation failed");
        exit(1);
      }
    }

    const size_t nsize_integration = 300 + 50 * (Ntable.high_def_integration);
    if (gslw != NULL)
    {
      gsl_integration_glfixed_table_free(gslw);    
    }
    gslw = gsl_integration_glfixed_table_alloc(nsize_integration);

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
    for (int i=0; i<tomo.shear_Nbin; i++) 
    {
      double ar[1] = {(double) i};
      gsl_function F;
      F.params = ar;
      F.function = int_for_zmean_source;
      
      table[i] = gsl_integration_glfixed(&F, tomo.shear_zmin[i], 
                                         tomo.shear_zmax[i], gslw);
    }
  }

  if (ni < 0 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  return table[ni];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double max_g_tomo(int zs) 
{
  int i = 0;
  double max = 0;
  
  for (double z=0.; z<tomo.shear_zmax[zs]; z += tomo.shear_zmax[zs] / 50.) 
  {
    if (i != 0)
    {
      if (!(fabs(z) > 0))
      {
        log_fatal("z = 0 with 1/z needed");
        exit(1);
      }
      const double a = 1. / (1. + z);
      const double g = g_tomo(a, zs) * (1. + z) * f_K(chi(a));
      if (g > max) 
      {
        max = g;
      }
    }
    else
    {
      i++;
    }
  }
  return max;
}

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// Bin-averaged lens efficiencies
// ------------------------------------------------------------------------
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
  
  if (ni < -1 || ni > tomo.shear_Nbin - 1)
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
  static gsl_integration_glfixed_table* gslw = NULL;

  const double amin = 1.0/(redshift.shear_zdistrpar_zmax + 1.0);
  const double amax = 0.999999;
  const double da = (amax - amin)/((double) Ntable.N_a - 1.0);

  if (table == NULL || recompute_table(numtable)) 
  {
    if (table == NULL)
    {
      free(table);
      table = NULL;
    }

    table = (double**) malloc(sizeof(double*)*(tomo.shear_Nbin + 1) + 
                              sizeof(double)*(tomo.shear_Nbin + 1)*Ntable.N_a);
    if (table == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }
    for (int i=0; i<(tomo.shear_Nbin + 1); i++)
    {
      table[i] = ((double*)(table + (tomo.shear_Nbin + 1)) + Ntable.N_a*i);
      for (int j=0; j<Ntable.N_a; j++)
      {
        table[i][j] = 0.0;
      }
    }

    if (gslw != NULL)
    {
      gsl_integration_glfixed_table_free(gslw);    
    }
    const size_t nsize_integration = 210 + 50 * (Ntable.high_def_integration);
    
    gslw = gsl_integration_glfixed_table_alloc(nsize_integration);
    if (gslw == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }
  }

  if (recompute_zphot_shear(N) || 
      recompute_cosmo3D(C) || 
      recompute_table(numtable)) 
  {
    { // init static variables
      double ar[2];
      ar[1] = chi(amin);
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-variable"
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
      {
        const int j = -1;
        ar[0] = (double) j; 
        double init = int_for_g_tomo(amin, (void*) ar);
      } 
      if (tomo.shear_Nbin > 0) 
      { 
        const int j = 0;      
        ar[0] = (double) j;
        double init = int_for_g_tomo(amin, (void*) ar);
      }
      #pragma GCC diagnostic pop
      #pragma GCC diagnostic pop
    }

    #pragma omp parallel for collapse(2)
    for (int j=-1; j<tomo.shear_Nbin; j++) 
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
        
        table[j+1][i] = gsl_integration_glfixed(&F, amin, a, gslw);
      }
    }

    update_cosmopara(&C);
    update_nuisance(&N);
    update_table(&numtable);
  }

  if(ni < -1 || ni > tomo.shear_Nbin - 1)
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
    res = interpol(table[ni+1], Ntable.N_a, amin, amax, da, ainput, 1.0, 1.0); 
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
  
  if (ni < -1 || ni > tomo.shear_Nbin - 1)
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
  static gsl_integration_glfixed_table* gslw = NULL;

  const double amin = 1.0/(redshift.shear_zdistrpar_zmax + 1.0);
  const double amax = 0.999999;
  const double da = (amax - amin)/((double) Ntable.N_a - 1.0);

  if (table == 0 || recompute_table(numtable)) 
  {
    if (table != NULL)
    {
      free(table);
      table = NULL;
    }

    const int NSIZE=(tomo.shear_Nbin + 1);
    table = (double**) malloc(sizeof(double*)*NSIZE + 
                              sizeof(double)*NSIZE*Ntable.N_a);
    if (table == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }
    for (int i=0; i<NSIZE; i++)
    {
      table[i] = ((double*)(table + NSIZE) + Ntable.N_a*i);
      for (int j=0; j<Ntable.N_a; j++)
      {
        table[i][j] = 0.0;
      }
    }

    if (gslw != NULL)
    {
      gsl_integration_glfixed_table_free(gslw);    
    }
    const size_t nsize_integration = 210 + 50 * (Ntable.high_def_integration);
    gslw = gsl_integration_glfixed_table_alloc(nsize_integration);
    if (gslw == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }
  }

  if (recompute_zphot_shear(N) || 
      recompute_cosmo3D(C) || 
      recompute_table(numtable)) 
  {
    { // init static variables
      double ar[2];
      ar[1] = chi(amin);
      {
        const int j = -1;
        ar[0] = (double) j; 
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-variable"
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
        {
          double init = int_for_g2_tomo(amin, (void*) ar);
        }
        #pragma GCC diagnostic pop
        #pragma GCC diagnostic pop
      } 
      if (tomo.shear_Nbin > 0) 
      { 
        const int j = 0;      
        ar[0] = (double) j;

        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-variable"
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
        {
          double init = int_for_g2_tomo(amin, (void*) ar);
        }
        #pragma GCC diagnostic pop
        #pragma GCC diagnostic pop
      }
    }

    #pragma omp parallel for collapse(2)
    for (int j=-1; j<tomo.shear_Nbin; j++) 
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

        table[j+1][i] = gsl_integration_glfixed(&F, amin, a, gslw);
      }
    }

    update_nuisance(&N);
    update_cosmopara(&C);
    update_table(&numtable);
  }

  if(ni < -1 || ni > tomo.shear_Nbin - 1)
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
    res = interpol(table[ni + 1], Ntable.N_a, amin, amax, da, a, 1.0, 1.0);
  }  
  return res;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double int_for_g_lens(double aprime, void* params) 
{
  double *ar = (double*)params;

  const int ni = (int) ar[0];
  if(ni < -1 || ni > tomo.clustering_Nbin - 1)
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
  static gsl_integration_glfixed_table* gslw = NULL;

  const double amin = 1.0/(redshift.clustering_zdistrpar_zmax + 1.0);
  const double amax = 0.999999;
  const double da = (amax - amin)/((double) Ntable.N_a - 1.0);
  const double amin_shear = 1. / (redshift.shear_zdistrpar_zmax + 1.);

  if (table == NULL || recompute_table(numtable)) 
  {
    if (table != NULL)
    {
      free(table);
      table = NULL;
    }

    const int NSIZE = (tomo.clustering_Nbin + 1);
    table = (double**) malloc(sizeof(double*)*NSIZE + 
                              sizeof(double)*NSIZE*Ntable.N_a);
    if (table == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }
    for (int i=0; i<NSIZE; i++)
    {
      table[i] = ((double*)(table + NSIZE) + Ntable.N_a*i);
      for (int j=0; j<Ntable.N_a; j++)
      {
        table[i][j] = 0.0;
      }
    }

    if (gslw != NULL)
    {
      gsl_integration_glfixed_table_free(gslw);    
    }
    const size_t nsize_integration = 250 + 50 * (Ntable.high_def_integration);
    gslw = gsl_integration_glfixed_table_alloc(nsize_integration);
    if (gslw == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }
  }

  if (recompute_zphot_clustering(N) || 
      recompute_cosmo3D(C) || 
      recompute_table(numtable)) 
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
    { // init static variables
      double ar[2];
      ar[1] = chi(amin);
      const int j = -1;
      ar[0] = (double) j; // j = -1, no tomography is being done 
      {
        double init = int_for_g_lens(amin_shear, (void*) ar);
      }
      if (tomo.clustering_Nbin > 0) 
      {
        const int j = 0;
        ar[0] = (double) j;
        {
          double init = int_for_g_lens(amin_shear, (void*) ar);
        }
        
      }
    }
    #pragma GCC diagnostic pop 
    #pragma GCC diagnostic pop

    #pragma omp parallel for collapse(2)
    for (int j=-1; j<tomo.clustering_Nbin; j++) 
    {
      for (int i=0; i<Ntable.N_a; i++) 
      {
        double ar[2];
        ar[0] = (double) j; // if j=-1, no tomography is being done

        const double a =  amin + i*da;
        ar[1] = chi(a);

        gsl_function F;
        F.params = ar;
        F.function = int_for_g_lens;

        table[j+1][i] = gsl_integration_glfixed(&F, amin_shear, a, gslw);
      }
    }

    update_table(&numtable);
    update_nuisance(&N);
    update_cosmopara(&C);
  }

  if(ni < -1 || ni > tomo.clustering_Nbin - 1)
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
    res = interpol(table[ni + 1], Ntable.N_a, amin, amax, da, a, 1.0, 1.0); 
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



























































// TODO: CLUSTER ROUTINES BELOW ARE IN ALPHA STAGE



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
    table = (double**) malloc(sizeof(double*)*tomo.cluster_Nbin);
    for (int i=0; i<tomo.cluster_Nbin; i++)
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
      for (int k=0; k<tomo.cluster_Nbin; k++)
      {
        fscanf(ein," %le", &table[k][i]);
      }
    }
    fclose(ein);

    dz = (z_v[p - 1] - z_v[0]) / ((double) p - 1.0);
    zhisto_max = z_v[p - 1] + dz;
    zhisto_min = z_v[0];
    
    for (int k=0; k<tomo.cluster_Nbin; k++)
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
    
    if (zhisto_max < tomo.cluster_zmax[tomo.cluster_Nbin - 1] || zhisto_min > tomo.cluster_zmin[0])
    {
      log_fatal("%e %e %e %e", zhisto_min,tomo.cluster_zmin[0], 
        zhisto_max,tomo.cluster_zmax[tomo.cluster_Nbin-1]);
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

    if (ni < 0 || ni > tomo.cluster_Nbin - 1 || nj < 0 || nj > njmax - 1)
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

  if (nz < -1 || nz > tomo.cluster_Nbin - 1)
  {
    log_fatal("invalid bin input nz = %d (max %d)", nz, tomo.cluster_Nbin);
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
 
    table = (double**) malloc(sizeof(double*)*(tomo.cluster_Nbin+1));
    for (int i=0; i<(tomo.cluster_Nbin+1); i++)
    {
      table[i] = (double*) malloc(sizeof(double)*zbins);
    }
    z_v = (double*) malloc(sizeof(double)*zbins);

    for (int i=0; i<tomo.clustering_Nbin+1; i++) 
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
      for (int k=0; k<tomo.cluster_Nbin; k++)
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
    for (int i=0; i<tomo.cluster_Nbin; i++)
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
    for (int i=0; i<tomo.cluster_Nbin; i++)
    { // calculate normalized overall redshift distribution (without bins), store in table[0][:]
      norm += NORM[i];
    }

    for (int k=0; k<zbins; k++)
    {
      table[0][k] = 0; 
      for (int i=0; i<tomo.cluster_Nbin; i++)
      {
        table[0][k] += table[i+1][k]*NORM[i]/norm;
      }
    }

    #pragma omp parallel for
    for (int i=-1; i<tomo.cluster_Nbin; i++) 
    {
      int status = gsl_spline_init(photoz_splines[i + 1], z_v, table[i + 1], zbins);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
    }
  }

  if (nz > tomo.cluster_Nbin - 1 || nz < -1)
  {
    log_fatal("pz_cluster(z, %d) outside tomo.cluster_Nbin range", nz);
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

  const int N_z = tomo.clustering_Nbin;
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
  
  if (ni < -1 || ni > tomo.cluster_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", ni, tomo.cluster_Nbin);
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

  const int N_z = tomo.cluster_Nbin;
  const int N_l = Cluster.N200_Nbin;
  if (nl < 0 || nl > N_l - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", nl, N_l);
    exit(1);
  } 
  
  const double amin = 1.0/(tomo.cluster_zmax[tomo.cluster_Nbin - 1] + 1.);
  const double amax = 0.999999;
  const double da = (amax - amin)/((double) Ntable.N_a - 1.0);
  
  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*(tomo.cluster_Nbin + 1));
    for (int i=0; i<tomo.cluster_Nbin+1; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*Ntable.N_a);
    }

    const size_t nsize_integration = 300 + 50 * (Ntable.high_def_integration);
    w = gsl_integration_glfixed_table_alloc(nsize_integration);
  }

  if (recompute_cosmo3D(C) /* there is no nuisance bias/sigma parameters yet*/)
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
          gsl_integration_glfixed(&F, 1.0/(redshift.shear_zdistrpar_zmax + 1.0), aa, w);
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

// ni = Cluster Nbin, nj = Galaxy Nbin
int N_CGCL(const int ni, const int nj)
{
  static int N[MAX_SIZE_ARRAYS][MAX_SIZE_ARRAYS] = {{-42}};
  if (N[0][0] < 0)
  {
    int n = 0;
    for (int i=0; i<tomo.cluster_Nbin; i++)
    {
      for (int j=0; j<tomo.clustering_Nbin; j++)
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
  if (ni < 0 || ni > tomo.cluster_Nbin - 1 || nj < 0 || nj > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input (ni (cluster nbin), nj (galaxy nbin)) = (%d, %d)", ni, nj);
    exit(1);
  }
  return N[ni][nj];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Cluster-Galaxy Lensing bins (redshift overlap tests)
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
    for (int i=0; i<tomo.cluster_Nbin; i++) 
    {
      for (int j=0; j<tomo.shear_Nbin; j++) 
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
    for (int i=0; i<tomo.cluster_Nbin; i++) 
    {
      for (int j=0; j<tomo.shear_Nbin; j++) 
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
    for (int i=0; i<tomo.cluster_Nbin; i++) 
    {
      for (int j=0; j<tomo.shear_Nbin; j++) 
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
  if (ni < 0 || ni > tomo.cluster_Nbin - 1 || nj < 0 || nj > tomo.shear_Nbin - 1)
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
    for (int i=0; i<tomo.cluster_Nbin; i++) 
    {
      for (int j=i; j<tomo.cluster_Nbin; j++) 
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
    for (int i=0; i<tomo.cluster_Nbin; i++) 
    {
      for (int j=i; j<tomo.cluster_Nbin; j++) 
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
    for (int i=0; i<tomo.cluster_Nbin; i++) 
    {
      for (int j=i; j<tomo.cluster_Nbin; j++) 
      {
        N[i][j] = n;
        N[j][i] = n;
        n++;
      }
    }
  }
  if (ni < 0 || ni > tomo.cluster_Nbin - 1 || nj < 0 || nj > tomo.cluster_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
  return N[ni][nj];
}