#include <assert.h>
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

#define Z_SPLINE_TYPE gsl_interp_cspline
static int GSL_WORKSPACE_SIZE = 1024;
static double SQRT2 = 1.41421356237;

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// integration boundary routines
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

double amin_source(int ni) 
{
  const int photoz = redshift.shear_photoz;
  
  if (ni == -1 || photoz == 1 || photoz == 2 || photoz == 4) 
  {
    return 1. / (redshift.shear_zdistrpar_zmax + 1.);
  }
  
  if(ni < -1 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  if (photoz == 0) 
  {
    return 1. / (1 + tomo.shear_zmax[ni]);
  }
  
  return 1. /(1 + fmin(tomo.shear_zmax[ni] + 5.0 * nuisance.sigma_zphot_shear[ni] +
    fabs(nuisance.bias_zphot_shear[ni]), redshift.shear_zdistrpar_zmax));
}

double amax_source(int i __attribute__((unused))) 
{
  return 1. / (1. + fmax(redshift.shear_zdistrpar_zmin, 0.001));
}

double amax_source_IA(int ni) 
{
  const int photoz = redshift.shear_photoz;
  
  if (ni == -1 || photoz == 1 || photoz == 2 || photoz == 4) 
  {
    return 1. / (1. + fmax(redshift.shear_zdistrpar_zmin, 0.001));
  }

  if(ni < -1 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  if (photoz == 0)
  {
    return 1. / (1. + fmax(tomo.shear_zmin[ni], 0.001));
  }
  
  return 1. /(1 + fmax(tomo.shear_zmin[ni] - 5. * nuisance.sigma_zphot_shear[ni] -
    fabs(nuisance.bias_zphot_shear[ni]), 0.001));
}

double amin_lens(int ni) 
{
  const int photoz = redshift.clustering_photoz;
  
  if (ni == -1 || photoz == 1 || photoz == 2) 
  {
    return 1. / (redshift.clustering_zdistrpar_zmax + 1.);
  }

  if(ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  if (photoz == 0) 
  {
    return 1. / (1 + tomo.clustering_zmax[ni]);
  }
  else if (photoz == 4) 
  {
    return 1. / (1 + tomo.clustering_zmax[ni] + 2. * fabs(nuisance.bias_zphot_clustering[ni]));
  }
  
  return 1. / (1 + fmin(tomo.clustering_zmax[ni] + 5. * nuisance.sigma_zphot_clustering[ni] +
    fabs(nuisance.bias_zphot_clustering[ni]), redshift.clustering_zdistrpar_zmax));
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
  
  const int photoz = redshift.clustering_photoz;

  if (ni == -1 || photoz == 1 || photoz == 2) 
  {
    return 1. / (1. + fmax(redshift.clustering_zdistrpar_zmin, 0.001));
  }
  
  if (photoz == 0) 
  {
    return 1. / (1. + fmax(tomo.clustering_zmin[ni], 0.001));
  }
  
  if (photoz == 4) 
  {
    return 1. / (1 + fmax(tomo.clustering_zmin[ni] -
      2. * fabs(nuisance.bias_zphot_clustering[ni]), 0.001));
  }
  
  return 1. / (1 + fmax(tomo.clustering_zmin[ni] - 5. * nuisance.sigma_zphot_clustering[ni] -
    fabs(nuisance.bias_zphot_clustering[ni]), 0.001));
}

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// redshift overlap tests, allowed tomography combinations
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

// test whether the (l, ni) bin is in the linear clustering regime
int test_kmax(double l, int ni) // return 1 if true, 0 otherwise
{
  static double chiref[10] = {-1.};
    
  if (chiref[0] < 0)
  {
    for (int i=0; i<tomo.clustering_Nbin; i++)
    {
      chiref[i] = chi(1.0/(1. + 0.5 * (tomo.clustering_zmin[i] + tomo.clustering_zmax[i])));
    }
  }
  
  if(ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  
  const double R_min = like.Rmin_bias; // set minimum scale to which
                                       // we trust our bias model, in Mpc/h
  const double kmax = constants.twopi / R_min * cosmology.coverH0;
  
  int res = 0.0;
  if ((l + 0.5) / chiref[ni] < kmax)
  {
    res = 1;
  }

  return res;
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// Galaxy-Galaxy Lensing bins (redshift overlap tests)
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

int test_zoverlap(int ni, int nj) // test whether source bin nj is behind lens bin ni
{ 
  if(ni < 0 || ni > tomo.clustering_Nbin - 1 || nj < 0 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  if (ggl_efficiency(ni, nj) > survey.ggl_overlap_cut) 
  {
    return 1;
  }
  
  const int lphotoz = redshift.clustering_photoz;
  const int sphotoz = redshift.shear_photoz;
  if (sphotoz < 4 && tomo.clustering_zmax[ni] <= tomo.shear_zmin[nj]) 
  {
    return 1;
  }
  
  const double zmeanlens = zmean(ni);
  const double zmeansource = zmean_source(nj);
  
  if (sphotoz == 4 && lphotoz != 4 && tomo.clustering_zmax[ni] < zmeansource) 
  {
    return 1;
  }
  if (sphotoz == 4 && lphotoz == 4 && zmeanlens + 0.1 < zmeansource) 
  {
    return 1;
  }
  
  return 0;
}

int ZL(int ni) 
{
  static int N[100] = {-42};
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
  if(ni < 0 || ni > tomo.ggl_Npowerspectra - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", ni, tomo.ggl_Npowerspectra);
    exit(1);
  }
  return N[ni];
}

int ZS(int nj) 
{
  static int N[100] = {-42};
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
  if(nj < 0 || nj > tomo.ggl_Npowerspectra - 1)
  {
    log_fatal("invalid bin input nj = %d (max %d)", nj, tomo.ggl_Npowerspectra);
    exit(1);
  }
  return N[nj];
}

int N_ggl(int ni, int nj) 
{ // ni = redshift bin of the lens, nj = redshift bin of the source
  static int N[10][10] = {{-42}};
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
  if(ni < 0 || ni > tomo.clustering_Nbin - 1 || nj < 0 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
  return N[ni][nj];
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// Cluster-Galaxy Lensing bins (redshift overlap tests)
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

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

int ZC(int ni) 
{
  static int N[100] = {-42};
  if (N[0] < -1) 
  {
    int n = 0;
    for (int i = 0; i < tomo.cluster_Nbin; i++) 
    {
      for (int j = 0; j < tomo.shear_Nbin; j++) 
      {
        if (test_zoverlap_c(i, j)) 
        {
          N[n] = i;
          n++;
        }
      }
    }
  }
  if(ni < 0 || ni > tomo.cgl_Npowerspectra - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", ni, tomo.cgl_Npowerspectra);
    exit(1);
  }
  return N[ni];
}

int ZSC(int nj) 
{
  static int N[100] = {-42};
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
  if(nj < 0 || nj > tomo.cgl_Npowerspectra - 1)
  {
    log_fatal("invalid bin input nj = %d (max %d)", nj, tomo.cgl_Npowerspectra);
    exit(1);
  }
  return N[nj];
}

int N_cgl(int ni, int nj) 
{
  static int N[10][10] = {{-42}};
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
  if(ni < 0 || ni > tomo.cluster_Nbin - 1 || nj < 0 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
  return N[ni][nj];
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// Cosmic Shear bins (redshift overlap tests)
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

// find z1 of tomography combination (z1, z2) constituting shear tomo bin Nbin
int Z1(int ni) 
{
  static int N[55] = {-42};
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
  if(ni < 0 || ni > tomo.shear_Npowerspectra - 1)
  {
    log_fatal("invalid bin input ni = %d (max = %d)", ni, tomo.shear_Npowerspectra);
    exit(1);
  }
  return N[ni];
}

int Z2(int nj) 
{ // find z2 of tomography combination (z1,z2) constituting shear tomo bin Nbin
  static int N[55] = {-42};
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
  if(nj < 0 || nj > tomo.shear_Npowerspectra - 1)
  {
    log_fatal("invalid bin input nj = %d (max = %d)", nj, tomo.shear_Npowerspectra);
    exit(1);
  }
  return N[nj];
}

int N_shear(int ni, int nj) 
{ // find shear tomography bin number N_shear of tomography combination (z1, z2)
  static int N[10][10] = {{-42}};
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
  if(ni < 0 || ni > tomo.shear_Nbin - 1 || nj < 0 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d) (max = %d)", ni, nj, tomo.shear_Nbin);
    exit(1);
  }
  return N[ni][nj];
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// Galaxy Clustering (redshift overlap tests)
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

int ZCL1(int ni) 
{ // find ZCL1 of tomography combination (zcl1,zcl2) constituting g-clustering tomo bin Nbin
  static int N[55] = {-42};
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
  if(ni < 0 || ni > tomo.clustering_Npowerspectra - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", ni, tomo.clustering_Npowerspectra);
    exit(1);
  }
  return N[ni];
}

int ZCL2(int nj) 
{ // find ZCL2 of tomography combination (zcl1, zcl2) constituting g-clustering tomo bin Nbin
  static int N[55] = {-42};
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
  if(nj < 0 || nj > tomo.clustering_Npowerspectra - 1)
  {
    log_fatal("invalid bin input nj = %d (max %d)", nj, tomo.clustering_Npowerspectra);
    exit(1);
  }
  return N[nj];
}

int N_CL(int ni, int nj) 
{
  static int N[10][10] = {{-42}};
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
  if(ni < 0 || ni > tomo.clustering_Nbin - 1 || nj < 0 || nj > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
  return N[ni][nj];
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// Cluster Clustering (redshift overlap tests)
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

int ZCCL1(int ni) 
{ // find ZCCL1 of tomography combination (zccl1, zccl2) constituting c-clustering tomo bin Nbin
  static int N[55] = {-42};
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
  if(ni < 0 || ni > tomo.cc_clustering_Npowerspectra - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", ni, tomo.cc_clustering_Npowerspectra);
    exit(1);
  }
  return N[ni];
}

int ZCCL2(int nj) 
{ // find ZCCL2 of tomography combination (zcl1, zcl2) constituting c-clustering tomo bin Nbin
  static int N[55] = {-42};
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
  if(nj < 0 || nj > tomo.cc_clustering_Npowerspectra - 1)
  {
    log_fatal("invalid bin input nj = %d (max %d)", nj, tomo.cc_clustering_Npowerspectra);
    exit(1);
  }
  return N[nj];
}

int N_CCL(int ni, int nj) 
{
  static int N[10][10] = {{-42}};
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
  if(ni < 0 || ni > tomo.cluster_Nbin - 1 || nj < 0 || nj > tomo.cluster_Nbin - 1)
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
  if(ni < 0 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  return (1.0 + z) * nuisance.bias_zphot_shear[ni];
}

double sigma_zphot_shear(double z, int ni) 
{
  if(ni < 0 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  return (1.0 + z) * nuisance.sigma_zphot_shear[ni];
}

double zdistr_histo_1(double z, void* params __attribute__((unused))) 
{ // return nz(z) based on redshift file with one redshift distribution
  static double *tab = 0;
  static double zhisto_max, zhisto_min, dz;
  static int nimax;

  if (tab == 0) 
  {
    const int zbins = line_count(redshift.shear_REDSHIFT_FILE);
    nimax = zbins;
    tab = create_double_vector(0, zbins - 1);
    double* z_v = create_double_vector(0, zbins - 1);
    
    FILE* ein = fopen(redshift.shear_REDSHIFT_FILE, "r");
    if (ein == NULL)
    {
      log_fatal("file not opened");
      exit(1);
    }
    
    int p = 0;
    for (int i=0; i<zbins; i++) 
    {
      double space1, space2;
      fscanf(ein, "%le %le %le %le\n", &z_v[i], &space1, &space2, &tab[i]);
      p++;
      if (i > 0 && z_v[i] < z_v[i - 1]) 
      {
        break;
      }
    }
    fclose(ein);
    
    dz = (z_v[p - 1] - z_v[0]) / (p - 1.);
    zhisto_max = z_v[p - 1] + dz;
    zhisto_min = z_v[0];
    
    redshift.shear_zdistrpar_zmin = zhisto_min;
    redshift.shear_zdistrpar_zmax = zhisto_max;

    free_double_vector(z_v, 0, zbins - 1);
  }

  double res = 0.0;
  if ((z >= zhisto_min) && (z < zhisto_max)) 
  {
    const int ni = (int) floor((z - zhisto_min) / dz);
    
    if(ni < 0 || ni > nimax - 1)
    {
      log_fatal("invalid bin input ni = %d", ni);
      exit(1);
    } 
    
    res =  tab[ni];
  }

  return res;
}

double zdistr_histo_n(double z, void* params)
{ // return nz(z,j) based on redshift file with structure z[i] nz[0][i] .. nz[tomo.shear_Nbin-1][i]
  static double** tab;
  static double zhisto_max, zhisto_min, dz;
  static int njmax;
  
  if (tab == 0) 
  {
    int zbins = line_count(redshift.shear_REDSHIFT_FILE);
    njmax = zbins;
    tab = create_double_matrix(0, tomo.shear_Nbin - 1, 0, zbins - 1);
    double* z_v = create_double_vector(0, zbins - 1);
    
    FILE* ein = fopen(redshift.shear_REDSHIFT_FILE, "r");
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
      if (i > 0 && z_v[i] < z_v[i - 1]) 
      {
        break;
      }
      
      for (int k=0; k<tomo.shear_Nbin; k++) 
      {
        fscanf(ein, "%le", &tab[k][i]);
      }
    }
    fclose(ein);
    
    dz = (z_v[p - 1] - z_v[0]) / (p - 1.);
    zhisto_max = z_v[p - 1] + dz;
    zhisto_min = z_v[0];
    zbins = p;
    
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

      log_info("tomo.shear_zmin[%d] = %.3f,tomo.shear_zmax[%d] = %.3f", k,
        tomo.shear_zmin[k], k, tomo.shear_zmax[k]);
    }
    
    free_double_vector(z_v, 0, zbins - 1);
    
    if (zhisto_max < tomo.shear_zmax[tomo.shear_Nbin - 1] || zhisto_min > tomo.shear_zmin[0]) 
    {
      log_fatal("zhisto_min = %e,zhisto_max = %e", zhisto_min, zhisto_max);
      log_fatal("tomo.shear_zmin[0] = %e, tomo.shear_zmax[N-1] = %e", tomo.shear_zmin[0], 
        tomo.shear_zmax[tomo.shear_Nbin - 1]);
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
    const int nj = (int) floor((z - zhisto_min) / dz);
    
    if(ni < 0 || ni > tomo.shear_Nbin - 1 || nj < 0 || nj > njmax - 1)
    {
      log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
      exit(1);
    } 
    
    res = tab[ni][nj];
  }

  return res;
}

double n_of_z(double z, int nz) 
{
  double sigma = 0.15;
  double x = 1.0 + 0.5 * nz - z;
  return 1. / sqrt(2. * M_PI * sigma * sigma) * exp(-x * x / (2. * sigma * sigma));
}

// returns n(ztrue | j), works only with binned; distributions; j>= 0 -> tomography bin j
double zdistr_photoz(double zz, int nj) 
{
  static double** table = 0;
  static double* z_v = 0;
  static nuisancepara N;
  static int zbins = -1;
  static gsl_spline* photoz_splines[MAX_SIZE_ARRAYS];

  if (redshift.shear_photoz == -1) 
  {
    return n_of_z(zz, nj);
  }
  
  if ((redshift.shear_photoz != 4 && recompute_zphot_shear(N)) || table == 0) 
  {
    update_nuisance(&N);
    
    if (table == 0) 
    {
      const int zbins1 = line_count(redshift.shear_REDSHIFT_FILE);
      
      if (redshift.shear_photoz != 4) 
      {
        zbins = zbins1 * 20;
      } 
      else 
      {
        zbins = zbins1;
      }
      
      table = create_double_matrix(0, tomo.shear_Nbin, 0, zbins - 1);
      z_v = create_double_vector(0, zbins - 1);
      
      for (int i=0; i<tomo.shear_Nbin+1; i++) 
      {
        photoz_splines[i] = gsl_spline_alloc(Z_SPLINE_TYPE, zbins);
        if (photoz_splines[i] == NULL)
        {
          log_fatal("fail allocation");
          exit(1);
        }
      }

      if (redshift.shear_photoz == 4) 
      { // if multihisto, force zmin, zmax, tomo bins to match supplied file
        
        FILE* ein = fopen(redshift.shear_REDSHIFT_FILE, "r");
        if (ein == NULL)
        {
          log_fatal("file not opened");
          exit(1);
        }

        int p=0;
        for (int i=0; i<zbins1; i++) 
        {
          p++;
          
          fscanf(ein, "%le", &z_v[i]);
          
          if (i > 0 && z_v[i] < z_v[i - 1]) 
          {
            break;
          }
          
          for (int k=0; k<tomo.shear_Nbin; k++) 
          {
            double space;
            fscanf(ein, "%le", &space);
          }
        }
        fclose(ein);
        
        redshift.shear_zdistrpar_zmin = fmax(z_v[0], 1.e-5);
        redshift.shear_zdistrpar_zmax = z_v[p - 1] + (z_v[p - 1] - z_v[0]) / (zbins1 - 1.);
      }
    }

    const double zhisto_min = redshift.shear_zdistrpar_zmin;
    const double zhisto_max = redshift.shear_zdistrpar_zmax;
    const double da = (zhisto_max - zhisto_min) / ((double) zbins);
    
    for (int i = 0; i < zbins; i++) 
    {
      z_v[i] = zhisto_min + (i + 0.5) * da;
    }
    
// the outlier fraction (outfrac) should be specified externally. This is a temporary hack.

    const double outfrac = 0.05;  

// eta = 0.68; // this is the approximate mixing parameter for the pseudo
// voigt for equal FWHM L and G distributions. Note than when using this
// formulation, the non-equivalent profile re-scalings must be accounted for below.

    const double eta = 0.736; 

// this is the approximate mixing parameter for the pseudo voigt
// for equally rescaled L and G distributions. Note that this leads to a slightly broader
// profile than when using equal FWHM L and G distributions (at a given gaussian sigma). Currently,
// the input sigma parameter has not been rescaled. It still gives the width for the Gaussian
// component of the psuedo Voigt profile. In general, for a given FWHM, the corresponding sigma
// is smaller for the Voigt profile. For equal FWHM in the component profiles, eta = 0.68 and the
// effective FWHM is 1.64 that of the Guassian case. For equivalent rescalings, [which is the 
// currently implemented case] eta=0.736 and the effective FWHM is 1.79 that of the Gaussian case.

    switch (redshift.shear_photoz) 
    { // initialization that allows openmp threading of the next loop
      case 0:
      {
        const int i = 0;
        double ar[4];
        ar[0] = tomo.shear_zmin[i];
        ar[1] = tomo.shear_zmax[i];
        zdistr_histo_1(0.0, (void*) ar);
        break;
      }
      case 1:
      {
        const int i = 0;
        double ar[4];
        ar[0] = tomo.shear_zmin[i];
        ar[1] = tomo.shear_zmax[i];
        zdistr_histo_1(0.0, (void*) ar);
        break;
      }
      case 2:
      {
        const int i = 0;
        double ar[4];
        ar[0] = tomo.shear_zmin[i];
        ar[1] = tomo.shear_zmax[i];
        zdistr_histo_1(0.0, (void*) ar);
        break;
      }
      case 3:
      {
        const int i = 0;
        double ar[4];
        ar[0] = tomo.shear_zmin[i];
        ar[1] = tomo.shear_zmax[i];
        zdistr_histo_1(0.0, (void*) ar);
        break;
      }
      case 4:
      {
        const int i = 0;
        double ar[4]; 
        ar[0] = (double) i;
        zdistr_histo_n(0.0, (void*) ar);
        break;
      }
    }

    double NORM[MAX_SIZE_ARRAYS];

    #pragma omp parallel for
    for (int i=0; i<tomo.shear_Nbin; i++) 
    {
      double ar[4];
      ar[0] = tomo.shear_zmin[i];
      ar[1] = tomo.shear_zmax[i];
      
      double norm = 0.0;

      switch (redshift.shear_photoz) 
      {
        case 0: // no photo-zs, split 'true' n(z) histogram in bins
        {  
          norm = int_gsl_integrate_medium_precision(zdistr_histo_1, (void*) ar, ar[0], ar[1], 
            NULL, GSL_WORKSPACE_SIZE);
          
          for (int k=0; k < zbins; k++) 
          {
            table[i + 1][k] = 0;
            if (z_v[k] >= ar[0] && z_v[k] <= ar[1]) 
            {
              table[i + 1][k] = zdistr_histo_1(z_v[k], (void*) ar) / norm;
            }
          }
          break;
        }
        case 1: // Pseudo Voigt (Lorentzian + Gaussian) photo-zs
        {
          norm = 0.;
          
          if (sigma_zphot_shear(0., i) == 0.) 
          {
            log_fatal("Source galaxy photo-z model underdetermined!\nredshift.shear_photoz = 1,"
              "but  nuisance.sigma_zphot_shear[%d] not set", i);
            exit(1);
          }
          
          for (int k=0; k < zbins; k++) 
          {
            const double x1 = (ar[0] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
                 (SQRT2 * sigma_zphot_shear(z_v[k], i));
            const double x2 = (ar[1] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
                 (SQRT2 * sigma_zphot_shear(z_v[k], i));
            
            table[i + 1][k] =
                zdistr_histo_1(z_v[k], (void*) ar) *
                (eta * 0.3183 * (atan(x2) - atan(x1)) +
                 ((1.0 - eta) * 0.5 * (gsl_sf_erf(x2) - gsl_sf_erf(x1))));
            // this creates a pseudo Voigt profile by adding a Gaussian and
            // Lorentian (convolved with a tophat) with the correct weights. See,
            // eg, the Wikipedia article on Voigt profiles for an explanation of
            // where the numbers come from.
            
            norm += table[i + 1][k] * da;
          }
          for (int k=0; k < zbins; k++) 
          {
            table[i + 1][k] /= norm;
          }
          break;
        }
        case 2: // Pseudo Voigt (Lorentzian + Gaussian) + outlier (currently flat) photo-zs
        {
          norm = 0.;

          if (sigma_zphot_shear(0., i) == 0.) 
          {
            log_fatal("Source galaxy photo-z model underdetermined!\nredshift.shear_photoz = 2,"
              " but nuisance.sigma_zphot_shear[%d] not set",i);
            exit(1);
          }
          
          for (int k=0; k < zbins; k++) 
          {
            const double x1 = (ar[0] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
                 (SQRT2 * sigma_zphot_shear(z_v[k], i));

            const double x2 = (ar[1] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
                 (SQRT2 * sigma_zphot_shear(z_v[k], i));

            table[i + 1][k] = zdistr_histo_1(z_v[k], (void*) ar) *
              (eta * 0.3183 * (atan(x2) - atan(x1)) +
              ((1.0 - eta) * 0.5 * (gsl_sf_erf(x2) - gsl_sf_erf(x1))));

            norm += table[i + 1][k] * da;
          }
          
          for (int k=0; k < zbins; k++) 
          {
            table[i + 1][k] = outfrac / da + (1. - outfrac) * table[i + 1][k] / norm;
          }
          break;
        }
        case 3: // Gaussian photo-zs
        {
          norm = 0.;
          
          if (sigma_zphot_shear(0., i) == 0.) 
          {
            log_fatal("Source galaxy photo-z model underdetermined!\nredshift.shear_photoz = 3,"
              " but nuisance.sigma_zphot_shear[%d] not set", i);
            exit(1);
          }
          
          for (int k=0; k < zbins; k++) 
          {
            const double x1 = (ar[0] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
                 (SQRT2 * sigma_zphot_shear(z_v[k], i));
            const double x2 = (ar[1] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
                 (SQRT2 * sigma_zphot_shear(z_v[k], i));
            
            table[i + 1][k] = 0.5 * zdistr_histo_1(z_v[k], (void*) ar) *
                              (gsl_sf_erf(x2) - gsl_sf_erf(x1));
            
            norm += table[i + 1][k] * da;
          }
          
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] /= norm;
          }
          break;
        }
        case 4: // histogram file contains n(z) estimates for each bin
        {  
          ar[0] = (double) i;        
          norm = 0.0;
          
          for (int k=0; k < zbins; k++) 
          {
            table[i + 1][k] = zdistr_histo_n(z_v[k], (void*) ar); // norm;
            norm += table[i + 1][k] * da;
          }
          
          for (int k=0; k < zbins; k++) 
          {
            table[i + 1][k] /= norm;
          }
          break;
        }
        default:
        {
          log_fatal("redshift.shear_photoz = %d not supported in this cosmolike "
                 "version",redshift.shear_photoz);
          exit(1);          
        }
      }

      NORM[i] = norm;
    }
    
    // calculate normalized overall redshift distribution (without bins), store in table[0][:]
    double norm = 0;
    for (int i = 0; i < tomo.shear_Nbin; i++) 
    {
      norm += NORM[i];
    }

    for (int k=0; k<zbins; k++) 
    {
      table[0][k] = 0;
      for (int i = 0; i < tomo.shear_Nbin; i++) 
      {
        table[0][k] += table[i + 1][k] * NORM[i] / norm;
      }
    }

    #pragma omp parallel for
    for (int i=-1; i<tomo.shear_Nbin; i++) 
    {
      int status = gsl_spline_init(photoz_splines[i + 1], z_v, table[i + 1], zbins);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
    }
  }
  
  if (nj < -1 || nj > tomo.shear_Nbin - 1) 
  {
    log_fatal("zdistr_photoz(z, %d) outside bin range", nj);
    exit(1);
  }
  
  if (redshift.shear_photoz == 4) 
  {
    zz = zz - nuisance.bias_zphot_shear[nj];
  }
  
  double res; 
  if (zz <= z_v[0] || zz >= z_v[zbins - 1])
  {
    res = 0.0;
  }
  else
  {
    res = gsl_spline_eval(photoz_splines[nj + 1], zz, NULL);
  }

  return res;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Clustering routines for redshift distributions, including photo-zs * (optional)
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

double bias_zphot_clustering(double z, int ni) 
{
  if(ni < 0 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return (1. + z) * nuisance.bias_zphot_clustering[ni];
}

double sigma_zphot_clustering(double z, int ni) 
{
  if(ni < 0 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return (1. + z) * nuisance.sigma_zphot_clustering[ni];
}

double pf_histo(double z, void* params __attribute__((unused))) 
{ // return pf(z) based on redshift file with one redshift distribution
  static double *tab = 0;
  static double zhisto_max, zhisto_min, dz;

  if (tab == 0) 
  {
    const int zbins = line_count(redshift.clustering_REDSHIFT_FILE);
    tab = create_double_vector(0, zbins - 1);
    double* z_v = create_double_vector(0, zbins - 1);
    
    FILE* ein = fopen(redshift.clustering_REDSHIFT_FILE, "r");
    if (ein == NULL)
    {
      log_fatal("file not opened");
      exit(1);
    }

    int p = 0;
    for (int i=0; i<zbins; i++) 
    {
      double space1, space2;
      fscanf(ein, "%le %le %le %le\n", &z_v[i], &space1, &space2, &tab[i]);
      
      p++;
      if (i > 0 && z_v[i] < z_v[i - 1]) 
      {
        break;
      }
    }
    fclose(ein);

    dz = (z_v[p - 1] - z_v[0]) / (p - 1.);
    zhisto_max = z_v[p - 1] + dz;
    zhisto_min = z_v[0];

    free_double_vector(z_v, 0, zbins - 1);
    
    if (zhisto_max < tomo.clustering_zmax[tomo.clustering_Nbin - 1] ||
        zhisto_min > tomo.clustering_zmin[0]) 
    {
      log_fatal("Error in redshift_spline.c:pf_histo.c: %s parameters "
             "incompatible with tomo.clustering bin choice", redshift.clustering_REDSHIFT_FILE);
      exit(1);
    }
  }

  double res = 0.0;
  if ((z >= zhisto_min) && (z < zhisto_max)) 
  {
    const int ni = (int) floor((z - zhisto_min) / dz);
    
    if(ni < 0 || ni > tomo.clustering_Nbin - 1)
    {
      log_fatal("invalid bin input ni = %d", ni);
      exit(1);
    } 
    
    res = tab[ni];
  }

  return res;
}

double pf_histo_n(double z, void* params) 
{ // return pf(z, nj) based on z file with structure z[i] nz[0][i] .. nz[tomo.clustering_Nbin-1][i]
  static double** tab;
  static double zhisto_max, zhisto_min, dz;
  static int njmax;

  if (tab == 0) 
  {
    const int zbins = line_count(redshift.clustering_REDSHIFT_FILE);
    njmax = zbins;
    tab = create_double_matrix(0, tomo.clustering_Nbin - 1, 0, zbins - 1);
    double* z_v = create_double_vector(0, zbins - 1);
    
    FILE* ein = fopen(redshift.clustering_REDSHIFT_FILE, "r");
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
      if (i > 0 && z_v[i] < z_v[i - 1]) 
      {
        break;
      }

      for (int k=0; k < tomo.clustering_Nbin; k++) 
      {
        fscanf(ein, " %le", &tab[k][i]);
      }
    }
    fclose(ein);
    
    dz = (z_v[p - 1] - z_v[0]) / (p - 1.);
    zhisto_max = z_v[p - 1] + dz;
    zhisto_min = z_v[0];
    
    // now, set tomography bin boundaries
    for (int k=0; k < tomo.clustering_Nbin; k++) 
    {
      double max = tab[k][0];
      for (int i = 1; i < zbins; i++) 
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
        tomo.clustering_zmin[k] = z_v[i];
      }
      {
        int i = zbins - 1;
        while (tab[k][i] < 1.e-8 * max && i > 0) 
        {
          i--;
        }
        tomo.clustering_zmax[k] = z_v[i]; 
      }

      log_info("tomo.clustering_zmin[%d] = %.3f,tomo.clustering_zmax[%d] = %.3f",
          k, tomo.clustering_zmin[k], k, tomo.clustering_zmax[k]);
    }

    free_double_vector(z_v, 0, zbins - 1);
    
    if (zhisto_max < tomo.clustering_zmax[tomo.clustering_Nbin - 1] ||
        zhisto_min > tomo.clustering_zmin[0]) 
    {
      log_fatal("%e %e   %e %e", zhisto_min, tomo.clustering_zmin[0], zhisto_max, 
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
    const int nj = (int) floor((z - zhisto_min) / dz);
    
    if(ni < 0 || ni > tomo.clustering_Nbin - 1 || nj < 0 || nj > njmax - 1)
    {
      log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
      exit(1);
    } 
    res = tab[ni][nj];
  }

  return res;
}

double pf_photoz(double zz, int nj) 
{ // works only with binned distributions; nj =-1 -> no tomo; nj>= 0 -> tomo bin nj
  static double** table = 0;
  static double* z_v = 0;
  static nuisancepara N;
  static int zbins = -1;
  static gsl_spline* photoz_splines[MAX_SIZE_ARRAYS];

  if (redshift.clustering_photoz == -1) 
  {
    return n_of_z(zz, nj);
  }

  if ((redshift.clustering_photoz != 4 && recompute_zphot_clustering(N)) || table == 0) 
  {
    update_nuisance(&N);

    if (table == 0) 
    {
      zbins = line_count(redshift.clustering_REDSHIFT_FILE);
      
      if (redshift.clustering_photoz != 4 && redshift.clustering_photoz != 0) 
      {
        pf_histo(0.5, NULL);
        zbins *= 20;
      } // upsample if convolving with analytic photo-z model
      
      table = create_double_matrix(0, tomo.clustering_Nbin, 0, zbins - 1);
      z_v = create_double_vector(0, zbins - 1);
      
      for (int i=0; i<tomo.clustering_Nbin+1; i++) 
      {
        photoz_splines[i] = gsl_spline_alloc(Z_SPLINE_TYPE, zbins);
        if (photoz_splines[i] == NULL)
        {
          log_fatal("fail allocation");
          exit(1);
        }
      }
      
      if (redshift.clustering_photoz == 4) 
      { // if multihisto, force zmin, zmax, tomo bins to match supplied file
        
        FILE* ein = fopen(redshift.clustering_REDSHIFT_FILE, "r");
        if (ein == NULL)
        {
          log_fatal("file not opened");
          exit(1);
        }

        int p = 0;
        for (int i=0; i<zbins; i++) 
        {
           p++;
          fscanf(ein, "%le", &z_v[i]);
          
          if (i > 0 && z_v[i] < z_v[i - 1]) 
          {
            break;
          }

          for (int k=0; k < tomo.clustering_Nbin; k++) 
          {
            double space;
            fscanf(ein, "%le", &space);
          }
        }
        fclose(ein);
        redshift.clustering_zdistrpar_zmin = fmax(z_v[0], 1.e-5);
        redshift.clustering_zdistrpar_zmax = z_v[p - 1] + (z_v[p - 1] - z_v[0]) / (zbins - 1.);
      }
    }

    const double zhisto_min = redshift.clustering_zdistrpar_zmin;
    const double zhisto_max = redshift.clustering_zdistrpar_zmax;
    const double da = (zhisto_max - zhisto_min) / ((double) zbins);
    
    for (int i=0; i<zbins; i++) 
    {
      z_v[i] = zhisto_min + (i + 0.5) * da;
    }
    
// the outlier fraction (outfrac) should be specified externally. This is a temporary hack.

    const double outfrac = 0.05;  

// eta = 0.68; // this is the approximate mixing parameter for the pseudo
// voigt for equal FWHM L and G distributions. Note than when using this
// formulation, the non-equivalent profile re-scalings must be accounted for below.

    const double eta = 0.736; 

// this is the approximate mixing parameter for the pseudo voigt
// for equally rescaled L and G distributions. Note that this leads to a slightly broader
// profile than when using equal FWHM L and G distributions (at a given gaussian sigma). Currently,
// the input sigma parameter has not been rescaled. It still gives the width for the Gaussian
// component of the psuedo Voigt profile. In general, for a given FWHM, the corresponding sigma
// is smaller for the Voigt profile. For equal FWHM in the component profiles, eta = 0.68 and the
// effective FWHM is 1.64 that of the Guassian case. For equivalent rescalings, [which is the 
// currently implemented case] eta=0.736 and the effective FWHM is 1.79 that of the Gaussian case.

    switch (redshift.clustering_photoz) 
    { // initialization that allows openmp threading of the next loop
      case 0:
      {
        const int i = 0;
        double ar[2];
        ar[0] = tomo.clustering_zmin[i];
        ar[1] = tomo.clustering_zmax[i];
        pf_histo(0.0, (void*) ar);
        
        break;
      }
      case 1:
      {
        const int i = 0;
        double ar[2];
        ar[0] = tomo.clustering_zmin[i];
        ar[1] = tomo.clustering_zmax[i];
        pf_histo(0.0, (void*) ar);
        
        break;
      }
      case 2:
      {
        const int i = 0;
        double ar[2];
        ar[0] = tomo.clustering_zmin[i];
        ar[1] = tomo.clustering_zmax[i];
        pf_histo(0.0, (void*) ar);
        
        break;
      }
      case 3:
      {
        const int i = 0;
        double ar[2];
        ar[0] = tomo.clustering_zmin[i];
        ar[1] = tomo.clustering_zmax[i];
        pf_histo(0.0, (void*) ar);
        
        break;
      }
      case 4:
      {
        const int i = 0;
        double ar[2];        
        ar[0] = (double) i;
        pf_histo_n(0., (void*) ar);
        
        break;
      }
    }

    double NORM[MAX_SIZE_ARRAYS];
    
    #pragma omp parallel for
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    {
      double ar[4];
      ar[0] = tomo.clustering_zmin[i];
      ar[1] = tomo.clustering_zmax[i];

      double norm = 0.0;

      switch (redshift.clustering_photoz) 
      {
        case 0: // no photo-zs, split 'true' n(z) histogram in bins
        {
          norm = int_gsl_integrate_medium_precision(pf_histo, (void*) ar, ar[0], 
            ar[1], NULL, GSL_WORKSPACE_SIZE);
          for (int k=0; k<zbins; k++) 
          {
            double zi = z_v[k];
            table[i + 1][k] = 0.;
            if (zi >= ar[0] && zi <= ar[1]) 
            {
              table[i + 1][k] = pf_histo(zi, (void*) ar) / norm;
            }
          }
          
          break;
        }
        case 1: // Pseudo Voigt (Lorentzian + Gaussian) photo-zs
        {
          norm = 0.;
          if (sigma_zphot_clustering(0., i) == 0.) 
          {
            log_fatal("Lens galaxy photo-z model underdetermined!\nredshift.clustering_photoz = 1, "
              "but nuisance.sigma_zphot_clustering[%d] not set", i);
            exit(1);
          }
          
          for (int k=0; k<zbins; k++) 
          {
            const double zi = z_v[k];
            
            const double x1 = 
              (ar[0] - zi + bias_zphot_clustering(zi, i))/(SQRT2*sigma_zphot_clustering(zi, i));
            const double x2 = 
              (ar[1] - zi + bias_zphot_clustering(zi, i))/(SQRT2*sigma_zphot_clustering(zi, i));
            
            table[i + 1][k] = pf_histo(zi, (void*) ar)*(eta * 0.3183 * (atan(x2) - atan(x1)) +
              ((1.0 - eta) * 0.5 * (gsl_sf_erf(x2) - gsl_sf_erf(x1))));
            
            // this creates a pseudo Voigt profile by adding a Gaussian and Lorentian (convolved  
            // with a tophat) with the correct weights. See, eg, the Wikipedia article on Voigt
            // profiles for an explanation of where the numbers come from.
            norm += table[i + 1][k] * da;
          }
          
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] /= norm;
          }
          
          break;
        }
        case 2: // Pseudo Voigt (Lorentzian + Gaussian) + outlier (currently flat) photo-zs
        {
          norm = 0.;

          if (sigma_zphot_clustering(0., i) == 0.) 
          {
            log_fatal("Lens galaxy photo-z model underdetermined!\n redshift.clustering_photoz = 2,"
              "but nuisance.sigma_zphot_clustering[%d] not set", i);
            exit(1);
          }
          
          for (int k=0; k<zbins; k++) 
          {
            const double zi = z_v[k];
            
            double x1 = (ar[0]-zi+bias_zphot_clustering(zi,i))/(SQRT2*sigma_zphot_clustering(zi,i));
            double x2 = (ar[1]-zi+bias_zphot_clustering(zi,i))/(SQRT2*sigma_zphot_clustering(zi,i));
            
            table[i + 1][k] = pf_histo(zi, (void*) ar)*(eta * 0.3183 * (atan(x2) - atan(x1)) +
              ((1.0 - eta) * 0.5 * (gsl_sf_erf(x2) - gsl_sf_erf(x1))));
            
            norm += table[i + 1][k] * da;
          }
          
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] = outfrac / da + (1. - outfrac) * table[i + 1][k] / norm;
          }
          
          break;
        }
        case 3: // Gaussian photo-zs
        {
          if (sigma_zphot_clustering(0., i) == 0.) 
          {
            log_fatal("Lens galaxy photo-z model underdetermined!\nredshift.clustering_photoz = 3,"
              "but nuisance.sigma_zphot_clustering[%d] not set", i);
            exit(1);
          }
          
          norm = 0.;
          for (int k=0; k<zbins; k++) 
          {
            const double zi = z_v[k];
            double x1 = (ar[0]-zi+bias_zphot_clustering(zi,i))/(SQRT2*sigma_zphot_clustering(zi,i));
            double x2 = (ar[1]-zi+bias_zphot_clustering(zi,i))/(SQRT2*sigma_zphot_clustering(zi,i));
            table[i + 1][k] = 0.5 * pf_histo(zi, (void*) ar)*(gsl_sf_erf(x2) - gsl_sf_erf(x1));
            norm += table[i + 1][k] * da;
          }
          
          for (int k=0; k < zbins; k++) 
          {
            table[i + 1][k] /= norm;
          }
          
          break;
        }
        case 4: // histogram file contains n(z) estimates for each bin
        { 
          ar[0] = (double) i;
          
          norm = int_gsl_integrate_medium_precision(pf_histo_n, (void*) ar, tomo.clustering_zmin[i],
            tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE);
          
          if (norm == 0) 
          {
            log_fatal("pf_photoz:norm(nz=%d)=0", i);
            exit(1);
          }
          
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] = pf_histo_n(z_v[k], (void*) ar)/norm;
          }
          
          break;
        }
        default:
        {
          log_fatal("%d not supported in this cosmolike version", redshift.clustering_photoz);
          exit(1);
        }
      }
      NORM[i] = norm;
    }
    
    // calculate normalized overall redshift distribution (without bins), store in table[0][:]
    double norm = 0;
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    {
      norm += NORM[i];
    }
    
    for (int k=0; k<zbins; k++) 
    {
      table[0][k] = 0;
      for (int i=0; i<tomo.clustering_Nbin; i++) 
      {
        table[0][k] += table[i + 1][k] * NORM[i] / norm;
      }
    }
    
    #pragma omp parallel for
    for (int i=-1; i<tomo.clustering_Nbin; i++) 
    {
      int status = gsl_spline_init(photoz_splines[i + 1], z_v, table[i + 1], zbins);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
    }
  }
  
  if (nj > tomo.clustering_Nbin - 1 || nj < -1)  
  {
    log_fatal("pf_photoz(z, %d) outside tomo.clustering_Nbin range", nj);
    exit(1);
  }
  
  if (redshift.clustering_photoz == 4) 
  {
    zz = zz - nuisance.bias_zphot_clustering[nj];
  }
  
  double res;
  if (zz <= z_v[0] || zz >= z_v[zbins - 1]) 
  {
    res = 0.0;
  }
  else
  {
    double result = 0.0;
    int status = gsl_spline_eval_e(photoz_splines[nj + 1], zz, NULL, &result);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    res = result;
  }

  return res;
}

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// routines calculating the number of source and lens galaxies per bin
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

double int_nsource(double z, void* param __attribute__((unused))) 
{
  return zdistr_photoz(z, -1);
}

// returns n_gal for shear tomography bin j, works only with binned
double nsource(int ni) // ni =-1 -> no tomography; ni>= 0 -> tomography bin ni
{
  static double** table = 0;

  if (table == 0 || table[0][0] != survey.n_gal) 
  {
    zdistr_photoz(0.0, 0.0);
    
    table = create_double_matrix(0, tomo.shear_Nbin, 0, 1);
    
    if (redshift.shear_photoz == 4) 
    {
      for (int i=0; i<tomo.shear_Nbin; i++) 
      {
        table[i + 1][0] = tomo.n_source[i];
        
        if (tomo.n_source[i] < 0.01 || tomo.n_source[i] > 100.) 
        {
          log_fatal("nsource(%d) = %e\n outside expected range", i, tomo.n_source[i]);
          exit(0);
        }
      }
    } 
    else 
    {
      const double norm = int_gsl_integrate_medium_precision(int_nsource, NULL, 
        redshift.shear_zdistrpar_zmin, redshift.shear_zdistrpar_zmax, NULL, GSL_WORKSPACE_SIZE);

      {
        const int i = 0;
        double ar[3] = {0.0, 0.0, 0.0}; // not used
        
        const double res = int_gsl_integrate_medium_precision(int_nsource, (void*) ar, 
          tomo.shear_zmin[i], tomo.shear_zmax[i], NULL, GSL_WORKSPACE_SIZE);
        
        table[i + 1][0] = res / norm * survey.n_gal;
        
        if (table[i + 1][0] < 0.01 || table[i + 1][0] > 100) 
        {
          log_fatal("nsource(%d) = %e\n outside expected range", i, table[i + 1][0]);
          exit(1);
        }
      }
      #pragma omp parallel for
      for (int i=1; i<tomo.shear_Nbin; i++) 
      {
        double ar[3] = {0.0, 0.0, 0.0}; // not used
        
        const double res = int_gsl_integrate_medium_precision(int_nsource, (void*) ar, 
          tomo.shear_zmin[i], tomo.shear_zmax[i], NULL, GSL_WORKSPACE_SIZE);
        
        table[i + 1][0] = res / norm * survey.n_gal;
        
        if (table[i + 1][0] < 0.01 || table[i + 1][0] > 100) 
        {
          log_fatal("nsource(%d) = %e\n outside expected range", i, table[i + 1][0]);
          exit(1);
        }
      }
    }

    table[0][0] = survey.n_gal;
  }
  
  if(ni < -1 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return table[ni + 1][0];
}

double int_nlens(double z, void *param __attribute__((unused))) 
{
  return pf_photoz(z, -1);
}

// returns n_gal for clustering tomography bin ni, works only with binned distributions
double nlens(int ni) // ni =-1 -> no tomography; ni>= 0 -> tomography bin ni
{
  static double** table = 0;

  if (table == 0 || table[0][0] != survey.n_lens) 
  {
    pf_photoz(0., 0);
    
    table = create_double_matrix(0, tomo.clustering_Nbin, 0, 1);
    
    if (redshift.clustering_photoz == 4) 
    {      
      for (int i=0; i<tomo.clustering_Nbin; i++) 
      {
        table[i + 1][0] = tomo.n_lens[i];
      }
    } 
    else 
    {
      double norm = 0.0;
      {
        double ar[2];
        ar[0] = tomo.clustering_zmin[0];
        ar[1] = tomo.clustering_zmax[tomo.clustering_Nbin - 1];
        
        norm = int_gsl_integrate_medium_precision(int_nlens, NULL, ar[0], ar[1], NULL, GSL_WORKSPACE_SIZE); 
      }

      {
        const int i = 0;
        double ar[2];
        ar[0] = tomo.clustering_zmin[i];
        ar[1] = tomo.clustering_zmax[i];
        
        const double res = int_gsl_integrate_medium_precision(int_nlens, (void*) ar, ar[0], 
          ar[1], NULL, GSL_WORKSPACE_SIZE);
        
        table[i + 1][0] = res / norm * survey.n_lens;
        
        if (table[i + 1][0] < 0.001 || table[i + 1][0] > 100) 
        {
          log_fatal("nlens(%d) = %e\n outside expected range", i, table[i + 1][0]);
          exit(1);
        }
      }     
      #pragma omp parallel for
      for (int i=1; i<tomo.clustering_Nbin; i++) 
      {
        double ar[2];
        ar[0] = tomo.clustering_zmin[i];
        ar[1] = tomo.clustering_zmax[i];
        
        const double res = int_gsl_integrate_medium_precision(int_nlens, (void*) ar, ar[0], 
          ar[1], NULL, GSL_WORKSPACE_SIZE);
        
        table[i + 1][0] = res / norm * survey.n_lens;
        
        if (table[i + 1][0] < 0.001 || table[i + 1][0] > 100) 
        {
          log_fatal("nlens(%d) = %e\n outside expected range", i, table[i + 1][0]);
          exit(1);
        }
      }
    }
    
    table[0][0] = survey.n_lens;
  }

  if(ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return table[ni + 1][0];
}

double int_for_zmean(double z, void* params) 
{
  double* ar = (double*) params;
  const int ni = (int) ar[0];
  if(ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return z * pf_photoz(z, ni);
}

double norm_for_zmean(double z, void* params) 
{
  double* ar = (double*) params;
  const int ni = (int) ar[0];
  if(ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return pf_photoz(z, ni);
}

double zmean(int ni) 
{ // mean true redshift of galaxies in tomography bin j
  static double** table = 0;
  if (table == 0) 
  {
    {
      double ar[1];
      ar[0] = pf_photoz(0., 0);
    }

    table = create_double_matrix(0, tomo.clustering_Nbin, 0, 1);

    {
      const int i = 0;
      double ar[1];
      ar[0] = (double) i;

      table[i][0] = int_gsl_integrate_low_precision(int_for_zmean, (void*) ar, 
        tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE) /
        int_gsl_integrate_low_precision(norm_for_zmean, (void*) ar, tomo.clustering_zmin[i],
        tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE);
    }
    #pragma omp parallel for
    for (int i=1; i<tomo.clustering_Nbin; i++) 
    {
      double ar[1];
      ar[0] = (double) i;
      table[i][0] = int_gsl_integrate_low_precision(int_for_zmean, (void*) ar, 
        tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE) /
        int_gsl_integrate_low_precision(norm_for_zmean, (void*) ar, tomo.clustering_zmin[i],
        tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE);
    }
  }
  if(ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return table[ni][0];
}

double int_for_zmean_source(double z, void* params) 
{
  double* ar = (double*) params;
  const int ni = (int) ar[0];
  if(ni < -1 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return z * zdistr_photoz(z, ni);
}

double zmean_source(int ni) 
{ // mean true redshift of source galaxies in tomography bin j
  static double** table = 0;
  
  if (table == 0) 
  {
    {
      double ar[1];
      ar[0] = zdistr_photoz(0., 0);
    }
    
    table = create_double_matrix(0, tomo.shear_Nbin - 1, 0, 1);

    {
      const int i = 0;
      double ar[1];
      ar[0] = (double) i;
      table[i][0] = int_gsl_integrate_medium_precision(int_for_zmean_source, (void*) ar, 
        tomo.shear_zmin[i], tomo.shear_zmax[i], NULL, GSL_WORKSPACE_SIZE);
    }
    #pragma omp parallel for
    for (int i=1; i<tomo.shear_Nbin; i++) 
    {
      double ar[1];
      ar[0] = (double) i;
      table[i][0] = int_gsl_integrate_medium_precision(int_for_zmean_source, (void*) ar, 
        tomo.shear_zmin[i], tomo.shear_zmax[i], NULL, GSL_WORKSPACE_SIZE);
    }
  }
  if(ni < 0 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return table[ni][0];
}

double max_g_tomo(int zs) 
{
  int i = 0;
  double max = 0;
  
  for (double z=0.; z<tomo.shear_zmax[zs]; z += tomo.shear_zmax[zs] / 50.) 
  {
    if(i != 0)
    {
      if(!fabs(z>0)) 
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
  if(!(aprime>0) || !(aprime<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double *ar = (double*)params;

  const int ni = (int) ar[0];
  if(ni < -1 || ni > tomo.shear_Nbin - 1)
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
  static double** table = 0;
  
  const double da = (0.999999 - 1.0/(redshift.shear_zdistrpar_zmax + 1.))/(Ntable.N_a - 1.0);
  const double z_shear_max = redshift.shear_zdistrpar_zmax;
  const double amin = 1. / (z_shear_max + 1.);

  if (table == 0 || recompute_zphot_shear(N) || recompute_cosmo3D(C)) 
  {
    if (table == 0) 
    {
      table = create_double_matrix(0, tomo.shear_Nbin, 0, Ntable.N_a - 1);
    }

    { // j=-1, no tomography is being done
      const int j = -1;
      {
        const int i = 0;
        const double a = amin + i*da;
        
        double ar[2];
        ar[0] = (double) j; // if j=-1, no tomography is being done
        ar[1] = chi(a);
        
        table[j + 1][i] = 
          int_gsl_integrate_medium_precision(int_for_g_tomo, (void*) ar, amin, a, NULL, 4000);
      }
      #pragma omp parallel for
      for (int i=1; i<Ntable.N_a; i++) 
      {
        const double a = amin + i*da;
        
        double ar[2];
        ar[0] = (double) j; // if j=-1, no tomography is being done
        ar[1] = chi(a);
        
        table[j + 1][i] = 
          int_gsl_integrate_medium_precision(int_for_g_tomo, (void*) ar, amin, a, NULL, 4000);
      }
    }
    
    if (tomo.shear_Nbin > 0) 
    {
      {
        const int j = 0;
        {
          const int i = 0;
          const double a = amin + i*da;
          
          double ar[2];
          ar[0] = (double) j;
          ar[1] = chi(a);
          
          table[j + 1][i] = 
            int_gsl_integrate_medium_precision(int_for_g_tomo, (void*) ar, amin, a, NULL, 4000);
        }
        #pragma omp parallel for
        for (int i=1; i<Ntable.N_a; i++)  
        {
          const double a = amin + i*da;
  
          double ar[2];
          ar[0] = (double) j;
          ar[1] = chi(a);
          
          table[j + 1][i] = 
            int_gsl_integrate_medium_precision(int_for_g_tomo, (void*) ar, amin, a, NULL, 4000);
        }
      }
      #pragma omp parallel for
      for (int j=1; j<tomo.shear_Nbin; j++) 
      {
        for (int i=0; i<Ntable.N_a; i++)  
        {
          const double a = amin + i*da;

          double ar[2];
          ar[0] = (double) j;
          ar[1] = chi(a);
          
          table[j + 1][i] = 
            int_gsl_integrate_medium_precision(int_for_g_tomo, (void*) ar, amin, a, NULL, 4000);
        }
      }
    }
    
    update_nuisance(&N);
    update_cosmopara(&C);
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
    res = interpol(table[ni + 1], Ntable.N_a, amin, 0.999999, da, ainput, 1.0, 1.0); 
  }

  return res;
}

double int_for_g2_tomo(double aprime, void* params) 
{ // \int n(z') W(z,z')^2 routines for source clustering
  if(!(aprime>0) || !(aprime<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double *ar = (double*) params;
  
  const int ni = (int) ar[0];
  if(ni < -1 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 

  const double chi1 = ar[1];
  const double chi_prime = chi(aprime);
  
  return zdistr_photoz(1./aprime-1., ni) * f_K(chi_prime - chi1)/
    f_K(chi_prime - chi1)/(f_K(chi_prime)*f_K(chi_prime))/(aprime*aprime);
}

double g2_tomo(double a, int ni) 
{ // for tomography bin ni
  static nuisancepara N;
  static cosmopara C;
  static double** table = 0;
  
  const double da = (0.999999-1./(redshift.shear_zdistrpar_zmax+1.))/(Ntable.N_a-1);

  if (table == 0 || recompute_zphot_shear(N) || recompute_cosmo3D(C)) 
  {
    if (table == 0) 
    {
      table = create_double_matrix(0, tomo.shear_Nbin, 0, Ntable.N_a-1);
    }

    { // j=-1, no tomography is being done
      const int j = -1;
      {
        const int i = 0;
        const double a = 1./(redshift.shear_zdistrpar_zmax + 1.) + i*da; 
        
        double ar[2];
        ar[0] = (double) j; 
        ar[1] = chi(a);
        
        table[j+1][i] = int_gsl_integrate_medium_precision(int_for_g2_tomo, (void*) ar,
          1./(redshift.shear_zdistrpar_zmax+1.), a, NULL, 4000);
      }      
      #pragma omp parallel for
      for (int i=1; i<Ntable.N_a; i++) 
      {
        const double a = 1./(redshift.shear_zdistrpar_zmax + 1.) + i*da; 
        
        double ar[2];
        ar[0] = (double) j; 
        ar[1] = chi(a);
        
        table[j+1][i] = int_gsl_integrate_medium_precision(int_for_g2_tomo, (void*) ar,
          1./(redshift.shear_zdistrpar_zmax+1.), a, NULL, 4000);
      }
    }
  	if (tomo.shear_Nbin > 0) 
    {
      {
        const int j = 0;
        {
          const int i = 0;
          const double a = 1./(redshift.shear_zdistrpar_zmax + 1.) + i*da;
          
          double ar[2];
          ar[0] = (double) j;
          ar[1] = chi(a);
          
          table[j+1][i] = int_gsl_integrate_medium_precision(int_for_g2_tomo,
            (void*) ar, 1./(redshift.shear_zdistrpar_zmax+1.), a, NULL, 4000);
        }
        #pragma omp parallel for
        for (int i=1; i<Ntable.N_a; i++) 
        {
          const double a = 1./(redshift.shear_zdistrpar_zmax + 1.) + i*da;
          
          double ar[2];
          ar[0] = (double) j;
          ar[1] = chi(a);
          
          table[j+1][i] = int_gsl_integrate_medium_precision(int_for_g2_tomo,
            (void*) ar, 1./(redshift.shear_zdistrpar_zmax+1.), a, NULL, 4000);
        }
      }
      #pragma omp parallel for
      for (int j=1; j<tomo.shear_Nbin; j++) 
      {
        for (int i=0; i<Ntable.N_a; i++) 
        {
          const double a = 1./(redshift.shear_zdistrpar_zmax+1.) + i*da;
          
          double ar[2];
          ar[0] = (double) j;          
          ar[1] = chi(a);
          
          table[j+1][i] = int_gsl_integrate_medium_precision(int_for_g2_tomo,
             (void*) ar, 1./(redshift.shear_zdistrpar_zmax+1.), a, NULL, 4000);
        }
      }
    }

    update_nuisance(&N);
    update_cosmopara(&C);
  }

  if(ni < -1 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 

  double res;
  if (a <= 1./(redshift.shear_zdistrpar_zmax+1.) || a > 1.0 - da) 
  {
    res = 0.0;
  }
  else 
  {
    const double amin =  1./(redshift.shear_zdistrpar_zmax + 1.);
    res = interpol(table[ni + 1], Ntable.N_a, amin, 0.999999, da, a, 1.0, 1.0);
  }
  
  return res;
}

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

double g_lens(double a, int ni) 
{ // for *lens* tomography bin ni
  static nuisancepara N;
  static cosmopara C;
  static double** table = 0;
  
  const double da = (0.999999 - 1. / (redshift.clustering_zdistrpar_zmax + 1.))/(Ntable.N_a - 1);

  if (table == 0 || recompute_zphot_clustering(N) || recompute_cosmo3D(C)) 
  {
    if (table == 0) 
    {
      table = create_double_matrix(0, tomo.clustering_Nbin, 0, Ntable.N_a - 1);
    }

    {
      const int j = -1;
      {
        const int i = 0;
        const double a = 1./(redshift.clustering_zdistrpar_zmax + 1.) + i*da;
        
        double ar[2];
        ar[0] = (double) j; // j=-1, no tomography is being done        
        ar[1] = chi(a);
        
        table[j + 1][i] = int_gsl_integrate_medium_precision(int_for_g_lens, (void*) ar,
          1. / (redshift.shear_zdistrpar_zmax + 1.), a, NULL, 4000);
      }
      #pragma omp parallel for
      for (int i=1; i<Ntable.N_a; i++) 
      {
        const double a = 1./(redshift.clustering_zdistrpar_zmax + 1.) + i*da;
        
        double ar[2];
        ar[0] = (double) j; // j=-1, no tomography is being done        
        ar[1] = chi(a);
        
        table[j + 1][i] = int_gsl_integrate_medium_precision(int_for_g_lens, (void*) ar,
          1. / (redshift.shear_zdistrpar_zmax + 1.), a, NULL, 4000);
      }
    }
    
    if (tomo.clustering_Nbin > 0) 
    {
      {
        const int j = 0;
        {
          const int i = 0;
          const double a = 1./(redshift.clustering_zdistrpar_zmax + 1.) + i*da;
          
          double ar[2];
          ar[0] = (double) j; 
          ar[1] = chi(a);
          
          table[j + 1][i] = int_gsl_integrate_medium_precision(int_for_g_lens,
            (void*) ar, 1./(redshift.shear_zdistrpar_zmax + 1.), a, NULL, 4000);
        }
        #pragma omp parallel for
        for (int i=1; i<Ntable.N_a; i++) 
        {
          const double a = 1./(redshift.clustering_zdistrpar_zmax + 1.) + i*da;
          
          double ar[2];
          ar[0] = (double) j; 
          ar[1] = chi(a);
          
          table[j + 1][i] = int_gsl_integrate_medium_precision(int_for_g_lens,
            (void*) ar, 1./(redshift.shear_zdistrpar_zmax + 1.), a, NULL, 4000);
        }
      }
      #pragma omp parallel for
      for (int j=1; j<tomo.clustering_Nbin; j++) 
      {
        double ar[2];
        ar[0] = (double) j; // if j=-1, no tomography is being done
        for (int i=0; i<Ntable.N_a; i++) 
        {
          const double a =  1./(redshift.clustering_zdistrpar_zmax + 1.) + i*da;
          
          ar[1] = chi(a);
          
          table[j + 1][i] = int_gsl_integrate_medium_precision(int_for_g_lens, (void*) ar,
            1. / (redshift.shear_zdistrpar_zmax + 1.), a, NULL, 4000);
        }
      }
    }

    update_nuisance(&N);
    update_cosmopara(&C);
  }

  if(ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  double res;
  if (a <= 1. / (redshift.clustering_zdistrpar_zmax + 1.) || a > 1.0 - da) 
  {
    res = 0.0;
  }
  else 
  {
    const double amin = 1. / (redshift.clustering_zdistrpar_zmax + 1.);
    res = interpol(table[ni + 1], Ntable.N_a, amin, 0.999999, da, a, 1.0, 1.0); 
  }
  
  return res;
}

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

  return f_K(chi_cmb-chi(a))/fchi_cmb;
}

double int_for_ggl_efficiency(double z, void* params) 
{
  if(!fabs(z>0)) 
  {
    log_fatal("z = 0 with 1/z needed");
    exit(1);
  }
  
  double* ar = (double*) params;
  const int ni = (int) ar[0];
  const int nj = (int) ar[1];
  if(ni < -1 || ni > tomo.clustering_Nbin - 1 || nj < -1 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  const double a = 1. / (1. + z);
  return pf_photoz(z, ni)*g_tomo(a, nj) * (1.0 + z)*f_K(chi(a));
}

double ggl_efficiency(int ni, int nj) 
{
  static double** table = 0;
  
  if (table == 0) 
  {
    {
      double init = pf_photoz(0, 0);
      init = zdistr_photoz(0, 0);
    }
    table = create_double_matrix(0, tomo.clustering_Nbin, 0, tomo.shear_Nbin);
    
    {
      const int i = 0;
      {
        const int j=0;

        double ar[2];
        ar[0] = (double) i;
        ar[1] = (double) j;
        table[i][j] = int_gsl_integrate_medium_precision(int_for_ggl_efficiency, (void*) ar,
          tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE)/max_g_tomo(j);
      }
      #pragma omp parallel for
      for (int j=1; j<tomo.shear_Nbin; j++) 
      {
        double ar[2];
        ar[0] = (double) i;
        ar[1] = (double) j;
        
        table[i][j] = int_gsl_integrate_medium_precision(int_for_ggl_efficiency, (void*) ar,
          tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE)/max_g_tomo(j);
      }
    }
    #pragma omp parallel for
    for (int i=1; i<tomo.clustering_Nbin; i++) 
    {
      for (int j=0; j<tomo.shear_Nbin; j++) 
      {
        double ar[2];
        ar[0] = (double) i;
        ar[1] = (double) j;
        table[i][j] = int_gsl_integrate_medium_precision(int_for_ggl_efficiency, (void*) ar,
          tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE)/max_g_tomo(j);
      }
    }
  }

  if(ni < 0 || ni > tomo.clustering_Nbin - 1 || nj < 0 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }  
  return table[ni][nj];
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// SELECTION FUNCTION (CLUSTER)
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double dV_cluster(double chi, double hoverh0)
{ 
  const double fK = f_K(chi);
  return fK*fK/hoverh0;
}

double dV_cluster2(double z, void* param __attribute__((unused)))
{ 
  const double a = 1./(1.+ z);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  return  dV_cluster(chidchi.chi, hoverh0);
}

double norm_z_cluster(int nz)
{
  static cosmopara C;
  static double **table;
  
  if (table == 0) 
  {
    table = create_double_matrix(0, tomo.cluster_Nbin - 1, 0, 1);
  }
  if(recompute_cosmo3D(C))
  {
    double ar[2] = {0.0, 0.0}; // not used
    {
      const int i = 0;
      table[i][0] = int_gsl_integrate_high_precision(dV_cluster2, (void*) ar, 
        tomo.cluster_zmin[i], tomo.cluster_zmax[i], NULL, GSL_WORKSPACE_SIZE);
    }
    #pragma omp parallel for
    for (int i=1; i<tomo.cluster_Nbin; i++)
    {
      table[i][0] = int_gsl_integrate_high_precision(dV_cluster2, (void*) ar, 
        tomo.cluster_zmin[i], tomo.cluster_zmax[i], NULL, GSL_WORKSPACE_SIZE);
    }
    update_cosmopara(&C);
  }
  return table[nz][0];
}

double zdistr_cluster(int nz, double z, double chi, double hoverh0) 
{ //simplfied selection function, disregards evolution of N-M relation+mass function within z bin
  if (z > tomo.cluster_zmax[nz] || z < tomo.cluster_zmin[nz]) 
  {
    return 0.0;
  }
  return dV_cluster(chi, hoverh0)/norm_z_cluster(nz);
}
