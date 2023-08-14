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


static int GSL_WORKSPACE_SIZE = 1024;
static double SQRT2 = 1.41421356237;

// In lens = source LSST systematic paper, Supranta made a case where photo-z = 4 
// with one important change: the test_overlap would be given by 
// tomo.clustering_zmax[ni] <= tomo.shear_zmin[nj]. We will call this IA = 7

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// integration boundary routines
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

double amin_source(int ni) 
{
  if(ni < -1 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  const int photoz = (redshift.shear_photoz == 7) ? 4 : redshift.shear_photoz;
  
  if (ni == -1 || photoz == 1 || photoz == 2 || photoz == 4) 
  {
    return 1. / (redshift.shear_zdistrpar_zmax + 1.);
  }
  else if (photoz == 0) 
  {
    return 1. / (1 + tomo.shear_zmax[ni]);
  }
  else if (photoz == 3)
  {
    return 1. /(1 + fmin(tomo.shear_zmax[ni] + 5.0 * nuisance.sigma_zphot_shear[ni] +
      fabs(nuisance.bias_zphot_shear[ni]), redshift.shear_zdistrpar_zmax));
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

  const int photoz = (redshift.shear_photoz == 7) ? 4 : redshift.shear_photoz;
  
  if (ni == -1 || photoz == 1 || photoz == 2 || photoz == 4) 
  {
    return 1. / (1. + fmax(redshift.shear_zdistrpar_zmin, 0.001));
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
  if(ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  const int photoz = (redshift.clustering_photoz == 7) ? 4 : redshift.clustering_photoz;
  
  if (ni == -1 || photoz == 1 || photoz == 2) 
  {
    return 1. / (redshift.clustering_zdistrpar_zmax + 1.);
  }

  if (photoz == 0) 
  {
    return 1. / (1 + tomo.clustering_zmax[ni]);
  }
  if (photoz == 4) 
  {
    return 1. / (1 + tomo.clustering_zmax[ni] + 2. * fabs(nuisance.bias_zphot_clustering[ni]));
  }
  if (photoz == 5) 
  {
    return 1.0/(redshift.clustering_zdistrpar_zmax + 1.0);
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
  
  const int photoz = (redshift.clustering_photoz == 7) ? 4 : redshift.clustering_photoz;

  if (ni == -1 || photoz == 1 || photoz == 2) 
  {
    return 1. / (1. + fmax(redshift.clustering_zdistrpar_zmin, 0.001));
  }
  else if (photoz == 0) 
  {
    return 1. / (1. + fmax(tomo.clustering_zmin[ni], 0.001));
  }
  else if (photoz == 4) 
  {
    return 1. / (1 + fmax(tomo.clustering_zmin[ni] -
      2. * fabs(nuisance.bias_zphot_clustering[ni]), 0.001));
  }
  
  return 1. / (1 + fmax(tomo.clustering_zmin[ni] - 5. * nuisance.sigma_zphot_clustering[ni] -
    fabs(nuisance.bias_zphot_clustering[ni]), 0.001));
}

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// redshift overlap tests, allowed tomography combinations
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

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
  const double kmax = 2.0*M_PI / R_min * cosmology.coverH0;
  
  int res = 0.0;
  if ((l + 0.5) / chiref[ni] < kmax)
  {
    res = 1;
  }

  return res;
}

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// Galaxy-Galaxy Lensing bins (redshift overlap tests)
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------


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

  if (sphotoz < 4 || sphotoz == 5 || sphotoz == 6 || sphotoz == 7) 
  {
    if (tomo.clustering_zmax[ni] <= tomo.shear_zmin[nj])
    {
      return 1;
    }
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
  
  if (ni < 0 || ni > tomo.clustering_Nbin - 1 || nj < 0 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
  return N[ni][nj];
}

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// Cluster-Galaxy Lensing bins (redshift overlap tests)
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// Cosmic Shear bins (redshift overlap tests)
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

// find z1 of tomography combination (z1, z2) constituting shear tomo bin Nbin
int Z1(int ni) 
{
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
    log_fatal("invalid bin input ni = %d (max = %d)", ni, tomo.shear_Npowerspectra);
    exit(1);
  }
  return N[ni];
}

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
    log_fatal("invalid bin input nj = %d (max = %d)", nj, tomo.shear_Npowerspectra);
    exit(1);
  }
  return N[nj];
}

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
  if (ni < 0 || ni > tomo.shear_Nbin - 1 || nj < 0 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d) (max = %d)", ni, nj, tomo.shear_Nbin);
    exit(1);
  }
  return N[ni][nj];
}

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// Galaxy Clustering (redshift overlap tests)
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

int ZCL1(int ni) 
{ // find ZCL1 of tomography combination (zcl1,zcl2) constituting g-clustering tomo bin Nbin
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
    log_fatal("invalid bin input ni = %d (max %d)", ni, tomo.clustering_Npowerspectra);
    exit(1);
  }
  return N[ni];
}

int ZCL2(int nj) 
{ // find ZCL2 of tomography combination (zcl1, zcl2) constituting g-clustering tomo bin Nbin
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
    log_fatal("invalid bin input nj = %d (max %d)", nj, tomo.clustering_Npowerspectra);
    exit(1);
  }
  return N[nj];
}

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
  if (ni < 0 || ni > tomo.clustering_Nbin - 1 || nj < 0 || nj > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
  return N[ni][nj];
}

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// Cluster Clustering (redshift overlap tests)
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// Cluster-Galaxy Cross Clustering
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

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

double sigma_zphot_shear(double z, int ni) 
{
  if (ni < 0 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  return (1.0 + z) * nuisance.sigma_zphot_shear[ni];
}

double zdistr_histo_1(double z, void* params __attribute__((unused))) 
{ // return nz(z) based on redshift file with one redshift distribution
  static double* tab = 0;
  static int zbins;
  static double zhisto_max; 
  static double zhisto_min;
  static double dz_histo;

  const int photoz = redshift.shear_photoz;

  if (tab == 0) 
  {
    zbins = line_count(redshift.shear_REDSHIFT_FILE);
    
    tab = (double*) malloc(sizeof(double)*zbins);
    if (tab == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
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
    
    if(photoz == 5)
    {
      for (int i=0; i<zbins; i++) 
      {
        double space1, space2, space3, space4;
        
        int count = fscanf(ein,"%le %le %le %le %le %le", &z_v[i], &space1, 
            &space2, &tab[i], &space3, &space4);
        if(count != 6)
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
        for (int j=0; j<tomo.shear_Nbin; j++) 
        {
          count = fscanf(ein,"%le", &space1);
          if(count != 1)
          {
            log_fatal("fscanf failed to read the file");
            exit(1);
          } 
        }
      } 
    }
    else
    {
      for (int i=0; i<zbins; i++) 
      {
        double space1, space2;
        
        int count = fscanf(ein, "%le %le %le %le\n", &z_v[i], &space1, 
            &space2, &tab[i]);
        if(count != 4)
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
      }
    }

    dz_histo = (z_v[zbins - 1] - z_v[0]) / ((double) zbins - 1.);
    zhisto_min = z_v[0];
    zhisto_max = z_v[zbins - 1] + dz_histo;

    redshift.shear_zdistrpar_zmin = zhisto_min;
    redshift.shear_zdistrpar_zmax = zhisto_max;

    fclose(ein);
    free(z_v);
  }

  double res = 0.0;
  if ((z >= zhisto_min) && (z < zhisto_max)) 
  {
    const int ni = (int) floor((z - zhisto_min) / dz_histo);  
    if(ni < 0 || ni > zbins - 1)
    {
      log_fatal("invalid bin input ni = %d", ni);
      exit(1);
    } 
    res = tab[ni];
  }
  return res;
}

double zdistr_histo_n(double z, void* params)
{ //VM: function is called if photo-z = 4
  // return nz(z, j) based on redshift file with structure 
  // z[i] nz[0][i] .. nz[tomo.shear_Nbin-1][i]
  static double** tab;
  static int zbins;
  static double zhisto_max; 
  static double zhisto_min; 
  static double dz_histo;
  
  if (tab == 0) 
  {
    zbins = line_count(redshift.shear_REDSHIFT_FILE);
    
    tab = (double**) malloc(sizeof(double*)*tomo.shear_Nbin);
    if (tab == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }
    for (int i=0; i<tomo.shear_Nbin; i++)
    {
      tab[i] = (double*) malloc(sizeof(double)*zbins);
      if (tab[i] == NULL)
      {
        log_fatal("array allocation failed");
        exit(1);
      }
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

      log_info("tomo.shear_zmin[%d] = %.3f,tomo.shear_zmax[%d] = %.3f", k,
        tomo.shear_zmin[k], k, tomo.shear_zmax[k]);
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
    
    if(ni < 0 || ni > tomo.shear_Nbin - 1 || nj < 0 || nj > zbins - 1)
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
  // PZ == 4 || PZ == 5
  static int zbins_file = -1;
  static double* z_v_file = 0;
  static double* nz_old = 0;
  static double* nz_diag = 0; 
  static double* nz_ext = 0;
  static double** nz_ext_bin = 0;

  const int photoz = (redshift.shear_photoz == 7) ? 4 : redshift.shear_photoz;

  if (photoz == -1) 
  {
    return n_of_z(zz, nj);
  } 
  else if ((photoz != 4 && recompute_zphot_shear(N)) || table == 0) 
  {
    update_nuisance(&N);
    
    if (table == 0) 
    {
      if (photoz == 4) 
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
        
        const double zmin_file = z_v_file[0];
        const double zmax_file = z_v_file[zbins_file - 1];
        const double dz_file = (zmax_file - zmin_file) / ((double) zbins_file - 1.0);

        redshift.shear_zdistrpar_zmin = fmax(zmin_file, 1.e-5);
        redshift.shear_zdistrpar_zmax = zmax_file + dz_file;

        fclose(ein);
      }
      else if (photoz == 5) 
      { // extreme outlier nz from simulation
        zbins_file = line_count(redshift.shear_REDSHIFT_FILE);
        z_v_file = (double*) malloc(sizeof(double)*zbins_file);
        if (z_v_file == NULL)
        {
          log_fatal("array allocation failed");
          exit(1);
        }
        nz_old = (double*) malloc(sizeof(double)*zbins_file);
        if (nz_old == NULL)
        {
          log_fatal("array allocation failed");
          exit(1);
        }
        nz_diag = (double*) malloc(sizeof(double)*zbins_file);
        if (nz_diag == NULL)
        {
          log_fatal("array allocation failed");
          exit(1);
        }
        nz_ext = (double*) malloc(sizeof(double)*zbins_file);
        if (nz_ext == NULL)
        {
          log_fatal("array allocation failed");
          exit(1);
        }
        nz_ext_bin = (double**) malloc(sizeof(double*)*tomo.shear_Nbin);
        if (nz_ext_bin == NULL)
        {
          log_fatal("array allocation failed");
          exit(1);
        }
        for (int i=0; i<tomo.shear_Nbin; i++)
        {
          nz_ext_bin[i] = (double*) malloc(sizeof(double)*zbins_file);
          if (nz_ext_bin[i] == NULL)
          {
            log_fatal("array allocation failed");
            exit(1);
          }
        }

        FILE* ein = fopen(redshift.shear_REDSHIFT_FILE, "r");
        if (ein == NULL)
        {
          log_fatal("file not opened");
          exit(1);
        }
        
        for (int i=0; i<zbins_file; i++) 
        {
          double space1, space2;
          int count = fscanf(ein,"%le %le %le %le %le %le", &z_v_file[i], &space1, 
            &space2, &nz_old[i], &nz_diag[i], &nz_ext[i]);
          if(count != 6)
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
          for(int k=0; k<tomo.shear_Nbin; k++) 
          {
            count = fscanf(ein, "%le", &nz_ext_bin[k][i]);
            if(count != 1)
            {
              log_fatal("fscanf failed to read the file");
              exit(1);
            }
          }
        }
        
        const double zmin_file = z_v_file[0];
        const double zmax_file = z_v_file[zbins_file - 1];
        const double dz_file = (zmax_file - zmin_file) / ((double) zbins_file - 1.0);

        redshift.shear_zdistrpar_zmin = zmin_file;
        redshift.shear_zdistrpar_zmax = zmax_file + dz_file;

        fclose(ein);
      }
      else
      {
        // VM: this seems not well designed for analytical photo-z
        zbins_file = line_count(redshift.shear_REDSHIFT_FILE);
      }

      zbins = zbins_file * Ntable.acc_boost_photoz_sampling; //VM: upsampling if needed
      if (photoz != 4 && photoz != 5) 
      {
        zbins *= 20;
      } // upsample if convolving with analytic photo-z model
      
      table = (double**) malloc(sizeof(double*)*(tomo.shear_Nbin+1));
      if (table == NULL)
      {
        log_fatal("array allocation failed");
        exit(1);
      }
      for (int i=0; i<(tomo.shear_Nbin+1); i++)
      {
        table[i] = (double*) malloc(sizeof(double)*zbins);
        if (table[i] == NULL)
        {
          log_fatal("array allocation failed");
          exit(1);
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
    }

    const double zhisto_min = redshift.shear_zdistrpar_zmin;
    const double zhisto_max = redshift.shear_zdistrpar_zmax;
    const double dz_histo = (zhisto_max - zhisto_min) / ((double) zbins);  
    for (int i=0; i<zbins; i++) 
    {
      z_v[i] = zhisto_min + (i + 0.5) * dz_histo;
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

    switch (photoz) 
    { // initialization that allows openmp threading of the next loop
      case 0:
      {
        zdistr_histo_1(0.0, NULL);
        break;
      }
      case 1:
      {
        zdistr_histo_1(0.0, NULL);
        break;
      }
      case 2:
      {
        zdistr_histo_1(0.0, NULL);
        break;
      }
      case 3:
      {
        zdistr_histo_1(0.0, NULL);
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
      case 6:
      {
        zdistr_histo_1(0.0, NULL);
        break;
      }
    }

    double NORM[MAX_SIZE_ARRAYS];
    #pragma omp parallel for
    for (int i=0; i<tomo.shear_Nbin; i++) 
    {
      double norm;

      switch (photoz) 
      {
        case 0: // no photo-zs, split 'true' n(z) histogram in bins
        {  
          double ar[2] = {tomo.shear_zmin[i], tomo.shear_zmax[i]};

          norm = like.high_def_integration > 0 ?
            int_gsl_integrate_high_precision(zdistr_histo_1, NULL, tomo.shear_zmin[i],
              tomo.shear_zmax[i], NULL, GSL_WORKSPACE_SIZE) :
            int_gsl_integrate_medium_precision(zdistr_histo_1, NULL, tomo.shear_zmin[i],
              tomo.shear_zmax[i], NULL, GSL_WORKSPACE_SIZE);
          
          for (int k=0; k<zbins; k++) 
          {
            table[i+1][k] = 0;
            if (z_v[k] >= ar[0] && z_v[k] <= ar[1]) 
            {
              table[i + 1][k] = zdistr_histo_1(z_v[k], NULL) / norm;
            }
          }
          break;
        }
        case 1: // Pseudo Voigt (Lorentzian + Gaussian) photo-zs
        {          
          if (sigma_zphot_shear(0., i) == 0.) 
          {
            log_fatal("Source galaxy photo-z model underdetermined! redshift.shear_photoz = 1,"
              "but nuisance.sigma_zphot_shear[%d] not set", i);
            exit(1);
          }

          norm = 0.0;
          for (int k=0; k<zbins; k++) 
          {
            const double x1 = (tomo.shear_zmin[i] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
                 (M_SQRT2 * sigma_zphot_shear(z_v[k], i));
            const double x2 = (tomo.shear_zmax[i] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
                 (M_SQRT2 * sigma_zphot_shear(z_v[k], i));
            
            
            gsl_sf_result result;
            {
              int status = gsl_sf_erf_e(x1, &result);
              if (status)
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
            }
            const double erf_sf_x1 = result.val;
            {
              int status = gsl_sf_erf_e(x2, &result);
              if (status)
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
            }
            const double erf_sf_x2 = result.val;

            table[i + 1][k] = 
              zdistr_histo_1(z_v[k], NULL) * (eta * 0.3183 * (atan(x2) - atan(x1)) +
              ((1.0 - eta) * 0.5 * (erf_sf_x2 - erf_sf_x1)));
            // this creates a pseudo Voigt profile by adding a Gaussian and Lorentian (convolved
            // with a tophat) with the correct weights. See, eg, the Wikipedia article on Voigt 
            // profiles for an explanation of where the numbers come from.
            
            norm += table[i + 1][k] * dz_histo;
          }
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] /= norm;
          }
          break;
        }
        case 2: // Pseudo Voigt (Lorentzian + Gaussian) + outlier (currently flat) photo-zs
        {
          if (sigma_zphot_shear(0., i) == 0.) 
          {
            log_fatal("Source galaxy photo-z model underdetermined!\nredshift.shear_photoz = 2,"
              " but nuisance.sigma_zphot_shear[%d] not set",i);
            exit(1);
          }

          norm = 0.0;
          for (int k=0; k<zbins; k++) 
          {
            const double x1 = (tomo.shear_zmin[i] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
                 (M_SQRT2 * sigma_zphot_shear(z_v[k], i));

            const double x2 = (tomo.shear_zmax[i] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
                 (M_SQRT2 * sigma_zphot_shear(z_v[k], i));

            gsl_sf_result result;
            {
              int status = gsl_sf_erf_e(x1, &result);
              if (status)
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
            }
            const double erf_sf_x1 = result.val;
            {
              int status = gsl_sf_erf_e(x2, &result);
              if (status)
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
            }
            const double erf_sf_x2 = result.val;

            table[i + 1][k] = zdistr_histo_1(z_v[k], NULL) * (eta * 0.3183 * (atan(x2) - atan(x1)) +
              ((1.0 - eta) * 0.5 * (erf_sf_x2 - erf_sf_x1)));

            norm += table[i + 1][k] * dz_histo;
          }
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] = outfrac / dz_histo + (1. - outfrac) * table[i + 1][k] / norm;
          }
          break;
        }
        case 3: // Gaussian photo-zs
        {
          if (sigma_zphot_shear(0., i) == 0.) 
          {
            log_fatal("Source galaxy photo-z model underdetermined!\n"
              "redshift.shear_photoz = 3, "
              "but nuisance.sigma_zphot_shear[%d] not set", i);
            exit(1);
          }

          norm = 0.0;
          for (int k=0; k<zbins; k++) 
          {
            const double x1 = (tomo.shear_zmin[i] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
                 (M_SQRT2 * sigma_zphot_shear(z_v[k], i));
            const double x2 = (tomo.shear_zmax[i] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
                 (M_SQRT2 * sigma_zphot_shear(z_v[k], i));
            
            gsl_sf_result result;
            {
              int status = gsl_sf_erf_e(x1, &result);
              if (status)
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
            }
            const double erf_sf_x1 = result.val;
            {
              int status = gsl_sf_erf_e(x2, &result);
              if (status)
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
            }
            const double erf_sf_x2 = result.val;

            table[i + 1][k] = 0.5 * zdistr_histo_1(z_v[k], NULL) * (erf_sf_x2 - erf_sf_x1);
            
            norm += table[i + 1][k] * dz_histo;
          }
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] /= norm;
          }
          break;
        }
        case 4: // histogram file contains n(z) estimates for each bin
        { 
          norm = 0.0;
          double ar[1] = {(double) i};

          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] = zdistr_histo_n(z_v[k], (void*) ar); 
            norm += table[i + 1][k] * dz_histo;
          }
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] /= norm;
          }
          break;
        }
        case 5: 
        { //Gaussian photo-zs + outlier probability matrix
          
          if (sigma_zphot_shear(0., i) == 0.) 
          {
            log_fatal("Source galaxy photo-z model underdetermined!\n"
              "redshift.shear_photoz = 5, "
              "but nuisance.sigma_zphot_shear[%d] not set", i);
            exit(1);
          }

          norm = 0.0;
          const double zmin_file = z_v_file[0];
          const double dz_file = 
            (z_v_file[zbins_file - 1] - zmin_file)/((double) zbins_file - 1.0);
          
          for (int k = 0; k<zbins; k++)
          {
            const double x1 = (tomo.shear_zmin[i] - z_v[k] + bias_zphot_shear(z_v[k], i)) / 
                  (M_SQRT2 * sigma_zphot_shear(z_v[k], i));
            const double x2 = (tomo.shear_zmax[i] - z_v[k] + bias_zphot_shear(z_v[k], i)) / 
                  (M_SQRT2 * sigma_zphot_shear(z_v[k], i));
            
            const int z_index = (int) floor((z_v[k] - zmin_file)/dz_file);
            
            if ((z_index >= 0) && z_index < zbins_file) 
            {
              gsl_sf_result result;
              {
                int status = gsl_sf_erf_e(x1, &result);
                if (status)
                {
                  log_fatal(gsl_strerror(status));
                  exit(1);
                }
              }
              const double erf_sf_x1 = result.val;
              {
                int status = gsl_sf_erf_e(x2, &result);
                if (status)
                {
                  log_fatal(gsl_strerror(status));
                  exit(1);
                }
              }
              const double erf_sf_x2 = result.val;

              table[i + 1][k] = nz_ext_bin[i][z_index] + 
                  0.5*nz_diag[z_index]*(erf_sf_x2 - erf_sf_x1);
            }
            else
            {
              log_fatal("bad index z_index");
              exit(1);
            }
            norm += table[i + 1][k]*dz_histo;
          }
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] /= norm;
          }
          break;
        }
        case 6:
        { // simple outlier model: 2 islands + Gaussian photo-zs
          if (sigma_zphot_shear(0., i) == 0.) 
          {
            log_fatal("Source galaxy photo-z model underdetermined!\n"
              "redshift.shear_photoz = 5, "
              "but nuisance.sigma_zphot_shear[%d] not set", i);
            exit(1);
          }

          double ar[2] = {tomo.shear_zmin[i], tomo.shear_zmax[i]};

          int zpart = 0;
          double zlow, zhigh;
          if (ar[0] < 0.5) 
          {
            zlow = ar[0];
            zhigh = ((ar[1] >= 0.5) ? 0.5 : ar[1]);
            zpart = 1; // denote low-z
          }
          else if (ar[1] > 2.0) 
          {
            zlow = ((ar[0] <= 2.0) ? 2.0 : ar[0]);
            zhigh = ar[1];
            zpart = 2; // denote high-z
          }

          double n_out = 0.;
          if(zpart > 0) 
          {
            for (int k = 0; k<zbins; k++)
            {
              if ((z_v[k] >= zlow) && (z_v[k] <= zhigh))
              {
                n_out += zdistr_histo_1(z_v[k], NULL);
              }
            }
            if (zpart == 1)
            {
              n_out *= (nuisance.frac_lowz * dz_histo);
            }
            else if (zpart == 2)
            {
              n_out *= (nuisance.frac_highz * dz_histo);
            }
          }

          norm = 0.0;
          for (int k = 0; k<zbins; k++)
          {
            const double x1 = (ar[0] - z_v[k] + bias_zphot_shear(z_v[k], i)) / 
                  (M_SQRT2 * sigma_zphot_shear(z_v[k], i));
            const double x2 = (ar[1] - z_v[k] + bias_zphot_shear(z_v[k], i)) / 
                  (M_SQRT2 * sigma_zphot_shear(z_v[k], i));

            gsl_sf_result result;
            {
              int status = gsl_sf_erf_e(x1, &result);
              if (status)
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
            }
            const double erf_sf_x1 = result.val;
            {
              int status = gsl_sf_erf_e(x2, &result);
              if (status)
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
            }
            const double erf_sf_x2 = result.val;            

            const double DELTA = (erf_sf_x2 - erf_sf_x1);
            const double ZDIST = zdistr_histo_1(z_v[k], NULL);
            
            if (zpart == 0) 
            {
              table[i+1][k] = 0.5*ZDIST*DELTA;
            }
            else if (zpart == 1) 
            {
              if (z_v[k] <= 0.5)
              {
                table[i+1][k] = 0.5*(1.0 - nuisance.frac_lowz)*ZDIST*DELTA;
              }
              else if (z_v[k] > 2.0)
              {
                table[i+1][k] = n_out/(3.5 - 2.0) + 0.5*(1.0 - nuisance.frac_highz)*ZDIST*DELTA;
              }
              else
              {
                table[i+1][k] = 0.5*ZDIST*DELTA;
              }
            }
            else if (zpart == 2) 
            {
              if(z_v[k] <= 0.5)
              {
                table[i + 1][k] = n_out/(0.5) + 0.5*(1.0 - nuisance.frac_lowz)*ZDIST*DELTA;
              }
              else if (z_v[k] > 2.0)
              {
                table[i + 1][k] = 0.5*(1.0 - nuisance.frac_highz)*ZDIST*DELTA;
              }
              else
              {
                table[i + 1][k] = 0.5*ZDIST*DELTA;
              }
            }

            norm += table[i+1][k]*dz_histo;
          }
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] /= norm;
          }
          break;
        }
        default:
        {
          log_fatal("redshift.shear_photoz = %d not supported in this cosmolike version", photoz);
          exit(1);          
        }
      }

      NORM[i] = norm;
    }
    
    // calculate normalized overall redshift distribution (without bins), store in table[0][:]
    double norm = 0;
    for (int i=0; i<tomo.shear_Nbin; i++) 
    {
      norm += NORM[i];
    }
    #pragma omp parallel for
    for (int k=0; k<zbins; k++) 
    {
      table[0][k] = 0;
      for (int i=0; i<tomo.shear_Nbin; i++) 
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
    log_fatal("nj = %d bin outside range (max = %d)", nj, tomo.shear_Nbin);
    exit(1);
  }
  
  if (photoz == 4) 
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
    int status = gsl_spline_eval_e(photoz_splines[nj + 1], zz, NULL, &res);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
  }

  return res;
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// Clustering routines for redshift distributions, including photo-zs * (optional)
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double bias_zphot_clustering(double z, int ni) 
{
  if (ni < 0 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return (1. + z) * nuisance.bias_zphot_clustering[ni];
}

double sigma_zphot_clustering(double z, int ni) 
{
  if (ni < 0 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return (1. + z) * nuisance.sigma_zphot_clustering[ni];
}

double pf_histo(double z, void* params __attribute__((unused))) 
{ // return pf(z) based on redshift file with one redshift distribution
  static double* tab = 0;
  static int zbins = 0;
  static double zhisto_min = 0; 
  static double zhisto_max = 0; 
  static double dz_histo = 0;

  const int photoz = (redshift.clustering_photoz == 7) ? 4 : redshift.clustering_photoz;

  if (tab == 0) 
  {
    zbins = line_count(redshift.clustering_REDSHIFT_FILE);
    
    tab = (double*) malloc(sizeof(double)*zbins);
    if (tab == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
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

    if (photoz == 5)
    {
      for (int i=0; i<zbins; i++) 
      {
        double space1, space2, space3, space4;

        int count = fscanf(ein, "%le %le %le %le %le %le", &z_v[i], &space1, 
          &space2, &tab[i], &space3, &space4);
        if(count != 6)
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
      }
    }
    else
    {
      for (int i=0; i<zbins; i++) 
      {
        double space1, space2;
        
        int count = fscanf(ein, "%le %le %le %le\n", &z_v[i], &space1, 
          &space2, &tab[i]);
        if(count != 4)
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
      }
    }

    dz_histo = (z_v[zbins - 1]  - z_v[0]) / ((double) zbins - 1.0);
    zhisto_min = z_v[0];
    zhisto_max = z_v[zbins - 1] + dz_histo;
    
    fclose(ein);
    free(z_v);

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
    const int ni = (int) floor((z - zhisto_min) / dz_histo);
    
    if (ni < 0 || ni > zbins - 1)
    {
      log_fatal("invalid bin input ni = %d", ni);
      exit(1);
    } 
    res = tab[ni];
  }
  return res;
}

// Function only called on photo-z = 4
double pf_histo_n(double z, void* params) 
{ // return pf(z, nj) based on z file with structure z[i] nz[0][i] .. nz[tomo.clustering_Nbin-1][i]
  static double** tab;
  static int zbins = 0; 
  static double zhisto_min = 0;
  static double zhisto_max = 0; 
  static double dz_histo = 0;

  if (tab == 0) 
  {
    zbins = line_count(redshift.clustering_REDSHIFT_FILE);
    
    tab = (double**) malloc(sizeof(double*)*tomo.clustering_Nbin);
    if (tab == NULL)
    {
      log_fatal("array allocation failed");
      exit(1);
    }
    for (int i=0; i<tomo.clustering_Nbin; i++)
    {
      tab[i] = (double*) malloc(sizeof(double)*zbins);
      if (tab[i] == NULL)
      {
        log_fatal("array allocation failed");
        exit(1);
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
        count = fscanf(ein, "%le", &tab[k][i]);
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
    
    // VM: added lines below (JULY 2003). They are already in `pf_photoz` function
    redshift.clustering_zdistrpar_zmin = fmax(zhisto_min, 1.e-5);
    redshift.clustering_zdistrpar_zmax = zhisto_max;

    // now, set tomography bin boundaries
    for (int k=0; k<tomo.clustering_Nbin; k++) 
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

    free(z_v);
    
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
    const int nj = (int) floor((z - zhisto_min) / dz_histo);
    
    if(ni < 0 || ni > tomo.clustering_Nbin - 1 || nj < 0 || nj > zbins - 1)
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
  static int zbins_file = -1;    // PZ == 4 || PZ == 5
  static double* z_v_file=0;     // PZ == 4 || PZ == 5
  static double* nz_old=0;       // PZ == 5
  static double* nz_diag=0;      // PZ == 5
  static double* nz_ext=0;       // PZ == 5
  static double** nz_ext_bin=0;  // PZ == 5

  const int photoz = (redshift.clustering_photoz == 7) ? 4 : redshift.clustering_photoz;

  if (photoz == -1) 
  {
    return n_of_z(zz, nj);
  }
  else if ((photoz != 4 && recompute_zphot_clustering(N)) || table == 0) 
  {
    update_nuisance(&N);

    if (table == 0) 
    {
      if (photoz == 4) 
      { 
        zbins_file = line_count(redshift.clustering_REDSHIFT_FILE);
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
            fscanf(ein, "%le", &space);
            if(count != 1)
            {
              log_fatal("fscanf failed to read the file");
              exit(1);
            }
          }
        }
        
        const double zmin_file = z_v_file[0];
        const double zmax_file = z_v_file[zbins_file - 1];
        const double dz_file = (zmax_file - zmin_file)/ ((double) zbins_file - 1.0);
        
        redshift.clustering_zdistrpar_zmin = fmax(zmin_file, 1.e-5);
        redshift.clustering_zdistrpar_zmax = zmax_file + dz_file;
        
        fclose(ein);
      }
      else if (photoz == 5)
      {
        zbins_file = line_count(redshift.clustering_REDSHIFT_FILE);
        z_v_file = (double*) malloc(sizeof(double)*zbins_file);
        if (z_v_file == NULL)
        {
          log_fatal("array allocation failed");
          exit(1);
        }
        nz_old  = (double*) malloc(sizeof(double)*zbins_file);
        if (nz_old == NULL)
        {
          log_fatal("array allocation failed");
          exit(1);
        }
        nz_diag = (double*) malloc(sizeof(double)*zbins_file);
        if (nz_diag == NULL)
        {
          log_fatal("array allocation failed");
          exit(1);
        }
        nz_ext  = (double*) malloc(sizeof(double)*zbins_file);
        if (nz_ext == NULL)
        {
          log_fatal("array allocation failed");
          exit(1);
        }
        nz_ext_bin = (double**) malloc(sizeof(double*)*tomo.clustering_Nbin);
        if (nz_ext_bin == NULL)
        {
          log_fatal("array allocation failed");
          exit(1);
        }
        for (int i=0; i<tomo.clustering_Nbin; i++)
        {
          nz_ext_bin[i] = (double*) malloc(sizeof(double)*zbins_file);
          if (nz_ext_bin[i] == NULL)
          {
            log_fatal("array allocation failed");
            exit(1);
          }
        }

        FILE* ein = fopen(redshift.clustering_REDSHIFT_FILE, "r");
        if (ein == NULL)
        {
          log_fatal("file not opened");
          exit(1);
        }

        for (int i=0; i<zbins_file; i++) 
        {
          double space1, space2;       
          int count = fscanf(ein,"%le %le %le %le %le %le", &z_v_file[i], &space1, 
              &space2, &nz_old[i], &nz_diag[i], &nz_ext[i]);
          if(count != 6)
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
          for(int k=0; k<tomo.clustering_Nbin; k++) 
          {
            count = fscanf(ein,"%le", &nz_ext_bin[k][i]);
            if(count != 1)
            {
              log_fatal("fscanf failed to read the file");
              exit(1);
            }
          }
        }

        const double zmin_file = z_v_file[0];
        const double zmax_file = z_v_file[zbins_file - 1];
        const double dz_file = (zmax_file - zmin_file)/ ((double) zbins_file - 1.0);

        redshift.clustering_zdistrpar_zmin = zmin_file;
        redshift.clustering_zdistrpar_zmax = zmax_file + dz_file;
        
        fclose(ein);
      }
      else
      {
        // VM: this seems not well designed for analytical photo-z
        zbins_file = line_count(redshift.clustering_REDSHIFT_FILE);
      }

      zbins = zbins_file * Ntable.acc_boost_photoz_sampling; //VM: upsampling if needed
      if (photoz != 0 && photoz != 4 && photoz != 5) 
      {
        pf_histo(0.5, NULL);
        zbins *= 20;
      } // upsample if convolving with analytic photo-z model
      
      table = (double**) malloc(sizeof(double*)*(tomo.clustering_Nbin+1));
      if (table == NULL)
      {
        log_fatal("array allocation failed");
        exit(1);
      }
      for (int i=0; i<(tomo.clustering_Nbin+1); i++)
      {
        table[i] = (double*) malloc(sizeof(double)*zbins);
        if (table[i] == NULL)
        {
          log_fatal("array allocation failed");
          exit(1);
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
    }

    const double zhisto_min = redshift.clustering_zdistrpar_zmin;
    const double zhisto_max = redshift.clustering_zdistrpar_zmax;
    const double dz_histo = (zhisto_max - zhisto_min) / ((double) zbins);
    for (int i=0; i<zbins; i++) 
    {
      z_v[i] = zhisto_min + (i + 0.5) * dz_histo;
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

    switch (photoz) 
    { // initialization that allows openmp threading of the next loop
      case 0:
      {
        pf_histo(0.0, NULL);
        break;
      }
      case 1:
      {
        pf_histo(0.0, NULL);
        break;
      }
      case 2:
      {
        pf_histo(0.0, NULL);
        break;
      }
      case 3:
      {
        pf_histo(0.0, NULL);
        break;
      }
      case 4:
      {
        const int i = 0;
        double ar[1];        
        ar[0] = (double) i;
        pf_histo_n(0., (void*) ar);
        
        break;
      }
      case 6:
      {
        pf_histo(0.0, NULL);
        break;
      }
    }

    double NORM[MAX_SIZE_ARRAYS]; 
    #pragma omp parallel for
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    {
      double norm = 0.0;

      switch (photoz) 
      {
        case 0: // no photo-zs, split 'true' n(z) histogram in bins
        {
          norm = like.high_def_integration > 0  ?
            int_gsl_integrate_high_precision(pf_histo, NULL, tomo.clustering_zmin[i],
              tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE) :
            int_gsl_integrate_medium_precision(pf_histo, NULL, tomo.clustering_zmin[i],
              tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE);

          for (int k=0; k<zbins; k++)
          {
            table[i + 1][k] = 0.;
            if (z_v[k] >= tomo.clustering_zmin[i] && z_v[k] <= tomo.clustering_zmax[i]) 
            {
              table[i + 1][k] = pf_histo(z_v[k], NULL) / norm;
            }
          } 
          break;
        }
        case 1: // Pseudo Voigt (Lorentzian + Gaussian) photo-zs
        {
          if (sigma_zphot_clustering(0., i) == 0.) 
          {
            log_fatal("Lens galaxy photo-z model underdetermined!\nredshift.clustering_photoz = 1, "
              "but nuisance.sigma_zphot_clustering[%d] not set", i);
            exit(1);
          }

          norm = 0.0;
          for (int k=0; k<zbins; k++) 
          {            
            const double x1 = (tomo.clustering_zmin[i] - z_v[k] + 
              bias_zphot_clustering(z_v[k], i))/(M_SQRT2*sigma_zphot_clustering(z_v[k], i));
            const double x2 = (tomo.clustering_zmax[i] - z_v[k] + 
                bias_zphot_clustering(z_v[k], i))/(M_SQRT2*sigma_zphot_clustering(z_v[k], i));
            
            table[i + 1][k] = 
              pf_histo(z_v[k], NULL)*(eta * 0.3183 * (atan(x2) - atan(x1)) +
              ((1.0 - eta) * 0.5 * (gsl_sf_erf(x2) - gsl_sf_erf(x1))));
            
            // this creates a pseudo Voigt profile by adding a Gaussian and Lorentian (convolved  
            // with a tophat) with the correct weights. See, eg, the Wikipedia article on Voigt
            // profiles for an explanation of where the numbers come from.
            norm += table[i + 1][k] * dz_histo;
          }
          
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] /= norm;
          }
          break;
        }
        case 2: // Pseudo Voigt (Lorentzian + Gaussian) + outlier (currently flat) photo-zs
        {
          if (sigma_zphot_clustering(0., i) == 0.) 
          {
            log_fatal("Lens galaxy photo-z model underdetermined!\n redshift.clustering_photoz = 2,"
              "but nuisance.sigma_zphot_clustering[%d] not set", i);
            exit(1);
          }
          
          norm = 0.0;
          for (int k=0; k<zbins; k++)
          {            
            const double x1 = (tomo.clustering_zmin[i] - z_v[k] + 
              bias_zphot_clustering(z_v[k], i))/(M_SQRT2*sigma_zphot_clustering(z_v[k], i));
            const double x2 = (tomo.clustering_zmax[i] - z_v[k] + 
              bias_zphot_clustering(z_v[k], i))/(M_SQRT2*sigma_zphot_clustering(z_v[k], i));
            
            table[i + 1][k] = pf_histo(z_v[k], NULL)*(eta * 0.3183 * (atan(x2) - atan(x1)) +
              ((1.0 - eta) * 0.5 * (gsl_sf_erf(x2) - gsl_sf_erf(x1))));
            
            norm += table[i + 1][k] * dz_histo;
          }
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] = outfrac / dz_histo + (1. - outfrac) * table[i + 1][k] / norm;
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
         
          norm = 0.0;
          for (int k=0; k<zbins; k++)
          {
            const double x1 = 
                (tomo.clustering_zmin[i] - z_v[k] + bias_zphot_clustering(z_v[k], i)) /
                (M_SQRT2*sigma_zphot_clustering(z_v[k], i));
            const double x2 = 
                (tomo.clustering_zmax[i] - z_v[k] + bias_zphot_clustering(z_v[k], i)) /
                (M_SQRT2*sigma_zphot_clustering(z_v[k], i));
            
            table[i + 1][k] = 0.5 * pf_histo(z_v[k], NULL)*(gsl_sf_erf(x2) - gsl_sf_erf(x1));
            
            norm += table[i + 1][k] * dz_histo;
          }
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] /= norm;
          } 
          break;
        }
        case 4: // histogram file contains n(z) estimates for each bin
        { 
          double ar[1] = {(double) i};

          norm = like.high_def_integration > 0 ?
            int_gsl_integrate_high_precision(pf_histo_n, (void*) ar,
              tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE) :
            int_gsl_integrate_medium_precision(pf_histo_n, (void*) ar,
              tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE);

          if (norm == 0) 
          {
            log_fatal("pf_photoz: norm(nz = %d) = 0", i);
            exit(1);
          }
          for (int k=0; k<zbins; k++) 
          {
            table[i + 1][k] = pf_histo_n(z_v[k], (void*) ar)/norm;
          }
          break;
        }
        case 5: //Gaussian photo-zs + outlier probability matrix
        {
          if (sigma_zphot_clustering(0., i) == 0.) 
          {
            log_fatal("Lens galaxy photo-z model underdetermined!\nredshift.clustering_photoz = 3,"
              "but nuisance.sigma_zphot_clustering[%d] not set", i);
            exit(1);
          }
          
          norm = 0.0;
          const double zmin_file = z_v_file[0];
          const double dz_file = 
            (z_v_file[zbins_file - 1] - zmin_file)/((double) zbins_file - 1.0);
  
          for (int k=0; k<zbins; k++)
          {
            const double x1 = (tomo.clustering_zmin[i] - z_v[k] + 
              bias_zphot_clustering(z_v[k], i))/(M_SQRT2*sigma_zphot_clustering(z_v[k], i));
            const double x2 = (tomo.clustering_zmax[i] - z_v[k] + 
              bias_zphot_clustering(z_v[k], i))/(M_SQRT2*sigma_zphot_clustering(z_v[k], i));
            
            const int z_index = (int) floor((z_v[k] - zmin_file)/dz_file);
            if ((z_index >= 0) && z_index < zbins_file) 
            {
              table[i+1][k] = nz_ext_bin[i][z_index] 
                + 0.5*nz_diag[z_index]*(gsl_sf_erf(x2) - gsl_sf_erf(x1));
            }
            else
            {
              log_fatal("bad index z_index");
              exit(1);
            }

            norm += table[i+1][k]*dz_file;
          }
          break;
        }
        case 6: // simple outlier model: 2 islands + Gaussian photo-zs
        { 
           if (sigma_zphot_clustering(0., i) == 0.) 
          {
            log_fatal("Lens galaxy photo-z model underdetermined!\nredshift.clustering_photoz = 3,"
              "but nuisance.sigma_zphot_clustering[%d] not set", i);
            exit(1);
          }

          double ar[2] = {tomo.clustering_zmin[i], tomo.clustering_zmax[i]};
          
          int zpart=0;
          double zlow, zhigh;
          if (ar[0] < 0.5) 
          {
            zlow = ar[0];
            zhigh = (ar[1] >= 0.5) ? 0.5 : ar[1];
            zpart = 1; // denote low-z
          }
          else if (ar[1] > 2.0) 
          {
            zhigh = ar[1];
            zlow = (ar[0] <= 2.0) ? 2.0 : ar[0];
            zpart = 2; // denote high-z
          }

          double n_out=0.;
          if (zpart > 0) 
          {
            for (int k = 0;k<zbins; k++)
            {
              if ((z_v[k] >= zlow) && (z_v[k] <= zhigh))
              {
                n_out += pf_histo(z_v[k], NULL);
              }
            }
            if (zpart == 1)
            {
              n_out *= (nuisance.frac_lowz * dz_histo);
            }
            else if (zpart == 2)
            {
              n_out *= (nuisance.frac_highz * dz_histo);
            }
          }

          norm = 0.;
          for (int k = 0;k<zbins; k++)
          {            
            const double x1 = 
                (tomo.clustering_zmin[i] - z_v[k] + bias_zphot_clustering(z_v[k], i)) /
                (M_SQRT2*sigma_zphot_clustering(z_v[k], i));
            const double x2 = 
                (tomo.clustering_zmax[i] - z_v[k] + bias_zphot_clustering(z_v[k], i)) /
                (M_SQRT2*sigma_zphot_clustering(z_v[k], i));
            
            gsl_sf_result result;
            {
              int status = gsl_sf_erf_e(x1, &result);
              if (status)
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
            }
            const double erf_sf_x1 = result.val;
            {
              int status = gsl_sf_erf_e(x2, &result);
              if (status)
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
            }
            const double erf_sf_x2 = result.val;

            const double DELTA = erf_sf_x2 - erf_sf_x1;
            const double ZDIST = pf_histo(z_v[k], NULL);

            if (zpart == 0) 
            {
              table[i+1][k] = 0.5*ZDIST*DELTA;
            }
            else if (zpart == 1) 
            {
              if (z_v[k] <= 0.5)
              {
                table[i+1][k] = 0.5*(1.0 - nuisance.frac_lowz)*ZDIST*DELTA;
              }
              else if(z_v[k] > 2.0)
              {
                table[i+1][k] = n_out/(3.5-2.0) + 0.5*(1.0 - nuisance.frac_highz)*ZDIST*DELTA;
              }
              else
              {
                table[i+1][k] = 0.5*ZDIST*DELTA;
              }
            }
            else if (zpart == 2) 
            {
              if (z_v[k] <= 0.5)
              {
                table[i+1][k] = n_out/(0.5) + 0.5*(1.-nuisance.frac_lowz)*ZDIST*DELTA;
              }
              else if (z_v[k] > 2.0)
              {
                table[i+1][k] = 0.5*(1.-nuisance.frac_highz)*ZDIST*DELTA;
              }
              else
              {
                table[i+1][k] = 0.5*ZDIST*DELTA;
              }
            }
            norm += table[i+1][k]*dz_histo;
          }
          for (int k = 0; k<zbins; k++)
          {
            table[i+1][k]/= norm;
          }
          break;
        }
        default:
        {
          log_fatal("%d not supported in this cosmolike version", photoz);
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
  
  if (photoz == 4) 
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
    int status = gsl_spline_eval_e(photoz_splines[nj + 1], zz, NULL, &res);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
  }
  return res;
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// cluster routines (for redshift distributions, including photo-zs * (optional))
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

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
    double NORM[11]; 
    #pragma omp parallel for
    for (int i=0; i<tomo.cluster_Nbin; i++)
    {
      double ar[1] = {(double) i};

      const double norm = int_gsl_integrate_medium_precision(pf_cluster_histo_n, (void*) ar, 1E-5, 
        tomo.cluster_zmax[i] + 1.0, NULL, GSL_WORKSPACE_SIZE) / 
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

  const int N_z = tomo.clustering_Nbin;
  const double zmin = tomo.cluster_zmin[nz];
  const double zmax = tomo.cluster_zmax[nz];
  
  if (table == 0) 
  {
    table = (double*) malloc(sizeof(double)*N_z);
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
      double params[1] = {(double) i};
      table[i] = int_gsl_integrate_high_precision(dV_cluster, (void*) params, zmin, zmax, NULL, 
        GSL_WORKSPACE_SIZE);
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

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// routines calculating the number of source and lens galaxies per bin
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double int_nsource(double z, void* param __attribute__((unused))) 
{
  return zdistr_photoz(z, -1);
}

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

    const int photoz = (redshift.shear_photoz == 7) ? 4 : redshift.shear_photoz;
    if (photoz == 4) 
    {
      for (int i=0; i<tomo.shear_Nbin; i++) 
      {
        table[i + 1] = tomo.n_source[i];
        
        if (tomo.n_source[i] < 0.01 || tomo.n_source[i] > 100.) 
        {
          log_fatal("nsource(%d) = %e\n outside expected range", i, tomo.n_source[i]);
          exit(0);
        }
      }
    } 
    else 
    {
      const double norm = like.high_def_integration > 0 ?
        int_gsl_integrate_high_precision(int_nsource, NULL, redshift.shear_zdistrpar_zmin,
          redshift.shear_zdistrpar_zmax, NULL, GSL_WORKSPACE_SIZE) :
        int_gsl_integrate_medium_precision(int_nsource, NULL, redshift.shear_zdistrpar_zmin,
          redshift.shear_zdistrpar_zmax, NULL, GSL_WORKSPACE_SIZE);

      #pragma omp parallel for
      for (int i=0; i<tomo.shear_Nbin; i++) 
      {        
        const double res = like.high_def_integration > 0 ?
          int_gsl_integrate_high_precision(int_nsource, NULL, 
            tomo.shear_zmin[i], tomo.shear_zmax[i], NULL, GSL_WORKSPACE_SIZE) :
          int_gsl_integrate_medium_precision(int_nsource, NULL, 
            tomo.shear_zmin[i], tomo.shear_zmax[i], NULL, GSL_WORKSPACE_SIZE);
        
        table[i + 1] = res / norm * survey.n_gal;
        
        if (table[i + 1] < 0.01 || table[i + 1] > 100) 
        {
          log_fatal("nsource(%d) = %e\n outside expected range", i, table[i + 1]);
          exit(1);
        }
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

double int_nlens(double z, void *param __attribute__((unused))) 
{
  return pf_photoz(z, -1);
}

// returns n_gal for clustering tomography bin ni, works only with binned distributions
double nlens(int ni) // ni =-1 -> no tomography; ni>= 0 -> tomography bin ni
{
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
    const int photoz = (redshift.clustering_photoz == 7) ? 4 : redshift.clustering_photoz;
    if (photoz == 4) 
    {      
      for (int i=0; i<tomo.clustering_Nbin; i++) 
      {
        table[i + 1] = tomo.n_lens[i];
      }
    } 
    else 
    {
      const double norm = like.high_def_integration > 0 ?
        int_gsl_integrate_high_precision(int_nlens, NULL, 
          tomo.clustering_zmin[0], tomo.clustering_zmax[tomo.clustering_Nbin - 1], NULL,
          GSL_WORKSPACE_SIZE) :
        int_gsl_integrate_medium_precision(int_nlens, NULL, 
          tomo.clustering_zmin[0], tomo.clustering_zmax[tomo.clustering_Nbin - 1], NULL,
          GSL_WORKSPACE_SIZE);

      #pragma omp parallel for 
      for (int i=0; i<tomo.clustering_Nbin; i++) 
      { // OMP: nlens depends on pf_photoz which is already initialized        
        const double res = like.high_def_integration > 0  ?
          int_gsl_integrate_high_precision(int_nlens, NULL, 
            tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE) :
          int_gsl_integrate_medium_precision(int_nlens, NULL, 
            tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE);

        table[i + 1] = res / norm * survey.n_lens;

        if (table[i + 1] < 0.001 || table[i + 1] > 100) 
        {
          log_fatal("nlens(%d) = %e\n outside expected range", i, table[i + 1]);
          exit(1);
        }
      }
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

double zmean(const int ni)
{ // mean true redshift of galaxies in tomography bin j
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
    #pragma omp parallel for
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    {
      double ar[1] = {(double) i};
      table[i] = 
        like.high_def_integration == 2 ?
        int_gsl_integrate_high_precision(int_for_zmean, (void*) ar, 
          tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE) /
        int_gsl_integrate_high_precision(norm_for_zmean, (void*) ar, tomo.clustering_zmin[i],
          tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE) :
        like.high_def_integration == 1 ?
        int_gsl_integrate_medium_precision(int_for_zmean, (void*) ar, 
          tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE) /
        int_gsl_integrate_medium_precision(norm_for_zmean, (void*) ar, tomo.clustering_zmin[i],
          tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE) :
        int_gsl_integrate_low_precision(int_for_zmean, (void*) ar, 
          tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE) /
        int_gsl_integrate_low_precision(norm_for_zmean, (void*) ar, tomo.clustering_zmin[i],
          tomo.clustering_zmax[i], NULL, GSL_WORKSPACE_SIZE);
    }
  }
  if (ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return table[ni];
}

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

double zmean_source(int ni) 
{ // mean true redshift of source galaxies in tomography bin j
  static double* table = 0;
  
  if (table == 0) 
  {
    table = (double*) malloc(sizeof(double)*tomo.shear_Nbin);
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = zdistr_photoz(0., 0);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for
    for (int i=0; i<tomo.shear_Nbin; i++) 
    {
      double ar[1] = {(double) i};
      table[i] = like.high_def_integration > 0 ? 
        int_gsl_integrate_high_precision(int_for_zmean_source, (void*) ar, 
          tomo.shear_zmin[i], tomo.shear_zmax[i], NULL, GSL_WORKSPACE_SIZE) :
        int_gsl_integrate_medium_precision(int_for_zmean_source, (void*) ar, 
          tomo.shear_zmin[i], tomo.shear_zmax[i], NULL, GSL_WORKSPACE_SIZE);
    }
  }
  if (ni < 0 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  } 
  return table[ni];
}

double max_g_tomo(int zs) 
{
  int i = 0;
  double max = 0;
  for (double z=0.; z<tomo.shear_zmax[zs]; z += tomo.shear_zmax[zs] / 50.) 
  {
    if (i != 0)
    {
      if (!fabs(z>0))
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
  static double** table = 0;
  
  const int N_z = tomo.shear_Nbin;
  const int N_a = Ntable.N_a; 
  const double amin = 1.0/(redshift.shear_zdistrpar_zmax + 1.0);
  const double amax = 0.999999;
  const double da = (amax - amin)/((double) N_a - 1.0);

  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*(N_z+1));
    for (int i=0; i<(N_z+1); i++)
    {
      table[i] = (double*) malloc(sizeof(double)*N_a);
    }
  }
  if (recompute_zphot_shear(N) || recompute_cosmo3D(C)) 
  {
    { // init static variables
      double ar[2];
      ar[1] = chi(amin);
      {
        const int j = -1;
        ar[0] = (double) j; 
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-variable"
        {
          double init = int_for_g_tomo(amin, (void*) ar);
        }
        #pragma GCC diagnostic pop
      } 
      if (tomo.shear_Nbin > 0) 
      { 
        const int j = 0;      
        ar[0] = (double) j;

        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-variable"
        {
          double init = int_for_g_tomo(amin, (void*) ar);
        }
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
        
        table[j+1][i] = 
          like.high_def_integration == 2 ?
          int_gsl_integrate_high_precision(int_for_g_tomo,
            (void*) ar, amin, a, NULL, GSL_WORKSPACE_SIZE) :          
          like.high_def_integration == 1 ? 
          int_gsl_integrate_medium_precision(int_for_g_tomo,
            (void*) ar, amin, a, NULL, GSL_WORKSPACE_SIZE) :
          int_gsl_integrate_low_precision(int_for_g_tomo,
            (void*) ar, amin, a, NULL, GSL_WORKSPACE_SIZE);
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
    res = interpol(table[ni + 1], Ntable.N_a, amin, amax, da, ainput, 1.0, 1.0); 
  }

  return res;
}

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

double g2_tomo(double a, int ni) 
{ // for tomography bin ni
  static nuisancepara N;
  static cosmopara C;
  static double** table = 0;

  const int N_z = tomo.shear_Nbin; 
  const int N_a = Ntable.N_a; 
  const double amin = 1.0/(redshift.shear_zdistrpar_zmax + 1.0);
  const double amax = 0.999999;
  const double da = (amax - amin)/((double) N_a - 1.0);

  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*(N_z+1));
    for (int i=0; i<(N_z+1); i++)
    {
      table[i] = (double*) malloc(sizeof(double)*N_a);
    }
  }
  if (recompute_zphot_shear(N) || recompute_cosmo3D(C)) 
  {
    { // init static variables
      double ar[2];
      ar[1] = chi(amin);
      {
        const int j = -1;
        ar[0] = (double) j; 
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-variable"
        {
          double init = int_for_g2_tomo(amin, (void*) ar);
        }
        #pragma GCC diagnostic pop
      } 
      if (tomo.shear_Nbin > 0) 
      { 
        const int j = 0;      
        ar[0] = (double) j;

        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-variable"
        {
          double init = int_for_g2_tomo(amin, (void*) ar);
        }
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

        table[j+1][i] = like.high_def_integration > 0 ? 
          int_gsl_integrate_high_precision(int_for_g2_tomo,
            (void*) ar, amin, a, NULL, GSL_WORKSPACE_SIZE) :
          int_gsl_integrate_medium_precision(int_for_g2_tomo,
            (void*) ar, amin, a, NULL, GSL_WORKSPACE_SIZE);
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
  
  const int N_z = tomo.clustering_Nbin; 
  const int N_a = Ntable.N_a;
  const double amin = 1.0/(redshift.clustering_zdistrpar_zmax + 1.0);
  const double amax = 0.999999;
  const double da = (amax - amin)/((double) N_a - 1.0);
  const double amin_shear = 1. / (redshift.shear_zdistrpar_zmax + 1.);

  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*(N_z+1));
    for (int i=0; i<(N_z+1); i++)
    {
      table[i] = (double*) malloc(sizeof(double)*N_a);
    }
  }
  if (recompute_zphot_clustering(N) || recompute_cosmo3D(C)) 
  {
    {
      double ar[2];
      ar[1] = chi(amin);
      { // init static variables
        const int j = -1;
        ar[0] = (double) j; // j = -1, no tomography is being done        
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-variable"
        {
          double init = int_for_g_lens(amin_shear, (void*) ar);
        }
        #pragma GCC diagnostic pop
        if (tomo.clustering_Nbin > 0) 
        {
          const int j = 0;
          ar[0] = (double) j;
          #pragma GCC diagnostic push
          #pragma GCC diagnostic ignored "-Wunused-variable"
          {
            double init = int_for_g_lens(amin_shear, (void*) ar);
          }
          #pragma GCC diagnostic pop 
        }
      }
    }
    #pragma omp parallel for collapse(2)
    for (int j=-1; j<tomo.clustering_Nbin; j++) 
    {
      for (int i=0; i<Ntable.N_a; i++) 
      {
        double ar[2];
        ar[0] = (double) j; // if j=-1, no tomography is being done

        const double a =  amin + i*da;
        ar[1] = chi(a);
        
        table[j + 1][i] = like.high_def_integration == 2 ?
          int_gsl_integrate_high_precision(int_for_g_lens,
            (void*) ar, amin_shear, a, NULL, GSL_WORKSPACE_SIZE) : 
          like.high_def_integration == 1 ?
          int_gsl_integrate_medium_precision(int_for_g_lens,
            (void*) ar, amin_shear, a, NULL, GSL_WORKSPACE_SIZE) :
          int_gsl_integrate_low_precision(int_for_g_lens,
            (void*) ar, amin_shear, a, NULL, GSL_WORKSPACE_SIZE);
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

  const int N_z = tomo.cluster_Nbin;
  const int N_l = Cluster.N200_Nbin;
  if (nl < 0 || nl > N_l - 1)
  {
    log_fatal("invalid bin input ni = %d (max %d)", nl, N_l);
    exit(1);
  } 
  const int na =  Ntable.N_a; 
  const double amin = 1.0/(tomo.cluster_zmax[tomo.cluster_Nbin - 1] + 1.);
  const double amax = 0.999999;
  const double da = (amax - amin)/((double) na - 1.0);
  
  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*(tomo.cluster_Nbin + 1));
    for (int i=0; i<tomo.cluster_Nbin+1; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*na);
    }
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
      for (int i=0; i<na; i++) 
      {
        const double aa = amin + i*da;
        double ar[3];
        ar[0] = (double) j; 
        ar[1] = aa;
        ar[2] = nl;
        table[j + 1][i] = int_gsl_integrate_medium_precision(int_for_g_lens_cl, (void*) ar, 
          1.0/(redshift.shear_zdistrpar_zmax + 1.0), aa, NULL, 4000); 
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
  return interpol(table[nz + 1], na, amin, amax, da, a, 1.0, 1.0); 
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// routines for association of a pair redshift bin numbers and power spectrum tomography bin
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double int_for_ggl_efficiency(double z, void* params) 
{
  if (!fabs(z>0)) 
  {
    log_fatal("z = 0 with 1/z needed");
    exit(1);
  }
  
  double* ar = (double*) params;
  const int ni = (int) ar[0];
  const int nj = (int) ar[1];
  if (ni < -1 || ni > tomo.clustering_Nbin - 1 || nj < -1 || nj > tomo.shear_Nbin - 1)
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
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = pf_photoz(0, 0);
      double init2 = zdistr_photoz(0, 0);
    }
    #pragma GCC diagnostic pop
    
    table = (double**) malloc(sizeof(double*)*(tomo.clustering_Nbin+1));
    for (int i=0; i<(tomo.clustering_Nbin+1); i++)
    {
      table[i] = (double*) malloc(sizeof(double)*(tomo.shear_Nbin+1));
    }
    {
      const int i = 0;
      const int j = 0;
      double ar[2];
      ar[0] = (double) i;
      ar[1] = (double) j;
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-variable"
      double init = int_for_ggl_efficiency(tomo.clustering_zmin[i], (void*) ar);
      double init2 = max_g_tomo(j);
      #pragma GCC diagnostic pop
    }
    #pragma omp parallel for collapse (2)
    for (int i=0; i<tomo.clustering_Nbin; i++) 
    {
      for (int j=0; j<tomo.shear_Nbin; j++) 
      {
        double ar[2];
        ar[0] = (double) i;
        ar[1] = (double) j;
        table[i][j] = like.high_def_integration > 0 ?
          int_gsl_integrate_high_precision(int_for_ggl_efficiency, (void*) ar,
            tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, 
            GSL_WORKSPACE_SIZE)/max_g_tomo(j) :
          int_gsl_integrate_medium_precision(int_for_ggl_efficiency, (void*) ar,
            tomo.clustering_zmin[i], tomo.clustering_zmax[i], NULL, 
            GSL_WORKSPACE_SIZE)/max_g_tomo(j);
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
