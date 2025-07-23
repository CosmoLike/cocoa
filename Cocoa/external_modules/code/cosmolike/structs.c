#include <stdio.h>
#include "basics.h"
#include "structs.h"

#include "log.c/src/log.h"

FPT FPTIA;
FPT FPTbias;
cosmopara cosmology;
tomopara tomo;
redshiftparams redshift;

nuisanceparams nuisance;
likepara like;
sur survey;
//clusterparams Cluster;

pdeltapara pdeltaparams =
{
  .runmode = "Halofit"
};

Cmb cmb =
{
  .fwhm = 0.0,
  .sensitivity = 0.0
};

lim limits = 
{
  .a_min = 1.0/(1.0 + 40.0),    // a_min (z = 40, needed for CMB lensing)
  .k_min_cH0 = 2.e-2,           // k_min_cH0
  .k_max_cH0 = 3.e+6,           // k_max_cH0
  .M_min = 1.0e+6,              // M_min
  .M_max = 1.0e+17,             // M_max
  .LMIN_tab = 20,               // LMIN_tab
  .LMAX_NOLIMBER = 150,         // LMAX_NOLIMBER
  .cluster_util_log_M_min = 12.0,
  .cluster_util_log_M_max = 15.9,
  .binned_P_lambda_obs_given_M_zmin_table = 0.20,
  .binned_P_lambda_obs_given_M_zmax_table = 0.80,
  .SDSS_P_lambda_obs_given_M_true_lambda_min = 3.0,
  .SDSS_P_lambda_obs_given_M_true_lambda_max = 160.0,
  .halo_exclusion_k_min = 1E-2,
  .halo_exclusion_k_max = 3E6,
  .halo_exclusion_k_min_hankel = 5.0E-4,
  .halo_exclusion_k_max_hankel = 1.0E8,
  .halo_exclusion_R_min = 0.0,
  .halo_exclusion_R_max = 15.0/2997.,
  .halo_uKS_cmin = 0.1,     // halo.c u_KS(double c, double k, double rv)
  .halo_uKS_cmax = 50.0,    // halo.c u_KS(double c, double k, double rv)
  .halo_uKS_xmin = 1e-10,   // halo.c u_KS(double c, double k, double rv)
  .halo_uKS_xmax = 5e3      // halo.c u_KS(double c, double k, double rv)
};

Ntab Ntable;

//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------
//  RESET STRUCT 
//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------

void reset_like_struct()
{
  like.bias = 0;
  like.Ncl = 0;
  like.Ncos = 0;
  like.Ndata = 0;
  like.lmin = 0;
  like.lmax = 0;
  if((like.ell != NULL) == 1)
  {
    free(like.ell);
    like.ell = NULL;
  }
  like.cosmax = 0;
  like.Rmin_bias = 0;
  like.Rmin_shear = 0;
  like.lmax_shear = 0;
  like.lmin_kappacmb = 0;
  like.lmax_kappacmb = 0;
  
  like.bias = 0;
  like.clusterN = 0;
  like.clusterWL = 0;
  like.clusterCG = 0;
  like.clusterCC = 0;
  like.shear_shear = 0;
  like.shear_pos = 0;
  like.pos_pos = 0;
  like.gk = 0;
  like.kk = 0;
  like.ks = 0;
  like.gy = 0;
  like.sy = 0;
  like.ky = 0;
  like.yy = 0;
  like.adopt_limber_gg = 0;
}

void reset_cosmology_struct()
{
  cosmology.Omega_nu = 0.;
  cosmology.coverH0 = 2997.92458;
  cosmology.rho_crit = 7.4775e+21;
  cosmology.MGSigma = 0.0;
  cosmology.MGmu = 0.0;
  cosmology.Omega_b = 0.0;
  cosmology.Omega_m = 0.0;
  cosmology.Omega_v = 0.0;
  cosmology.h0 = 0.0;
  cosmology.Omega_nu = 0.0;
  cosmology.random = 0.0;
  cosmology.sigma_8 = 0.0;
  cosmology.lnP_nk = 0;
  cosmology.lnP_nz = 0;
  cosmology.lnP = NULL;
  cosmology.lnPL_nk = 0;
  cosmology.lnPL_nz = 0;
  cosmology.lnPL = NULL;
  cosmology.chi_nz = 0;
  cosmology.chi = NULL;
  cosmology.G_nz = 0;
  cosmology.G = NULL;
}

void reset_tomo_struct()
{
  tomo.shear_Npowerspectra = 0;
  tomo.clustering_Npowerspectra = 0;
  tomo.ggl_Npowerspectra = 0;
  if (tomo.ggl_exclude != NULL) {
    free(tomo.ggl_exclude);
    tomo.ggl_exclude = NULL;
  }
  tomo.N_ggl_exclude = 0;
}

void reset_redshift_struct()
{
  redshift.random_shear = 0.0;
  redshift.random_clustering = 0.0;

  redshift.shear_nbin = 0;
  redshift.shear_photoz = 0;
  if (redshift.shear_zdist_table != NULL) {
    free(redshift.shear_zdist_table);
    redshift.shear_zdist_table = NULL;
  }
  redshift.shear_nzbins = 0;
  redshift.shear_zdist_zmin_all = 0.0;
  redshift.shear_zdist_zmax_all = 0.0;

  redshift.clustering_nbin = 0;
  redshift.clustering_nzbins = 0;
  if (redshift.clustering_zdist_table != NULL) {
    free(redshift.clustering_zdist_table);
    redshift.clustering_zdist_table = NULL;
  }
  redshift.clustering_photoz = 0;
  redshift.clustering_zdist_zmin_all = 0.0;
  redshift.clustering_zdist_zmax_all = 0.0;

  for (int i=0; i<MAX_SIZE_ARRAYS; i++) {
    redshift.shear_zdist_zmin[i] = 0.0;
    redshift.shear_zdist_zmax[i] = 0.0; 
    redshift.clustering_zdist_zmin[i] = 0.0;
    redshift.clustering_zdist_zmax[i] = 0.0;
    redshift.clustering_zdist_zmean[i] = 0.0;
  }
}

void reset_survey_struct()
{
  survey.area = 0.0;
  survey.n_gal = 0.0;
  survey.sigma_e = 0.0;
  survey.area_conversion_factor =
    60.0 * 60.0 * 2.90888208665721580e-4 * 2.90888208665721580e-4;
  survey.n_gal_conversion_factor =
    1.0 / 2.90888208665721580e-4 / 2.90888208665721580e-4;
  survey.n_lens = 0.0;
  survey.m_lim = 0.0;
  sprintf(survey.name, "%s", "");
}

void reset_pdeltaparams_struct()
{
  sprintf(pdeltaparams.runmode, "%s", "Halofit");
}

void reset_nuisance_struct()
{
  nuisance.random_ia = 0.0;
  nuisance.random_photoz_shear = 0.0;
  nuisance.random_photoz_clustering = 0.0;

  for (int i=0; i<MAX_SIZE_ARRAYS; i++) {
    nuisance.shear_calibration_m[i] = 0.0;
    nuisance.gc[i] = 0.0;
    nuisance.gas[i] = 0.0;
    for (int j=0; j<MAX_SIZE_ARRAYS; j++) {
      nuisance.ia[i][j] = 0.0;
      nuisance.ia[i][j] = 0.0;
      nuisance.ia[i][j] = 0.0;
      nuisance.gb[i][j] = 0.0;
      nuisance.hod[i][j] = 0.0;
      for (int k=0; k<MAX_SIZE_ARRAYS; k++) {
        if (j==1) {
          nuisance.photoz[i][j][k] = 1.0;
        } // photo-z stretch params
        else {
          nuisance.photoz[i][j][k] = 0.0;
        }
      }
    }
  }
  nuisance.oneplusz0_ia = 0.0;
  nuisance.c1rhocrit_ia = 0.01389;
  nuisance.IA = 0;
  nuisance.IA_MODEL = 0;
}

void reset_cmb_struct()
{
  sprintf(cmb.name, "%s", "");
  cmb.fwhm = 0.0;
  cmb.sensitivity = 0.0;
  sprintf(cmb.pathLensRecNoise, "%s", "");
}

void reset_Ntable_struct()
{
  Ntable.LMAX     = 50000;
  Ntable.random   = 0.0;
  Ntable.N_a      = 200;                   // N_a       
  Ntable.N_k_lin  = 500;                   // N_k_lin
  Ntable.N_k_nlin = 500;                   // N_k_nlin
  Ntable.N_ell    = 250;                   // N_ell      
  Ntable.Ntheta   = 225;                   // N_theta  
  Ntable.N_M      = 1000;                  // N_M, M = mass (Halo Model)
  Ntable.NL_Nchi  = 500;                   // Cosmo2D - NL = NonLimber (NL_Nchi)
  Ntable.photoz_interpolation_type = 0;
  Ntable.high_def_integration = 0;
  Ntable.FPTboost=0;
}

/*
void reset_cluster_struct()
{
  Cluster.N200_min = 0.0;
  Cluster.N200_max = 0.0;
  Cluster.N200_Nbin = 0;
  for(int i=0; i<MAX_SIZE_ARRAYS; i++)
  {
    Cluster.N_min[i] = 0.0;
    Cluster.N_max[i] = 0.0;
  }
  sprintf(Cluster.model, "%s", "default");
}
*/

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------