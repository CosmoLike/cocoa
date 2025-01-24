#include <stdio.h>
#include "basics.h"
#include "structs.h"

#include "log.c/src/log.h"

FPT FPTIA =
{
  .k_min = 1.e-5,
  .k_max = 1.e+3,
  .N = 800,
  .tab = NULL
};

FPT FPTbias =
{
  .k_min = 1.e-5,
  .k_max = 1.e+3,
  .N = 800,
  .tab = NULL
};

likepara like =
{
  .IA = 0.,
  .bias = 0,
  .ell = NULL,
  .theta = NULL,
  .shear_shear = 0,
  .shear_pos = 0,
  .pos_pos = 0,
  .gk = 0,
  .kk = 0,
  .ks = 0,
  .gy = 0,
  .sy = 0,
  .ky = 0,
  .yy = 0,
  .clusterN = 0,
  .clusterWL = 0,
  .clusterCG = 0,
  .clusterCC = 0,
  .adopt_limber_gg = 0,
  .adopt_limber_gammat = 1,
  .use_ggl_efficiency_zoverlap = 1
};

cosmopara cosmology =
{
  .random = 0.0,
  .Omega_nu = 0.,
  .coverH0 = 2997.92458,
  .rho_crit = 7.4775e+21,
  .MGSigma = 0.0,
  .MGmu = 0.0,
  .sigma_8 = 0.0,
  .Omega_b = 0.0,
  .Omega_v = 0.0,
  .h0 = 0.0,
  .lnP_nk = 0,
  .lnP_nz = 0,
  .lnP = NULL,
  .lnPL_nk = 0,
  .lnPL_nz = 0,
  .lnPL = NULL,
  .chi_nz = 0,
  .chi = NULL,
  .G_nz = 0,
  .G = NULL
};

tomopara tomo = 
{
  .shear_Npowerspectra = 0,
  .ggl_Npowerspectra = 0,
  .clustering_Npowerspectra = 0
};

redshiftparams redshift =
{
  .random_shear = 0.0,
  .random_clustering = 0.0,
  .shear_nbin = 0,
  .shear_photoz = 0,
  .shear_nzbins = 0,
  .shear_zdist_zmin_all = 0.0,
  .shear_zdist_zmax_all = 0.0,
  .shear_zdist_zmin = {0},
  .shear_zdist_zmax = {0},
  .shear_zdist_table = NULL,
  .clustering_nbin = 0,
  .clustering_photoz = 0,
  .clustering_nzbins = 0,
  .clustering_zdist_zmin_all = 0.0,
  .clustering_zdist_zmax_all = 0.0,
  .clustering_zdist_zmin = {0},
  .clustering_zdist_zmax = {0},
  .clustering_zdist_table = NULL
};

clusterparams Cluster =
{
  .interpolate_survey_area = 1,
  .bias_model = 0,
  .hmf_model = 0,
  .nonlinear_bias = 0,
  .N200_min = 0.0,
  .N200_max = 0.0,
  .N200_Nbin = 0,
  .N_min = {0},
  .N_max = {0},
  .model = "default"
};

pdeltapara pdeltaparams =
{
  .runmode = "Halofit"
};

nuisanceparams nuisance =
{
  .random_photoz_shear = 0.0,
  .random_photoz_clustering = 0.0,
  .c1rhocrit_ia = 0.01389,
  .ia = {0},
  .shear_calibration_m = {0},
  .gb = {0},
  .gas = {0},
  .gc = {0},
  .hod = {0}
};

Cmb cmb =
{
  .fwhm = 0.0,
  .sensitivity = 0.0
};

sur survey =
{
  .area_conversion_factor  = 60.0 * 60.0 * 2.90888208665721580e-4 * 2.90888208665721580e-4,
  .n_gal_conversion_factor = 1.0 / 2.90888208665721580e-4 / 2.90888208665721580e-4
};

lim limits = 
{
  .a_min = 1.0/(1.0 + 10.0),    // a_min (z = 10, needed for CMB lensing)
  .k_min_cH0 = 2.e-2,           // k_min_cH0
  .k_max_cH0 = 3.e+6,           // k_max_cH0
  .M_min = 1.0e+6,              // M_min
  .M_max = 1.0e+17,             // M_max
  .LMIN_tab = 20,               // LMIN_tab
  .LMAX_NOLIMBER = 250,         // LMAX_NOLIMBER
  .LMAX = 100000,
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

Ntab Ntable = 
{
  .random = 0.0,
  .N_a = 350,                         // N_a        (modified by COCOA from 100)
  .N_k_lin = 500,                     // N_k_lin
  .N_k_nlin = 500,                    // N_k_nlin
  .N_ell  = 350,                      // N_ell      (modified by COCOA from 200)
  .Ntheta = 250,                      // N_theta    (modified by COCOA from 200)
  .N_M  = 1000,                       // N_M, M = mass (Halo Model)
  .NL_Nell_block = 50,                // Cosmo2D - NL = NonLimber (NL_Nell_block)
  .NL_Nchi = 500,                     // Cosmo2D - NL = NonLimber (NL_Nchi)
  .photoz_interpolation_type = 0,
  .high_def_integration = 0,
  .FPTboost=0,
  // ---------------------------------------------------
  // CLUSTER ROUTINES (ALPHA STAGE)
  // ---------------------------------------------------
  .N_a_halo_exclusion = 100,          // N_a for binned_p_cc_incl_halo_exclusion (cluster_util.c)
  .N_k_halo_exclusion = 100,          // N_k for binned_p_cc_incl_halo_exclusion (cluster_util.c)
  .N_k_hankel_halo_exclusion = 3192,  // N for 3D Hankel Transform (pk_to_xi and xi_to_pk) 
  .N_R_halo_exclusion = 64,
  .binned_P_lambda_obs_given_M_size_z_table = 10,
  .binned_P_lambda_obs_given_M_size_M_table = 50,
  .binned_p_cm_size_a_table = 30,
  .halo_uKS_nc = 20,
  .halo_uks_nx = 200,
};

//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------
//  RESET STRUCT 
//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------

void reset_like_struct()
{
  like.IA = 0.;
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
  if((like.theta != NULL) == 1)
  {
    free(like.theta);
    like.theta = NULL;
  }
  like.cosmax = 0;
  like.Rmin_bias = 0;
  like.Rmin_shear = 0;
  like.lmax_shear = 0;
  like.lmin_kappacmb = 0;
  like.lmax_kappacmb = 0;
  like.IA = 0;
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
}

void reset_tomo_struct()
{
  tomo.shear_Npowerspectra = 0;
  tomo.clustering_Npowerspectra = 0;
  tomo.ggl_Npowerspectra = 0; 
/*
  tomo.cg_clustering_Npowerspectra = 0;
  tomo.cc_clustering_Npowerspectra = 0;
  tomo.cgl_Npowerspectra = 0;
*/
}

void reset_redshift_struct()
{
  redshift.random_shear = 0.0;
  redshift.random_clustering = 0.0;

  redshift.shear_nbin = 0;
  redshift.shear_photoz = 0;
  if (redshift.shear_zdist_table != NULL)
  {
    free(redshift.shear_zdist_table);
    redshift.shear_zdist_table = NULL;
  }
  redshift.shear_nzbins = 0;
  redshift.shear_zdist_zmin_all = 0.0;
  redshift.shear_zdist_zmax_all = 0.0;

  redshift.clustering_nbin = 0;
  redshift.clustering_nzbins = 0;
  if (redshift.clustering_zdist_table != NULL)
  {
    free(redshift.clustering_zdist_table);
    redshift.clustering_zdist_table = NULL;
  }
  redshift.clustering_photoz = 0;
  redshift.clustering_zdist_zmin_all = 0.0;
  redshift.clustering_zdist_zmax_all = 0.0;

  for(int i=0; i<MAX_SIZE_ARRAYS; i++)
  {
    redshift.shear_zdist_zmin[i] = 0.0;
    redshift.shear_zdist_zmax[i] = 0.0; 
    redshift.clustering_zdist_zmin[i] = 0.0;
    redshift.clustering_zdist_zmax[i] = 0.0;
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

void reset_pdeltaparams_struct()
{
  sprintf(pdeltaparams.runmode, "%s", "Halofit");
}

void reset_nuisance_struct()
{
  nuisance.random_ia = 0.0;
  nuisance.random_photoz_shear = 0.0;
  nuisance.random_photoz_clustering = 0.0;

  for(int i=0; i<MAX_SIZE_ARRAYS; i++)
  {
    nuisance.shear_calibration_m[i] = 0.0;
    nuisance.gc[i] = 0.0;
    nuisance.gas[i] = 0.0;

    for(int j=0; j<MAX_SIZE_ARRAYS; j++)
    {
      nuisance.ia[i][j] = 0.0;
      nuisance.ia[i][j] = 0.0;
      nuisance.ia[i][j] = 0.0;

      nuisance.gb[i][j] = 0.0;
      nuisance.hod[i][j] = 0.0;
      
      for(int k=0; k<MAX_SIZE_ARRAYS; k++)
      {
        nuisance.photoz[i][j][k] = 0.0;
      }
    }
        
    /*
    nuisance.cluster_MOR[i] = 0.0;
    nuisance.cluster_selection[i] = 0.0;
    */
  }
  nuisance.oneplusz0_ia = 0.0;
  nuisance.c1rhocrit_ia = 0.01389;
}

void reset_cmb_struct()
{
  sprintf(cmb.name, "%s", "");
  cmb.fwhm = 0.0;
  cmb.sensitivity = 0.0;
  sprintf(cmb.pathLensRecNoise, "%s", "");
}

// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// UPDATE
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

void update_cosmopara(cosmopara *C) 
{
  C->Omega_b = cosmology.Omega_b;
  C->Omega_m = cosmology.Omega_m;
  C->Omega_v = cosmology.Omega_v;
  C->Omega_nu = cosmology.Omega_nu;
  C->h0 = cosmology.h0;
  C->MGSigma = cosmology.MGSigma;
  C->MGmu = cosmology.MGmu;
  C->random = cosmology.random;
  //C->w0 = cosmology.w0;    // cluster cosmology need that information
  //C->A_s = cosmology.A_s;  // cluster cosmology need that information
  //C->n_s = cosmology.n_s;  // cluster cosmology need that information
  C->random = cosmology.random;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------