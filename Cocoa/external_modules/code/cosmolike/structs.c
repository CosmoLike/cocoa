#include <stdio.h>
#include "basics.h"
#include "structs.h"

#include "log.c/src/log.h"

//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------
//  STRUCT INITIALIZATION  
//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------

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
  .use_full_sky_shear = 1,
  .use_full_sky_ggl = 1,
  .use_full_sky_clustering = 1,
  .use_full_sky_gk = 1,
  .use_full_sky_sk = 1,
  .adopt_limber_gg = 0,
  .adopt_limber_gammat = 1,
  .use_ggl_efficiency_zoverlap = 1
};

cosmopara cosmology =
{
  .Omega_nu = 0.,
  .coverH0 = 2997.92458,
  .rho_crit = 7.4775e+21,
  .MGSigma = 0.0,
  .MGmu = 0.0,
  .sigma_8 = 0.0,
  .Omega_b = 0.0,
  .Omega_v = 0.0,
  .h0 = 0.0,
  //.w0 = -1,
  //.A_s =  2e-9,
  //.n_s = 0.96,
};


tomopara tomo = 
{
  .shear_Npowerspectra = 0,
  .ggl_Npowerspectra = 0,
  .clustering_Npowerspectra = 0
};

redshiftparams redshift =
{
  .random = 0.0,

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

  //.cluster_Nbin = 0
};

galpara gbias =
{
  .b2 = {0},
  .bs2 = {0},
  .b1_function = &b1_per_bin,
  .b_mag = {0},
  .cg = {0},
  .hod = {0}
}; // default: point to old bgal_z routin

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

nuisancepara nuisance =
{
  .c1rhocrit_ia = 0.01389,
  .A_z = {0},
  .A2_z = {0},
  .b_ta_z = {0},
  .shear_calibration_m = {0},
  .sigma_zphot_shear = {0},
  .bias_zphot_shear = {0},
  .sigma_zphot_clustering = {0},
  .bias_zphot_clustering = {0},
};

ynuisancepara ynuisance =
{
  .gas_Gamma_KS = 0,
  .gas_beta = 0,
  .gas_lgM0 = 0, 
  .gas_eps1 = 0,
  .gas_eps2 = 0,
  .gas_alpha = 0,
  .gas_A_star = 0,
  .gas_lgM_star = 0,
  .gas_sigma_star = 0,
  .gas_lgT_w = 0,
  .gas_f_H = 0
};

barypara bary =
{
  .is_Pk_bary = 0,
  .Na_bins = 0,
  .Nk_bins = 0,
  .a_bins = NULL,
  .logk_bins = NULL,
  .log_PkR = NULL,
  .interp2d = NULL,
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

FPTpara FPT =
{
  .k_min = 1.e-5,
  .k_max = 1.e+3,
  .N = 800,
  .N_per_dec = 100,
  .N_AB = 7,
  .N_IA = 10
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
  .LMIN_hankel = 0.0001,
  .LMAX_hankel = 5.0e6,
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
  .N_ell  = 325,                      // N_ell      (modified by COCOA from 200)
  .Ntheta = 250,                      // N_theta    (modified by COCOA from 200)
  .N_thetaH = 2048,                   // N_theta for Hankel
  .N_M  = 1000,                       // N_M, M = mass (Halo Model)
  .N_ell_TATT = 200,                  // N_ell_TATT (modified by COCOA from 60)
  .NL_Nell_block = 50,                // Cosmo2D - NL = NonLimber (NL_Nell_block)
  .NL_Nchi = 500,                     // Cosmo2D - NL = NonLimber (NL_Nchi)
  .N_a_halo_exclusion = 100,          // N_a for binned_p_cc_incl_halo_exclusion (cluster_util.c)
  .N_k_halo_exclusion = 100,          // N_k for binned_p_cc_incl_halo_exclusion (cluster_util.c)
  .N_k_hankel_halo_exclusion = 3192,  // N for 3D Hankel Transform (pk_to_xi and xi_to_pk) 
  .N_R_halo_exclusion = 64,
  .binned_P_lambda_obs_given_M_size_z_table = 10,
  .binned_P_lambda_obs_given_M_size_M_table = 50,
  .binned_p_cm_size_a_table = 30,
  .halo_uKS_nc = 20,
  .halo_uks_nx = 200,
  .photoz_interpolation_type = 0,
  .high_def_integration = 0,
};

pre precision = 
{
  .medium = 5e-4, 
  .high   = 1e-5, 
  .insane = 1e-7    
};

TinkerEmuParameters tinkerEmuParam =
{
  .tinker_bias_ncosmo = 7,                // do not change
  .tinker_bias_nparam = 4,                // do not change
  .tinker_bias_nsamp = 40, 
  .tinker_bias_nparam_redshift = 6,       // number of parameters of equation 7 in 1907.13167.;
  .tinker_hmf_ncosmo = 7,                 // do not change
  .tinker_hmf_nparam = 6,                 // do not change
  .tinker_hmf_nsamp = 40, 
  .tinker_hmf_nparam_redshift = 4,        // number of parameters of eq 2 in 1804.05866;
  .tinker_bias_extrapolation_cut_in = 40,
  .tinker_bias_extrapolation_cut_out = 45,
  .tinker_hmf_extrapolation_cut_in = 40,
  .tinker_hmf_extrapolation_cut_out = 45,
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
  like.use_full_sky_shear = 1;
  like.use_full_sky_ggl = 1;
  like.use_full_sky_clustering = 1;
  like.use_full_sky_gk = 1;
  like.use_full_sky_sk = 1;
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
  redshift.random = 0.0;

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

void reset_gbias_struct()
{
  for(int i=0; i<10; i++)
  {
    gbias.b[i] = 0.0;
    gbias.b2[i] = 0.0;
    gbias.bs2[i] = 0.0;
    gbias.b_mag[i] = 0.0;
  }
  gbias.b1_function = &b1_per_bin;
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

void reset_bary_struct()
{
  bary.is_Pk_bary = 0;
  bary.T = NULL;
  bary.Na_bins = 0;
  bary.Nk_bins = 0;

  if((bary.a_bins != NULL))
  {
    free(bary.a_bins);
    bary.a_bins = NULL;
  }

  if((bary.logk_bins != NULL))
  {
    free(bary.logk_bins);
    bary.logk_bins = NULL;
  }
  
  if((bary.log_PkR != NULL))
  {
    free(bary.log_PkR);
    bary.log_PkR = NULL;
  }
  
  if((bary.interp2d != NULL))
  {
    gsl_interp2d_free(bary.interp2d);
    bary.interp2d = NULL;
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void reset_pdeltaparams_struct()
{
  sprintf(pdeltaparams.runmode, "%s", "Halofit");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void reset_nuisance_struct()
{
  for(int i=0; i<MAX_SIZE_ARRAYS; i++)
  {
    nuisance.A_z[i] = 0.0;
    nuisance.A2_z[i] = 0.0;
    nuisance.b_ta_z[i] = 0.0;
    nuisance.fred[i] = 0.0;
    
    nuisance.shear_calibration_m[i] = 0.0;
    
    nuisance.sigma_zphot_shear[i] = 0.0;
    nuisance.bias_zphot_shear[i] = 0.0;
    
    nuisance.sigma_zphot_clustering[i] = 0.0;
    nuisance.bias_zphot_clustering[i] = 0.0;
    
    nuisance.cluster_MOR[i] = 0.0;
    nuisance.cluster_selection[i] = 0.0;
  }

  nuisance.A_ia = 0.0;
  nuisance.A2_ia = 0.0;
  nuisance.beta_ia = 0.0;
  nuisance.eta_ia = 0.0;
  nuisance.eta_ia_tt = 0.0;
  nuisance.eta_ia_highz = 0.0;
  nuisance.oneplusz0_ia = 0.0;
  nuisance.c1rhocrit_ia = 0.01389;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void reset_ynuisance_struct()
{ // Compton-Y related variables  
  ynuisance.gas_Gamma_KS =  0.0;
  ynuisance.gas_beta = 0.0;     
  ynuisance.gas_lgM0 = 0.0;     
  ynuisance.gas_eps1 = 0.0;
  ynuisance.gas_eps2 = 0.0;
  ynuisance.gas_alpha = 0.0;
  ynuisance.gas_A_star = 0.0;
  ynuisance.gas_lgM_star = 0.0;
  ynuisance.gas_sigma_star = 0.0;
  ynuisance.gas_lgT_w = 0.0;
  ynuisance.gas_f_H = 0.0;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

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

void update_galpara(galpara *G) 
{
  for (int i=0; i<MAX_SIZE_ARRAYS; i++) 
  {
    G->b[i] = gbias.b[i];
    G->b2[i] = gbias.b2[i];
    G->bs2[i] = gbias.bs2[i];  
    for(int j=0; j<MAX_SIZE_ARRAYS; j++)
    {
      G->hod[i][j] = gbias.hod[i][j];
    }
    G->cg[i]= gbias.cg[i];
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void update_nuisance(nuisancepara* N)
{
  N->A_ia          = nuisance.A_ia;
  N->beta_ia       = nuisance.beta_ia;
  N->eta_ia        = nuisance.eta_ia;
  N->eta_ia_highz  = nuisance.eta_ia_highz;
  N->A2_ia         = nuisance.A2_ia;
  N->eta_ia_tt     = nuisance.eta_ia_tt;
  
  for (int i=0; i<MAX_SIZE_ARRAYS; i++) 
  {
    N->sigma_zphot_clustering[i] = nuisance.sigma_zphot_clustering[i];
    N->bias_zphot_clustering[i]  = nuisance.bias_zphot_clustering[i];
  }
  
  for (int i=0; i<MAX_SIZE_ARRAYS; i++) 
  {
    N->sigma_zphot_shear[i] = nuisance.sigma_zphot_shear[i];
    N->bias_zphot_shear[i]  = nuisance.bias_zphot_shear[i];
    
    N->A_z[i]    = nuisance.A_z[i];
    N->A2_z[i]   = nuisance.A2_z[i];
    N->b_ta_z[i] = nuisance.b_ta_z[i];

    N->cluster_MOR[i] = nuisance.cluster_MOR[i];
    N->cluster_selection[i] = nuisance.cluster_selection[i];
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void update_ynuisance(ynuisancepara* N)
{ // Compton-Y related variables
  N->gas_Gamma_KS = ynuisance.gas_Gamma_KS;
  N->gas_beta = ynuisance.gas_beta;     
  N->gas_lgM0 = ynuisance.gas_lgM0;     
  N->gas_eps1 = ynuisance.gas_eps1;
  N->gas_eps2 = ynuisance.gas_eps2;
  N->gas_alpha = ynuisance.gas_alpha;
  N->gas_A_star = ynuisance.gas_A_star;
  N->gas_lgM_star = ynuisance.gas_lgM_star;
  N->gas_sigma_star = ynuisance.gas_sigma_star;
  N->gas_lgT_w = ynuisance.gas_lgT_w;
  N->gas_f_H = ynuisance.gas_f_H;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void update_table(Ntab* N)
{
  N->random = Ntable.random;

  N->N_a = Ntable.N_a;
  N->N_ell = Ntable.N_ell;
  N->Ntheta = Ntable.Ntheta;
  N->N_ell_TATT = Ntable.N_ell_TATT;
  N->high_def_integration = Ntable.high_def_integration;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
void update_redshift(redshiftparams* N)
{
  N->random = redshift.random;

  if (redshift.clustering_zdist_table != NULL)
  {
    N->clustering_nbin = redshift.clustering_nbin;
    N->clustering_nzbins = redshift.clustering_nzbins;
    N->clustering_zdist_zmin_all = redshift.clustering_zdist_zmin_all;
    N->clustering_zdist_zmax_all = redshift.clustering_zdist_zmax_all;

    for(int i=0; i<MAX_SIZE_ARRAYS; i++)
    {
      N->clustering_zdist_zmin[i] = redshift.clustering_zdist_zmin[i];
      N->clustering_zdist_zmax[i] = redshift.clustering_zdist_zmax[i];
    }

    const int ntomo  = redshift.clustering_nbin;
    const int nzbins = redshift.clustering_nzbins;
    
    if (N->clustering_zdist_table != NULL) free(N->clustering_zdist_table);
    N->clustering_zdist_table = (double**) malloc2d(ntomo + 1, nzbins);

    for (int i=0; i<ntomo+1; i++)
    {
      for (int k=0; k<nzbins; k++) 
      {
        N->clustering_zdist_table[i,k] = redshift.clustering_zdist_table[i,k];
      }
    }
  }

  if (redshift.shear_zdist_table != NULL)
  {
    N->shear_nbin = redshift.shear_nbin;
    N->shear_nzbins = redshift.shear_nzbins;
    N->shear_zdist_zmin_all = redshift.shear_zdist_zmin_all;
    N->shear_zdist_zmax_all = redshift.shear_zdist_zmax_all;

    for(int i=0; i<MAX_SIZE_ARRAYS; i++)
    {
      N->shear_zdist_zmin[i] = redshift.shear_zdist_zmin[i];
      N->shear_zdist_zmax[i] = redshift.shear_zdist_zmax[i];
    }

    const int nzbins = redshift.shear_nzbins;
    const int ntomo  = redshift.shear_nbin;

    if (N->shear_zdist_table != NULL) free(N->shear_zdist_table);
    N->shear_zdist_table = (double**) malloc2d(ntomo + 1, nzbins);

    #pragma omp parallel for
    for (int i=0; i<ntomo+1; i++)
    {
      for (int k=0; k<nzbins; k++) 
      {
        N->shear_zdist_table[i,k] = redshift.shear_zdist_table[i,k];
      }
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------