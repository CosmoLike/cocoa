#include <stdio.h>
#include "structs.h"

likepara like =
{
  .IA = 0.,
  .bias = 0,
  .clusterMobs = 0,
  .ell = NULL,
  .theta = NULL,
};

void reset_like_struct()
{
  like.IA = 0.;
  like.bias = 0;
  like.clusterMobs = 0;
  like.Ncl = 0;
  like.Ntheta = 0;
  like.Ncos = 0;
  like.Ndata = 0;
  like.lmin = 0;
  like.lmax = 0;
  like.vtmax = 0;
  like.vtmin = 0;
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
  like.lmax_kappacmb = 0;
  like.IA = 0;
  like.bias = 0;
  like.clusterMobs = 0;
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
}

cosmopara cosmology =
{
  .Omega_nu = 0.,
  .coverH0 = 2997.92458,
  .rho_crit = 7.4775e+21,
  .MGSigma = 0.0,
  .MGmu = 0.0,
  .is_cached = 0
};

void reset_cosmology_struct()
{
  cosmology.Omega_nu = 0.;
  cosmology.coverH0 = 2997.92458;
  cosmology.rho_crit = 7.4775e+21;
  cosmology.MGSigma = 0.0;
  cosmology.MGmu = 0.0;

  cosmology.is_cached = 0;

  cosmology.Omega_m = 0.0;
  cosmology.Omega_v = 0.0;
  cosmology.h0 = 0.0;
  cosmology.Omega_nu = 0.0;
  cosmology.random = 0.0;
}

tomopara tomo = {
  .n_source = {0},
  .n_lens = {0}
};

void reset_tomo_struct()
{
  for(int i=0; i<MAX_SIZE_ARRAYS; i++)
  {
    tomo.shear_zmax[i] = 0.0;
    tomo.shear_zmin[i] = 0.0;
    tomo.n_source[i] = 0.0;

    tomo.clustering_zmax[i] = 0.0;
    tomo.clustering_zmin[i] = 0.0;
    tomo.n_lens[i] = 0.0;

    tomo.cluster_zmax[i] = 0.0;
    tomo.cluster_zmin[i] = 0.0;

    tomo.magnification_zmax[i] = 0.0;
    tomo.magnification_zmin[i] = 0.0;
  }

  tomo.shear_Nbin = 0;
  tomo.shear_Npowerspectra = 0;
  tomo.clustering_Nbin = 0;
  tomo.clustering_Npowerspectra = 0;
  tomo.cluster_Nbin = 0;
  tomo.cluster_cg_Npowerspectra = 0;
  tomo.cgl_Npowerspectra = 0;
  tomo.ggl_Npowerspectra = 0;
  tomo.magnification_Nbin = 0;
  tomo.magnification_Npowerspectra = 0;
}

redshiftpara redshift;

void reset_redshift_struct()
{
  redshift.shear_photoz = 0;
  redshift.shear_zdistrpar_zmin = 0.0;
  redshift.shear_zdistrpar_zmax = 0.0;
  redshift.shear_histogram_zbins = 0;
  sprintf(redshift.shear_REDSHIFT_FILE,"");

  redshift.clustering_photoz = 0;
  redshift.clustering_zdistrpar_zmin = 0.0;
  redshift.clustering_zdistrpar_zmax = 0.0;
  redshift.clustering_histogram_zbins = 0;
  sprintf(redshift.clustering_REDSHIFT_FILE,"");

  redshift.magnification_photoz = 0;
  redshift.magnification_zdistrpar_zmin = 0.0;
  redshift.magnification_zdistrpar_zmax = 0.0;
  redshift.magnification_histogram_zbins = 0;
  sprintf(redshift.magnification_REDSHIFT_FILE,"");
}

sur survey =
{
  .area_conversion_factor =
                  60.0 * 60.0 * 2.90888208665721580e-4 * 2.90888208665721580e-4,
  .n_gal_conversion_factor =
                  1.0 / 2.90888208665721580e-4 / 2.90888208665721580e-4,
  .ggl_overlap_cut = 1.
};

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
  sprintf(survey.name,"");
  survey.ggl_overlap_cut = 1.0;
}

galpara gbias =
{
  .b2 = {0},
  .bs2 = {0},
  .b1_function = &b1_per_bin,
  .b_mag = {0},
}; // default: point to old bgal_z routin

void reset_gbias_struct()
{
  for(int i=0; i<10; i++)
  {
    gbias.b[i] = 0.0;
    gbias.b2[i] = 0.0;
    gbias.bs2[i] = 0.0;
    gbias.rcorr[i] = 0.0;
    gbias.cg[i] = 0.0;
    gbias.n_hod[i] = 0.0;
    gbias.b_mag[i] = 0.0;
  }
  gbias.b1_function = &b1_per_bin;
}

clusterpara Cluster =
{
  .model = "default"
};

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
  Cluster.lbin = 0;
  Cluster.l_min = 0.0;
  Cluster.l_min = 0.0;
  sprintf(Cluster.model,"default");
}


pdeltapara pdeltaparams =
{
  .runmode = "Halofit"
};

void reset_pdeltaparams_struct()
{
  sprintf(pdeltaparams.runmode, "Halofit");
}

FPTpara FPT =
{
  .k_min = 1.e-5,
  .k_max = 1.e+3,
  .N = 800,
  .N_per_dec = 100,
  .N_AB = 7,
  .N_IA = 10
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

void reset_nuisance_struct()
{
  for(int i=0; i<MAX_SIZE_ARRAYS; i++)
  {
    nuisance.A_z[i];
    nuisance.A2_z[i];
    nuisance.b_ta_z[i];
    nuisance.fred[i];
    nuisance.shear_calibration_m[i];
    nuisance.sigma_zphot_shear[i];
    nuisance.bias_zphot_shear[i];
    nuisance.sigma_zphot_clustering[i];
    nuisance.sigma_zphot_magnification[i];
    nuisance.bias_zphot_magnification[i];
    nuisance.cluster_completeness[i];
    nuisance.cluster_MOR[i];
    nuisance.cluster_selection[i];
  }

  nuisance.A_ia = 0.0;
  nuisance.A2_ia = 0.0;
  nuisance.beta_ia = 0.0;
  nuisance.eta_ia = 0.0;
  nuisance.eta_ia_tt = 0.0;
  nuisance.eta_ia_highz = 0.0;
  nuisance.oneplusz0_ia = 0.0;
  nuisance.c1rhocrit_ia = 0.01389;

  nuisance.cluster_Mobs_lgM0 = 0.0;
  nuisance.cluster_Mobs_sigma = 0.0;
  nuisance.cluster_Mobs_alpha = 0.0;
  nuisance.cluster_Mobs_beta = 0.0;
  nuisance.cluster_Mobs_N_pivot = 0.0;
  nuisance.cluster_Mobs_lgN0 = 0.0;
  nuisance.cluster_Mobs_sigma0 = 0.0;
  nuisance.cluster_Mobs_sigma_qm = 0.0;
  nuisance.cluster_Mobs_sigma_qz = 0.0;
  nuisance.cluster_centering_f0 = 0.0;
  nuisance.cluster_centering_alpha = 0.0;
  nuisance.cluster_centering_sigma = 0.0;
  nuisance.cluster_centering_M_pivot = 0.0;

  nuisance.N_cluster_MOR = 0;
  nuisance.N_cluster_selection = 0;
}

barypara bary = {
  .is_Pk_bary = 0,
  .Na_bins = 0,
  .Nk_bins = 0,
  .a_bins = NULL,
  .logk_bins = NULL,
  .log_PkR = NULL,
  .interp2d = NULL,
};

Cmb cmb;

void reset_cmb_struct()
{
  sprintf(cmb.name,"");
  cmb.fwhm = 0.0;
  cmb.sensitivity = 0.0;
  sprintf(cmb.pathLensRecNoise,"");
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
