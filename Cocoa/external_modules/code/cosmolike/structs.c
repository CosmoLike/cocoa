#include "structs.h"

likepara like = {.baryons = 0,
                 .IA = 0.,
                 .bias = 0,
                 .wlphotoz = 0,
                 .clphotoz = 0,
                 .shearcalib = 0,
                 .clusterMobs = 0,
                 .BAO = 0,
                 .SN_WFIRST = 0,
                 .GRS = 0,
                 .SRD = 0,
                 .Planck15_BAO_H070p6_JLA_w0wa = 0,
                 .Planck18_BAO_Riess18_Pantheon_w0wa = 0,
                 .Planck18_BAO_w0wa = 0,
                 .Planck18_w0 = 0,
                 .theta_s = 0};

cosmopara cosmology = {.Omega_nu = 0.,
                       .coverH0 = 2997.92458,
                       .rho_crit = 7.4775e+21,
                       .MGSigma = 0.0,
                       .MGmu = 0.0};

tomopara tomo = {.n_source = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
                 .n_lens = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}};

redshiftpara redshift;

sur survey = {.area_conversion_factor =
                  60.0 * 60.0 * 2.90888208665721580e-4 * 2.90888208665721580e-4,
              .n_gal_conversion_factor =
                  1.0 / 2.90888208665721580e-4 / 2.90888208665721580e-4,
              .ggl_overlap_cut = 1.};

Cmb cmb;

galpara gbias = {.b2 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 .bs2 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 .b1_function = &b1_per_bin,
                 .b_mag = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0}}; // default: point to old bgal_z routin

clusterpara Cluster = {.model = "default"};

pdeltapara pdeltaparams = {.runmode = "Halofit"};

FPTpara FPT = {.k_min = 1.e-5,
               .k_max = 1.e+3,
               .N = 800,
               .N_per_dec = 100,
               .N_AB = 7,
               .N_IA = 10};

nuisancepara nuisance = {
    .c1rhocrit_ia = 0.01389,
    .A_z = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    .A2_z = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    .b_ta_z = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    .shear_calibration_m = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    .sigma_zphot_shear = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    .bias_zphot_shear = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    .sigma_zphot_clustering = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    .bias_zphot_clustering = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}};


barypara bary = {.isPkbary = 0};

fft_optimize fft_int;