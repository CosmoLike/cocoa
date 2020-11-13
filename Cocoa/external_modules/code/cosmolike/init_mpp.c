#include "init_mpp.h"
#include "basics.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void init_probes_5x2pt(char *probes) {
  printf("\n");
  printf("------------------------------\n");
  printf("Initializing 5x2pt Probes\n");
  printf("------------------------------\n");

  sprintf(like.probes, "%s", probes);
  like.shear_shear = 0;
  like.shear_pos = 0;
  like.pos_pos = 0;
  like.gk = 0;
  like.ks = 0;
  if (strcmp(probes, "all_2pt") == 0) {
    printf("init_mpp.c: called init_probes_5x2pt with outdated argument "
           "probes=%s\n",
           probes);
    printf("Update your code to use either \"3x2pt\" or \"3x2pt\"\n");
    exit(1);
  }
  char *type1 = "xi";
  if (strcmp(probes, "shear_shear") == 0 || strstr(probes, type1) != NULL) {
    like.shear_shear = 1;
    printf("Shear-Shear computation initialized\n");
  }

  char *type2 = "wtheta";
  if (strcmp(probes, "pos_pos") == 0 || strstr(probes, type2) != NULL) {
    like.pos_pos = 1;
    printf("Position-Position computation initialized\n");
  }

  char *type3 = "ggl";
  char *type4 = "gammat";
  if (strcmp(probes, "shear_pos") == 0 || strstr(probes, type3) != NULL ||
      strstr(probes, type4) != NULL) {
    if (tomo.ggl_Npowerspectra == 0) {
      void init_ggl_tomo();
    }
    like.shear_pos = 1;
    printf("Position-Shear computation initialized\n");
  }
  char *type5 = "skappa";
  if (strstr(probes, type5) != NULL) {
    like.ks = 1;
    printf("Kappa-Shear computation initialized\n");
  }
  char *type6 = "gkappa";
  if (strstr(probes, type6) != NULL) {
    like.gk = 1;
    printf("Kappa-galaxy computation initialized\n");
  }
  if (strcmp(probes, "3x2pt") == 0) {
    if (tomo.ggl_Npowerspectra == 0) {
      void init_ggl_tomo();
    }
    like.shear_shear = 1;
    like.shear_pos = 1;
    like.pos_pos = 1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  }

  if (strcmp(probes, "5x2pt") == 0) {
    if (tomo.ggl_Npowerspectra == 0) {
      void init_ggl_tomo();
    }
    like.shear_shear = 1;
    like.shear_pos = 1;
    like.pos_pos = 1;
    like.ks = 1;
    like.gk = 1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
    printf("Kappa-Shear computation initialized\n");
    printf("Kappa-galaxy computation initialized\n");
  }
  if (strcmp(probes, "ggl_cl") == 0) {
    if (tomo.ggl_Npowerspectra == 0) {
      void init_ggl_tomo();
    }
    like.shear_pos = 1;
    like.pos_pos = 1;
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  }
  printf("like.pos_pos = %d, like.shear_pos = %d,like.shear_shear = %d,like.ks "
         "= %d,like.gl = %d\n\n",
         like.pos_pos, like.shear_pos, like.shear_shear, like.ks, like.gk);
  printf("like.Ntheta = %d\n", like.Ntheta);
  printf("tomo.shear_Npowerspectra = %d\n", tomo.shear_Npowerspectra);
  printf("tomo.ggl_Npowerspectra = %d\n", tomo.ggl_Npowerspectra);
  printf("tomo.clustering_Npowerspectra = %d\n", tomo.clustering_Npowerspectra);

  like.Ndata =
      like.Ntheta *
      (2 * tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra +
       tomo.clustering_Npowerspectra + tomo.shear_Nbin + tomo.clustering_Nbin);

  printf("Total number of data points like.Ndata=%d\n", like.Ndata);
}

void init_probes_real_mpp(char *probes) {
#ifdef DEBUG 
  printf("\n");
  printf("------------------------------\n");
  printf("Initializing 3x2pt Probes\n");
  printf("------------------------------\n");
#endif 
  sprintf(like.probes, "%s", probes);
  like.shear_shear = 0;
  like.shear_pos = 0;
  like.pos_pos = 0;
  like.gk = 0;
  like.ks = 0;

  if (strcmp(probes, "all_2pt") == 0) {   
    printf("init_mpp.c: called init_probes_real_mpp with outdated argument "
           "probes=%s\n",
           probes);
    printf("Update your code to use \"3x2pt\" instead\n");   
    exit(1);
  }

  char *type1 = "xi";
  if (strcmp(probes, "shear_shear") == 0 || strstr(probes, type1) != NULL) {
    like.shear_shear = 1;
#ifdef DEBUG     
    printf("Shear-Shear computation initialized\n");
#endif     
  }

  char *type2 = "wtheta";
  if (strcmp(probes, "pos_pos") == 0 || strstr(probes, type2) != NULL) {
    like.pos_pos = 1;
#ifdef DEBUG    
    printf("Position-Position computation initialized\n");
#endif     
  }

  char *type3 = "ggl";
  char *type4 = "gammat";
  if (strcmp(probes, "shear_pos") == 0 || strstr(probes, type3) != NULL ||
      strstr(probes, type4) != NULL) {
    if (tomo.ggl_Npowerspectra == 0) {
      void init_ggl_tomo();
    }
    like.shear_pos = 1;
#ifdef DEBUG     
    printf("Position-Shear computation initialized\n");
#endif     
  }

  //if (strcmp(probes, "all_2pt") == 0 || strcmp(probes, "3x2pt") == 0) {
  if (strcmp(probes, "3x2pt") == 0) {  
    if (tomo.ggl_Npowerspectra == 0) {
      void init_ggl_tomo();
    }
    like.shear_shear = 1;
    like.shear_pos = 1;
    like.pos_pos = 1;
#ifdef DEBUG    
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
#endif    
  }

  if (strcmp(probes, "ggl_cl") == 0) {
    if (tomo.ggl_Npowerspectra == 0) {
      void init_ggl_tomo();
    }
    like.shear_pos = 1;
    like.pos_pos = 1;
#ifdef DEBUG    
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
#endif    
  }
#ifdef DEBUG
  printf("like.pos_pos = %d, like.shear_pos = %d,like.shear_shear = %d\n\n",
         like.pos_pos, like.shear_pos, like.shear_shear);
  printf("like.Ntheta = %d\n", like.Ntheta);
  printf("tomo.shear_Npowerspectra = %d\n", tomo.shear_Npowerspectra);
  printf("tomo.ggl_Npowerspectra = %d\n", tomo.ggl_Npowerspectra);
  printf("tomo.clustering_Npowerspectra = %d\n", tomo.clustering_Npowerspectra);
#endif
  like.Ndata =
      like.Ntheta * (2 * tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra +
                     tomo.clustering_Npowerspectra);
#ifdef DEBUG
  printf("Total number of data points like.Ndata=%d\n", like.Ndata);
#endif  
}

void init_survey_mpp(char *surveyname, double area, double sigma_e) {
  sprintf(survey.name, "%s", surveyname);
  survey.area = area;
  survey.sigma_e = sigma_e;
}

void init_source_sample_mpp(char *multihisto_file, int Ntomo) {
  sprintf(redshift.shear_REDSHIFT_FILE, "%s", multihisto_file);
  redshift.shear_photoz = 4;
  tomo.shear_Nbin = Ntomo;
  tomo.shear_Npowerspectra = tomo.shear_Nbin * (tomo.shear_Nbin + 1) / 2;
  printf(
      "Source redshifts: multi-histo file %s, containing %d tomography bins\n",
      multihisto_file, tomo.shear_Nbin);
  for (int i = 0; i < tomo.shear_Nbin; i++) {
    printf("bin %d: <z_s>=%f\n", i, zmean_source(i));
    // tomo.n_source[i]= n_source[i];
    nuisance.bias_zphot_shear[i] = 0.0;
  }
}

void init_ggl_tomo() {
#ifdef DEBUG  
  if (tomo.clustering_Nbin == 0) {
    printf("WARNING! init_mpp.c: init_ggl_tomo called while "
           "tomo.clustering_Nbin =0\n");
  }
  if (tomo.shear_Nbin == 0) {
    printf(
        "WARNING! init_mpp.c: init_ggl_tomo called while tomo.shear_Nbin =0\n");
  }
#endif  
  int n = 0;
  for (int i = 0; i < tomo.clustering_Nbin; i++) {
    for (int j = 0; j < tomo.shear_Nbin; j++) {
      n += test_zoverlap(i, j);
#ifdef DEBUG      
      // printf("GGL combinations zl=%d zs=%d accept=%d; <z_l> = %.3f, <z_s> =
      // %.3f\n",i,j,test_zoverlap(i,j), zmean(i),zmean_source(j));
#endif       
    }
  }
  tomo.ggl_Npowerspectra = n;
#ifdef DEBUG  
  printf("%d GGL Powerspectra\n", tomo.ggl_Npowerspectra);
#endif   
}

void init_lens_sample_mpp(char *multihisto_file, int Ntomo, double *b1,
                          double *b2, double ggl_cut) {
  sprintf(redshift.clustering_REDSHIFT_FILE, "%s", multihisto_file);
  redshift.clustering_photoz = 4;
  tomo.clustering_Nbin = Ntomo;
  tomo.clustering_Npowerspectra = tomo.clustering_Nbin;
  if (ggl_cut > 0) {
    survey.ggl_overlap_cut = ggl_cut;
  }
  printf("Lens redshifts: multi-histo file %s, containing %d tomography bins\n",
         multihisto_file, tomo.clustering_Nbin);
  pf_photoz(0.1, 0);
  for (int i = 0; i < tomo.clustering_Nbin; i++) {
    gbias.b1_function = &b1_per_bin;
    //   tomo.n_lens[i]= n_lens[i];
    gbias.b[i] = b1[i];
    gbias.b2[i] = b2[i];
    nuisance.bias_zphot_clustering[i] = 0.0;
    //   printf("bin %d: <z_l>=%.3f, b_1=%.3f,
    //   b_2=%.3f\n",i,zmean(i),gbias.b[i],gbias.b2[i]);
  }
  init_ggl_tomo();
}

void init_binning_mpp(int Ntheta, double theta_min_arcmin,
                      double theta_max_arcmin) {
  like.Ntheta = Ntheta;
  like.vtmin = theta_min_arcmin * constants.arcmin;
  like.vtmax = theta_max_arcmin * constants.arcmin;
  if (tomo.shear_Npowerspectra + tomo.clustering_Nbin +
      tomo.ggl_Npowerspectra == 0) {
    printf("init_mpp.c: init_binning_mpp called with Npowerspectra = 0.\n "
           "EXIT!\n");
    exit(1);
  }
}

void init_IA_mpp(int N) {
  if (N == 3) {
    like.IA = N;
#ifdef DEBUG    
    printf("Set like.IA =3; supply one A_IA parameter per source tomo bin\n");
#endif    
  } else if (N == 4) {
    like.IA = N;
#ifdef DEBUG    
    printf("Set like.IA =4; supply nuisance.A_ia and nuisance.eta_ia\n");
#endif    
  } else {
    printf("like.IA = %d not supported in des_mpp\nEXIT\n", N);
    exit(1);
  }
}
