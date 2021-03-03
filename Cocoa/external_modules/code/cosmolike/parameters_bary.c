#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>

#include "basics.h"
#include "parameters_bary.h"
#include "structs.h"

void init_bary(char *scenario) {
  if (strcmp(scenario, "dmo") == 0)
    bary.isPkbary = 0;
  if (strcmp(scenario, "mb2") == 0)
    set_bary_parameters_to_mb2(); // go to check parameters_bary.c
  if (strcmp(scenario, "illustris") == 0)
    set_bary_parameters_to_illustris();
  if (strcmp(scenario, "eagle") == 0)
    set_bary_parameters_to_eagle();
  if (strcmp(scenario, "HzAGN") == 0)
    set_bary_parameters_to_HzAGN();
  if (strcmp(scenario, "TNG100") == 0)
    set_bary_parameters_to_TNG100();
  if (strcmp(scenario, "owls_AGN") == 0)
    set_bary_parameters_to_owls_AGN();
  if (strcmp(scenario, "owls_DBLIMFV1618") == 0)
    set_bary_parameters_to_owls_DBLIMFV1618();
  if (strcmp(scenario, "owls_NOSN") == 0)
    set_bary_parameters_to_owls_NOSN();
  if (strcmp(scenario, "owls_NOSN_NOZCOOL") == 0)
    set_bary_parameters_to_owls_NOSN_NOZCOOL();
  if (strcmp(scenario, "owls_NOZCOOL") == 0)
    set_bary_parameters_to_owls_NOZCOOL();
  if (strcmp(scenario, "owls_REF") == 0)
    set_bary_parameters_to_owls_REF();
  if (strcmp(scenario, "owls_WDENS") == 0)
    set_bary_parameters_to_owls_WDENS();
  if (strcmp(scenario, "owls_WML1V848") == 0)
    set_bary_parameters_to_owls_WML1V848();
  if (strcmp(scenario, "owls_WML4") == 0)
    set_bary_parameters_to_owls_WML4();
}

void set_bary_parameters_to_TNG100() {
  sprintf(bary.FILE_logPkR,
          "../cosmolike_core/logPkRatio/logPkRatio_TNG100.dat");
  sprintf(bary.scenario, "TNG100");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 13;
  double z[13] = {3.71, 3.49, 3.28, 2.90, 2.44, 2.1, 1.74,
                  1.41, 1.04, 0.7,  0.35, 0.18, 0.0};

  for (int i=0; i<bary.Nabins; i++)
  {
    bary.z_bins[i] = z[i];
  }
}

void set_bary_parameters_to_HzAGN() {
  sprintf(bary.FILE_logPkR,
          "../cosmolike_core/logPkRatio/logPkRatio_HzAGN.dat");
  sprintf(bary.scenario, "HzAGN");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 11;
  double z[11] = {4.9285,  4.249,    3.7384,  3.33445,  3.00295, 1.96615,
                  1.02715, 0.519195, 0.22878, 0.017865, 0.0};

  for (int i = 0; i < bary.Nabins; i++)
  {
    bary.z_bins[i] = z[i];
  }
}

void set_bary_parameters_to_mb2() {
  sprintf(bary.FILE_logPkR, "../cosmolike_core/logPkRatio/logPkRatio_mb2.dat");
  sprintf(bary.scenario, "mb2");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 21;
  double z[21] = {3.5, 3.25, 2.8, 2.45, 2.1, 2.0, 1.8,  1.7, 1.6,    1.4, 1.2,
                  1.1, 1.0,  0.8, 0.7,  0.6, 0.4, 0.35, 0.2, 0.0625, 0.0};

  for (int i = 0; i < bary.Nabins; i++)
  {
    bary.z_bins[i] = z[i];
  }
}

void set_bary_parameters_to_illustris() {
  sprintf(bary.FILE_logPkR, "../cosmolike_core/logPkRatio/logPkRatio_ill1.dat");
  sprintf(bary.scenario, "illustris");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 23;
  double z[23] = {3.5,  3.49, 3.28, 3.08, 2.90, 2.73, 2.44, 2.1,
                  2.0,  1.82, 1.74, 1.6,  1.41, 1.21, 1.04, 1.0,
                  0.79, 0.7,  0.6,  0.4,  0.35, 0.2,  0.0};

  for (int i = 0; i < bary.Nabins; i++)
  {
    bary.z_bins[i] = z[i];
  }
}

void set_bary_parameters_to_eagle() {
  sprintf(bary.FILE_logPkR,
          "../cosmolike_core/logPkRatio/logPkRatio_eagle.dat");
  sprintf(bary.scenario, "eagle");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 13;
  double z[13] = {3.53, 3.02, 2.48, 2.24, 2.01, 1.74, 1.49,
                  1.26, 1.0,  0.74, 0.5,  0.27, 0.0};

  for (int i = 0; i < bary.Nabins; i++)
  {
    bary.z_bins[i] = z[i];
  }
}

void set_bary_parameters_to_owls_AGN() {
  sprintf(bary.FILE_logPkR, "../cosmolike_core/logPkRatio/logPkRatio_owls_AGN.dat");
  sprintf(bary.scenario, "owls_AGN");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 16;
  double z[16] = {3.5,  3.25, 3.0,  2.75, 2.25,  2.00, 1.75,  1.50,
                  1.25, 1.00, 0.75, 0.50, 0.375, 0.25, 0.125, 0.0};

  for (int i = 0; i < bary.Nabins; i++)
    bary.z_bins[i] = z[i];
}

void set_bary_parameters_to_owls_DBLIMFV1618() {
  sprintf(bary.FILE_logPkR, "../cosmolike_core/logPkRatio/logPkRatio_owls_DBLIMFV1618.dat");
  sprintf(bary.scenario, "owls_DBLIMFV1618");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 16;
  double z[16] = {3.5,  3.25, 3.0,  2.75, 2.25,  2.00, 1.75,  1.50,
                  1.25, 1.00, 0.75, 0.50, 0.375, 0.25, 0.125, 0.0};

  for (int i = 0; i < bary.Nabins; i++)
    bary.z_bins[i] = z[i];
}

void set_bary_parameters_to_owls_NOSN() {
  sprintf(bary.FILE_logPkR, "../cosmolike_core/logPkRatio/logPkRatio_owls_NOSN.dat");
  sprintf(bary.scenario, "owls_NOSN");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 15;
  double z[15] = {
      3.5,  3.25, 3.0,  2.75, 2.25,  2.00, 1.75, 1.50,
      1.25, 1.00, 0.75, 0.50, 0.375, 0.25, 0.0}; // no z=0.125 col for NOSN

  for (int i = 0; i < bary.Nabins; i++)
    bary.z_bins[i] = z[i];
}

void set_bary_parameters_to_owls_NOSN_NOZCOOL() {
  sprintf(bary.FILE_logPkR, "../cosmolike_core/logPkRatio/logPkRatio_owls_NOSN_NOZCOOL.dat");
  sprintf(bary.scenario, "owls_NOSN_NOZCOOL");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 16;
  double z[16] = {3.5,  3.25, 3.0,  2.75, 2.25,  2.00, 1.75,  1.50,
                  1.25, 1.00, 0.75, 0.50, 0.375, 0.25, 0.125, 0.0};

  for (int i = 0; i < bary.Nabins; i++)
    bary.z_bins[i] = z[i];
}

void set_bary_parameters_to_owls_NOZCOOL() {
  sprintf(bary.FILE_logPkR,
          "../cosmolike_core/logPkRatio/logPkRatio_owls_NOZCOOL.dat");
  sprintf(bary.scenario, "owls_NOZCOOL");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 16;
  double z[16] = {3.5,  3.25, 3.0,  2.75, 2.25,  2.00, 1.75,  1.50,
                  1.25, 1.00, 0.75, 0.50, 0.375, 0.25, 0.125, 0.0};

  for (int i = 0; i < bary.Nabins; i++)
    bary.z_bins[i] = z[i];
}

void set_bary_parameters_to_owls_REF() {
  sprintf(bary.FILE_logPkR,
          "../cosmolike_core/logPkRatio/logPkRatio_owls_REF.dat");
  sprintf(bary.scenario, "owls_REF");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 16;
  double z[16] = {3.5,  3.25, 3.0,  2.75, 2.25,  2.00, 1.75,  1.50,
                  1.25, 1.00, 0.75, 0.50, 0.375, 0.25, 0.125, 0.0};

  for (int i = 0; i < bary.Nabins; i++)
    bary.z_bins[i] = z[i];
}

void set_bary_parameters_to_owls_WDENS() {
  sprintf(bary.FILE_logPkR,
          "../cosmolike_core/logPkRatio/logPkRatio_owls_WDENS.dat");
  sprintf(bary.scenario, "owls_WDENS");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 16;
  double z[16] = {3.5,  3.25, 3.0,  2.75, 2.25,  2.00, 1.75,  1.50,
                  1.25, 1.00, 0.75, 0.50, 0.375, 0.25, 0.125, 0.0};

  for (int i = 0; i < bary.Nabins; i++)
    bary.z_bins[i] = z[i];
}

void set_bary_parameters_to_owls_WML1V848() {
  sprintf(bary.FILE_logPkR,
          "../cosmolike_core/logPkRatio/logPkRatio_owls_WML1V848.dat");
  sprintf(bary.scenario, "owls_WML1V848");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 16;
  double z[16] = {3.5,  3.25, 3.0,  2.75, 2.25,  2.00, 1.75,  1.50,
                  1.25, 1.00, 0.75, 0.50, 0.375, 0.25, 0.125, 0.0};

  for (int i = 0; i < bary.Nabins; i++)
    bary.z_bins[i] = z[i];
}

void set_bary_parameters_to_owls_WML4() {
  sprintf(bary.FILE_logPkR,
          "../cosmolike_core/logPkRatio/logPkRatio_owls_WML4.dat");
  sprintf(bary.scenario, "owls_WML4");
  bary.isPkbary = 1;
  bary.Nkbins = line_count(bary.FILE_logPkR) - 1;
  bary.Nabins = 16;
  double z[16] = {3.5,  3.25, 3.0,  2.75, 2.25,  2.00, 1.75,  1.50,
                  1.25, 1.00, 0.75, 0.50, 0.375, 0.25, 0.125, 0.0};

  for (int i = 0; i < bary.Nabins; i++)
    bary.z_bins[i] = z[i];
}