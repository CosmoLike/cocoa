#include "init_baryon.h"
#include "parameters_bary.h"
#include "structs.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void init_bary(char *scenario) {

  printf("\n");
  printf("-------------------------------\n");
  printf("Initializing baryonic Parameters\n");
  printf("-------------------------------\n");

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
