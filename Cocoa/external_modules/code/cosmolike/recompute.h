#include "structs.h"

#ifndef __COSMOLIKE_RECOMPUTE_H
#define __COSMOLIKE_RECOMPUTE_H
#ifdef __cplusplus
extern "C" {
#endif

void update_cosmopara(cosmopara *C);

void update_galpara(galpara *G);

void update_nuisance(nuisancepara *N);

int recompute_cosmo3D(cosmopara C);

int recompute_zphot_shear(nuisancepara N);

int recompute_zphot_clustering(nuisancepara N);

int recompute_zphot_magnification(nuisancepara N);

// for shear 2-pt statics
int recompute_shear(cosmopara C, nuisancepara N);

// for shear 2-pt statics
int recompute_ii(cosmopara C, nuisancepara N);

// for gg-lensing statistics
int recompute_ggl(cosmopara C, galpara G, nuisancepara N, int i);

// clustering
int recompute_clustering(cosmopara C, galpara G, nuisancepara N, int i, int j);

// recompute criteria
int recompute_clusters(cosmopara C, nuisancepara N);

int recompute_PkRatio(barypara B);

void update_PkRatio(barypara *B);

// recompute criteria
int recompute_DESclusters(cosmopara C, nuisancepara N);

int recompute_galaxies(galpara G, int i);

int recompute_ks(cosmopara C, galpara G, nuisancepara N, int i);

int recompute_gk(cosmopara C, galpara G, nuisancepara N, int i);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD