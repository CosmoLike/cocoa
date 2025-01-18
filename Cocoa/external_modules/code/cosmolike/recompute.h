#include "structs.h"

#ifndef __COSMOLIKE_RECOMPUTE_H
#define __COSMOLIKE_RECOMPUTE_H
#ifdef __cplusplus
extern "C" {
#endif

// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// RECOMPUTE
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

int recompute_cosmo3D(cosmopara C);

/*

int recompute_clusters(cosmopara C, nuisancepara N);

int recompute_cc(cosmopara C, nuisancepara N);

int recompute_cg(cosmopara C, galpara G, nuisancepara N);

int recompute_cs(cosmopara C, galpara G, nuisancepara N);
*/

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD