#ifndef __COSMOLIKE_PARAMETERS_BARYONS_H
#define __COSMOLIKE_PARAMETERS_BARYONS_H
#ifdef __cplusplus
extern "C" {
#endif

void init_baryons(const char* sim);

void init_baryons_from_hdf5_file(const char* sim, int tag, const char* allsims);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD