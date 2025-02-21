#ifndef __COSMOLIKE_PARAMETERS_BARYONS_H
#define __COSMOLIKE_PARAMETERS_BARYONS_H
#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  int is_Pk_bary;
  int Na_bins;
  int Nk_bins;
  double* a_bins;
  double* logk_bins;
  double* log_PkR;
  gsl_interp2d_type* T;
  gsl_interp2d* interp2d;
} barypara;

extern barypara bary;

void reset_bary_struct();

void init_baryons(const char* sim);

void init_baryons_from_hdf5_file(const char* sim, int tag, const char* allsims);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD