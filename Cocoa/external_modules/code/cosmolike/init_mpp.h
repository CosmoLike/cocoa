#ifndef __COSMOLIKE_INIT_MPP_H
#define __COSMOLIKE_INIT_MPP_H
#ifdef __cplusplus
extern "C" {
#endif

// these routines need to exposed to python
void init_probes_5x2pt(char *probes);

void init_probes_real_mpp(char *probes);

void init_source_sample_mpp(char *multihisto_file, int Ntomo);

void init_lens_sample_mpp(char *multihisto_file, int Ntomo, double *b1,
                          double *b2, double ggl_cut);

void init_binning_mpp(int Ntheta, double theta_min_arcmin,
                      double theta_max_arcmin);

void init_IA_mpp(int N);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD