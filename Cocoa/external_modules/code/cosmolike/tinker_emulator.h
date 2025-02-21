#ifndef __TINKER_EMULATOR_H
#define __TINKER_EMULATOR_H
#ifdef __cplusplus
extern "C" {
#endif

typedef struct 
{
  int tinker_bias_ncosmo;
  int tinker_bias_nparam;
  int tinker_bias_nsamp;
  int tinker_bias_nparam_redshift;
  int tinker_bias_extrapolation_cut_in;
  int tinker_bias_extrapolation_cut_out;
  int tinker_hmf_ncosmo;
  int tinker_hmf_nparam;
  int tinker_hmf_nsamp;
  int tinker_hmf_nparam_redshift;
  int tinker_hmf_extrapolation_cut_in;
  int tinker_hmf_extrapolation_cut_out;
} TinkerEmuParameters;

extern TinkerEmuParameters tinkerEmuParam;

double tinker_emulator_massfunc(const double m, const double a); // Neff hardcoded

double tinker_emulator_B1(const double m, const double a);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD