#ifndef __COSMOLIKE_IA_H
#define __COSMOLIKE_IA_H
#ifdef __cplusplus
extern "C" {
#endif

void set_LF_GAMA(void);

void set_LF_DEEP2(void);

double A_IA_Joachimi(const double a);

void IA_A1_Z1Z2(const double a, 
  const double growfac_a, 
  const int n1, const int n2, double res[2]);

double IA_A1_Z1(const double a, const double growfac_a, const int n1);

void IA_A2_Z1Z2(const double a, 
  const double growfac_a, 
  const int n1, const int n2, double res[2]);

double IA_A2_Z1(const double a, const double growfac_a, const int n1);

void IA_BTA_Z1Z2(const double a __attribute__((unused)), 
  const double growfac_a __attribute__((unused)), 
  const int n1, const int n2, double res[2]);

double IA_BTA_Z1(const double a, const double growfac_a, const int n1);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD