#ifndef __COSMOLIKE_PT_C_FASTPT_H
#define __COSMOLIKE_PT_C_FASTPT_H
#ifdef __cplusplus
extern "C" {
#endif

void get_FPT_IA(void);

double PT_d1d2(double k_coverH0);

double PT_d1d3(double k_coverH0);

double PT_d2d2(double k_coverH0);

double PT_d1s2(double k_coverH0);

double PT_d2s2(double k_coverH0);

double PT_s2s2(double k_coverH0);

double PT_sigma4(double k_coverH0);

double TATT_II_EE(double k_coverH0, double a, double C1, double C2, double b_ta,
	double C1_2, double C2_2, double b_ta_2, double growfac_a, double pdelta_ak);

double TATT_II_BB(double k_coverH0, double a, double C1, double C2, double b_ta,
  double C1_2, double C2_2, double b_ta_2, double growfac_a);

double TATT_GI_E(double k_coverH0, double a, double C1, double C2, double b_ta, 
	double growfac_a, double pdelta_ak);

void FPT_input(double k[], double P[]);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD