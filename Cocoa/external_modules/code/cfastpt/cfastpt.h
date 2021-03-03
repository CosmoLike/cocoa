#ifndef __CFASTPT_CFASTPT_H
#define __CFASTPT_CFASTPT_H
#ifdef __cplusplus
extern "C" {
#endif

typedef struct fastpt_config {
	double nu; // only used in scalar; in tensor, nu1,nu2 are computed by alpha,beta
	double c_window_width;
	long N_pad;
	long N_extrap_low;
	long N_extrap_high;
} fastpt_config;

void fastpt_scalar(int *alpha_ar, int *beta_ar, int *ell_ar, int *isP13type_ar, double *coeff_A_ar, int Nterms, double *Pout, double *k, double *Pin, int Nk);

void J_abl_ar(double *x, double *fx, long N, int *alpha, int *beta, int *ell, int *isP13type, int Nterms, fastpt_config *config, double **Fy);

void J_abl(double *x, double *fx, int alpha, int beta, long N, fastpt_config *config, int ell, double *Fy);

void fastpt_tensor(int *alpha_ar, int *beta_ar, int *J1_ar, int *J2_ar, int *Jk_ar, double *coeff_AB_ar, int Nterms, double *Pout, double *k, double *Pin, int Nk);

void J_abJ1J2Jk_ar(double *x, double *fx, long N, int *alpha, int *beta, int *J1, int *J2, int *Jk, int Nterms, fastpt_config *config, double **Fy);

typedef struct fastpt_todo {
	int isScalar;
	double *alpha;
	double *beta;
	double *ell;
	int *isP13type;
	double *coeff_ar;
	int Nterms;
} fastpt_todo;

typedef struct fastpt_todolist {
	fastpt_todo *fastpt_todo;
	int N_todo;
} fastpt_todolist;


void Pd1d2(double *k, double *Pin, long Nk, double *Pout);
void Pd2d2(double *k, double *Pin, long Nk, double *Pout);
void Pd1s2(double *k, double *Pin, long Nk, double *Pout);
void Pd2s2(double *k, double *Pin, long Nk, double *Pout);
void Ps2s2(double *k, double *Pin, long Nk, double *Pout);

void IA_tt(double *k, double *Pin, long Nk, double *P_E, double *P_B);
void IA_ta(double *k, double *Pin, long Nk, double *P_dE1, double *P_dE2, double *P_0E0E, double *P_0B0B);
void IA_mix(double *k, double *Pin, long Nk, double *P_A, double *P_B, double *P_DEE, double *P_DBB);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD