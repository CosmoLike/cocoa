#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "utils_cfastpt.h"

void extrap_log_linear(double *fk, int N_origin, int N_extra, double *large_fk) {
	double dln_left, dln_right;
	int i;

	dln_left = log(fk[1]/fk[0]);
	if(fk[0]<=0.) {
		for(i=0; i<N_extra; i++) {
			large_fk[i] = 0.;
		}
	}
	else{
		for(i=0; i<N_extra; i++) {
			large_fk[i] = exp(log(fk[0]) + (i - N_extra) * dln_left);
		}
	}

	for(i=N_extra; i< N_extra+N_origin; i++) {
		large_fk[i] = fk[i - N_extra];
	}

	dln_right = log(fk[N_origin-1]/fk[N_origin-2]);
	if(fk[N_origin-1]<=0.) {
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i] = 0.;
		}
	}
	else {
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i] = exp(log(fk[N_origin-1]) + (i - N_extra - N_origin +1) * dln_right);
		}
	}
}

int J_table(int *alpha_ar, int *beta_ar, int *l1_ar, int *l2_ar, int *l_ar, 
double *coeff_A_ar, int Nterms, int *alpha_ar_new, int *beta_ar_new, 
int *J1_ar, int *J2_ar, int *Jk_ar, double *coeff_AB_ar) {
	int row = 0;
	for(int i=0;i<Nterms;i++){
		const int l1 = l1_ar[i]; 
		const int l2 = l2_ar[i]; 
		const int l = l_ar[i];
		const int l1_m_l2 = abs(l1-l2); 
		const int l1_m_l = abs(l1-l); 
		const int l_m_l2 = abs(l-l2);
		const int l1_p_l2 = l1+l2; 
		const int l1_p_l = l1+l; 
		const int l_p_l2 = l+l2;		
		for(int J1=l_m_l2; J1<=l_p_l2; J1++){
			for(int J2=l1_m_l; J2<=l1_p_l; J2++){
				for(int Jk=l1_m_l2; Jk<=l1_p_l2; Jk++){
					const double B = coeff_B(l1,l2,l,J1,J2,Jk);
					if(B!=0) {
						alpha_ar_new[row] = alpha_ar[i];
						beta_ar_new[row] = beta_ar[i];
						J1_ar[row] = J1;
						J2_ar[row] = J2;
						Jk_ar[row] = Jk;
						coeff_AB_ar[row] = coeff_A_ar[i] * B;						
						row++;
					}
				}
			}
		}
	}
	/*
	int row = 0;
	for(int i=0;i<Nterms;i++){
		const int l1 = l1_ar[i]; 
		const int l2 = l2_ar[i]; 
		const int l = l_ar[i];
		const int l1_m_l2 = abs(l1-l2); 
		const int l1_m_l = abs(l1-l); 
		const int l_m_l2 = abs(l-l2);
		const int l1_p_l2 = l1+l2; 
		const int l1_p_l = l1+l; 
		const int l_p_l2 = l+l2;				
		const int sz1 = (l_p_l2-l_m_l2+1);
		const int sz2 = (l1_p_l-l1_m_l+1);
		const int sz3 = (l1_p_l2-l1_m_l2+1);
		double tmp[sz1*sz2*sz3];
		#pragma omp parallel for
		for(int J1=l_m_l2; J1<=l_p_l2; J1++){
			for(int J2=l1_m_l; J2<=l1_p_l; J2++){			
				for(int Jk=l1_m_l2; Jk<=l1_p_l2; Jk++){
					const int tmp_j = J1-l_m_l2;
					const int tmp_l = J2-l1_m_l;
					const int tmp_k = Jk-l1_m_l2;
					tmp[tmp_j*(sz2*sz3)+tmp_l*sz3+tmp_k] = coeff_B(l1,l2,l,J1,J2,Jk);			
				}
			}
		}
		for(int J1=l_m_l2; J1<=l_p_l2; J1++){
			for(int J2=l1_m_l; J2<=l1_p_l; J2++){			
				for(int Jk=l1_m_l2; Jk<=l1_p_l2; Jk++){
					const int tmp_j = J1-l_m_l2;
					const int tmp_l = J2-l1_m_l;
					const int tmp_k = Jk-l1_m_l2;
					const double B = tmp[tmp_j*(sz2*sz3)+tmp_l*sz3+tmp_k];					
					if(B!=0) {
						alpha_ar_new[row] = alpha_ar[i];
						beta_ar_new[row] = beta_ar[i];
						J1_ar[row] = J1;
						J2_ar[row] = J2;
						Jk_ar[row] = Jk;
						coeff_AB_ar[row] = coeff_A_ar[i] * B;						
						row++;
					}
				}
			}
		}
	}
*/
	if(row==0) {
		printf("J_table empty! Check input coefficients!\n");exit(1);
	}
	return row;
}

double coeff_B(int l1, int l2, int l, int J1, int J2, int Jk){
	double pf;
	int sign;
	double B;

	if(((J1+l2+l)%2==0)&&((l1+J2+l)%2==0)&&((l1+l2+Jk)%2==0)&&((J1+J2+Jk)%2==0)){
		sign = ((l+(J1+J2+Jk)/2)%2? -1:1);
		pf = sign * (2*J1+1)*(2*J2+1)*(2*Jk+1) / (M_PI*M_PI*M_PI);
		B = pf *wigner_3j_jjj_000(J1,l2,l)*wigner_3j_jjj_000(l1,J2,l)*wigner_3j_jjj_000(l1,l2,Jk) \
			   *wigner_3j_jjj_000(J1,J2,Jk)*wigner_6j(J1,J2,Jk,l1,l2,l);
		return B;
	}else{
		return 0.;
	}
}

long factorial(int n){
	if(n<0){printf("factorial(n): n=%d \n",n);exit(1);}
	static long FACTORIAL_LIST[] = {1,1,2,6,24,120,\
									 720,5040,40320,362880,\
									 39916800,479001600,6227020800,87178291200};
	if(n>14){printf("factorial(n): n=%d too large, need to extend FACTORIAL_LIST!\n",n);exit(1);}
	return FACTORIAL_LIST[n];
}


// special case of wigner 3j symbol
double wigner_3j_jjj_000(int j1, int j2, int j3){
	int J = j1+j2+j3;
	int halfJ, sign;
	double wigner3j;
	double pf;
	if(J%2==0){
		halfJ = J/2;
		sign = (halfJ%2? -1:1);
		pf = Delta_coef(j1,j2,j3);
		if(pf==0) return 0.;
		wigner3j = sign* pf * (double)factorial(halfJ) / ((double)factorial(halfJ-j1)*(double)factorial(halfJ-j2)*(double)factorial(halfJ-j3));
		return wigner3j;
	}else{
		return 0.;
	}
}

// special case: only integer angular momenta
double Delta_coef(int a, int b, int c){
	int ab_c, ac_b, bc_a;
	ab_c = a+b-c; ac_b = a+c-b; bc_a = b+c-a;
	if(ab_c<0) return 0;
	if(ac_b<0) return 0;
	if(bc_a<0) return 0;

	return sqrt((double)factorial(ab_c)*(double)factorial(ac_b)*(double)factorial(bc_a)/(double)factorial(a+b+c+1));
}

double Racah(int a, int b, int c, int d, int e, int f){
	double pf=Delta_coef(a,b,e)*Delta_coef(c,d,e)*Delta_coef(a,c,f)*Delta_coef(b,d,f);
	if(pf==0.) return 0.;

	int imin, imax;
	imin = ((a+b+e)>(c+d+e)? (a+b+e):(c+d+e));
	if(a+c+f>imin) imin=a+c+f;
	if(b+d+f>imin) imin=b+d+f;

	imax = ((a+b+c+d)<(a+d+e+f)? (a+b+c+d):(a+d+e+f));
	if(b+c+e+f<imax) imax=b+c+e+f;

	double sum = 0.0;
	long denom;
	int i;
	int sign;
	for(i=imin;i<=imax;i++){
		denom=factorial(i-a-b-e)*factorial(i-c-d-e)*factorial(i-a-c-f)*factorial(i-b-d-f)\
			  *factorial(a+b+c+d-i)*factorial(a+d+e+f-i)*factorial(b+c+e+f-i);
		sign = (i%2? -1:1);
		sum += (sign*factorial(i+1))/(double)denom;
	}

	sign = ((a+b+c+d)%2? -1:1);
	return sign*pf*sum;
}

double wigner_6j(int j1, int j2, int j3, int j4, int j5, int j6){
	int sign = ((j1+j2+j4+j5)%2? -1:1);
	return sign * Racah(j1, j2, j5, j4, j3, j6);
}