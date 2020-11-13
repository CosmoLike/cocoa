#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#include <time.h>

#include <fftw3.h>

#include "utils.h"
#include "utils_complex.h"
#include "cfastpt.h"

int main(int argc, char const *argv[])
{
	FILE *finput;
	finput = fopen("Pk_test_cl", "r");
	long Nk = 70;
	long line_num;
	double k[Nk], Pin[Nk];
	double dummy; 
	if(finput == NULL) 
	{printf("File not found\n");}
	else
	{
		line_num = 0;
		while(!feof(finput)) {
			fscanf(finput, "%lg %lg %lg %lg", &(k[line_num]), &(Pin[line_num]), &dummy, &dummy );
			line_num++;
		}
		fclose(finput);
	}

	// int alpha_ar[] = {0,0,0,2,1,1,2};
	// int beta_ar[]  = {0,0,0,-2,-1,-1,-2};
	// int ell_ar[]   = {0,2,4,2,1,3,0};


	double Pout[Nk];
	double IA_tt_EE[Nk],IA_tt_BB[Nk];
	double IA_ta_dE1[Nk], IA_ta_dE2[Nk], IA_ta_0E0E[Nk], IA_ta_0B0B[Nk];
	double IA_mix_A[Nk], IA_mix_B[Nk], IA_mix_DEE[Nk], IA_mix_DBB[Nk];

	clock_t t1, t2;
	t1 = clock();
	// Pd1d2(k, Pin, Nk, Pout);
	// Pd2d2(k, Pin, Nk, Pout);
	// Pd1s2(k, Pin, Nk, Pout);
	// Pd2s2(k, Pin, Nk, Pout);
	// Ps2s2(k, Pin, Nk, Pout);

	IA_tt(k, Pin, Nk, IA_tt_EE, IA_tt_BB);
	printf("tt finished\n");
	IA_ta(k, Pin, Nk, IA_ta_dE1, IA_ta_dE2, IA_ta_0E0E, IA_ta_0B0B);
	printf("ta finished\n");
	IA_mix(k,Pin, Nk, IA_mix_A, IA_mix_B, IA_mix_DEE, IA_mix_DBB);
	printf("mix finished\n");
	for (int i = 0; i < Nk; i++){
		printf("%e %e   %e %e %e %e %e %e %e %e %e %e\n",k[i],Pin[i], IA_tt_EE[i],IA_tt_BB[i], IA_ta_dE1[i],IA_ta_dE2[i],IA_ta_0E0E[i],IA_ta_0B0B[i],IA_mix_A[i], IA_mix_B[i], IA_mix_DEE[i],IA_mix_DBB[i]);
	}
	t2 = clock();
	printf("time: %lg\n", (double)(t2 - t1) / CLOCKS_PER_SEC);



	FILE *fout;
	fout = fopen("out.txt", "w");
	for(line_num=0; line_num<Nk; line_num++){
		// fprintf(fout, "%lg %lg %lg\n", k[line_num], Pin[line_num], Pout[line_num]);
		fprintf(fout, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", k[line_num], Pin[line_num], IA_tt_EE[line_num], IA_tt_BB[line_num], \
														   IA_ta_dE1[line_num], IA_ta_dE2[line_num], IA_ta_0E0E[line_num], IA_ta_0B0B[line_num],\
														   IA_mix_A[line_num], IA_mix_B[line_num], IA_mix_DEE[line_num], IA_mix_DBB[line_num]);
	}
	fclose(fout);
	return 0;
}