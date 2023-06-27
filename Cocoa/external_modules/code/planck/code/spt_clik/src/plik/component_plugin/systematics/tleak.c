#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* 
   ====================================
   set of low level routines necessary to compute the coupling of power spectra due to beam mismatch

   v0.1: EH, IAP, 2014-06-25
   ====================================
*/

/* -------------------------------------- */
double polynomial(long, double*, double);
int compute_cl_t2pmatrix(double, double, long, double*);
/* -------------------------------------- */

double polynomial(long ncoeff,   /* number of coefficients */
      double *coeff, /* list of coefficients */
      double x)      /* variable */
     /* return coeff[0] + coeff[1]*x + coeff[2]*x^2 + .... */
{
  double p = coeff[ncoeff-1];
  long d;

  for (d = ncoeff-2 ; d >= 0 ; d--){
    p *= x ;
    p += coeff[d];
  }
  return p;
}

int compute_cl_t2pmatrix(double e, /* value of epsilon polynomial in l */
       double b, /* value of beta    polynomial in l */
       long nk,  /* linear size of matrix */
       double *matrix)
     /* fill matrix coupling TT, EE, ... spectra because of beam mismatch
       cf Eq 11 of  http://hfilfi.planck.fr/index.php/WG/CorrectingPolClFromBeamMismatch */
{
  long i,j;
  double zero=(double)0;
  double one=(double)1;

  if (nk != 6) {
    printf("Error: matrix must be 6x6\n");
    printf("current size: %ldx%ld \n",nk,nk );
    return -1;
  }

  /* identity matrix */
  for (i=0; i<nk ; i++){
    for (j=0; j<nk ; j++){ matrix[i*nk+j] = zero;}
    matrix[i*nk+i] = one;
  }

  /* spectra are in order TT, EE, BB, TE, TB, EB */
  matrix[1*nk + 0] = e*e ; matrix[1*nk + 3] = 2*e ;
  matrix[2*nk + 0] = b*b ;                        matrix[2*nk + 4] = 2*b ;
  matrix[3*nk + 0] = e ;
  matrix[4*nk + 0] = b ;
  matrix[5*nk + 0] = e*b ; matrix[5*nk + 3] = b ; matrix[5*nk + 4] = e;

  return 0;

}

/* ------------------------------------------------------------------- */
/* interface with likelihood code */
int t2p_beam_leakage( int nt,
          int np,
          int *teb,
          long ncoeff,
          double *ecoeff,
          double *bcoeff,
          long multipole,
          double *matrix)
{
  int flag=0;
  long nk=6;
  long i;
  double xl=(double)multipole;
  double *matrix1;

  /* treat one pair of detector at a time */
  /* generate polynomials in l from their coefficients */
  double e = polynomial(ncoeff, ecoeff, xl);
  double b = polynomial(ncoeff, bcoeff, xl);

  /* build matrix */
  matrix1 = (double *)malloc(nk*nk*sizeof(double));
  compute_cl_t2pmatrix(e, b, nk, matrix1);

  /* copy to final matrix */
  for (i=0; i<nk*nk; i++){ 
    matrix[i] = matrix1[i];
  }

  return flag;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

/* test program */
int main() {

  int nt=1;
  int np=1;
  int *teb;
  long ncoeff=5;
  double *ecoeff;
  double *bcoeff;
  long l;
  double *matrix;

  long nk=6;
  long i,j;

  teb = (int *)malloc(3*sizeof(int));
  teb[0]=1;
  teb[1]=1;
  teb[2]=1;
  ecoeff = (double *)malloc(ncoeff*sizeof(double));
  bcoeff = (double *)malloc(ncoeff*sizeof(double));
  ecoeff[0] = 0; ecoeff[1] = 0 ; ecoeff[2] = 1. ; ecoeff[3]=0 ; ecoeff[4]=0.01 ;
  bcoeff[0] = 0; bcoeff[1] = 0 ; bcoeff[2] = -2. ; bcoeff[3]=0 ; bcoeff[4]=0 ;
  l = 3 ; 
    

  matrix = (double *)malloc(nk*nk*sizeof(double));

  t2p_beam_leakage(nt, np, teb, ncoeff, ecoeff, bcoeff, l, matrix);

  for (i=0; i<nk; i++){
    for (j=0; j<nk; j++){
      printf("%10.5f ",matrix[i*nk+j]);
    }
    printf("\n");
  }

} /* end of main */
