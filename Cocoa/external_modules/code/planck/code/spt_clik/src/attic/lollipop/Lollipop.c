#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


typedef struct {
  double l, cl, na, nb;
} dataset;

int Lollipop_Init( char *InFileName, int Lmin, int Lmax, int dl, dataset *data);
double Lollipop_H2( double y);
double Lollipop_H3( double y);
double Lollipop_EdgeWorthSeries( int dl, double fsky, double cl0, dataset data);
double Lollipop_computeLikelihood(const unsigned int *l,
				  double *cltt,
				  double *clte,
				  double *clee,
				  double *clbb,
				  dataset *data, 
				  int Lmin, int Lmax, double fsky, int dl);



int main()
{
  int i, nel, Lmin, Lmax, dl;
  double fsky, lnL;
  dataset *data;
  
  Lmin=2;
  Lmax=13;
  dl=2;
  fsky=0.77521502;

  char InFileName[256] = {0};
  sprintf( InFileName, "cross_100x143.dat");

  nel = (Lmax-Lmin+1)/dl;
  data = (dataset *) calloc( nel, sizeof(dataset));

  //Init datasets
  int status = Lollipop_Init( InFileName, Lmin, Lmax, dl, data);
  if( status !=0) {
    printf( "Error in Lollipop_Init\n");
    return(-1);
  }

  //read cl
  double tt, ee, bb, te;
  FILE *fp = fopen( "planck_beam_corbeam.bestfit_cl_zre8_r0.2.data", "r");
  unsigned int *ell  = (unsigned int *) malloc( (Lmax+1)*sizeof(unsigned int));
  double       *clee = (double       *) malloc( (Lmax+1)*sizeof(double      ));
  for( i=0; i<=Lmax; i++) {
    fscanf(fp,"%lf  %lf  %lf  %lf", &tt, &ee, &bb, &te);
    ell[i]  = i;
    clee[i] = ee*1e12;  //in muK2
  }
  fclose( fp);

  //compute Likelihood
  lnL = Lollipop_computeLikelihood( ell, NULL, NULL, clee, NULL, data, Lmin, Lmax, fsky, dl);
  printf( "Likelihood = %e\n", 2*lnL);

  return( 0);
}




/* Init routine               */
/* fill the structure dataset */
int Lollipop_Init( char *InFileName, int Lmin, int Lmax, int dl, dataset *data)
{
  int nel, n;
  double *l, *cl, *na, *nb;
  double var1, var2, var3, var4;
  FILE *fp=NULL;

  nel = (Lmax-Lmin+1)/dl;

  fp = fopen( InFileName, "r");
  if( fp == NULL) {
    printf( "File not found : %s\n", InFileName);
    return( -1);
  }

  n=0;
  while( fscanf(fp,"%lf  %lf  %lf  %lf", &var1, &var2, &var3, &var4) == 4){
    if( var1 < Lmin ||var1 > Lmax) continue;
    data[n/dl].l  += var1/dl;
    data[n/dl].cl += var2/dl;
    data[n/dl].na += var3/dl;
    data[n/dl].nb += var4/dl;
    n++;
  }
  fclose( fp);

  printf( "n=%d\n", n);

  for( n=0; n<nel; n++)
    printf( "%f %e %e %e\n", data[n].l, data[n].cl, data[n].na, data[n].nb);

  return( 0);
}



/* Compute Likelihood                                    */
/* for a given dataset and cls (in muK2) return the -lnL */
double Lollipop_computeLikelihood(const unsigned int *l,
				  double *cltt,
				  double *clte,
				  double *clee,
				  double *clbb,
				  dataset *data, 
				  int Lmin, int Lmax, double fsky, int dl)
//cl in muK2
{
  int i, j, b, nel;
  double bin, model;

  double lnL=0;

  nel = (Lmax-Lmin+1)/dl;

  //Likelihood EE
  for( i=0; i<Lmax; i+=dl) {
    bin = 0.;
    model = 0.;
    
    if( l[i] < Lmin || l[i] > Lmax) continue;
    
    for( j=i; j<i+dl; j++) {
      bin   += (double)(l[j])/dl;
      model += clee[j]/dl;
    }
    
    for( b=0; b<nel; b++) {
      if( data[b].l == bin)
	lnL += log( Lollipop_EdgeWorthSeries( dl, fsky, model, data[b]));
      
    }
  }

  //return -ln(L)
  return( -lnL);
  
}





double Lollipop_EdgeWorthSeries( int dl, double fsky, double cl0, dataset data)
{
  double na = data.na;
  double nb = data.nb;

  //get cumulants
  double K1 = cl0;
  double K2 =     ( 2.*cl0*cl0 +  2.*cl0*(na+nb) + 1.2*na*nb) / ((2.*data.l+1)*fsky*(double)(dl));
  double K3 = cl0*( 8.*cl0*cl0 + 15.*cl0*(na+nb) + 6.0*na*nb) / ((2.*data.l+1)*(2*data.l+1)*fsky*fsky*fsky*(double)(dl*dl));
/*   printf( "%e %e %e   fsky=%e  dl=%d\t", K1, K2, K3, fsky, dl); */

  K3 /= sqrt(K2)*sqrt(K2)*sqrt(K2);

  //write Edgeworth Serie Expansion
  double y = (data.cl - K1) / sqrt(K2);
  double g = 1./sqrt( 2.*M_PI*K2) * exp( -0.5*y*y );
//   printf( "y=%e g=%e \t", y, g);
  
  double f = g * ( 1. + K3/6. * Lollipop_H3(y));

  if( f < 0) f = 1e-3;

//   printf( "L[%f]=%e\n", data.l, f);

  return(f);

}

double Lollipop_H2( double y)
{
  return( y*y-1.);
}

double Lollipop_H3( double y)
{
  return( y*y*y - 3.*y);
}
