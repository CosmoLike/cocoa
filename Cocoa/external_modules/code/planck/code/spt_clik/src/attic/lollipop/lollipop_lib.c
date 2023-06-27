#include "lollipop_lib.h"

#include <unistd.h>



/* Init routine               */
/* fill the structure dataset */
int Lollipop_Init( char *DataFile, char *FiducialFile, int Lmin, int Lmax, dataset *data)
{
  int nel, n, i, j;
  double *l, *cl, *na, *nb;
  double var1, var2, var_l, var_tt, var_ee, var_bb, var_te;
  FILE *fp=NULL;

  nel = (Lmax-Lmin+1);

  //read data file
  fp = fopen( DataFile, "r");
  if( fp == NULL) {
    printf( "File not found : %s\n", DataFile);
    return( -1);
  }

  n=0;
  while( fscanf(fp,"%lf  %lf  %lf", &var_l, &var1, &var2) == 3){
    if( var_l < Lmin ||var_l > Lmax) continue;
    data[n].cldata = var1;
    data[n].offset = var2;
    n++;
  }
  fclose( fp);
  if( n != nel) {
    printf( "Wrong number of l in datafile (%d)\n", n);
    exit(-1);
  }

  //read fiducial file
  fp = fopen( FiducialFile, "r");
  if( fp == NULL) {
    printf( "File not found : %s\n", FiducialFile);
    return( -1);
  }

  n=0;
  while( fscanf(fp,"%lf  %lf  %lf  %lf  %lf", &var_l, &var_tt, &var_ee, &var_bb, &var_te) == 5){
    if( var_l < Lmin ||var_l > Lmax) continue;
    data[n].fiducial = var_ee;
    n++;
  }
  fclose( fp);
  if( n != nel) {
    printf( "Wrong number of l in fiducialfile (%d)\n", n);
    exit(-1);
  }


/*   for( n=0; n<nel; n++) */
/*     printf( "%d %lf %lf %lf\n", Lmin+n, data[n].cldata, data[n].offset, data[n].fiducial); */


  return( 0);
}


/* Init routine     */
/* fill the cov mat */
int Lollipop_Cov( char *CovFile, int nel, double *invcov)
{
  int i, j;
  FILE *fp=NULL;

  fp = fopen( CovFile, "r");
  if( fp == NULL) {
    printf( "File not found : %s\n", CovFile);
    return( -1);
  }

  for(i = 0; i < nel; i++) {
    for(j = 0; j < nel; j++) {
      if (!fscanf(fp, "%lf", &invcov[i*nel+j])) break;
/*       printf("%lf\t", invcov[i*nel+j]); */
    }
/*     printf( "\n"); fflush(stdout); */

  }
  fclose( fp);

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
				  double *invcov,
				  int Lmin, int Lmax)
//cl in muK2
{
  int nel, i, c=0;
  double chi2=0;
  double *model=NULL;

  nel = Lmax-Lmin+1;

  model = (double *) malloc( nel*sizeof(double));
  if( model == NULL) exit(0);

  for( i=0; i<=Lmax; i++) {
    if( l[i] < Lmin || l[i] > Lmax) continue;
    model[c] = clee[i];
/*     printf( "%d l=%d : %lf\n", c, l[i], clee[i]); */
    c++;
  }
  if( c != nel) exit(-1);

  //Likelihood EE
  chi2 = Lollipop_oHL( nel, model, data, invcov);

  //free
  free( model);

  return( chi2/2.);
  
}





double Lollipop_oHL( int nel, double *model, dataset *data, double *invcov)
{
  int i, j, n;
  double g;
  double *X=NULL;

  //alloc
  X = (double *) calloc( nel, sizeof(double));
  if( X == NULL) {
    printf( "LOLLIPOP: Allocation failed\n");
    exit( -1);
  }

  //change variable
  for( i=0; i<nel; i++) {
    if( (data[i].cldata+data[i].offset) > 0.) {
      g = Lollipop_ghl( (data[i].cldata+data[i].offset) / (model[i]+data[i].offset));
      
      X[i] = sqrt(data[i].fiducial+data[i].offset) * g * sqrt(data[i].fiducial+data[i].offset);
    }
  }
  
  //compute chi2
  double f = 0.;
  for( i=0; i<nel; i++)
    for( j=0; j<nel; j++)
      f += X[i] * invcov[i*nel+j] * X[j];

  //free
  free( X);
  
  return(f);

}


double Lollipop_ghl( double x)
{
  
  return( sqrt( 2*(x-log(x)-1)));

}


