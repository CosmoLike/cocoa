#include "lollipop_lib.h"


int main()
{
  int i, nel, Lmin, Lmax, dl, l;
  double fsky, lnL;
  dataset *data;
  double *invcov;

  Lmin=4;
  Lmax=13;
  nel = Lmax-Lmin+1;

  char DataFile[256]     = {0};
  char FiducialFile[256] = {0};
  char CovFile[256]      = {0};
  sprintf( DataFile,     "data_70GHzx100GHz_HL.dat");
  sprintf( FiducialFile, "fiducial_HL.dat");
  sprintf( CovFile,      "invcov_70GHzx100GHz_HL.dat");

  //alloc dataset
  data   = (dataset *) malloc(     nel*sizeof(dataset));
  invcov = (double  *) malloc( nel*nel*sizeof(double ));

  //Init datasets
  Lollipop_Init( DataFile, FiducialFile, Lmin, Lmax, data);
  Lollipop_Cov( CovFile, nel, invcov);
  
  //read cl
  double tt, ee, bb, te;
  FILE *fp = fopen( "fiducial_HL.dat", "r");
  unsigned int *ell  = (unsigned int *) malloc( (Lmax+1)*sizeof(unsigned int));
  double       *clee = (double       *) malloc( (Lmax+1)*sizeof(double      ));
  fprintf(stderr,"lmax %d\n",Lmax);
  for( i=0; i<=Lmax; i++) {
    fscanf(fp,"%d  %lf  %lf  %lf  %lf", &l, &tt, &ee, &bb, &te);
    ell[i]  = l;
    clee[i] = ee;  //in muK2
    fprintf(stderr,"%d %d %d\n",i,l,ell[i]);
  }
  fclose( fp);

  //compute Likelihood
  lnL = Lollipop_computeLikelihood( ell, NULL, NULL, clee, NULL, data, invcov, Lmin, Lmax);
  printf( "Likelihood = %lf (-66.513262 70GHzx100GHz or  -53.458413 100GHzx143GHz)\n", 2*lnL);

  //free
  free( data);
  free( invcov);

  return( 0);
}

