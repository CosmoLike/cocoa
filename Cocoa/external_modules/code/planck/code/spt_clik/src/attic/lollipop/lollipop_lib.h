#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


typedef struct {
  double cldata;
  double offset;
  double fiducial;
} dataset;


int Alloc_Lollipop_Data( int nel, dataset *data);
int Dealloc_Lollipop_Data( dataset *data);
int Lollipop_Init( char *DataFile, char *FiducialFile, int Lmin, int Lmax, dataset *data);
int Lollipop_Cov( char *CovFile, int nel, double *invcov);
double Lollipop_computeLikelihood(const unsigned int *l,
				  double *cltt,
				  double *clte,
				  double *clee,
				  double *clbb,
				  dataset *data, double *invcov,
				  int Lmin, int Lmax);
double Lollipop_oHL( int nel, double* model, dataset *data, double* invcov);
double Lollipop_ghl( double x);
