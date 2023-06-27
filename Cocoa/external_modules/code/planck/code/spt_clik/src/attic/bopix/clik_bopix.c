#include "clik.h"
#include "clik_helper.h"
#include <errno.h>
#include <string.h>

typedef struct {
  char tmpdir[800];
  } bopix;


void free_bopix(void **none) {
  bopix_extra_free_();
}

double bopix_lkl(void* none, double* pars, error **err) {
  double lkl;
  
  bopix_extra_lkl_(&lkl,pars);
  return lkl;
}

cmblkl* clik_bopix_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096];
  char *bopix_nthreads_env;
  int bopix_nthreads;
  int bok;
  cmblkl *cing;
  int mlmax;
  
  bopix_extra_only_one_(&bok);
  testErrorRet(bok!=0,-100,"Bopix already initialized",*err,__LINE__,NULL);
  
  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  // call bopix_init
  bopix_extra_parameter_init_();
  
  cldf_external_cleanup(directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  bopix_nthreads = clik_getenviron_numthread("BOPIX",1,err);
  forwardError(*err,__LINE__,NULL);
  bopix_extra_set_bopix_nthreads_(&bopix_nthreads);
  
  mlmax = ell[nell-1];
  bopix_extra_set_lmax_(&mlmax);

  cing = init_cmblkl(NULL, &bopix_lkl, 
                     &free_bopix,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);
  return cing;
}
