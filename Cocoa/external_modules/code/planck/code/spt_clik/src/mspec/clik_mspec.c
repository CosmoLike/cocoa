
#include "clik.h"
#include "clik_helper.h"
#include "cldf/cfrd.h"


typedef struct {
  int pipin,pipout;
  int ntot;
  } mspec;


double mspec_lkl(void* pbic, double* pars, error **err) {
  double lkl;
  mspec *bic;
  char cmd[500];
  bic = pbic;
  int i;


  for(i=0;i<bic->ntot;i++) {
    sprintf(cmd,"%g\n",pars[i]);
    cfrd_send(bic->pipin,cmd,-1,err);
    forwardError(*err,__LINE__,0);
  }
  while(1) {
    cfrd_read_err(bic->pipout,cmd,err);
    forwardError(*err,__LINE__,0);
    if (strcmp(cmd,"READY")==0) {
      break;
    }
  }
  cfrd_read_err(bic->pipout,cmd,err);
  forwardError(*err,__LINE__,0);
  
  sscanf(cmd,"%lg",&lkl);
  return lkl;

}


void free_mspec(void **pbic) {
  mspec *bic;
  error *_err,**err;

  _err = NULL;
  err = &_err;
  bic = *pbic;
  cfrd_send(bic->pipin,"stop\n",-1,err);
  forwardError(*err,__LINE__,);
  close(bic->pipin);
  close(bic->pipout);
}


cmblkl* clik_mspec_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096],pwd2[4096],cmd[4096];
  cmblkl *cing;
  int pipin[2],pipout[2];
  int childpid;
  mspec * bic;
  int xdim;
  char *xnames_tot[300];
  char *nuisance;
  int i;
  
  bic = malloc_err(sizeof(mspec),err);
  forwardError(*err,__LINE__,NULL);

  childpid = cfrd_startup(pipin,pipout,"mspec_main",NULL,err);
  forwardError(*err,__LINE__,NULL);
    
  bic->pipin = pipin[1];
  bic->pipout = pipout[0];

  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
 
  testErrorRetVA(getcwd(pwd2,4096)==NULL,-101010,"can't get cwd name (cause = '%s')",*err,__LINE__,NULL,strerror(errno));
  sprintf(cmd,"%s/mspec_data\n",pwd2);
  cfrd_send(bic->pipin,cmd,-1,err);
  forwardError(*err,__LINE__,NULL);
  
  cldf_external_cleanup(directory_name,pwd,err);  
  forwardError(*err,__LINE__,NULL);
  
  xdim = cldf_readint(df,"n_nuisance",err);
  forwardError(*err,__LINE__,NULL);
  

  nuisance = cldf_readstr(df,"nuisance",NULL,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<xdim;i++) {
    xnames_tot[i] = &(nuisance[i*256]);
  }  
  
  cing = init_cmblkl(bic, &mspec_lkl, 
                     &free_mspec,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);

  while(1) {
    //_DEBUGHERE_("waiting for the child","");
    cfrd_read_err(bic->pipout,cmd,err);
    forwardError(*err,__LINE__,0);
    if (strcmp(cmd,"READY")==0) {
      break;
    }
  }
  //_DEBUGHERE_("OK","");
  cmblkl_set_names(cing, xnames_tot,err);
  forwardError(*err,__LINE__,NULL);
  
  free(nuisance);

  bic->ntot = xdim + (ell[nell-1]+1)*has_cl[0] + (ell[nell-1]+1)*has_cl[1] + (ell[nell-1]+1)*has_cl[2] + (ell[nell-1]+1)*has_cl[3] + (ell[nell-1]+1)*has_cl[4] + (ell[nell-1]+1)*has_cl[5];
  
  return cing;
}


