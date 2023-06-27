
#include "clik.h"
#include "clik_helper.h"
#include "cldf/cfrd.h"


typedef struct {
  int pipin,pipout;
  int ntot;
  } momento;


double momento_lkl(void* pbic, double* pars, error **err) {
  double lkl;
  momento *bic;
  char cmd[500];
  bic = pbic;
  int i;


  for(i=0;i<bic->ntot;i++) {
    sprintf(cmd,"val: %.10g\n",pars[i]);
    //_DEBUGHERE_("send %s",cmd);
    cfrd_send(bic->pipin,cmd,-1,err);
    forwardError(*err,__LINE__,0);
  }
  while(1) {
    cfrd_read_err(bic->pipout,cmd,err);
    forwardError(*err,__LINE__,0);
    if (strcmp(cmd,"rep: READY")==0) {
      break;
    } else {
      fprintf(stderr,"From child '%s'\n",cmd);
    }
  }
  cfrd_read_err(bic->pipout,cmd,err);
  forwardError(*err,__LINE__,0);
  
  sscanf(cmd,"%lg",&lkl);
  return lkl;

}


void free_momento(void **pbic) {
  momento *bic;
  error *_err,**err;

  _err = NULL;
  err = &_err;
  bic = *pbic;
  cfrd_send(bic->pipin,"stop\n",-1,err);
  forwardError(*err,__LINE__,);
  close(bic->pipin);
  close(bic->pipout);
}


cmblkl* clik_momento_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096],pwd2[4096],cmd[4096];
  cmblkl *cing;
  int pipin[2],pipout[2];
  int childpid;
  momento * bic;
  int xdim;
  char *xnames_tot[300];
  char *nuisance;
  int i;
  
  bic = malloc_err(sizeof(momento),err);
  forwardError(*err,__LINE__,NULL);

  childpid = cfrd_startup(pipin,pipout,"momento_main",NULL,err);
  forwardError(*err,__LINE__,NULL);
  
  bic->pipin = pipin[1];
  bic->pipout = pipout[0];

  sprintf(cmd,"%s\n",df->root);

  cfrd_send(bic->pipin,cmd,-1,err);
  forwardError(*err,__LINE__,NULL);
  
  //_DEBUGHERE_("%s",df->root);
  cing = init_cmblkl(bic, &momento_lkl, 
                     &free_momento,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);

  while(1) {
    //_DEBUGHERE_("waiting for the child","");
    cfrd_read_err(bic->pipout,cmd,err);
    forwardError(*err,__LINE__,0);
    if (strcmp(cmd,"rep: READY")==0) {
      break;
    }
  }
  //_DEBUGHERE_("OK","");
  
  
  bic->ntot = (ell[nell-1]+1-ell[0])*has_cl[0] + (ell[nell-1]+1-ell[0])*has_cl[1] + (ell[nell-1]+1-ell[0])*has_cl[2] + (ell[nell-1]+1-ell[0])*has_cl[3] + (ell[nell-1]+1-ell[0])*has_cl[4] + (ell[nell-1]+1-ell[0])*has_cl[5];
  //_DEBUGHERE_("%d %d %d %d %d %d %d %d",ell[nell-1],bic->ntot,has_cl[0],has_cl[1],has_cl[2],has_cl[3],has_cl[4],has_cl[5]);
  return cing;
}


