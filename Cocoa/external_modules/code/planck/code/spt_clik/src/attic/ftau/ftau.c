#include "pmc.h"
#include "cldf/cldf.h"
#include "clik_helper.h"

typedef struct {
  int nt;
  double *tv, *Pv;
  int vps;
} ftau_data;

void free_ftau(void **pft) {
  ftau_data *ft;

  ft = *pft;
  free(ft->tv);
  free(ft->Pv);
  free(ft);
}

double ftau_lkl(void* ot, double *pars, error **err) {
  ftau_data *ft;
  double res;
  double w,d,l,theta;
  int i;

  ft = ot;
  res = -1e10;
  w = pars[ft->vps];

  if (w<ft->tv[0] || w>ft->tv[ft->nt-1]) {
    return res;
  }
  for(i=1;i<ft->nt;i++) {
    if (w<=ft->tv[i]) {
      break;
    }
  }
  l = ft->tv[i]-ft->tv[i-1];
  d = w-ft->tv[i-1];
  theta = d/l;
  res = ft->Pv[i-1] * (1-theta) + ft->Pv[i] * theta;
  return res;
}


cmblkl* clik_ftau_init(cldf * df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
	cmblkl *cing;  
  ftau_data *ft;
  char* tauname;

  ft = malloc_err(sizeof(ftau_data),err);
  forwardError(*err,__LINE__,NULL);

  ft->nt = cldf_readint(df,"n_tab",err);
  forwardError(*err,__LINE__,NULL);

  ft->tv = cldf_readfloatarray(df,"tau",&(ft->nt),err);
  forwardError(*err,__LINE__,NULL);  

  ft->Pv = cldf_readfloatarray(df,"Ptau",&(ft->nt),err);
  forwardError(*err,__LINE__,NULL);  

  ft->vps = (ell[nell-1]+1)*has_cl[0] + (ell[nell-1]+1)*has_cl[1] + (ell[nell-1]+1)*has_cl[2] + (ell[nell-1]+1)*has_cl[3] + (ell[nell-1]+1)*has_cl[4] + (ell[nell-1]+1)*has_cl[5];

  cing = init_cmblkl(ft, &ftau_lkl, 
                     &free_ftau,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,1,err);
  forwardError(*err,__LINE__,NULL);
  
  tauname = cldf_readstr(df,"tau_name",NULL,err);
  cmblkl_set_names(cing, &tauname,err);
  forwardError(*err,__LINE__,NULL);
  
  free(tauname);
  return cing;  
}

