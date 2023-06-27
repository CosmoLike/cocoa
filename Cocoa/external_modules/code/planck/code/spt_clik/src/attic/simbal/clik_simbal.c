#include "pmc.h"
#include "cldf/cldf.h"
#include "clik_helper.h"

typedef struct {
  int lmin,nell,nstep;
  int clmin,clmax;
  double *tab;
  double step;
} simlow_data;

void free_simlow(void **pft) {
  simlow_data *ft;

  ft = *pft;
  free(ft->tab);
  free(ft);
}

double simlow_lkl(void* ot, double *pars, error **err) {
  simlow_data *ft;
  double res;
  int ell,xdfid;

  ft = ot;
  res = 0;
  for(ell=ft->clmin;ell<ft->clmax+1;ell++) {
    testErrorRetVA(pars[ell-ft->clmin]<0,-1233,"multipole %d <0 (got %g )",*err,__LINE__,-1e10,ell,pars[ell-ft->lmin]); 
    xdfid = (int) (pars[ell-ft->clmin]*ell*(ell+1)/2/M_PI/ft->step);
    testErrorRetVA(xdfid>ft->nstep,-1233,"multipole %d too large (got %g expected <%g)",*err,__LINE__,-1e10,ell,pars[ell-ft->lmin],ft->step*ft->nstep*2*M_PI/ell/(ell+1)); 
    res += ft->tab[xdfid+(ell-ft->lmin)*ft->nstep];
  }

  return res;
}



cmblkl* clik_simlow_init(cldf * df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  cmblkl *cing;  
  simlow_data *ft;
  int ndat;

  ft = malloc_err(sizeof(simlow_data),err);
  forwardError(*err,__LINE__,NULL);

  ft->lmin = cldf_readint(df,"lmin",err);
  forwardError(*err,__LINE__,NULL);

  ft->nell = cldf_readint(df,"nell",err);
  forwardError(*err,__LINE__,NULL);
  
  ft->nstep = cldf_readint(df,"nstep",err);
  forwardError(*err,__LINE__,NULL);

  ft->step = cldf_readfloat(df,"step",err);
  forwardError(*err,__LINE__,NULL);
  
  ft->clmin = ell[0];
  ft->clmax = ell[nell-1];
  testErrorRetVA(ft->clmax>ft->lmin+ft->nell+1,-1233,"lmax too large (got %d expected %d at most)",*err,__LINE__,NULL,ft->clmax,ft->lmin+ft->nell); 
  
  ndat = ft->nell*ft->nstep;  
  ft->tab = cldf_readfloatarray(df,"tab",&(ndat),err);
  forwardError(*err,__LINE__,NULL);  

  cing = init_cmblkl(ft, &simlow_lkl, 
                     &free_simlow,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);
  
  return cing;  
}


