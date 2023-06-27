#include "clik_helper.h"
#include "smica_ext.h"

typedef SmicaComp * clik_smica_comp_init_func(cldf * df,int nb, int mT,int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err);

#if 0 
cmblkl* clik_lowly_init(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  powly *self;
  cmblkl *cing;
  herr_t hstat;
  int neff,nnn,cli;
  
  self = malloc_err(sizeof(powly), err);
  forwardError(*err,__LINE__,NULL);
  SET_PRINT_STAT(self);
  
  
  hstat = H5LTget_attribute_int( group_id, ".", "neff",  &neff);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read neff in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  self->neff = neff;
  
  self->nell = nell;
  self->ell=malloc_err(sizeof(int)*nell,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(self->ell,ell,sizeof(int)*nell);
  self->ncl_fid = lowly_get_offset_cl(has_cl,self->offset_cl,self->nell);  
  self->nmode = lowly_nmode(self->nell,self->ell);
  
  self->offset_mode[0]=(-1,-1,-1);
  self->tot_mode = 0;
  if ((self->offset_cl[0]!=-1) || (self->offset_cl[3]!=-1) || (self->offset_cl[4]!=-1)) {
    self->offset_mode[0]=self->tot_mode; //T
    self->tot_mode += self->nmode;
  }
  if ((self->offset_cl[1]!=-1) || (self->offset_cl[3]!=-1) || (self->offset_cl[5]!=-1)) {
    self->offset_mode[1]=self->tot_mode; //E
    self->tot_mode += self->nmode;
  }
  if ((self->offset_cl[2]!=-1) || (self->offset_cl[4]!=-1) || (self->offset_cl[5]!=-1)) {
    self->offset_mode[1]=self->tot_mode; //B
    self->tot_mode += self->nmode;
  }
  
  self->cl_fid = hdf5_double_datarray(group_id, cur_lkl,"cl_fid",&self->ncl_fid,err);
  forwardError(*err,__LINE__,NULL);

  self->buffer = malloc_err(sizeof(double)*(_SZT_(self->tot_mode)*_SZT_(self->neff) + _SZT_(self->neff)*_SZT_(self->neff) + self->neff*2),err);
  forwardError(*err,__LINE__,NULL);
  self->H = self->buffer;
  self->a_bar = self->H + self->tot_mode*self->neff;
  self->tmp = self->a_bar + self->neff;
  self->Re = self->tmp + self->neff;
  
  nnn = self->tot_mode*self->neff;
  self->H = hdf5_double_datarray(group_id, cur_lkl,"H",&nnn,err);
  forwardError(*err,__LINE__,NULL);
  
  nnn = self->neff;
  self->a_bar = hdf5_double_datarray(group_id, cur_lkl,"a_bar",&nnn,err);
  forwardError(*err,__LINE__,NULL);
  
  self->time_build = 0;
  self->time_tot = 0;
  self->time_chol = 0;
  self->n_tests = 0;
  
  self->eMat = NULL;
  self->eHt = NULL;
  
  self->n_sym = 0;
  for(cli=0;cli<3;cli++) {
    self->n_sym += self->offset_cl[cli]==-1 ? 0:1;
  }
  self->n_cross = 0;
  for(cli=0;cli<3;cli++) {
    self->n_cross += self->offset_cl[cli+3]==-1 ? 0:1;
  }
  self ->nstore[0]=0;
  self ->nstore[1]=0;
  self ->tot_store=0;
  self->storage = NULL;
  
  cing = init_cmblkl(self, &powly_lkl, 
                     &free_powly,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);
  return cing;  
}

void al_bl_cor(hid_t group_id, char* cur_lkl, int *_n, double **_al, double **_bl, double **_nl, double **_cor,int *_isdiag,error **err) {
  double *al, *bl, *corr,*nl;
  int n,isdiag;
  herr_t hstat;
  
  n = *_n;
  // get al
  al = hdf5_double_datarray(group_id, cur_lkl,"al",&n,err);
  forwardError(*err,__LINE__,);

  // get bl
  bl = hdf5_double_datarray(group_id, cur_lkl,"bl",&n,err);
  forwardError(*err,__LINE__,);

  // get nl
  nl = NULL;
  hstat = H5LTfind_dataset(group_id, "nl");
  if (hstat==1) {
    nl = hdf5_double_datarray(group_id, cur_lkl,"nl",&n,err);
    forwardError(*err,__LINE__,);
  }

  // get cor
  corr = NULL;
  hstat = H5LTfind_dataset(group_id, "cor");
  if (hstat==1) {
    int ncor;
    ncor=-1;
    corr = hdf5_double_datarray(group_id, cur_lkl,"cor",&ncor,err);
    forwardError(*err,__LINE__,);
    testErrorRetVA((ncor!=n) && (ncor!=n*n),hdf5_base,"Bad size for %s in %s (got %d expected %d or %d)",*err,__LINE__,,"cor",cur_lkl,ncor,n,n*n);
    isdiag = 0;
    if (ncor==n) {
      isdiag = 1;
    } 
  }
  
  *_al = al;
  *_bl = bl;
  *_nl = nl;
  *_cor = corr;
  *_isdiag = isdiag;
  
}

cmblkl* clik_ivg_init(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  tease *ing;
  cmblkl *cing;
  int i,isdiag;
  double *al, *bl, *nl,*corr;
  int nlcst;
  double zero;
  int n;
  int *Ml;
  hsize_t ndum;
  H5T_class_t dum;
  size_t ddum;
  herr_t hstat;
  
  zero = 0;
  n=nell;
  Ml = ell;
  if (nbins!=0) {
    n=nbins;
    Ml = NULL;
  }
  
  al_bl_cor(group_id,cur_lkl, &n, &al, &bl, &nl, &corr, &isdiag,err);
  forwardError(*err,__LINE__,NULL);
  
  nlcst = 0;
  if (nl == NULL) {
    nl = &zero;
    nlcst = 1;    
  }
    
  ing = tease_init(n,Ml,al, bl, nl, nlcst, corr, isdiag, err);
  forwardError(*err,__LINE__,NULL);
  
  free(al);
  free(bl);
  if (nl!=&zero) {
    free(nl);
  }
  if (corr!=NULL) {
    free(corr);
  }
  
  cing = init_cmblkl(ing, &tease_log_pdf, 
                     &tease_free,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);
  return cing;  
}


cmblkl* clik_gauss_init(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  gausslkl *ing;
  cmblkl *cing;
  int i,isdiag;
  double *al, *bl, *nl,*corr;
  int nlcst;
  double zero;
  int n;
  int *Ml;
  hsize_t ndum;
  H5T_class_t dum;
  size_t ddum;  
  herr_t hstat;
  
  zero = 0;
  n=nell;
  Ml = ell;
  if (nbins!=0) {
    n=nbins;
    Ml = NULL;
  }
  
  al_bl_cor(group_id,cur_lkl, &n, &al, &bl, &nl, &corr, &isdiag,err);
  forwardError(*err,__LINE__,NULL);
  
  nlcst = 0;
  if (nl == NULL) {
    nl = &zero;
    nlcst = 1;    
  }
  
  ing = gausslkl_init(n,Ml,al, bl, nl, nlcst, corr, isdiag, err);
  forwardError(*err,__LINE__,NULL);
  
  free(al);
  free(bl);
  if (nl!=&zero) {
    free(nl);
  }
  if (corr!=NULL) {
    free(corr);
  }
  
  cing = init_cmblkl(ing, &gausslkl_log_pdf, 
                     &tease_free,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);
  return cing;  
}
#endif

cmblkl* clik_smica_init(cldf * df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  Smica *ing;
  cmblkl *cing;  
  double zero;
  int n;
  int *Ml;
  double *wq,*rq_hat,*rq_0;
  int m;
  SmicaComp **SCs;
  int ic;
  int mT,mP,nc;
  double *A_cmb;
  Smica *smic;
  int ncl,icl,nb;
  size_t ddum;
  int cnt,xdim;
  char **xnames;
  parname *xnames_buf;
  int hk;
  int mT_plus_mP;

  zero = 0;

  ncl = 0;
  for(icl=0;icl<6;icl++) {
    if (has_cl[icl]==1) {
      ncl++;
    }
  }
  
  nb = nell;
  if (nbins!=0) {
    nb = nbins/ncl;  
  }  
    
  // try to read the bin weights
  wq = NULL;
  hk = cldf_haskey(df,"wq",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    wq = cldf_readfloatarray(df,"wq",&nb,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  
  // read the number of channels
  mT = cldf_readint(df,"m_channel_T",err);
  forwardError(*err,__LINE__,NULL);
  mP = cldf_readint(df,"m_channel_P",err);
  forwardError(*err,__LINE__,NULL);
  
  m = mtot(mT, mP, has_cl);
  mT_plus_mP = mT+mP;

  // read rq_hat
  int nrq;
  
  nrq = nb*m*m;
  rq_hat = cldf_readfloatarray(df,"Rq_hat",&nrq,err);
  forwardError(*err,__LINE__,NULL);  
  
  // try to read rq_0
  rq_0 = NULL;
  hk = cldf_haskey(df,"Rq_0",err);
  forwardError(*err,__LINE__,NULL);  
  if (hk==1) {
    rq_0 = cldf_readfloatarray(df,"Rq_0",&nrq,err);
    forwardError(*err,__LINE__,NULL);  
  }
  
  // how many components ?
  nc = cldf_readint(df,"n_component",err);
  forwardError(*err,__LINE__,NULL);
  
  SCs = malloc_err(sizeof(SmicaComp*) * nc,err);
  forwardError(*err,__LINE__,NULL);
  
  // now deal with the CMB component
  // read A_cmb
  A_cmb = cldf_readfloatarray(df,"A_cmb",&mT_plus_mP,err);
  forwardError(*err,__LINE__,NULL);    

  // init cmb comp
  SCs[0] = comp_CMB_init(nb, mT,mP, has_cl, A_cmb, err);
  forwardError(*err,__LINE__,NULL);    
  SC_set_compname(SCs[0],"CMB");
  
  free(A_cmb);
  
  // deal with other components
  xdim = 0;
  for(ic=1;ic<nc;ic++) {
    char cur_cmp[256];
    char cur_cmp_tot[256];
    clik_smica_comp_init_func *smica_dl_init;
    void* dlhandle;
    char init_func_name[256];
    parname comp_type;
    cldf *comp_df;

#ifdef HAS_RTLD_DEFAULT 
    dlhandle = RTLD_DEFAULT;
#else
    dlhandle = NULL;
#endif
    
    SCs[ic] = NULL;
    
    sprintf(cur_cmp,"component_%d",ic);
    comp_df = cldf_openchild(df,cur_cmp,err);
    forwardError(*err,__LINE__,NULL);

    // get type
    char * cmt;

    cmt = cldf_readstr(comp_df,"component_type",NULL,err);
    forwardError(*err,__LINE__,NULL);
    sprintf(comp_type,"%s",cmt);
    free(cmt);

    sprintf(init_func_name,"clik_smica_comp_%s_init",comp_type);
    smica_dl_init = dlsym(dlhandle,init_func_name);
    testErrorRetVA(smica_dl_init==NULL,-1111,"Cannot initialize smica component type %s from %s dl error : %s",*err,__LINE__,NULL,comp_type,df->root,dlerror()); 
 
    SCs[ic] = smica_dl_init(comp_df,nb,mT,mP, nell, ell, has_cl, unit, wl, bins,nb,err);
 
    forwardError(*err,__LINE__,NULL);
    sprintf(comp_type,"%s_%d",comp_type,ic);
    SC_set_compname(SCs[ic],comp_type);
      
    cldf_close(&comp_df);
    
    xdim += SCs[ic]->ndim;
 
  
  }
  
  // deal with names and xdims
  if (xdim!=0) {
    xnames = malloc_err(sizeof(char*)*xdim,err);
    forwardError(*err,__LINE__,NULL);
  
    xnames_buf = malloc_err(sizeof(parname)*xdim,err);
    forwardError(*err,__LINE__,NULL);
    cnt = 0;
    for(ic=1;ic<nc;ic++) {
      int ix;
      for(ix=0;ix<SCs[ic]->ndim;ix++) {
        if (SCs[ic]->names!=NULL) {
          xnames[cnt]=(char*)&(SCs[ic]->names[ix]);
        } else {
          sprintf(xnames_buf[cnt],"SMICA_COMP_%d_%d",ic,ix);
          xnames[cnt]=(char*)&(xnames_buf[cnt]);
        }
        cnt++;
      }
    }
  }
    
  smic = Smica_init(nb, wq, m, rq_hat, rq_0, nc, SCs,err);
  forwardError(*err,__LINE__,NULL);
    
  // deal with criterion
  hk = cldf_haskey(df,"criterion",err);
  forwardError(*err,__LINE__,NULL);
  
  if (hk == 1) {
    char *crit_name;
    
    crit_name = cldf_readstr(df,"criterion",NULL,err);
    forwardError(*err,__LINE__,NULL);
    if(strcmp(crit_name,"classic")==0) {
      //nothing to do here, we are fine
    } else if(strcmp(crit_name,"mgauss")==0) {
      int select_size,sz2;
      double *sigma;
      int *ordering,*qmin,*qmax;
      int hk,lsz;

      select_size = cldf_readint(df,"mgauss/select_size",err);
      forwardError(*err,__LINE__,NULL);
      
      sz2 = select_size * select_size;
      sigma = cldf_readfloatarray(df,"mgauss/sigma",&sz2,err);
      forwardError(*err,__LINE__,NULL);
      
      hk = cldf_haskey(df,"mgauss/ordering",err);
      forwardError(*err,__LINE__,NULL);
      if (hk == 1) { 
        int i;
        lsz = m*(m+1);
        ordering = cldf_readintarray(df,"mgauss/ordering",&lsz,err);
        forwardError(*err,__LINE__,NULL);
        //_DEBUGHERE_("","");
      } else {
        int iv,jv,nv;
        ordering = malloc_err(sizeof(int)*m*(m+1),err);
        forwardError(*err,__LINE__,NULL);
        for(iv=0;iv<m;iv++) {
          for(jv=iv;jv<m;jv++) {
            ordering[nv] = iv;
            ordering[nv+1] = jv;
            nv+=2;
          }
        }
      }

      lsz = m*m;
      qmin = cldf_readintarray(df,"mgauss/qmin",&lsz,err);
      forwardError(*err,__LINE__,NULL);
      qmax = cldf_readintarray(df,"mgauss/qmax",&lsz,err);
      forwardError(*err,__LINE__,NULL);

      smica_set_crit_mgauss(smic,select_size,sigma,ordering,qmin,qmax,err);
      forwardError(*err,__LINE__,NULL);      
      free(sigma);
      free(ordering);
      free(qmin);
      free(qmax);

    } else if(strcmp(crit_name,"gauss")==0) {
      double *quad_crit;
      int *mask;
      int *ordering;
      int nqu;

      nqu = -1;
      quad_crit = cldf_readfloatarray(df,"criterion_gauss_mat",&nqu,err);
      forwardError(*err,__LINE__,NULL);
      mask = NULL;
      hk = cldf_haskey(df,"criterion_gauss_mask",err);
      forwardError(*err,__LINE__,NULL);
      if (hk == 1) { 
        int i;
        nqu = nb*m*m;
        mask = cldf_readintarray(df,"criterion_gauss_mask",&nqu,err);
        forwardError(*err,__LINE__,NULL);
      }
      ordering = NULL;
      hk = cldf_haskey(df,"criterion_gauss_ordering",err);
      forwardError(*err,__LINE__,NULL);
      if (hk == 1) { 
        int i;
        nqu = m*(m+1);
        ordering = cldf_readintarray(df,"criterion_gauss_ordering",&nqu,err);
        forwardError(*err,__LINE__,NULL);
        //_DEBUGHERE_("","");
      }

      smica_set_crit_gauss(smic, quad_crit,mask,ordering,err);
      forwardError(*err,__LINE__,NULL);
      
      if (mask!=NULL) {
        free(mask);  
      }
      if (ordering!=NULL) {
        free(ordering);  
      }
      free(quad_crit);
    } else if(strcmp(crit_name,"eig")==0) {
      double *nrm;
    
      nrm = NULL;
      //hstat = H5LTfind_dataset(group_id, "criterion_eig_norm");
      hk = cldf_haskey(df,"criterion_eig_norm",err);
      forwardError(*err,__LINE__,NULL);
      if (hk == 1) { 
        //nrm = hdf5_double_datarray(group_id, cur_lkl,"criterion_eig_norm",&nb,err);
        nrm = cldf_readfloatarray(df,"criterion_eig_norm",&nb,err);
        forwardError(*err,__LINE__,NULL);
      }
      smica_set_crit_eig(smic, nrm, err);
      forwardError(*err,__LINE__,NULL);
      if (nrm!=NULL) {
        free(nrm);    
      }
    } else if(strcmp(crit_name,"quad")==0) {
      double *fid;
      int nqu;
      int *mask;
      fid = NULL;
      hk = cldf_haskey(df,"criterion_quad_mat",err);
      forwardError(*err,__LINE__,NULL);
      //hstat = H5LTfind_dataset(group_id, "criterion_quad_mat");
      if (hk == 1) { 
        nqu = -1;
        //fid = hdf5_double_datarray(group_id, cur_lkl,"criterion_quad_mat",&nqu,err);
        fid = cldf_readfloatarray(df,"criterion_quad_mat",&nqu,err);
        forwardError(*err,__LINE__,NULL);
      }
      mask = NULL;
      //hstat = H5LTfind_dataset(group_id, "criterion_quad_mask");
      hk = cldf_haskey(df,"criterion_quad_mask",err);
      forwardError(*err,__LINE__,NULL);
      if (hk == 1) { 
        int i;
        nqu = m*m;
        //mask = hdf5_int_datarray(group_id, cur_lkl,"criterion_quad_mask",&nqu,err);
        mask = cldf_readintarray(df,"criterion_quad_mask",&nqu,err);
        forwardError(*err,__LINE__,NULL);
      }
      smica_set_crit_quad(smic, fid,mask,err);
      forwardError(*err,__LINE__,NULL);
      if (fid!=NULL) {
        free(fid);  
      }
      
      if (mask!=NULL) {
        free(mask);  
      }
    } else {
      testErrorRetVA(1==1,hdf5_base,"does not understand criterion '%s' in %s",*err,__LINE__,NULL,crit_name,df->root);
    }
    free(crit_name);
  }
  
  free(rq_hat);
  
  if (rq_0!=NULL) {
    free(rq_0);
  }    
  
  if (wq!=NULL) {
    free(wq);
  }  
  
  cing = init_cmblkl(smic, &Smica_lkl, 
                     &free_Smica,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);
  if (xdim!=0) {
    cmblkl_set_names(cing, xnames,err);
    forwardError(*err,__LINE__,NULL);
  
    free(xnames);
    free(xnames_buf);
  }
  
  
  
  return cing;  
}


SmicaComp * clik_smica_comp_diag_init(cldf *df,int nb, int mT, int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {
  SmicaComp *SC; 
  int m;

  m = mtot(mT,mP,has_cl);
  // try to read A

  SC = comp_diag_init(nb, m, err);
  forwardError(*err,__LINE__,NULL);    
  return SC;
}

SmicaComp * clik_smica_comp_cst_init(cldf *df,int nb, int mT,int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {
  SmicaComp *SC; 
  double *rq_0;
  int tt;
  int m;

  m = mtot(mT,mP,has_cl);
  
  tt = -1;
  //rq_0 =  hdf5_double_datarray(comp_id,cur_lkl,"Rq_0",&tt,err);
  rq_0 =  cldf_readfloatarray(df,"Rq_0",&tt,err);
  forwardError(*err,__LINE__,NULL);    
  testErrorRetVA(tt != m*m*nb,-22345,"%s:cst component does not have the correct number of data (expected %d got %d)",*err,__LINE__,NULL,df->root,m*m*nb,tt)

  SC = comp_cst_init(nb, m, rq_0, err);
  forwardError(*err,__LINE__,NULL);    

  free(rq_0);
  return SC;
}


SmicaComp * clik_smica_comp_gcal2_init(cldf *df,int nb, int mT,int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {
  SmicaComp *SC; 
  int npar,i;
  int *im,*jm;
  double *tpl;
  int tot_tpl;
  char *bnames, **xnames;
  int m;
  int dsz;

  m = mtot(mT,mP,has_cl);
  
  npar = -1;
  im =  cldf_readintarray(df,"im",&npar,err);
  forwardError(*err,__LINE__,NULL);    
  jm =  cldf_readintarray(df,"jm",&npar,err);
  forwardError(*err,__LINE__,NULL);    

  tot_tpl = npar * nb;
  tpl = cldf_readfloatarray(df,"tpl",&tot_tpl,err);
  forwardError(*err,__LINE__,NULL);    

  SC = comp_gcal2_init(nb, m, npar, im, jm,  tpl, err);
  forwardError(*err,__LINE__,NULL);    

  free(im);
  free(jm);
  free(tpl);

  dsz= -1;
  bnames = cldf_readstr(df,"names",&dsz, err);
  forwardError(*err,__LINE__,NULL); 
  xnames = malloc_err(sizeof(char*)*npar,err);
  for(i=0;i<npar;i++) {
    xnames[i] =&(bnames[i*256]);
  } 
  SC_setnames(SC, xnames, err);
  forwardError(*err,__LINE__,NULL);
  free(xnames); 
  free(bnames); 
  
  return SC;
}



SmicaComp * clik_smica_comp_gcal_log_init(cldf *df,int nb, int mT,int mP,int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  int *ngcal;
  double *gcaltpl;
  SmicaComp *SC;
  int mm,im,tt,ig,to;
  char **xnames,*bnames;
  int binned;
  int hk;
  int m;

  m = mtot(mT,mP,has_cl);
  
  mm = m;
  //ngcal = hdf5_int_attarray(comp_id,cur_lkl,"ngcal",&mm,err);
  ngcal = cldf_readintarray(df,"ngcal",&mm,err);
  forwardError(*err,__LINE__,NULL);    

  tt = 0;
  for(im=0;im<m;im++) {
    testErrorRetVA(ngcal[im]<0,hdf5_base,"%s: ngcal[%d] does make any sence (got %d)",*err,__LINE__,NULL,df->root,im,ngcal[im]);
    tt += ngcal[im];
  }
  
  //hstat = H5LTget_attribute_int(comp_id, ".", "binned",  &binned);
  //testErrorRetVA(hstat<0,hdf5_base,"cannot read binned in component %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  binned = cldf_readint(df,"binned",err);
  forwardError(*err,__LINE__,NULL);

  if (binned!=0) {
    binned = 1;
  }

  to = tt * (binned*nbins + (1-binned)*nell);
  gcaltpl =  cldf_readfloatarray(df,"gcaltpl",&to,err);
  //gcaltpl =  hdf5_double_datarray(comp_id,cur_lkl,"gcaltpl",&to,err);
  forwardError(*err,__LINE__,NULL);    
 
  SC = comp_gcal_log_init(nb,m, ngcal, gcaltpl,nell*(1-binned),bins,err);
  forwardError(*err,__LINE__,NULL);    

  //hstat = H5LTfind_attribute(comp_id, "names");
  hk = cldf_haskey(df,"names",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    int dz;
    dz =-1;
    //bnames = hdf5_char_attarray(comp_id,cur_lkl,"names",&dz, err);
    bnames = cldf_readstr(df,"names",&dz, err);
    forwardError(*err,__LINE__,NULL); 
  } else {
    int ii;
    bnames = malloc_err(sizeof(char)*256*tt,err);
    forwardError(*err,__LINE__,NULL); 
    ii=0;
    for(im=0;im<m;im++) {
      for(ig=0;ig<ngcal[im];ig++) {
        sprintf(&(bnames[ii*256]),"gcal_%d_%d",im,ig);
        ii++;
      }
    }
  }
  
  xnames = malloc_err(sizeof(char*)*tt,err);
  for(im=0;im<tt;im++) {
    xnames[im] =&(bnames[im*256]);
  } 
  SC_setnames(SC, xnames, err);
  forwardError(*err,__LINE__,NULL);
  free(xnames); 
  free(bnames); 
  free(ngcal);
  free(gcaltpl);
  
  return SC;
}

SmicaComp * clik_smica_comp_gcal_lin_init(cldf *df,int nb, int mT,int mP,int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  int *ngcal;
  double *gcaltpl;
  SmicaComp *SC;
  int mm,im,tt,ig,to;
  char **xnames,*bnames;
  int binned;
  int hk;
  int m;

  m = mtot(mT,mP,has_cl);
  
  //ngcal = hdf5_int_attarray(comp_id,cur_lkl,"ngcal",&mm,err);
  ngcal = cldf_readintarray(df,"ngcal",&mm,err);
  forwardError(*err,__LINE__,NULL);    

  tt = 0;
  for(im=0;im<m;im++) {
    testErrorRetVA(ngcal[im]<0,hdf5_base,"%s: ngcal[%d] does make any sence (got %d)",*err,__LINE__,NULL,df->root,im,ngcal[im]);
    tt += ngcal[im];
  }
  
  //hstat = H5LTget_attribute_int(comp_id, ".", "binned",  &binned);
  //testErrorRetVA(hstat<0,hdf5_base,"cannot read binned in component %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  binned = cldf_readint(df,"binned",err);
  forwardError(*err,__LINE__,NULL);

  if (binned!=0) {
    binned = 1;
  }

  to = tt * (binned*nbins + (1-binned)*nell);
  gcaltpl =  cldf_readfloatarray(df,"gcaltpl",&to,err);
  //gcaltpl =  hdf5_double_datarray(comp_id,cur_lkl,"gcaltpl",&to,err);
  forwardError(*err,__LINE__,NULL);    
 
  SC = comp_gcal_lin_init(nb,m, ngcal, gcaltpl,nell*(1-binned),bins,err);
  forwardError(*err,__LINE__,NULL);    

  //hstat = H5LTfind_attribute(comp_id, "names");
  hk = cldf_haskey(df,"names",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    int dz;
    dz =-1;
    //bnames = hdf5_char_attarray(comp_id,cur_lkl,"names",&dz, err);
    bnames = cldf_readstr(df,"names",&dz, err);
    forwardError(*err,__LINE__,NULL); 
  } else {
    int ii;
    bnames = malloc_err(sizeof(char)*256*tt,err);
    forwardError(*err,__LINE__,NULL); 
    ii=0;
    for(im=0;im<m;im++) {
      for(ig=0;ig<ngcal[im];ig++) {
        sprintf(&(bnames[ii*256]),"gcal_%d_%d",im,ig);
        ii++;
      }
    }
  }

  xnames = malloc_err(sizeof(char*)*tt,err);
  for(im=0;im<tt;im++) {
    xnames[im] =&(bnames[im*256]);
  } 
  SC_setnames(SC, xnames, err);
  forwardError(*err,__LINE__,NULL);
  
  free(xnames); 
  free(bnames); 
  free(ngcal);
  free(gcaltpl);

  return SC;
}

SmicaComp * clik_smica_comp_amp_diag_init(cldf *df,int nb, int mT,int mP,int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  int zz;
  double *tmpl;
  SmicaComp *SC;
  char *bnames;
  char **xnames;
  int im,ii;
  int hk;
  int m;

  m = mtot(mT,mP,has_cl);
  
  zz = m*nb;
  tmpl =  cldf_readfloatarray(df,"diag",&zz,err); 
  //tmpl =  hdf5_double_datarray(comp_id,cur_lkl,"diag",&zz,err); 
  forwardError(*err,__LINE__,NULL);  

  SC = amp_diag_init(nb, m, tmpl, err);

  hk = cldf_haskey(df,"names",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    int dz;
    dz =-1;
    //bnames = hdf5_char_attarray(comp_id,cur_lkl,"names",&dz, err);
    bnames = cldf_readstr(df,"names",&dz, err);
    forwardError(*err,__LINE__,NULL); 
  } else {
    bnames = malloc_err(sizeof(char)*256*m,err);
    forwardError(*err,__LINE__,NULL); 
    for(im=0;im<m;im++) {
      sprintf(&(bnames[ii*256]),"adiag_%d",im);
    }
  }

  xnames = malloc_err(sizeof(char*)*m,err);
  for(im=0;im<m;im++) {
    xnames[im] =&(bnames[im*256]);
  } 
  SC_setnames(SC, xnames, err);
  forwardError(*err,__LINE__,NULL);
  
  free(xnames); 
  free(bnames); 
  free(tmpl);
  return SC;
}

