
#include "clik_parametric.h"
#include "clik_parametric_addon.h"

double dust_spectrum(double nu, double T_dust, double beta_dust, double nu0);
#define PRM_NU0 143.


typedef struct {
  int n;
  int *m1,*m2;
  double *nrm;
  char **nrm_names;
  int off1,off2;
} pw_XX_payload;

#define EE_KIND 1
#define BB_KIND 2
#define TE_KIND 3
#define TB_KIND 4
#define EB_KIND 5

pw_XX_payload*  init_pw_XX_payload(int kind,int nT,int nP,int* has_TEB,error **err) {
  pw_XX_payload* payload;
  int mx,i,j;
  int lim1,lim2,off1,off2,mul1;

  payload  = malloc_err(sizeof(pw_XX_payload),err);  
  forwardError(*err,__LINE__,NULL);
  
  if (kind==EE_KIND) {
    lim1 = nP;
    lim2 = nP;
    off1 = nT;
    off2 = nT;
    mul1 = 1;
    payload->n = (nP*(nP+1))/2;
  } else if (kind==BB_KIND) {
    lim1 = nP;
    lim2 = nP;
    off1 = nT+nP*has_TEB[1];
    off2 = nT+nP*has_TEB[1];
    mul1 = 1;
    payload->n = (nP*(nP+1))/2;
  } else if (kind==TE_KIND) {
    lim1 = nT;
    lim2 = nP;
    off1 = 0;
    off2 = nT;
    mul1 = 0;
    payload->n = nT * nP;
   } else if (kind==TB_KIND) {
    lim1 = nT;
    lim2 = nP;
    off1 = 0;
    off2 = nT+nP*has_TEB[1];
    mul1 = 0;
    payload->n = nT * nP;
  } else if (kind==EB_KIND) {
    lim1 = nP;
    lim2 = nP;
    off1 = nT;
    off2 = nT+nP*has_TEB[1];
    mul1 = 0;
    payload->n = nP * nP;
  } else {
    testErrorRetVA(1==1,-130,"invalid kind (%d)",*err,__LINE__,NULL,kind)
  }
    
  payload->off1 = off1;
  payload->off2 = off2;

  payload->m1 = malloc_err(sizeof(int)*payload->n,err);
  forwardError(*err,__LINE__,NULL);
  payload->m2 = malloc_err(sizeof(int)*payload->n,err);
  forwardError(*err,__LINE__,NULL);
  
  mx = 0;

  for(i=0;i<lim1;i++) {
    for(j=0 + i*mul1;j<lim2;j++) {
      payload->m1[mx] = i + off1;
      payload->m2[mx] = j + off2;
      mx++;
    }
  }

  payload->nrm_names = NULL;

  //mx = nT + nP;
  //
  //if (mx < payload->n) {
  //  mx = payload->n;
  //} 

  payload->nrm = malloc_err(sizeof(double)*(nT+nP*2)*(nT+nP*2),err);
  forwardError(*err,__LINE__,NULL);
  
  return payload;
}

void pw_XX_free(void **PP) {
  pw_XX_payload *P;
  P = *PP;
  free(P->nrm);
  free(P->m1);
  free(P->m2);
  if (P->nrm_names!=NULL) {
    free(P->nrm_names[0]);
    free(P->nrm_names);  
  }
  free(P);
  *PP = NULL;
}

void powerlaw_free_emissivity_XX_compute(parametric* egl, double *Rq,  error **err);

parametric *powerlaw_free_emissivity_XX_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int mx,i,j,kind;
  pw_XX_payload *payload;
  int tl, mm, lbf,l,delta,m1,m2;
  pfchar name;
  
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_free_emissivity_XX_compute;
  egl->eg_free = &pw_XX_free;

  kind = parametric_get_value(egl,"pwfe_XX_kind",err);
  forwardError(*err,__LINE__,NULL);
  
  egl->payload  = init_pw_XX_payload(kind,egl->nfreq_T,egl->nfreq_P,has_TEB,err);
  forwardError(*err,__LINE__,NULL);
  
  payload = egl->payload;

  parametric_set_default(egl,"pwfe_XX_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_XX_l2_norm",1,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_XX_index",0,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<payload->n;i++) {
    m1 = payload->m1[i];
    m2 = payload->m2[i];
    if (m1-payload->off1 > m2-payload->off2) {
      continue;
    }
    if (m1-payload->off1==m2-payload->off2) {
      sprintf(name,"pwfe_XX_A_%d",(int)egl->freqlist[m1]);
    } else {
      sprintf(name,"pwfe_XX_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
    }
    parametric_set_default(egl,name,0,err);
    forwardError(*err,__LINE__,NULL);
  }

  return egl;
}

void powerlaw_free_emissivity_XX_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,i;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;
  pw_XX_payload *payload;
  double nrmit;
  int l2norm;

  l_pivot = parametric_get_value(egl,"pwfe_XX_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  l2norm = parametric_get_value(egl,"pwfe_XX_l2_norm",err);
  forwardError(*err,__LINE__,);

  nrmit = 1;
  if (l2norm==1) {
    nrmit = l_pivot*(l_pivot+1.)/2./M_PI;
  }
  
  index = parametric_get_value(egl,"pwfe_XX_index",err);
  forwardError(*err,__LINE__,);

  payload = egl->payload;

  A = payload->nrm;
  nfreq = egl->nfreq;

  for(i=0;i<payload->n;i++) {
    m1 = payload->m1[i];
    m2 = payload->m2[i];
    if (m1-payload->off1 > m2-payload->off2) {
      int mtmp;
      mtmp = m2-payload->off2+payload->off1;
      m2 = m1-payload->off1+payload->off2;
      m1 = mtmp;
    }
    if (m1-payload->off1==m2-payload->off2) {
      sprintf(name,"pwfe_XX_A_%d",(int)egl->freqlist[m1]);
    } else {
      sprintf(name,"pwfe_XX_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
    }
    v = 1;
    v = parametric_get_value(egl,name,err);
    A[i] = v/nrmit;
  }

  #pragma omp parallel for private(ell,v, mell,i,m1,m2,lA)
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/l_pivot,index);
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(i=0;i<payload->n;i++) {
      m1 = payload->m1[i];
      m2 = payload->m2[i];
      lA = A[i];
      Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
      Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
    }  
  }

  return;
}


CREATE_PARAMETRIC_POL_FILE_INIT(powerlaw_free_emissivity_XX,powerlaw_free_emissivity_XX_init);


void pw_TE_compute(parametric* exg, double *Rq, error **err);
parametric *pw_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int mx,i,j;
  pw_XX_payload *payload;
  int hk;
  char **chn, *bf,*cbf,*pch;
  char *channeln;
  int tl, mm, lbf,l,delta;

  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &pw_TE_compute;
  egl->eg_free = &pw_XX_free;

  egl->payload  = init_pw_XX_payload(TE_KIND,egl->nfreq_T,egl->nfreq_P,has_TEB,err);
  forwardError(*err,__LINE__,NULL);
    
  payload = egl->payload;

  mx = egl->nfreq_T + egl->nfreq_P;

  parametric_set_default(egl,"pw_TE_index",-1,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"pw_TE_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);
  
  bf = malloc_err(sizeof(char)*mx*256,err);
  forwardError(*err,__LINE__,NULL);
  cbf = malloc_err(sizeof(char)*mx*256,err);
  forwardError(*err,__LINE__,NULL);
  chn = malloc_err(sizeof(char*)*mx,err);
  forwardError(*err,__LINE__,NULL);

  for(i=0;i<egl->nfreq_T;i++) {
    sprintf(cbf,"%s,%d_T",cbf,i);
  }
  for(i=0;i<egl->nfreq_P;i++) {
    sprintf(cbf,"%s,%d_E",cbf,i);
  }

  hk = cdic_key_index(egl->pf,"pw_TE_channel_names",err);
  forwardError(*err,__LINE__,NULL);
  
  channeln = cbf;

  
  if (hk!=-1) {
    channeln = cdic_get(egl->pf,"pw_TE_channel_names",NULL,err);
    forwardError(*err,__LINE__,NULL);  
  }
  
  tl = strlen(channeln);
  mm = 0;
  lbf = 0;
  l=0;
  while(l<tl) {
    chn[mm] = bf + lbf;
    while(channeln[l]==' ' || channeln[l]==',') {
      l++;
    }
    while(l<tl && channeln[l]!=',') {
      bf[lbf] = channeln[l];
      l++;
      lbf++;
    }
    bf[lbf] = '\0';
    lbf++;
    //_DEBUGHERE_("%d |%s|",mm,chn[mm]);
    mm++;
  }
  testErrorRet(mm!=mx,-97664,"bad pw_TE_channel_names value",*err,__LINE__,NULL);

  pch = malloc_err(sizeof(char)*mx*256,err);
  forwardError(*err,__LINE__,NULL);  

  payload->nrm_names = malloc_err(sizeof(char*)*mx,err);
  forwardError(*err,__LINE__,NULL);  

  delta = 0;

  for(i=0;i<mx;i++) {
    payload->nrm_names[i] = pch + delta;
    sprintf(payload->nrm_names[i],"pw_TE_A_%s",chn[i]);
    parametric_set_default(egl,payload->nrm_names[i],1.,err);
    forwardError(*err,__LINE__,NULL); 
   
    delta = strlen(payload->nrm_names[i]) +1;
    pch = payload->nrm_names[i];
  }      

  free(bf);
  free(cbf);
  return egl;
}

void pw_TE_compute(parametric* egl, double *Rq, error **err) {
  pw_XX_payload *payload;
  int i,mx;
  int ell,m1,m2;
  double *A;
  double l_pivot, index,v;

  payload = egl->payload;
  mx = egl->nfreq_T + egl->nfreq_P;
  

  for(i=0;i<mx;i++) {
    payload->nrm[i] = parametric_get_value(egl,payload->nrm_names[i],err);
    forwardError(*err,__LINE__,);
  }
  A = payload->nrm;

  l_pivot = parametric_get_value(egl,"pw_TE_l_pivot",err);
  forwardError(*err,__LINE__,);

  index = parametric_get_value(egl,"pw_TE_index",err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index);
    for(i=0;i<payload->n;i++) {
      m1 = payload->m1[i];
      m2 = payload->m2[i];
      Rq[IDX_R(egl,ell,m1,m2)] = A[m1] * A[m2] * v;
      Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
    }
  }
}

void pw_EE_compute(parametric* exg, double *Rq, error **err);

parametric *pw_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int mx,i,j;
  pw_XX_payload *payload;
  int hk;
  char **chn, *bf,*cbf,*pch;
  char *channeln;
  int tl, mm, lbf,l,delta;

  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &pw_EE_compute;
  egl->eg_free = &pw_XX_free;

  egl->payload  = init_pw_XX_payload(TE_KIND,egl->nfreq_T,egl->nfreq_P,has_TEB,err);
  forwardError(*err,__LINE__,NULL);
  payload = egl->payload;

  
  mx = egl->nfreq_P;

  parametric_set_default(egl,"pw_EE_index",-1,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"pw_EE_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);
  
  bf = malloc_err(sizeof(char)*mx*256,err);
  forwardError(*err,__LINE__,NULL);
  cbf = malloc_err(sizeof(char)*mx*256,err);
  forwardError(*err,__LINE__,NULL);
  chn = malloc_err(sizeof(char*)*mx,err);
  forwardError(*err,__LINE__,NULL);

  for(i=0;i<egl->nfreq_P;i++) {
    sprintf(cbf,"%s,%d_E",cbf,i);
  }

  hk = cdic_key_index(egl->pf,"pw_EE_channel_names",err);
  forwardError(*err,__LINE__,NULL);
  
  channeln = cbf;

  
  if (hk!=-1) {
    channeln = cdic_get(egl->pf,"pw_EE_channel_names",NULL,err);
    forwardError(*err,__LINE__,NULL);  
  }
  
  tl = strlen(channeln);
  mm = 0;
  lbf = 0;
  l=0;
  while(l<tl) {
    chn[mm] = bf + lbf;
    while(channeln[l]==' ' || channeln[l]==',') {
      l++;
    }
    while(l<tl && channeln[l]!=',') {
      bf[lbf] = channeln[l];
      l++;
      lbf++;
    }
    bf[lbf] = '\0';
    lbf++;
    //_DEBUGHERE_("%d |%s|",mm,chn[mm]);
    mm++;
  }
  testErrorRet(mm!=mx,-97664,"bad pw_EE_channel_names value",*err,__LINE__,NULL);

  pch = malloc_err(sizeof(char)*mx*256,err);
  forwardError(*err,__LINE__,NULL);  

  payload->nrm_names = malloc_err(sizeof(char*)*mx,err);
  forwardError(*err,__LINE__,NULL);  

  delta = 0;

  for(i=0;i<mx;i++) {
    payload->nrm_names[i] = pch + delta;
    sprintf(payload->nrm_names[i],"pw_EE_A_%s",chn[i]);
    parametric_set_default(egl,payload->nrm_names[i],1.,err);
    forwardError(*err,__LINE__,NULL); 
   
    delta = strlen(payload->nrm_names[i]) +1;
    pch = payload->nrm_names[i];
  }      

  free(bf);
  free(cbf);
  return egl;
}

void pw_EE_compute(parametric* egl, double *Rq, error **err) {
  pw_XX_payload *payload;
  int i,mx;
  int ell,m1,m2;
  double *A;
  double l_pivot, index,v;

  payload = egl->payload;
  mx = egl->nfreq_P;
  

  for(i=0;i<mx;i++) {
    payload->nrm[i] = parametric_get_value(egl,payload->nrm_names[i],err);
    forwardError(*err,__LINE__,);
  }
  A = payload->nrm;

  l_pivot = parametric_get_value(egl,"pw_EE_l_pivot",err);
  forwardError(*err,__LINE__,);

  index = parametric_get_value(egl,"pw_EE_index",err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index);
    for(i=0;i<payload->n;i++) {
      m1 = payload->m1[i];
      m2 = payload->m2[i];
      Rq[IDX_R(egl,ell,m1,m2)] = A[m1 - egl->nfreq_T] * A[m2- egl->nfreq_T] * v;
      Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
    }
  }
}


CREATE_PARAMETRIC_POL_FILE_INIT(pw_TE,pw_TE_init);
CREATE_PARAMETRIC_POL_FILE_INIT(pw_EE,pw_EE_init);

void powerlaw_free_emissivity_TE_compute(parametric* egl, double *Rq,  error **err);

parametric *powerlaw_free_emissivity_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int mx,i,j;
  pw_XX_payload *payload;
  int tl, mm, lbf,l,delta,m1,m2;
  pfchar name;
  
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_free_emissivity_TE_compute;
  egl->eg_free = &pw_XX_free;

  egl->payload  = init_pw_XX_payload(TE_KIND,egl->nfreq_T,egl->nfreq_P,has_TEB,err);
  forwardError(*err,__LINE__,NULL);
  payload = egl->payload;

  parametric_set_default(egl,"pwfe_TE_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_TE_l2_norm",1,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_TE_index",0,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<payload->n;i++) {
    m1 = payload->m1[i];
    m2 = payload->m2[i];
    if (m1 > m2 - ndet_T) {
      continue;
    }
    if (m1==m2 - ndet_T) {
      sprintf(name,"pwfe_TE_A_%d",(int)egl->freqlist[m1]);
    } else {
      sprintf(name,"pwfe_TE_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
    }
    parametric_set_default(egl,name,0,err);
    forwardError(*err,__LINE__,NULL);
  }

  return egl;
}

void powerlaw_free_emissivity_TE_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,i;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;
  pw_XX_payload *payload;
  double nrmit;
  int l2norm;

  l_pivot = parametric_get_value(egl,"pwfe_TE_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  l2norm = parametric_get_value(egl,"pwfe_TE_l2_norm",err);
  forwardError(*err,__LINE__,);

  nrmit = 1;
  if (l2norm==1) {
    nrmit = l_pivot*(l_pivot+1.)/2./M_PI;
  }
  
  index = parametric_get_value(egl,"pwfe_TE_index",err);
  forwardError(*err,__LINE__,);

  payload = egl->payload;

  A = payload->nrm;
  nfreq = egl->nfreq;

  for(i=0;i<payload->n;i++) {
    m1 = payload->m1[i];
    m2 = payload->m2[i];
    if (m1 > m2 - egl->ndet_T) {
      int mtmp;
      mtmp = m2-egl->ndet_T;
      m2 = m1+egl->ndet_T;
      m1 = mtmp;
    }
    if (m1==m2-egl->ndet_T) {
      sprintf(name,"pwfe_TE_A_%d",(int)egl->freqlist[m1]);
    } else {
      sprintf(name,"pwfe_TE_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
    }
    v = 1;
    v = parametric_get_value(egl,name,err);
    A[i] = v/nrmit;
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/l_pivot,index);
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(i=0;i<payload->n;i++) {
      m1 = payload->m1[i];
      m2 = payload->m2[i];
      lA = A[i];
      Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
      Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
    }  
  }

  return;
}

void powerlaw_free_emissivity_EE_compute(parametric* egl, double *Rq,  error **err);

parametric *powerlaw_free_emissivity_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int mx,i,j;
  pw_XX_payload *payload;
  int tl, mm, lbf,l,delta,m1,m2;
  pfchar name;
  
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_free_emissivity_EE_compute;
  egl->eg_free = &pw_XX_free;

  egl->payload  = init_pw_XX_payload(EE_KIND,egl->nfreq_T,egl->nfreq_P,has_TEB,err);
  forwardError(*err,__LINE__,NULL);
  payload = egl->payload;

  parametric_set_default(egl,"pwfe_EE_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_EE_l2_norm",1,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_EE_index",0,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"pwfe_EE_check_defpo",0,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<payload->n;i++) {
    m1 = payload->m1[i];
    m2 = payload->m2[i];
    
    if (m1==m2) {
      sprintf(name,"pwfe_EE_A_%d",(int)egl->freqlist[m1]);
    } else {
      sprintf(name,"pwfe_EE_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
    }
    parametric_set_default(egl,name,0,err);
    forwardError(*err,__LINE__,NULL);
  }

  return egl;
}

void powerlaw_free_emissivity_EE_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;
  pw_XX_payload *payload;
  double nrmit;
  int l2norm,i;
  int defpo;

  l_pivot = parametric_get_value(egl,"pwfe_EE_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  l2norm = parametric_get_value(egl,"pwfe_EE_l2_norm",err);
  forwardError(*err,__LINE__,);

  nrmit = 1;
  if (l2norm==1) {
    nrmit = l_pivot*(l_pivot+1.)/2./M_PI;
  }
  
  index = parametric_get_value(egl,"pwfe_EE_index",err);
  forwardError(*err,__LINE__,);

  payload = egl->payload;

  A = payload->nrm;
  nfreq = egl->nfreq;

  for(i=0;i<payload->n;i++) {
    m1 = payload->m1[i];
    m2 = payload->m2[i];
    if (m1==m2) {
      sprintf(name,"pwfe_EE_A_%d",(int)egl->freqlist[m1]);
    } else {
      sprintf(name,"pwfe_EE_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
    }
    v = 1;
    v = parametric_get_value(egl,name,err);
    A[m1*nfreq+m2] = v/nrmit;
    A[m2*nfreq+m1] = A[m1*nfreq+m2];
  }

  defpo = parametric_get_value(egl,"pwfe_EE_check_defpo",err);
  forwardError(*err,__LINE__,);

  if (defpo!=0) {
    for(m1=0;m1<egl->nfreq;m1++) {
      if (A[m1*nfreq+m1]==0) {
        continue;
      }
      for(m2=m1+1;m2<egl->nfreq;m2++) {
      if (A[m2*nfreq+m2]==0) {
        continue;
      }
      testErrorRetVA(A[m1*nfreq+m2]>sqrt(A[m1*nfreq+m1] * A[m2*nfreq+m2]),-130,"invalid dust amplitude (%d %d)",*err,__LINE__,,m1,m2)
      }
    }
  }
  
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/l_pivot,index);
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(i=0;i<payload->n;i++) {
      m1 = payload->m1[i];
      m2 = payload->m2[i];
      lA = A[m1*nfreq+m2];
      Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
      Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
    }  
  }

  return;
}

CREATE_PARAMETRIC_POL_FILE_INIT(powerlaw_free_emissivity_EE,powerlaw_free_emissivity_EE_init);
CREATE_PARAMETRIC_POL_FILE_INIT(powerlaw_free_emissivity_TE,powerlaw_free_emissivity_TE_init);

