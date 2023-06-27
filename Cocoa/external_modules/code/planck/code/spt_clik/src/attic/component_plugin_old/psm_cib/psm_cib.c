#include "clik_parametric.h"
#include "clik_parametric_addon.h"

////////////////////////////////////////////////////////////////////////////////////////
// PSM CIB rigid model
////////////////////////////////////////////////////////////////////////////////////////
void psm_cib_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell,mell_in;
  double psm_cib_norm;
  double *A;
  int nfreqs_hfi = 4;
  int lmax_in = 3000;
  int *ind_freq;
  int ind1, ind2;
  template_payload *payload;

  payload = egl->payload;
  A = payload->template; // Stores input C_l(nu,nu') for all HFI frequencies and ell in [0,3000]
  nfreq = egl->nfreq;
  ind_freq = payload->ind_freq;

  psm_cib_norm = parametric_get_value(egl,"psm_cib_norm",err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell_in=ell*nfreqs_hfi*nfreqs_hfi;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        ind1 = ind_freq[m1];
        ind2 = ind_freq[m2];
        Rq[IDX_R(egl,ell,m1,m2)] = psm_cib_norm *1e12* A[mell_in + ind1*nfreqs_hfi + ind2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;
}

parametric *psm_cib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  pfchar type;
  char *pt;
  int lmax_in=3000;
  double hfi_freqlist[4] = {100,143,217,353};
  int nfreqs_hfi = 4;
  int m1,m2;

  testErrorRetVA(lmax>lmax_in,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,NULL,lmax,lmax_in);

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_template_payload_init(egl,template, (nfreqs_hfi*nfreqs_hfi*(lmax_in+1)),hfi_freqlist,nfreqs_hfi,err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &psm_cib_compute;
  egl->eg_free = &parametric_template_payload_free;

  parametric_set_default(egl,"psm_cib_norm",1.0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"psm_cib_norm",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}


////////////////////////////////////////////////////////////////////////////////////////
// PSM CIB free model
////////////////////////////////////////////////////////////////////////////////////////
void psm_cib_free_compute(parametric* egl, double *Rq, double *dRq,error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell,mell_in;
  double psm_cib_norm;
  double *A,*N;
  int lmax_in = 3000;
  pfchar name;

  N = egl->payload;
  A = N + egl->nfreq; // Stores input C_l
  nfreq = egl->nfreq;
  

  for(m1=0;m1<nfreq;m1++) {    
    sprintf(name,"psm_cib_free_A_%d",(int)egl->freqlist[m1]);
    N[m1] = parametric_get_value(egl,name,err);
    forwardError(*err,__LINE__,);
  }

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = N[m1]*N[m2] * A[ell];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  if (dRq!=NULL) {
    double *ldRq;

    for (iv=0;iv<egl->nvar;iv++) {
      int m,done;
      
      ldRq = PTR_DER(egl,iv,dRq);

      done = 0;
      for(m=0;m<nfreq;m++) {
        sprintf(name,"psm_cib_free_A_%d",(int)egl->freqlist[m]);
        if (strcmp(egl->varkey[iv],name)==0) {
          memset(ldRq,0,sizeof(double)*(egl->lmax-egl->lmin+1)*nfreq*nfreq);
          for (ell=egl->lmin;ell<=egl->lmax;ell++) {
            for (m2=0;m2<nfreq;m2++) {
              if (m2==m) {
                ldRq[IDX_R(egl,ell,m,m)] = 2*Rq[mell+m*nfreq+m]/N[m];
                continue;    
              }
              ldRq[IDX_R(egl,ell,m1,m2)] = N[m2]*A[ell];
              ldRq[IDX_R(egl,ell,m2,m1)] = ldRq[IDX_R(egl,ell,m1,m2)];
            }
          }
          done=1;
          break;
        }
      }
      if (done==1) {
        continue;
      }
      // error return
      parametric_end_derivative_loop(egl,&(dRq[mv]),egl->varkey[iv],err);
      forwardError(*err,__LINE__,);
    }
  }

  return;
}

parametric *psm_cib_free_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  pfchar type;
  char *pt;
  int lmax_in=3000;
  int m1;
  pfchar name;

  testErrorRetVA(lmax>lmax_in,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,NULL,lmax,lmax_in);
   
  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
   
  egl->eg_compute_and_deriv = &psm_cib_free_compute;
  egl->eg_free = &parametric_simple_payload_free;
   
  egl->payload = malloc_err(sizeof(double)*((lmax_in+1)+egl->nfreq),err); //6 frequencies (HFI), lmax_in=3000
  forwardError(*err,__LINE__,NULL);
  memcpy(egl->payload+egl->nfreq*sizeof(double),template,(lmax_in+1)*sizeof(double));
  
  for(m1=0;m1<egl->nfreq;m1++) {
    sprintf(name,"psm_cib_free_A_%d",(int)egl->freqlist[m1]);
    pt = name;

    parametric_set_default(egl,pt,1,err);
    forwardError(*err,__LINE__,NULL);
  }

  return egl;
}


CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(psm_cib,psm_cib_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(psm_cib_free,psm_cib_free_init);
