#include "clik_parametric.h"
#include "clik_parametric_addon.h"

////////////////////////////////////////////////////////////////////////////////////////
// PSM foreground
////////////////////////////////////////////////////////////////////////////////////////
void ffp6_foreground_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell,mell_in;
  double v;
  double *A;     // covariance template
  double C[100]; // recalibration factor
  int nfreqs_hfi = 3;
  int lmax_in = 2508;
  int *ind_freq;
  int ind1, ind2;
  template_payload *payload;
  pfchar name;

  payload = egl->payload;
  A = payload->template; // Stores input C_l(nu,nu') for all HFI frequencies and ell in [0,3000]
  nfreq = egl->nfreq;
  ind_freq = payload->ind_freq;

  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"ffp6_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      forwardError(*err,__LINE__,);
      v = 1;
      v = parametric_get_value(egl,name,err);
      C[m1*nfreq+m2] = v;
      //C[m2*nfreq+m1] = v;
    }
  }

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell_in=ell*nfreqs_hfi*nfreqs_hfi;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        ind1 = ind_freq[m1];
        ind2 = ind_freq[m2];
        Rq[IDX_R(egl,ell,m1,m2)] = C[m1*nfreq+m2] * A[mell_in + ind1*nfreqs_hfi + ind2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;
}


void ffp6_foreground_derivative(parametric* egl, int iv,double *Rq, double *dRq, error **err) {
  int ell,m1,m2,mell,nfreq,mv,lell,mell_in;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;
  int nfreqs_hfi = 3;
  int lmax_in = 2508;
  int ind1,ind2;
  int* ind_freq;
  template_payload *payload;

  payload = egl->payload;
  A = payload->template; // Stores input C_l(nu,nu') for all HFI frequencies and ell in [0,3000]
  nfreq = egl->nfreq;
  ind_freq = payload->ind_freq;

  stop = 0;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"ffp6_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      if (strcmp(egl->varkey[iv],name)==0) {
        stop=1;
        memset(dRq,0,sizeof(double)*(egl->lmax+1-egl->lmin)*nfreq*nfreq);
	
	ind1 = ind_freq[m1];
	ind2 = ind_freq[m2];

        for(ell=egl->lmin;ell<=egl->lmax;ell++) {
	  mell_in=ell*nfreqs_hfi*nfreqs_hfi;
          v = A[mell_in + ind1*nfreqs_hfi + ind2];
          dRq[IDX_R(egl,ell,m1,m2)] = v;
          dRq[IDX_R(egl,ell,m2,m1)] = v;
        }
        break;
      }
    }
    if (stop==1) {
      return;

    }
  }      
  // error return
  parametric_end_derivative_loop(egl,&(dRq[mv]),egl->varkey[iv],err);
  forwardError(*err,__LINE__,);

  return;
}


parametric *ffp6_thermaldust_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  pfchar type;
  char *pt;
  int lmax_in=2508;
  double hfi_freqlist[3] = {100,143,217};
  int nfreqs_hfi = 3;
  int m1,m2;
  pfchar name;

  testErrorRetVA(lmax>lmax_in,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,NULL,lmax,lmax_in);

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_template_payload_init(egl,template, (nfreqs_hfi*nfreqs_hfi*(lmax_in+1)),hfi_freqlist,nfreqs_hfi,err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &ffp6_foreground_compute;
  egl->eg_free = &parametric_template_payload_free;

  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {  
      sprintf(name,"ffp6_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      parametric_set_default(egl,name,1,err); 
      forwardError(*err,__LINE__,NULL);       
    }
  }
  parametric_add_derivative_function(egl,"any",&ffp6_foreground_derivative,err);
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}



parametric *ffp6_firb_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  pfchar type;
  char *pt;
  int lmax_in=2508;
  double hfi_freqlist[3] = {100,143,217};
  int nfreqs_hfi = 3;
  int m1,m2;
  pfchar name;

  testErrorRetVA(lmax>lmax_in,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,NULL,lmax,lmax_in);

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_template_payload_init(egl,template, (nfreqs_hfi*nfreqs_hfi*(lmax_in+1)),hfi_freqlist,nfreqs_hfi,err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &ffp6_foreground_compute;
  egl->eg_free = &parametric_template_payload_free;

  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {  
      sprintf(name,"ffp6_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      parametric_set_default(egl,name,1,err); 
      forwardError(*err,__LINE__,NULL);       
    }
  }
  parametric_add_derivative_function(egl,"any",&ffp6_foreground_derivative,err);
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}

CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(ffp6_thermaldust,ffp6_thermaldust_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(ffp6_firb,ffp6_firb_init);  
