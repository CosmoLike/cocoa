#include "clik_parametric.h"
#include "clik_parametric_addon.h"

////////////////////////////////////////////////////////////////////////////////////////
// PSM diffuse emission rigid model
////////////////////////////////////////////////////////////////////////////////////////
void psm_diffuse_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell,mell_in;
  double psm_diffuse_norm;
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

  psm_diffuse_norm = parametric_get_value(egl,"psm_diffuse_norm",err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell_in=ell*nfreqs_hfi*nfreqs_hfi;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        ind1 = ind_freq[m1];
        ind2 = ind_freq[m2];
        Rq[IDX_R(egl,ell,m1,m2)] = psm_diffuse_norm *1e12* A[mell_in + ind1*nfreqs_hfi + ind2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;
}

parametric *psm_diffuse_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
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
  
  egl->eg_compute = &psm_diffuse_compute;
  egl->eg_free = &parametric_template_payload_free;

  parametric_set_default(egl,"psm_diffuse_norm",1.0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"psm_diffuse_norm",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}

CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(psm_diffuse,psm_diffuse_init);
