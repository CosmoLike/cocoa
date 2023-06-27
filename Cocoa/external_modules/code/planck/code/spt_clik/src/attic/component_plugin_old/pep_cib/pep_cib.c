#include "clik_parametric.h"
#include "clik_parametric_addon.h"

////////////////////////////////////////////////////////////////////////////////////////
// PEP CIB model
////////////////////////////////////////////////////////////////////////////////////////

void ir_clustered_pep_compute(parametric* egl, double *Rq,error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,mell_in;
  double ir_clustered_pep_norm;
  double *A;
  int nfreqs_hfi = 6;
  int lmax_in = 3000;
  int *ind_freq;
  int ind1, ind2;
  template_payload *payload;

  payload = egl->payload;
  A = payload->template; // Stores input C_l(nu,nu') for all HFI frequencies and ell in [0,3000]
  nfreq = egl->nfreq;
  ind_freq = payload->ind_freq;

  ir_clustered_pep_norm = parametric_get_value(egl,"ir_clustered_pep_norm",err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell_in=ell*nfreqs_hfi*nfreqs_hfi;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        ind1 = ind_freq[m1];
        ind2 = ind_freq[m2];
        Rq[IDX_R(egl,ell,m1,m2)] = ir_clustered_pep_norm * A[mell_in + ind1*nfreqs_hfi + ind2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  return;
}

void ir_clustered_pep2_compute(parametric* egl, double *Rq,error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,mell_in;
  double ir_clustered_pep2_norm;
  double *A;
  int nfreqs_hfi = 6;
  int lmax_in = 3000;
  int *ind_freq;
  int ind1, ind2;
  template_payload *payload;

  payload = egl->payload;
  A = payload->template; // Stores input C_l(nu,nu') for all HFI frequencies and ell in [0,3000]
  nfreq = egl->nfreq;
  ind_freq = payload->ind_freq;

  ir_clustered_pep2_norm = parametric_get_value(egl,"ir_clustered_pep2_norm",err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell_in=ell*nfreqs_hfi*nfreqs_hfi;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        ind1 = ind_freq[m1];
        ind2 = ind_freq[m2];
        Rq[IDX_R(egl,ell,m1,m2)] = ir_clustered_pep2_norm * A[mell_in + ind1*nfreqs_hfi + ind2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  return;
}


parametric *ir_clustered_pep_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  pfchar type;
  char *pt;
  int lmax_in=3000;
  int nfreqs_hfi = 6;
  double hfi_freqlist[6] = {100.,143.,217.,353.,545.,857.};
  
  testErrorRetVA(lmax>lmax_in,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,NULL,lmax,lmax_in);

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);

  parametric_template_payload_init(egl,template, (nfreqs_hfi*nfreqs_hfi*(lmax_in+1)),hfi_freqlist,nfreqs_hfi,err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &ir_clustered_pep_compute;
  egl->eg_free = &parametric_template_payload_free;

  parametric_set_default(egl,"ir_clustered_pep_norm",1.0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"ir_clustered_pep_norm",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}


parametric *ir_clustered_pep2_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  pfchar type;
  char *pt;
  int lmax_in=3000;
  int nfreqs_hfi = 6;
  double hfi_freqlist[6] = {100.,143.,217.,353.,545.,857.};
  
  testErrorRetVA(lmax>lmax_in,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,NULL,lmax,lmax_in);

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);

  parametric_template_payload_init(egl,template, (nfreqs_hfi*nfreqs_hfi*(lmax_in+1)),hfi_freqlist,nfreqs_hfi,err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &ir_clustered_pep2_compute;
  egl->eg_free = &parametric_template_payload_free;

  parametric_set_default(egl,"ir_clustered_pep2_norm",1.0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"ir_clustered_pep2_norm",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}

void ir_poisson_pep_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double ir_poisson_pep_norm;
  double *A;
  int *ind_freq;
  int ind1, ind2;
  template_payload *payload;
  int nfreqs_hfi = 6;

  payload = egl->payload;
  A = payload->template; // Stores input C_l(nu,nu') for all HFI frequencies
  nfreq = egl->nfreq;
  ind_freq = payload->ind_freq;

  ir_poisson_pep_norm = parametric_get_value(egl,"ir_poisson_pep_norm",err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        ind1 = ind_freq[m1];
        ind2 = ind_freq[m2];
        Rq[IDX_R(egl,ell,m1,m2)] = ir_poisson_pep_norm * A[ind1*nfreqs_hfi + ind2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  return;
}

parametric *ir_poisson_pep_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  double hfi_freqlist[6] = {100.,143.,217.,353.,545.,857.};
  int nfreqs_hfi = 6;
  
  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_template_payload_init(egl,template, (nfreqs_hfi*nfreqs_hfi),hfi_freqlist,nfreqs_hfi,err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &ir_poisson_pep_compute;
  egl->eg_free = &parametric_template_payload_free;

  parametric_set_default(egl,"ir_poisson_pep_norm",1.0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"ir_poisson_pep_norm",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void ir_poisson_guilaine_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double ir_poisson_guilaine_norm;
  double *A;
  int *ind_freq;
  int ind1, ind2;
  template_payload *payload;
  int nfreqs_hfi = 6;

  payload = egl->payload;
  A = payload->template; // Stores input C_l(nu,nu') for all HFI frequencies
  nfreq = egl->nfreq;
  ind_freq = payload->ind_freq;

  ir_poisson_guilaine_norm = parametric_get_value(egl,"ir_poisson_guilaine_norm",err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        ind1 = ind_freq[m1];
        ind2 = ind_freq[m2];
        Rq[IDX_R(egl,ell,m1,m2)] = ir_poisson_guilaine_norm * A[ind1*nfreqs_hfi + ind2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  return;
}

parametric *ir_poisson_guilaine_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  double hfi_freqlist[6] = {100.,143.,217.,353.,545.,857.};
  int nfreqs_hfi = 6;
  
  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_template_payload_init(egl,template, (nfreqs_hfi*nfreqs_hfi),hfi_freqlist,nfreqs_hfi,err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &ir_poisson_guilaine_compute;
  egl->eg_free = &parametric_template_payload_free;

  parametric_set_default(egl,"ir_poisson_guilaine_norm",1.0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"ir_poisson_guilaine_norm",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(ir_clustered_pep,ir_clustered_pep_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(ir_clustered_pep2,ir_clustered_pep2_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(ir_poisson_pep,ir_poisson_pep_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(ir_poisson_guilaine,ir_poisson_guilaine_init);
