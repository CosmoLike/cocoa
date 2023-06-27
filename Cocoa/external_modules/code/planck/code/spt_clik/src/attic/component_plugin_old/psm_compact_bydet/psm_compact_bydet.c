#include "clik_parametric.h"
#include "clik_parametric_addon.h"

////////////////////////////////////////////////////////////////////////////////////////
// PSM compact emission rigid model, by detector/detset
////////////////////////////////////////////////////////////////////////////////////////
void psm_compact_bydet_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,ndet,iv,mv,lell,mell_in;
  double psm_compact_bydet_norm;
  double *A;
  int ndets_hfi = 19;
  int lmax_in = 3000;

  A = egl->payload; // Stores input C_l(det,det') for det,det' in all HFI detsets/detectors and ell in [0,3000]
  ndet = egl->ndet;

  psm_compact_bydet_norm = parametric_get_value(egl,"psm_compact_bydet_norm",err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell_in=ell*ndets_hfi*ndets_hfi;
    for (m1=0;m1<ndet;m1++) {
      for (m2=m1;m2<ndet;m2++) {
        // Assumes all 19 detectors are in the parameter file ...
        Rq[IDX_R(egl,ell,m1,m2)] = psm_compact_bydet_norm *1e12* A[mell_in + m1*ndets_hfi + m2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;
}

parametric *psm_compact_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  pfchar type;
  char *pt;
  int lmax_in=3000;
  //char hfi_detnames[19][7] = {"100_ds1","100_ds2","143_ds1","143_ds2","143_5","143_6","143_7","217_1","217_2","217_3","217_4","217_ds1","217_ds2","353_1","353_2","353_ds1","353_ds2","353_7","353_8"};
  int ndets_hfi = 19;
  int m1,m2;

  testErrorRetVA(lmax>lmax_in,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,NULL,lmax,lmax_in);

  egl = parametric_bydet_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);

  egl->payload = malloc_err(sizeof(double)*(ndets_hfi*ndets_hfi*(lmax_in+1)),err);
  forwardError(*err,__LINE__,NULL);
  memcpy(egl->payload,template,sizeof(double)*(ndets_hfi*ndets_hfi*(lmax_in+1)));
  
  egl->eg_compute = &psm_compact_bydet_compute;
  egl->eg_free = &parametric_simple_payload_free;

  parametric_set_default(egl,"psm_compact_bydet_norm",1.0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"psm_compact_bydet_norm",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}

CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(psm_compact_bydet,psm_compact_bydet_init);
