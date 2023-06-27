#include "clik_parametric.h"
#include "clik_parametric_addon.h"

#define PRM_NU0 143.

// SZ+CIB+CROSS

#define N_FREQ_CIB 3
#define lmin_sz_template 2
#define lmax_sz_template  10000                               
#define lmin_corr_template  2
#define lmax_corr_template  9999

typedef struct {
  int *ind_freq;
  double *template;
  double *sz_template;
  double *corr_template;
  double a_cib[N_FREQ_CIB];
  double r_cib[N_FREQ_CIB*N_FREQ_CIB];
  double *fnu;
} sz_cib_payload;

void sz_cib_A_sz_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_A_cib_100_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_A_cib_143_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_A_cib_217_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_r_cib_100_143_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_r_cib_100_217_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_r_cib_143_217_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_xi_sz_cib_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_cib_index_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void parametric_sz_cib_payload_free(void **pp);
void parametric_sz_cib_payload_init(parametric *egl, double *template, error **err);
void sz_cib_fill_cib_struct(parametric *egl,error **err);

void parametric_zero_rq(parametric *egl, double* rrq);

#define SZ_CIB_DEFS                          \
  double a_sz,  xi_sz_cib;                   \
  int lnorm = 3000;                          \
  double *sz_template, *corr_template, *fnu; \
  double *a_cib,*r_cib;                      \
  double cib_corr, cib_index;                \
  sz_cib_payload *payload;                   \
  int do_sz,do_cib,do_szxcib;                \
  int ell,m1,m2,mell,nfreq;                  \
  int *ind_freq;                             \
  int ind1, ind2;                            \
  double d3000, dell;
  

#define SZ_CIB_INITS                                                  \
  do_sz = pflist_get_int_value(egl->pf,"do_sz",NULL,err);             \
  forwardError(*err,__LINE__,);                                       \
  do_cib = pflist_get_int_value(egl->pf,"do_cib",NULL,err);           \
  forwardError(*err,__LINE__,);                                       \
  do_szxcib = pflist_get_int_value(egl->pf,"do_szxcib",NULL,err);     \
  forwardError(*err,__LINE__,);                                       \
  nfreq = egl->nfreq;                                                 \
  payload = egl->payload;                                             \
  ind_freq = payload->ind_freq;                                       \
  d3000 = (3000.0*3001.)/2.0/M_PI;                                    \
  if (do_cib == 1 || do_szxcib == 1) {                                \
    a_cib = payload->a_cib;                                           \
    r_cib = payload->r_cib;                                           \
    sz_cib_fill_cib_struct(egl,err);                                  \
    forwardError(*err,__LINE__,);                                     \
    cib_index = parametric_get_value(egl,"cib_index",err);            \
    forwardError(*err,__LINE__,);                                     \
  }                                                                   \
  if (do_szxcib == 1)  {                                              \
    xi_sz_cib = parametric_get_value(egl,"xi_sz_cib",err);            \
    forwardError(*err,__LINE__,);                                     \
    sz_template = payload->sz_template;                               \
    corr_template = payload->corr_template;                           \
  }                                                                   \
  if (do_sz == 1 || do_szxcib == 1)  {                                \
    a_sz = parametric_get_value(egl,"A_sz",err);                      \
    forwardError(*err,__LINE__,);                                     \
    sz_template = payload->sz_template;                               \
    fnu = payload->fnu;                                               \
  }                                                                    
void cib_A_derivative(parametric* egl, int iv,double *Rq, double *dRq, error **err);

void fill_offset_freq(int idreq,double *dreq, parametric *egl,int *mv,int def,error **err);