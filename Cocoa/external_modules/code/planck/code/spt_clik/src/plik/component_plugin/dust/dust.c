// group the dust stuff here
#include "dust.h"



// Dust emissivity normalized to 1 at nu0
double dust_spectrum(double nu, double T_dust, double beta_dust, double nu0) {

  double h_over_kT, x, x0, ex, ex0, res, res0;

  h_over_kT = 1./(T_dust * 20.836); // Assumes frequencies in GHz
  x0 = nu0 * h_over_kT;
  ex0 = exp(x0);
  res0 = pow(nu0,3.0+beta_dust)/(ex0 - 1.); // Overall normalization will go away
  x = nu * h_over_kT;
  ex = exp(x);
  res = pow(nu,3.0+beta_dust)/(ex - 1.);

  // Return dust emissivity normalized to one at nu0, in dT (CMB)
  return (res/res0/dBdT(nu,nu0));
}



// Non-thermal emissivity normalized to 1 at nu0 (e.g. synchrotron, free-free, etc.)
// The spectral index is called alpha here, to distinguish from the gray body case
// Expressed in dT (CMB)
double non_thermal_spectrum(double nu, double alpha_non_thermal, double nu0) {

  //return (pow(nu/nu0,alpha_non_thermal)/dBdT(nu,nu0));
  return (exp(alpha_non_thermal*log(nu/nu0))/dBdT(nu,nu0));
}

// TT

void hgal_compute(parametric *egl, double *Rq, error **err) {
  double *A;
  int nfreq,m1,m2,ell;
  double a,b,v,l_pivot,lA;
  pfchar name;
  int mell;

  A = egl->payload;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"hgal_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"hgal_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*nfreq+m2] = v;
      A[m2*nfreq+m1] = v;
    }
  }

  b = parametric_get_value(egl,"hgal_beta",err);
  forwardError(*err,__LINE__,);
  l_pivot = parametric_get_value(egl,"hgal_l_pivot",err);
  forwardError(*err,__LINE__,);
  
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = exp(b*(ell-l_pivot))/l_pivot/(l_pivot+1)*2*M_PI;
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        lA = A[m1*nfreq+m2];
        Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
        Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
      }  
    }
  }
}

parametric *hgal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int m1,m2;
  pfchar name;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &hgal_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"hgal_l_pivot",250,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"hgal_beta",-0.00849,err);
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"hgal_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"hgal_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      parametric_set_default(egl,name,50,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  
  
  return egl;
}
CREATE_PARAMETRIC_FILE_INIT(hgal,hgal_init);

void kgal_compute(parametric *egl, double *Rq, error **err) {
  double *A;
  int nfreq,m1,m2,ell;
  double a,b,v,l_pivot,lA,ns,nl;
  pfchar name;
  int mell;

  A = egl->payload;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"kgal_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"kgal_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*nfreq+m2] = v;
      A[m2*nfreq+m1] = v;
    }
  }

  nl = parametric_get_value(egl,"kgal_n_large",err);
  forwardError(*err,__LINE__,);
  ns = parametric_get_value(egl,"kgal_n_small",err);
  forwardError(*err,__LINE__,);
  l_pivot = parametric_get_value(egl,"kgal_l_pivot",err);
  forwardError(*err,__LINE__,);
  
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    //_DEBUGHERE_("%g %g %g %g %g %d %d",ell/l_pivot,nl,ns,pow(ell/l_pivot,nl),pow(ell/l_pivot,ns),(ell<l_pivot),(ell>=l_pivot))
    v = pow(ell/l_pivot,nl) * (ell<l_pivot) + pow(ell/l_pivot,ns) * (ell>=l_pivot);
    v = v/l_pivot/(l_pivot+1.)*2*M_PI;
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        lA = A[m1*nfreq+m2];
        Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
        Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
      }  
    }
  }
}

parametric *kgal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int m1,m2;
  pfchar name;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &kgal_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"kgal_l_pivot",250,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"kgal_n_large",-2.2,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"kgal_n_small",-2.8,err);
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"kgal_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"kgal_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  
  
  return egl;
}
CREATE_PARAMETRIC_FILE_INIT(kgal,kgal_init);

void galactic_component_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double norm,l_pivot,index,alpha_non_thermal,beta_dust,T_dust;
  double v,lA;
  double *A, *a, *B, *b;
  pfchar type;
  char* pt;
  int isdust;
  double prm_nu0;
  
  // dust or non_thermal

  norm = parametric_get_value(egl,"gal_norm",err);
  forwardError(*err,__LINE__,);
  
  l_pivot = parametric_get_value(egl,"gal_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"gal_index",err);
  forwardError(*err,__LINE__,);

  sprintf(type,"dust");
  isdust = 1;
  pt = (char*)type;
  pt = pflist_get_value(egl->pf,"gal_type",pt,err);
  forwardError(*err,__LINE__,);

  prm_nu0 = parametric_get_value(egl,"gal_nu0",err);
  forwardError(*err,__LINE__,);

  if (strcmp(pt,"dust")==0) {
    
    beta_dust = parametric_get_value(egl,"gal_beta_dust",err);
    forwardError(*err,__LINE__,);

    T_dust = parametric_get_value(egl,"gal_T_dust",err);
    forwardError(*err,__LINE__,);

  } else if (strcmp(pt,"non_thermal")==0) {

    isdust = 0;
    alpha_non_thermal = parametric_get_value(egl,"gal_alpha_non_thermal",err);
    forwardError(*err,__LINE__,);

  } else {
    testErrorRetVA(1==1,-1234,"Unknown Galactic component type '%s'",*err,__LINE__,,type);
  }
    

  A = egl->payload;            // Rank-one, tensorial emissivity
  nfreq = egl->nfreq;
  a = &(A[nfreq*nfreq]);       // Emissivity vector

  B = &(A[nfreq*(nfreq+1)]);   // Used for derivatives
  b = &(A[nfreq*(2*nfreq+1)]); // Idem

  for (m1=0;m1<nfreq;m1++) {
    if (isdust) {
      a[m1] = dust_spectrum((double)egl->freqlist[m1],T_dust,beta_dust,prm_nu0);
    } else {
      a[m1] = non_thermal_spectrum((double)egl->freqlist[m1],alpha_non_thermal,prm_nu0);
    }
  }

  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      v = a[m1]*a[m2];
      A[m1*nfreq+m2] = v;
      A[m2*nfreq+m1] = v;
    }
  }

  // R_q

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/l_pivot,index) * norm;
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        lA = A[m1*nfreq+m2];
        Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }

  return;
}

parametric *galactic_component_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  pfchar type;
  char* pt;
  int isdust;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  egl->eg_compute = &galactic_component_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*2*egl->nfreq*(egl->nfreq+1),err);
  forwardError(*err,__LINE__,NULL);
  
  // uK^2 at l=500, nu=143 GHz;
  parametric_set_default(egl,"gal_norm",1,err);
  forwardError(*err,__LINE__,NULL);
  parametric_declare_mandatory(egl,"gal_norm",err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"gal_norm",&parametric_norm_derivative,err);  
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"gal_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_index",0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_declare_mandatory(egl,"gal_index",err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"gal_index",&parametric_index_derivative,err);  
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_nu0",PRM_NU0,err);
  forwardError(*err,__LINE__,NULL);

  sprintf(type,"dust");
  isdust = 1;
  pt = (char*)type;
  pt = pflist_get_value(egl->pf,"gal_type",pt,err);
  forwardError(*err,__LINE__,NULL);

  if (strcmp(pt,"dust")==0) {
    parametric_set_default(egl,"gal_beta_dust",1.8,err);
    forwardError(*err,__LINE__,NULL);
    parametric_add_derivative_function(egl,"gal_beta_dust",&gal_beta_dust_derivative,err);  
    forwardError(*err,__LINE__,NULL);

    parametric_set_default(egl,"gal_T_dust",18,err);
    forwardError(*err,__LINE__,NULL);
    parametric_add_derivative_function(egl,"gal_T_dust",&gal_T_dust_derivative,err);  
    forwardError(*err,__LINE__,NULL);

  } else if (strcmp(pt,"non_thermal")==0) {
    //Intensity, = -3.0 in RJ
    parametric_set_default(egl,"gal_alpha_non_thermal",-1,err);
    forwardError(*err,__LINE__,NULL);
    parametric_add_derivative_function(egl,"gal_alpha_non_thermal",&gal_alpha_non_thermal_derivative,err);  
    forwardError(*err,__LINE__,NULL);
  } else {
    testErrorRetVA(1==1,-1234,"Unknown Galactic component type '%s'",*err,__LINE__,NULL,type);
    // return ?
  }

  return egl;
}

// P definitions
typedef struct {
  int n;
  int *m1,*m2;
  double *nrm;
  char **nrm_names;
} pw_XX_payload;

#define TE_KIND 1
#define EE_KIND 2
pw_XX_payload*  init_pw_XX_payload(int kind,int nT,int nP,error **err);

void pw_XX_free(void **PP);

//TE
void gal_TE_compute(parametric* exg, double *Rq, error **err);

parametric *gal_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int mx,i,j;
  pw_XX_payload *payload;
  int tl, mm, lbf,l,delta;

  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &gal_TE_compute;
  egl->eg_free = &pw_XX_free;

  egl->payload  = init_pw_XX_payload(TE_KIND,egl->nfreq_T,egl->nfreq_P,err);
  forwardError(*err,__LINE__,NULL);


  parametric_set_default(egl,"gal_TE_norm",1,err);
  forwardError(*err,__LINE__,NULL);
    
  parametric_set_default(egl,"gal_TE_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_TE_index",0,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"gal_TE_nu0",PRM_NU0,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_TE_beta_dust",1.8,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_TE_T_dust",18,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void gal_TE_compute(parametric* egl, double *Rq, error **err) {
  pw_XX_payload *payload;
  int i,mx;
  int ell,m1,m2;
  double *A;
  double l_pivot, index,v,norm;
  double prm_nu0,beta_dust,T_dust;

  payload = egl->payload;
  mx = egl->nfreq_T + egl->nfreq_P;
  

  prm_nu0 = parametric_get_value(egl,"gal_TE_nu0",err);
  forwardError(*err,__LINE__,);

  beta_dust = parametric_get_value(egl,"gal_TE_beta_dust",err);
  forwardError(*err,__LINE__,);

  T_dust = parametric_get_value(egl,"gal_TE_T_dust",err);
  forwardError(*err,__LINE__,);

  for(i=0;i<mx;i++) {
    payload->nrm[i] = dust_spectrum((double)egl->freqlist[i],T_dust,beta_dust,prm_nu0);
    forwardError(*err,__LINE__,);
  }
  A = payload->nrm;

  l_pivot = parametric_get_value(egl,"gal_TE_l_pivot",err);
  forwardError(*err,__LINE__,);

  index = parametric_get_value(egl,"gal_TE_index",err);
  forwardError(*err,__LINE__,);

  norm = parametric_get_value(egl,"gal_TE_norm",err);
  forwardError(*err,__LINE__,);
  

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index) * norm;
    for(i=0;i<payload->n;i++) {
      m1 = payload->m1[i];
      m2 = payload->m2[i];
      Rq[IDX_R(egl,ell,m1,m2)] = A[m1] * A[m2] * v;
      Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
    }
  }
}

//EE
void gal_EE_compute(parametric* exg, double *Rq, error **err);

parametric *gal_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int mx,i,j;
  pw_XX_payload *payload;
  int tl, mm, lbf,l,delta;
  double prm_nu0,beta_dust,T_dust;

  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &gal_EE_compute;
  egl->eg_free = &pw_XX_free;

  egl->payload  = init_pw_XX_payload(EE_KIND,egl->nfreq_T,egl->nfreq_P,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_EE_norm",1,err);
  forwardError(*err,__LINE__,NULL);
    
  parametric_set_default(egl,"gal_EE_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_EE_index",0,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"gal_EE_nu0",PRM_NU0,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_EE_beta_dust",1.8,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_EE_T_dust",18,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void gal_EE_compute(parametric* egl, double *Rq, error **err) {
  pw_XX_payload *payload;
  int i,mx;
  int ell,m1,m2;
  double *A;
  double l_pivot, index,v,norm;
  double prm_nu0,beta_dust,T_dust;

  payload = egl->payload;
  mx = egl->nfreq_P;
  

  prm_nu0 = parametric_get_value(egl,"gal_EE_nu0",err);
  forwardError(*err,__LINE__,);

  beta_dust = parametric_get_value(egl,"gal_EE_beta_dust",err);
  forwardError(*err,__LINE__,);

  T_dust = parametric_get_value(egl,"gal_EE_T_dust",err);
  forwardError(*err,__LINE__,);

  for(i=0;i<mx;i++) {
    payload->nrm[i] = dust_spectrum((double)egl->freqlist[i+egl->nfreq_T],T_dust,beta_dust,prm_nu0);
    forwardError(*err,__LINE__,);
  }
  A = payload->nrm;

  l_pivot = parametric_get_value(egl,"gal_EE_l_pivot",err);
  forwardError(*err,__LINE__,);

  index = parametric_get_value(egl,"gal_EE_index",err);
  forwardError(*err,__LINE__,);

  norm = parametric_get_value(egl,"gal_EE_norm",err);
  forwardError(*err,__LINE__,);
  

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index) * norm;
    for(i=0;i<payload->n;i++) {
      m1 = payload->m1[i];
      m2 = payload->m2[i];
      Rq[IDX_R(egl,ell,m1,m2)] = A[m1-egl->nfreq_T] * A[m2-egl->nfreq_T] * v;
      Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
    }
  }
}

// definitions
CREATE_PARAMETRIC_FILE_INIT(galametric,galactic_component_init);
CREATE_PARAMETRIC_POL_FILE_INIT(gal_TE,gal_TE_init);
CREATE_PARAMETRIC_POL_FILE_INIT(gal_EE,gal_EE_init);

#define GPE_DUST_DEFS double gpe_dust_freqlist[3] = {100.,143.,217.}; \
  int nfreqs_gpe_dust = 3; \
  double A=28.2846, alpha=0.538389, B=657.4222, beta=2.96484, dellc=2029.09, gamma=1.68974, delta=42.1039; \
  double r_gpe_dust[3][3]; \
  double color = 3.16; \
  r_gpe_dust[0][0] = 0.; \
  r_gpe_dust[1][0] = 0.; \
  r_gpe_dust[0][1] = 0.; \
  r_gpe_dust[1][1] = 0.; \
  r_gpe_dust[2][0] = 0.; \
  r_gpe_dust[0][2] = 0.; \
  r_gpe_dust[2][1] = color; \
  r_gpe_dust[1][2] = color; \
  r_gpe_dust[2][2] = color*color;


// Dust component 'a la GPE'
void gpe_dust_compute(parametric *egl, double *Rq, error **err);

parametric *gpe_dust_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  double fac;
  int ell, m1, m2;
  double dell;
  template_payload *payload;
  GPE_DUST_DEFS;

  egl = parametric_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);

  // Declare payload, allocate it and fill it

  parametric_check_freq(egl, gpe_dust_freqlist, nfreqs_gpe_dust,err);
  forwardError(*err,__LINE__,NULL);


  egl->payload = malloc_err(sizeof(template_payload),err);
  forwardError(*err,__LINE__,NULL);
  payload = egl->payload;

  // Get mapping between input frequencies and {143,217}
  payload->ind_freq = malloc_err(sizeof(int)*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);
  for (m1=0;m1<egl->nfreq;m1++) {
    payload->ind_freq[m1]=-1; //Init 
    for (m2=0;m2<nfreqs_gpe_dust;m2++) {
      if (fabs(egl->freqlist[m1]-gpe_dust_freqlist[m2])<1e-6) {
        payload->ind_freq[m1]=m2;
      }
    }
  }
  

  // Allocate template and precompute it (value at 143 GHz)
  payload->template = malloc_err(sizeof(double)*(lmax-lmin+1),err);
  forwardError(*err,__LINE__,NULL);
  for (ell=lmin;ell<=lmax;ell++) {
    dell = (double) ell;
    payload->template[ell-lmin] = 1./(color*color-1.) * ( A*pow(100./dell,alpha) + B* pow(dell/1000.,beta)/pow(1.+pow(dell/dellc,gamma),delta) );
  //payload->template[ell-lmin] = 1./(color*color-1.) * ( A*exp(alpha*log(100./dell)) + B* exp(beta*log(dell/1000.))/exp(delta*log(1.+exp(gamma*log(dell/dellc)))) );
    payload->template[ell-lmin] *= 2.0*M_PI/dell/(dell+1.0); // Converts to C_ell
  }
  
  egl->eg_compute = &gpe_dust_compute;
  egl->eg_free = &parametric_template_payload_free;
  
  parametric_set_default(egl,"gpe_dust_norm",1.0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"gpe_dust_norm",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void gpe_dust_compute(parametric *egl, double *Rq, error **err) {

  int ell, m1, m2, m2ll, nfreq;
  int *ind_freq;
  int ind1, ind2;
  double *template;
  template_payload *payload;
  double dell;
  double gpe_dust_norm;
  GPE_DUST_DEFS;

  nfreq = egl->nfreq;
  gpe_dust_norm = parametric_get_value(egl,"gpe_dust_norm",err);
  forwardError(*err,__LINE__,);
  payload = egl->payload;
  template = payload->template;

  ind_freq = payload->ind_freq;

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<nfreq;m1++) {
      ind1 = ind_freq[m1];
      for (m2=m1;m2<nfreq;m2++) {
        ind2 = ind_freq[m2];
        if ((ind1 >=0) && (ind2 >=0)) { //100, 143 or 217
          Rq[IDX_R(egl,ell,m1,m2)] = gpe_dust_norm * r_gpe_dust[ind1][ind2] * template[ell-egl->lmin];
          Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
        } else {
          Rq[IDX_R(egl,ell,m1,m2)] = 0.0;
          Rq[IDX_R(egl,ell,m2,m1)] = 0.0;
        }
      }
    }
  }
  return;
}

CREATE_PARAMETRIC_FILE_INIT(gpe_dust,gpe_dust_init);

void gpe_gal_compute(parametric *egl, double *Rq, error **err);

parametric *gpegal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  double fac;
  int ell, m1, m2;
  double dell;
  template_payload *payload;
  pfchar name;
  
  egl = parametric_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);

  
  parametric_set_default(egl,"gpegal_A",9.17742,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gpegal_alpha",0.650544,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gpegal_B",310.596,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gpegal_beta",2.6311,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gpegal_ell_c",1000,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gpegal_gamma",2.30962,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gpegal_delta",27.2319,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gpegal_l_pivot",200,err);
  forwardError(*err,__LINE__,NULL);

  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"gpegal_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"gpegal_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  
  

  // Declare payload, allocate it and fill it  

  egl->eg_compute = &gpe_gal_compute;
  egl->eg_free = &parametric_simple_payload_free;
  
  egl->payload = malloc_err(sizeof(double)*sizeof(double)*egl->nfreq*egl->nfreq,err);

  return egl;
}

void gpe_gal_compute(parametric *egl, double *Rq, error **err) {

  int ell, m1, m2, m2ll, nfreq;
  int *ind_freq;
  int ind1, ind2;
  double *template;
  template_payload *payload;
  double dell;
  double gpe_dust_norm;
  double A,alpha,B,beta,ell_c,gamma,delta,l_pivot;
  double *AA;
  double nrm,v;
  pfchar name;
  
  nfreq = egl->nfreq;
  
  l_pivot = parametric_get_value(egl,"gpegal_l_pivot",err);
  forwardError(*err,__LINE__,);

  A = parametric_get_value(egl,"gpegal_A",err);
  forwardError(*err,__LINE__,);
  alpha = parametric_get_value(egl,"gpegal_alpha",err);
  forwardError(*err,__LINE__,);
  B = parametric_get_value(egl,"gpegal_B",err);
  forwardError(*err,__LINE__,);
  beta = parametric_get_value(egl,"gpegal_beta",err);
  forwardError(*err,__LINE__,);
  ell_c = parametric_get_value(egl,"gpegal_ell_c",err);
  forwardError(*err,__LINE__,);
  gamma = parametric_get_value(egl,"gpegal_gamma",err);
  forwardError(*err,__LINE__,);
  delta = parametric_get_value(egl,"gpegal_delta",err);
  forwardError(*err,__LINE__,);
  
  AA = (double*) egl->payload;
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"gpegal_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"gpegal_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      AA[m1*nfreq+m2] = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      AA[m2*nfreq+m1] = AA[m1*nfreq+m2];      
    }
  }

  nrm = A * pow(100./l_pivot,alpha)+ B * pow(l_pivot/1000.,beta) * pow(1.+pow(l_pivot/ell_c,gamma),-delta);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = A * pow(100./ell,alpha)+ B * pow(ell/1000.,beta) * pow(1.+pow(ell/ell_c,gamma),-delta);
    v /= (ell*(ell+1)/2./M_PI);
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
      Rq[IDX_R(egl,ell,m1,m2)] = AA[m1*nfreq+m2] * v/nrm;
      Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;
}

CREATE_PARAMETRIC_FILE_INIT(gpegal,gpegal_init);

void gal545_compute(parametric *egl, double *Rq, error **err);

parametric *gal545_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  double fac;
  int ell, m1, m2;
  double dell;
  template_payload *payload;
  pfchar name;
  
  egl = parametric_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);

  
  parametric_set_default(egl,"gal545_h",2.3e-11,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gal545_k",5.05,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gal545_t",56,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gal545_index",-2.63,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal545_l_pivot",200,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal545_check_defpo",0,err);
  forwardError(*err,__LINE__,NULL);

  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"gal545_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"gal545_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  
  

  // Declare payload, allocate it and fill it  

  egl->eg_compute = &gal545_compute;
  egl->eg_free = &parametric_simple_payload_free;
  
  egl->payload = malloc_err(sizeof(double)*sizeof(double)*egl->nfreq*egl->nfreq,err);

  return egl;
}

void gal545_compute(parametric *egl, double *Rq, error **err) {

  int ell, m1, m2, m2ll, nfreq;
  int *ind_freq;
  int ind1, ind2;
  double *template;
  template_payload *payload;
  double dell;
  double gpe_dust_norm;
  double h,k,t,n,l_pivot;
  double *AA;
  double nrm,v;
  pfchar name;
  int defpo;

  nfreq = egl->nfreq;
  
  l_pivot = parametric_get_value(egl,"gal545_l_pivot",err);
  forwardError(*err,__LINE__,); 

  h = parametric_get_value(egl,"gal545_h",err);
  forwardError(*err,__LINE__,);
  
  k = parametric_get_value(egl,"gal545_k",err);
  forwardError(*err,__LINE__,);

  t = parametric_get_value(egl,"gal545_t",err);
  forwardError(*err,__LINE__,);

  n = parametric_get_value(egl,"gal545_index",err);
  forwardError(*err,__LINE__,);

  AA = (double*) egl->payload;
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"gal545_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"gal545_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      AA[m1*nfreq+m2] = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      AA[m2*nfreq+m1] = AA[m1*nfreq+m2];      
    }
  }
  
  defpo = parametric_get_value(egl,"gal545_check_defpo",err);
  forwardError(*err,__LINE__,);

  if (defpo!=0) {
    for(m1=0;m1<egl->nfreq;m1++) {
      if (AA[m1*nfreq+m1]==0) {
        continue;
      }
      for(m2=m1+1;m2<egl->nfreq;m2++) {
      if (AA[m2*nfreq+m2]==0) {
        continue;
      }
      testErrorRetVA(AA[m1*nfreq+m2]>sqrt(AA[m1*nfreq+m1] * AA[m2*nfreq+m2]),-130,"invalid dust amplitude (%d %d)",*err,__LINE__,,m1,m2)
      }
    }
  }

  nrm = (h * pow(l_pivot,k) * exp(-l_pivot/t) + 1) * 200*201/2./M_PI;

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = (h * pow(ell,k) * exp(-ell/t) + 1) * pow(ell/l_pivot,n);
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
      Rq[IDX_R(egl,ell,m1,m2)] = AA[m1*nfreq+m2] * v/nrm;
      Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;
}

CREATE_PARAMETRIC_FILE_INIT(gal545,gal545_init);

void gal545_80pc_compute(parametric *egl, double *Rq, error **err);

parametric *gal545_80pc_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  double fac;
  int ell, m1, m2;
  double dell;
  template_payload *payload;
  pfchar name;
  
  egl = parametric_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);


  parametric_set_default(egl,"gal545_80pc_h",3.03e-10,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gal545_80pc_k",4.87,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gal545_80pc_t",47.98,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gal545_80pc_index",-2.68,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal545_80pc_l_pivot",200,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal545_check_defpo",0,err);
  forwardError(*err,__LINE__,NULL);

  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"gal545_80pc_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"gal545_80pc_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  
  

  // Declare payload, allocate it and fill it  

  egl->eg_compute = &gal545_80pc_compute;
  egl->eg_free = &parametric_simple_payload_free;
  
  egl->payload = malloc_err(sizeof(double)*sizeof(double)*egl->nfreq*egl->nfreq,err);

  return egl;
}

void gal545_80pc_compute(parametric *egl, double *Rq, error **err) {

  int ell, m1, m2, m2ll, nfreq;
  int *ind_freq;
  int ind1, ind2;
  double *template;
  template_payload *payload;
  double dell;
  double gpe_dust_norm;
  double h,k,t,n,l_pivot;
  double *AA;
  double nrm,v;
  pfchar name;
  int defpo;

  nfreq = egl->nfreq;
  
  l_pivot = parametric_get_value(egl,"gal545_80pc_l_pivot",err);
  forwardError(*err,__LINE__,); 

  h = parametric_get_value(egl,"gal545_80pc_h",err);
  forwardError(*err,__LINE__,);
  
  k = parametric_get_value(egl,"gal545_80pc_k",err);
  forwardError(*err,__LINE__,);

  t = parametric_get_value(egl,"gal545_80pc_t",err);
  forwardError(*err,__LINE__,);

  n = parametric_get_value(egl,"gal545_80pc_index",err);
  forwardError(*err,__LINE__,);

  
  AA = (double*) egl->payload;
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"gal545_80pc_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"gal545_80pc_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      AA[m1*nfreq+m2] = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      AA[m2*nfreq+m1] = AA[m1*nfreq+m2];      
    }
  }

  defpo = parametric_get_value(egl,"gal545_check_defpo",err);
  forwardError(*err,__LINE__,);

  if (defpo!=0) {
    for(m1=0;m1<egl->nfreq;m1++) {
      if (AA[m1*nfreq+m1]==0) {
        continue;
      }
      for(m2=m1+1;m2<egl->nfreq;m2++) {
      if (AA[m2*nfreq+m2]==0) {
        continue;
      }
      testErrorRetVA(AA[m1*nfreq+m2]>sqrt(AA[m1*nfreq+m1] * AA[m2*nfreq+m2]),-130,"invalid dust amplitude (%d %d)",*err,__LINE__,,m1,m2)
      }
    }
  }
  
  nrm = (h * pow(l_pivot,k) * exp(-l_pivot/t) + 1) * 200*201/2./M_PI;

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = (h * pow(ell,k) * exp(-ell/t) + 1) * pow(ell/l_pivot,n);
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
      Rq[IDX_R(egl,ell,m1,m2)] = AA[m1*nfreq+m2] * v/nrm;
      Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;
}

CREATE_PARAMETRIC_FILE_INIT(gal545_80pc,gal545_80pc_init);


void t1gal_compute(parametric* egl, double *Rq, error **err);

parametric *t1gal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int lmin_t1gal_template= 0; //CHECK
  int lmax_t1gal_template= 3499; // CHECK
  int m1,m2;
  int l;
  pfchar name;
  testErrorRet(template==NULL,-11111,"template must be defined",*err,__LINE__,NULL);

  egl = parametric_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);

  egl->eg_compute = &t1gal_compute;
  egl->eg_free = &parametric_simple_payload_free;
  
  egl->payload = malloc_err(sizeof(double)*(lmax_t1gal_template-lmin_t1gal_template+1)+sizeof(double)*egl->nfreq*egl->nfreq,err);

  memcpy(egl->payload+sizeof(double)*egl->nfreq*egl->nfreq,template,sizeof(double)*(lmax_t1gal_template-lmin_t1gal_template+1));
  
  // Declare payload, allocate it and fill it
  parametric_set_default(egl,"t1gal_l_pivot",250,err);
  forwardError(*err,__LINE__,NULL);

  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"t1gal_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"t1gal_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  
  
  return egl;
}


void t1gal_compute(parametric* egl, double *Rq, error **err) {
  int nfreq,m1,m2,ell;
  double a,b,v,l_pivot,lA,ns,nl;
  pfchar name;
  int mell;
  int lmin_t1gal_template = 0; // CHECK
  int lmax_t1gal_template = 3499; // CHECK
  double *cl,*A;
  double nrm;
  
  nfreq = egl->nfreq;
  A = (double*) egl->payload;
  cl = A + egl->nfreq*egl->nfreq;
  
  l_pivot = parametric_get_value(egl,"t1gal_l_pivot",err);
  forwardError(*err,__LINE__,);
  
  nrm = cl[(int)l_pivot-lmin_t1gal_template]*l_pivot*(l_pivot+1.)/2./M_PI;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"t1gal_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"t1gal_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*nfreq+m2] = v/nrm;
      A[m2*nfreq+m1] = v/nrm;
    }
  }
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell = (ell-lmin_t1gal_template);
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
      Rq[IDX_R(egl,ell,m1,m2)] = A[m1*nfreq+m2] * cl[mell];
      Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;

}

CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(t1gal,t1gal_init);
