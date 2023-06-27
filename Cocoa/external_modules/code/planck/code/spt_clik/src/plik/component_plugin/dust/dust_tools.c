#include "dust.h"

// Group all unused derivatives stuff here
// Derivative of dust_spectrum with respect to T_dust
double d_dust_spectrum_d_T_dust(double nu, double T_dust, double beta_dust, double nu0) {

  double h_over_kT, x, x0, ex, ex0, res;

  h_over_kT = 1./(T_dust * 20.836); // Assumes frequencies in GHz
  x0 = nu0 * h_over_kT;
  ex0 = exp(x0);
  x = nu * h_over_kT;
  ex = exp(x);
  res = 1./T_dust * pow(nu/nu0,3.0+beta_dust) * (x*ex*(ex0-1.)-x0*ex0*(ex-1.)) / ((ex-1.)*(ex-1));
  res /= dBdT(nu,nu0);

  return(res);
  
} 

// Derivative of dust_spectrum with respect to beta_dust
double d_dust_spectrum_d_beta_dust(double nu, double T_dust, double beta_dust, double nu0) {
  
  return (log(nu/nu0)*dust_spectrum(nu,T_dust,beta_dust,nu0));

}

// Derivative of non_thermal_spectrum with respect to alpha
double d_non_thermal_spectrum_d_alpha_non_thermal(double nu, double alpha_non_thermal, double nu0) {

  return (log(nu/nu0)*non_thermal_spectrum(nu,alpha_non_thermal,nu0));
}

void gal_beta_dust_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  int ell,m1,m2;
  double *A, *B, *a, *b;
  int nfreq;
  double T_dust, beta_dust;
  double lA,v;
  double norm,l_pivot,index;
  double prm_nu0;

  A = egl->payload;            // Rank-one, tensorial emissivity
  nfreq = egl->nfreq;
  a = &(A[nfreq*nfreq]);       // Emissivity vector

  B = &(A[nfreq*(nfreq+1)]);   // Used for derivatives
  b = &(A[nfreq*(2*nfreq+1)]); // Idem

  beta_dust = parametric_get_value(egl,"gal_beta_dust",err);
  forwardError(*err,__LINE__,);

  T_dust = parametric_get_value(egl,"gal_T_dust",err);
  forwardError(*err,__LINE__,);

  norm = parametric_get_value(egl,"gal_norm",err);
  forwardError(*err,__LINE__,);
  
  l_pivot = parametric_get_value(egl,"gal_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"gal_index",err);
  forwardError(*err,__LINE__,);

  prm_nu0 = parametric_get_value(egl,"gal_nu0",err);
  forwardError(*err,__LINE__,);

  // dR/dbeta_dust
  // Get vector emissivity derivative
  for (m1=0;m1<nfreq;m1++) {
    b[m1] = d_dust_spectrum_d_beta_dust((double)egl->freqlist[m1],T_dust,beta_dust,prm_nu0);
  }
  for (m1=0;m1<nfreq;m1++) {
    for (m2=m1;m2<nfreq;m2++) {
      B[m1*nfreq+m2] = b[m1]*a[m2] + b[m2]*a[m1];
      B[m2*nfreq+m1] = B[m1*nfreq+m2];
    }
  }
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/egl->l_pivot,index) * norm;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        lA = B[m1*nfreq+m2];
        dRq[IDX_R(egl,ell,m1,m2)] = lA*v;
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  return;        
}

void gal_T_dust_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  int ell,m1,m2;
  double *A, *B, *a, *b;
  int nfreq;
  double T_dust, beta_dust;
  double lA,v;
  double norm,l_pivot,index;
  double prm_nu0;

  A = egl->payload;            // Rank-one, tensorial emissivity
  nfreq = egl->nfreq;
  a = &(A[nfreq*nfreq]);       // Emissivity vector

  B = &(A[nfreq*(nfreq+1)]);   // Used for derivatives
  b = &(A[nfreq*(2*nfreq+1)]); // Idem

  beta_dust = parametric_get_value(egl,"gal_beta_dust",err);
  forwardError(*err,__LINE__,);

  T_dust = parametric_get_value(egl,"gal_T_dust",err);
  forwardError(*err,__LINE__,);

  norm = parametric_get_value(egl,"gal_norm",err);
  forwardError(*err,__LINE__,);
  
  l_pivot = parametric_get_value(egl,"gal_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"gal_index",err);
  forwardError(*err,__LINE__,);

  prm_nu0 = parametric_get_value(egl,"gal_nu0",err);
  forwardError(*err,__LINE__,);

  // dR/dT_dust
  // Get vector emissivity derivative
  for (m1=0;m1<nfreq;m1++) {
    b[m1] = d_dust_spectrum_d_T_dust((double)egl->freqlist[m1],T_dust,beta_dust,prm_nu0);
  }
  for (m1=0;m1<nfreq;m1++) {
    for (m2=m1;m2<nfreq;m2++) {
      B[m1*nfreq+m2] = b[m1]*a[m2] + b[m2]*a[m1];
      B[m2*nfreq+m1] = B[m1*nfreq+m2];
    }
  }
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/egl->l_pivot,index) * norm;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        lA = B[m1*nfreq+m2];
        dRq[IDX_R(egl,ell,m1,m2)] = lA*v;
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
}

void gal_alpha_non_thermal_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  int ell,m1,m2;
  double *A, *B, *a, *b;
  int nfreq;
  double alpha_non_thermal;
  double lA,v;
  double norm,l_pivot,index;
  double prm_nu0;

  A = egl->payload;            // Rank-one, tensorial emissivity
  nfreq = egl->nfreq;
  a = &(A[nfreq*nfreq]);       // Emissivity vector

  B = &(A[nfreq*(nfreq+1)]);   // Used for derivatives
  b = &(A[nfreq*(2*nfreq+1)]); // Idem

  alpha_non_thermal = parametric_get_value(egl,"gal_alpha_non_thermal",err);
  forwardError(*err,__LINE__,);
  
  norm = parametric_get_value(egl,"gal_norm",err);
  forwardError(*err,__LINE__,);
  
  l_pivot = parametric_get_value(egl,"gal_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"gal_index",err);
  forwardError(*err,__LINE__,);

  prm_nu0 = parametric_get_value(egl,"gal_nu0",err);
  forwardError(*err,__LINE__,);

  // dR/dalpha_non_thermal
  // Get vector emissivity derivative
  for (m1=0;m1<nfreq;m1++) {
    b[m1] = d_non_thermal_spectrum_d_alpha_non_thermal((double)egl->freqlist[m1],alpha_non_thermal,prm_nu0);
  }
  for (m1=0;m1<nfreq;m1++) {
    for (m2=m1;m2<nfreq;m2++) {
      B[m1*nfreq+m2] = b[m1]*a[m2] + b[m2]*a[m1];
      B[m2*nfreq+m1] = B[m1*nfreq+m2];
    }
  }
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/egl->l_pivot,index) * norm;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        lA = B[m1*nfreq+m2];
        dRq[IDX_R(egl,ell,m1,m2)] = lA*v;
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }  
}
