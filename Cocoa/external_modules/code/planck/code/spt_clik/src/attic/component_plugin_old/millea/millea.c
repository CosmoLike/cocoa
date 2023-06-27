// group the old millea et al stuff here.
// that will probably be discarded in the near future...

#include "clik_parametric.h"
#include "clik_parametric_addon.h"

#define PRM_NU0 143.

void basic_compute_AB(parametric *egl, error **err) {
  int m1,m2,nfreq;
  double x, lnx, ln2x;
  double *A,*B,*vec;

  A = egl->payload;

  nfreq = egl->nfreq;
  B = A + nfreq*nfreq;
  vec = B + nfreq*nfreq;
  for (m1=0;m1<nfreq;m1++) {
    vec[m1] = dBdT((double)egl->freqlist[m1],PRM_NU0);
    for(m2=m1;m2<nfreq;m2++) {
      x = (double)egl->freqlist[m1]*(double)egl->freqlist[m2]/(PRM_NU0*PRM_NU0);
      lnx = log(x);
      ln2x = lnx*lnx;
      A[m1*nfreq+m2] = lnx;
      A[m2*nfreq+m1] = lnx;
      B[m1*nfreq+m2] = ln2x;
      B[m2*nfreq+m1] = ln2x;
    }
  }
  return;
}


void basic_alpha_derivative(parametric *egl,int iv, double *Rq, double *dRq, error **err) {
  double *A;
  int ell,m1,m2;

  A = egl->payload;
  // dR/dalpha_rg = log(nu1*nu2/nu0^2) * R
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<egl->nfreq;m1++) {
      for (m2=m1;m2<egl->nfreq;m2++) {
        dRq[IDX_R(egl,ell,m1,m2)] = A[m1*egl->nfreq+m2] * Rq[IDX_R(egl,ell,m1,m2)];
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
}

void basic_sigma_derivative(parametric *egl,int iv, double *Rq, double *dRq, error **err) {
  double *A,*B;
  double sigma;
  int ell,m1,m2;
  
  A = egl->payload;
  B = A + egl->nfreq*egl->nfreq;
  
  sigma = parametric_get_value(egl,egl->varkey[iv],err);
  forwardError(*err,__LINE__,);
  // dR/dsigma_rg = log(nu1*nu2/nu0^2)^2 * R
  
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<egl->nfreq;m1++) {
      for (m2=m1;m2<egl->nfreq;m2++) {
        dRq[IDX_R(egl,ell,m1,m2)] = sigma * B[m1*egl->nfreq+m2] * Rq[IDX_R(egl,ell,m1,m2)];
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  return;
}



// Millea et all reimplementation

// Millea et al. radio galaxies Poisson contribution. Based on a single population with given mean spectral index, as well as dispersion

void radiogal_compute(parametric* exg, double *Rq, error **err);

parametric *radiogal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  egl->eg_compute = &radiogal_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*(2*egl->nfreq*egl->nfreq + egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"radiogal_norm",78.5,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "radiogal_norm", &parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"radiogal_alpha",-0.36,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "radiogal_alpha", &basic_alpha_derivative,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"radiogal_sigma",0.64,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "radiogal_sigma", &basic_sigma_derivative,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void radiogal_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double norm_rg, alpha_rg, sigma_rg;
  double d3000;
  double x, lnx, ln2x;
  double *A,*B,*vec;

  d3000 = 3000.*3001./2./M_PI;
  
  A = egl->payload;
  
  nfreq = egl->nfreq;
  B = A + nfreq*nfreq;
  vec = B + nfreq*nfreq;
  
  norm_rg = parametric_get_value(egl,"radiogal_norm",err);
  forwardError(*err,__LINE__,);

  alpha_rg = parametric_get_value(egl,"radiogal_alpha",err);
  forwardError(*err,__LINE__,);

  sigma_rg = parametric_get_value(egl,"radiogal_sigma",err);
  forwardError(*err,__LINE__,);

  basic_compute_AB(egl, err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = norm_rg/d3000 * exp(alpha_rg*A[m1*nfreq+m2] +
              sigma_rg*sigma_rg/2.0 * B[m1*nfreq+m2])/(vec[m1]*vec[m2]);
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////
// CIB model. Millea et al. for now
////////////////////////////////////////////////////////////////////////////////////////////

void ir_poisson_compute(parametric* exg, double *Rq, error **err);

parametric *ir_poisson_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  egl->eg_compute = &ir_poisson_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*(2*egl->nfreq*egl->nfreq + egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ir_poisson_norm",5.9,err); // Millea et al. ref value
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "ir_poisson_norm", &parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ir_poisson_alpha",3.8,err); // Millea et al. ref value
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "ir_poisson_alpha", &basic_alpha_derivative,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ir_poisson_sigma",0.4,err); // Millea et al. ref value
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "ir_poisson_sigma", &basic_sigma_derivative,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void ir_poisson_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double ir_poisson_norm, ir_poisson_alpha, ir_poisson_sigma;
  double d3000;
  double x, lnx, ln2x;
  double *A,*B,*vec;

  basic_compute_AB(egl,err);
  forwardError(*err,__LINE__,);

  A = egl->payload;
  nfreq = egl->nfreq;
  B = A + nfreq*nfreq;
  vec = B + nfreq*nfreq;
  d3000 = 3000.*3001./2./M_PI;

  ir_poisson_norm = parametric_get_value(egl,"ir_poisson_norm",err);
  forwardError(*err,__LINE__,);

  ir_poisson_alpha = parametric_get_value(egl,"ir_poisson_alpha",err);
  forwardError(*err,__LINE__,);

  ir_poisson_sigma = parametric_get_value(egl,"ir_poisson_sigma",err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = ir_poisson_norm/d3000 * exp(ir_poisson_alpha*A[m1*nfreq+m2] +
              ir_poisson_sigma*ir_poisson_sigma/2.0 * B[m1*nfreq+m2])/(vec[m1]*vec[m2]);
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////////////
// CIB clustered part. Millea et al. for now.
////////////////////////////////////////////////////////////////////////////////////////

void ir_clustered_compute(parametric* exg, double *Rq, error **err);
void ir_clustered_step_derivative(parametric* egl, int iv, double *Rq, double *dRq, error **err);

parametric *ir_clustered_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  pfchar type;
  char *pt;


  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &ir_clustered_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)*(2*egl->nfreq*egl->nfreq + egl->nfreq + (lmax-lmin+1)),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ir_clustered_norm",3.9,err); // Millea et al. ref value
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "ir_clustered_norm", &parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ir_clustered_alpha",3.8,err); // Millea et al. ref value
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "ir_clustered_alpha", &basic_alpha_derivative,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ir_clustered_index",-1.25,err); // Karim's favorite
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "ir_clustered_index", &parametric_index_derivative,err);
  forwardError(*err,__LINE__,NULL);

  sprintf(type,"step");
  pt = (char*)type;
  pt = pflist_get_value(egl->pf,"ir_clustered_correlation",pt,err);
  forwardError(*err,__LINE__,NULL);

  //_DEBUGHERE_("%s",pt);
  
  if (strcmp(pt,"step")==0) {
    parametric_set_default(egl,"ir_clustered_correlation_step",0.5,err);
    forwardError(*err,__LINE__,NULL);
    parametric_add_derivative_function(egl, "ir_clustered_correlation_step", &ir_clustered_step_derivative,err);
    forwardError(*err,__LINE__,NULL);
  } else if (strcmp(pt,"matrix")==0) {
    int m1,m2;
    double vl;
    pfchar name;
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        sprintf(name,"ir_clustered_correlation_M_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
        parametric_set_default(egl,name,1,err);
        forwardError(*err,__LINE__,NULL);
      }
    }  
  } else {
    testErrorRetVA(1==1,-1234,"Unknown ir_clustered_correlation type '%s'",*err,__LINE__,,type);
    // return ?
  }

  return egl;
}

double ir_clustered_step_index(parametric* egl, int m1,int m2) {
  return fabs(log(egl->freqlist[m1])-log(egl->freqlist[m2]))*2.33;
}

void ir_clustered_step_derivative(parametric* egl, int iv, double *Rq, double *dRq, error **err) {
  int ell,m1,m2; 
  double step ;
  char *key;

  key = egl->varkey[iv];
  step = parametric_get_value(egl,key,err);
  forwardError(*err,__LINE__,);
  
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<egl->nfreq;m1++) {
      dRq[IDX_R(egl,ell,m1,m1)] = 0; //diagonal unchanged
      for (m2=m1+1;m2<egl->nfreq;m2++) {
        dRq[IDX_R(egl,ell,m1,m2)] = Rq[IDX_R(egl,ell,m1,m2)]/step * ir_clustered_step_index(egl,m1,m2);
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  
  return;
}

void ir_clustered_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double ir_clustered_norm, ir_clustered_alpha,ir_clustered_index;
  double d3000,t3000;
  double x, lnx, ln2x;
  double *A,*vec,*template,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;

  A = egl->payload; // Will store log(nu/nu0)
  d3000 = 3000.*3001./2./M_PI;
  
  egl->l_pivot = 3000;

  ir_clustered_norm = parametric_get_value(egl,"ir_clustered_norm",err);
  forwardError(*err,__LINE__,);

  ir_clustered_alpha = parametric_get_value(egl,"ir_clustered_alpha",err);
  forwardError(*err,__LINE__,);

  ir_clustered_index = parametric_get_value(egl,"ir_clustered_index",err);
  forwardError(*err,__LINE__,);

  nfreq = egl->nfreq;
  vec = A + nfreq*nfreq; // Will store dB/dT(nu,nu0)
  template = vec + nfreq; // Will store Cl spectrum (unnormalized)

  for (m1=0;m1<nfreq;m1++) {
    vec[m1] = dBdT((double)egl->freqlist[m1],PRM_NU0);
    for(m2=m1;m2<nfreq;m2++) {
      x = (double)egl->freqlist[m1]*(double)egl->freqlist[m2]/(PRM_NU0*PRM_NU0);
      lnx = log(x);
      A[m1*nfreq+m2] = lnx;
      A[m2*nfreq+m1] = lnx;
    }
  }

  // Create Cl template
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell=ell-egl->lmin;
    template[mell] = pow((double)ell,ir_clustered_index);
  }
  t3000 = pow(3000.,ir_clustered_index);

  //create correlation matrix
  dcm = template + egl->lmax-egl->lmin+1;

  sprintf(type,"step");
  pt = (char*)type;
  pt = pflist_get_value(egl->pf,"ir_clustered_correlation",pt,err);
  forwardError(*err,__LINE__,);
  //_DEBUGHERE_("%s",pt);

  isstep=0;
  if (strcmp(pt,"step")==0) {  
    isstep=1;

    step = parametric_get_value(egl,"ir_clustered_correlation_step",err);
    forwardError(*err,__LINE__,);
    for(m1=0;m1<egl->nfreq;m1++) {
      dcm[m1*egl->nfreq+m1] = 1.;
      for(m2=0;m2<m1;m2++) {
        dcm[m1*egl->nfreq+m2] = pow(step,ir_clustered_step_index(egl,m1,m2));
        dcm[m2*egl->nfreq+m1] = dcm[m1*egl->nfreq+m2];
        //_DEBUGHERE_("%d %d %g %g",m1,m2,dcm[m1*egl->nfreq+m2],fabs(log(egl->freqlist[m1])-log(egl->freqlist[m2]))*2.33)
      }
    }
  } else {
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        pfchar name;
        sprintf(name,"ir_clustered_correlation_M_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
        dcm[m1*egl->nfreq+m2] = parametric_get_value(egl,name,err);
        forwardError(*err,__LINE__,);
        dcm[m2*egl->nfreq+m1] = dcm[m1*egl->nfreq+m2];
        //_DEBUGHERE_("%d %d %g",m1,m2,dcm[m1*egl->nfreq+m2]);
      }
    }
  }

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell=(ell-egl->lmin)*nfreq*nfreq;
    lell = ell - egl->lmin;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = ir_clustered_norm * 
          template[lell]/(d3000*t3000) * 
          exp(ir_clustered_alpha*A[m1*nfreq+m2])/(vec[m1]*vec[m2]) * dcm[m1*egl->nfreq+m2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;
}

CREATE_PARAMETRIC_FILE_INIT(radiogal,radiogal_init);
CREATE_PARAMETRIC_FILE_INIT(ir_clustered,ir_clustered_init);
CREATE_PARAMETRIC_FILE_INIT(ir_poisson,ir_poisson_init);

