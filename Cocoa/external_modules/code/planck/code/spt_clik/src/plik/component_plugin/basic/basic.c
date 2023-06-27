#include "basic.h"


void poisson_tensor_bydet_compute(parametric* exg, double *Rq, error **err);

parametric *poisson_tensor_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;
  
  egl = parametric_bydet_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  egl->eg_compute = &poisson_tensor_bydet_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*(egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);

  sprintf(egl->tensor_norm_template,"%s","poisson_tensor_bydet_A_");
  egl->tensor_norm_template_len = strlen(egl->tensor_norm_template);
  
  for(ic=0;ic<egl->nfreq;ic++) {
    double nrm,frq;
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic);    
    frq = egl->freqlist[ic];
    nrm = sqrt(80/(3000.*3001./2./M_PI)*(1.98984227 * frq*frq*frq -27.70516974 *frq*frq + 125.04493152 *frq -181.36576083));
    parametric_set_default(egl,Ac,nrm,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  egl->eg_deriv_any = &parametric_tensor_norm_derivative;

  return egl;
}

void poisson_tensor_bydet_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double *A;
  int ic;
  pfchar Ac;

  A = egl->payload;

  parametric_tensor_fill(egl,A,err);
  forwardError(*err,__LINE__,);  
  
  nfreq = egl->nfreq;

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = A[m1] * A[m2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  
  return;
}

void powerlaw_tensor_bydet_compute(parametric* exg, double *Rq, error **err);

parametric *powerlaw_tensor_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;

  egl = parametric_bydet_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_tensor_bydet_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)*(egl->nfreq*egl->nfreq + egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"powerlaw_tensor_bydet_index",-1.25,err); // Karim's favorite
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "powerlaw_tensor_bydet_index", &parametric_index_derivative,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"powerlaw_tensor_bydet_correlation_step",1,err);
  forwardError(*err,__LINE__,NULL);
  //parametric_add_derivative_function(egl, "powerlaw_tensor_bydet_correlation_step", &ir_clustered_step_derivative,err);
  //forwardError(*err,__LINE__,NULL);

  sprintf(egl->tensor_norm_template,"%s","powerlaw_tensor_bydet_A_");
  egl->tensor_norm_template_len = strlen(egl->tensor_norm_template);
  
  for(ic=0;ic<egl->nfreq;ic++) {
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic);    
    parametric_set_default(egl,Ac,1,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  parametric_set_default(egl,"powerlaw_tensor_bydet_l_pivot",1000,err);
  forwardError(*err,__LINE__,NULL);
    
  egl->eg_deriv_any = &parametric_tensor_norm_derivative;
  return egl;
}

void powerlaw_tensor_bydet_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double l_pivot,index,v;
  double x, lnx, ln2x;
  double *A,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;
  
  A = egl->payload;
  parametric_tensor_fill(egl,A,err);
  forwardError(*err,__LINE__,);  

  l_pivot = parametric_get_value(egl,"powerlaw_tensor_bydet_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"powerlaw_tensor_bydet_index",err);
  forwardError(*err,__LINE__,);

  nfreq = egl->nfreq;
  step = parametric_get_value(egl,"powerlaw_tensor_bydet_correlation_step",err);  
  forwardError(*err,__LINE__,);

  dcm = A + egl->nfreq;

  for(m1=0;m1<egl->nfreq;m1++) {
    dcm[m1*egl->nfreq+m1] = 1.;
    for(m2=0;m2<m1;m2++) {
      dcm[m1*egl->nfreq+m2] = pow(step,powerlaw_tensor_step_index(egl,m1,m2));
      dcm[m2*egl->nfreq+m1] = dcm[m1*egl->nfreq+m2];
      //_DEBUGHERE_("%d %d %g %g",m1,m2,dcm[m1*egl->nfreq+m2],fabs(log(egl->freqlist[m1])-log(egl->freqlist[m2]))*2.33)
    }
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index);
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = A[m1]*A[m2] * dcm[m1*egl->nfreq+m2] * v;
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
  return;
}

parametric *poisson_tensor_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;
  
  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &poisson_tensor_bydet_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*(egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);

  sprintf(egl->tensor_norm_template,"%s","poisson_tensor_A_");
  egl->tensor_norm_template_len = strlen(egl->tensor_norm_template);
  
  for(ic=0;ic<egl->nfreq;ic++) {
    double nrm,frq;
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic);    
    frq = egl->freqlist[ic];
    nrm = sqrt(80/(3000.*3001./2./M_PI)*(1.98984227 * frq*frq*frq -27.70516974 *frq*frq + 125.04493152 *frq -181.36576083));
    parametric_set_default(egl,Ac,nrm,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  egl->eg_deriv_any = &parametric_tensor_norm_derivative;

  return egl;
}



void powerlaw_tensor_compute(parametric* exg, double *Rq, error **err);

parametric *powerlaw_tensor_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_tensor_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)*(egl->nfreq*egl->nfreq + egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"powerlaw_tensor_index",-1.25,err); // Karim's favorite
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "powerlaw_tensor_index", &parametric_index_derivative,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"powerlaw_tensor_correlation_step",1,err);
  forwardError(*err,__LINE__,NULL);
  //parametric_add_derivative_function(egl, "powerlaw_tensor_correlation_step", &ir_clustered_step_derivative,err);
  //forwardError(*err,__LINE__,NULL);

  sprintf(egl->tensor_norm_template,"%s","powerlaw_tensor_A_");
  egl->tensor_norm_template_len = strlen(egl->tensor_norm_template);
  
  for(ic=0;ic<egl->nfreq;ic++) {
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic);    
    parametric_set_default(egl,Ac,1,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  parametric_set_default(egl,"powerlaw_tensor_l_pivot",1000,err);
  forwardError(*err,__LINE__,NULL);
    
  egl->eg_deriv_any = &parametric_tensor_norm_derivative;
  return egl;
}



void powerlaw_tensor_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double l_pivot,index,v;
  double x, lnx, ln2x;
  double *A,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;
  
  A = egl->payload;
  parametric_tensor_fill(egl,A,err);
  forwardError(*err,__LINE__,);  

  l_pivot = parametric_get_value(egl,"powerlaw_tensor_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"powerlaw_tensor_index",err);
  forwardError(*err,__LINE__,);

  nfreq = egl->nfreq;
  step = parametric_get_value(egl,"powerlaw_tensor_correlation_step",err);  
  forwardError(*err,__LINE__,);

  dcm = A + egl->nfreq;

  for(m1=0;m1<egl->nfreq;m1++) {
    dcm[m1*egl->nfreq+m1] = 1.;
    for(m2=0;m2<m1;m2++) {
      dcm[m1*egl->nfreq+m2] = pow(step,powerlaw_tensor_step_index(egl,m1,m2));
      dcm[m2*egl->nfreq+m1] = dcm[m1*egl->nfreq+m2];
      //_DEBUGHERE_("%d %d %g %g",m1,m2,dcm[m1*egl->nfreq+m2],fabs(log(egl->freqlist[m1])-log(egl->freqlist[m2]))*2.33)
    }
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index);
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = A[m1]*A[m2] * dcm[m1*egl->nfreq+m2] * v;
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
  return;
}

void powerlaw_triangle_compute(parametric* exg, double *Rq, error **err);


parametric *powerlaw_triangle_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_triangle_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)*((egl->nfreq*egl->nfreq)*2),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"powerlaw_triangle_index",-1.25,err); // Karim's favorite
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "powerlaw_triangle_index", &parametric_index_derivative,err);
  forwardError(*err,__LINE__,NULL);

  
  sprintf(egl->tensor_norm_template,"%s","powerlaw_triangle_T_");
  egl->tensor_norm_template_len = strlen(egl->tensor_norm_template);
  
  for(ic=0;ic<egl->nfreq;ic++) {
    int jc;
    for(jc=ic;jc<egl->nfreq;jc++) {
      sprintf(Ac,"%s%d_%d",egl->tensor_norm_template,ic,jc);    
      parametric_set_default(egl,Ac,1,err);
      forwardError(*err,__LINE__,NULL);
    }
  }
  
  parametric_set_default(egl,"powerlaw_triangle_l_pivot",1000,err);
  forwardError(*err,__LINE__,NULL);
    
  egl->eg_deriv_any = &powerlaw_triangle_norm_derivative;
  return egl;
}

void powerlaw_triangle_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double l_pivot,index,v;
  double x, lnx, ln2x;
  double *A,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;
  
  A = egl->payload;
  parametric_triangle_fill(egl,A,err);
  forwardError(*err,__LINE__,);  

  l_pivot = parametric_get_value(egl,"powerlaw_triangle_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"powerlaw_triangle_index",err);
  forwardError(*err,__LINE__,);

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index);
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = A[m1*egl->nfreq+m2] * v;
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
  return;
}



void powerlaw_tanh_compute(parametric* exg, double *Rq, error **err);


parametric *powerlaw_tanh_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_tanh_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)*((egl->nfreq*egl->nfreq)*2),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"powerlaw_tanh_index",-1.25,err); // Karim's favorite
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "powerlaw_tanh_index", &parametric_index_derivative,err);
  forwardError(*err,__LINE__,NULL);

  
  sprintf(egl->tensor_norm_template,"%s","powerlaw_tanh_A_");
  egl->tensor_norm_template_len = strlen(egl->tensor_norm_template);
  
  for(ic=0;ic<egl->nfreq;ic++) {
    int jc;
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic);    
    parametric_set_default(egl,Ac,1,err);
    forwardError(*err,__LINE__,NULL);
    for(jc=ic+1;jc<egl->nfreq;jc++) {
      sprintf(Ac,"%sT_%d_%d",egl->tensor_norm_template,ic,jc);
      parametric_set_default(egl,Ac,0,err);
      forwardError(*err,__LINE__,NULL);
    }
  }
  
  parametric_set_default(egl,"powerlaw_tanh_l_pivot",1000,err);
  forwardError(*err,__LINE__,NULL);
    
  egl->eg_deriv_any = &powerlaw_tanh_norm_derivative;
  return egl;
}

void powerlaw_tanh_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double l_pivot,index,v;
  double x, lnx, ln2x;
  double *A,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;
  
  A = egl->payload;
  parametric_tanh_fill(egl,A,err);
  forwardError(*err,__LINE__,);  

  l_pivot = parametric_get_value(egl,"powerlaw_tanh_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"powerlaw_tanh_index",err);
  forwardError(*err,__LINE__,);

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index);
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = A[m1*egl->nfreq+m2] * v;
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
  return;
}










// PS a la GPE+Dunkley

void pointsource_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double l_pivot,index,v,lA,nrm;
  double *A;
  pfchar name;
  int stop;

  l_pivot = parametric_get_value(egl,"ps_l_pivot",err);
  forwardError(*err,__LINE__,);

  nrm = l_pivot*(l_pivot+1)/2/M_PI;
  
  A = egl->payload;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      v = 1;
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*nfreq+m2] = v/nrm;
      A[m2*nfreq+m1] = v/nrm;
    }
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = 1;
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        lA = A[m1*nfreq+m2];
        Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
        Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
      }  
    }
  }

  return;
}

void pointsource_A_derivative(parametric* egl, int iv,double *Rq, double *dRq, error **err);

parametric *pointsource_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int m1,m2;
  pfchar name;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  egl->eg_compute = &pointsource_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ps_l_pivot",3000,err);
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      parametric_set_default(egl,name,1,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  
  parametric_add_derivative_function(egl,"any",&pointsource_A_derivative,err);  
  forwardError(*err,__LINE__,NULL);

  return egl;
}

// PS a la GPE+Dunkley

void pointsource_bydet_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double l_pivot,index,v,lA,nrm;
  double *A;
  pfchar name;
  int stop;

  l_pivot = parametric_get_value(egl,"ps_l_pivot",err);
  forwardError(*err,__LINE__,);

  nrm = l_pivot*(l_pivot+1)/2/M_PI;
  
  A = egl->payload;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)m1,(int)m2);
      v = 1;
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*nfreq+m2] = v/nrm;
      A[m2*nfreq+m1] = v/nrm;
    }
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = 1;
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        lA = A[m1*nfreq+m2];
        Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
        Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
      }  
    }
  }

  return;
}

void pointsource_bydet_A_derivative(parametric* egl, int iv,double *Rq, double *dRq, error **err) {
  int ell,m1,m2,mell,nfreq,mv;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;
  double nrm;

  l_pivot = parametric_get_value(egl,"ps_l_pivot",err);
  forwardError(*err,__LINE__,);
  //
  //A = egl->payload;
  nfreq = egl->nfreq;
  //for(m1=0;m1<nfreq;m1++) {
  //  for(m2=m1;m2<nfreq;m2++) {
  //    sprintf(name,"ps_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
  //    v = 1;
  //    v = parametric_get_value(egl,name,err);
  //    A[m1*nfreq+m2] = v;
  //    A[m2*nfreq+m1] = v;
  //  }
  //}

  nrm = l_pivot*(l_pivot+1)/2/M_PI;
  stop = 0;
  v = 1/nrm;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)m1,m2);
      if (strcmp(egl->varkey[iv],name)==0) {
        stop=1;
        memset(dRq,0,sizeof(double)*(egl->lmax+1-egl->lmin)*nfreq*nfreq);
        for(ell=egl->lmin;ell<=egl->lmax;ell++) {
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
  parametric_end_derivative_loop(egl,&(dRq[iv]),egl->varkey[iv],err);
  forwardError(*err,__LINE__,);

  return;
}

parametric *pointsource_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int m1,m2;
  pfchar name;

  egl = parametric_bydet_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &pointsource_bydet_compute;
  forwardError(*err,__LINE__,NULL);
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ps_l_pivot",3000,err);
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)m1,(int)m2);
      parametric_set_default(egl,name,1,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  
  parametric_add_derivative_function(egl,"any",&pointsource_bydet_A_derivative,err);  
  forwardError(*err,__LINE__,NULL);

  return egl;
}





CREATE_PARAMETRIC_FILE_INIT(poisson_tensor_bydet,poisson_tensor_bydet_init);
CREATE_PARAMETRIC_FILE_INIT(powerlaw_tensor_bydet,powerlaw_tensor_bydet_init);
CREATE_PARAMETRIC_FILE_INIT(poisson_tensor,poisson_tensor_init);
CREATE_PARAMETRIC_FILE_INIT(powerlaw_tensor,powerlaw_tensor_init);
CREATE_PARAMETRIC_FILE_INIT(powerlaw_triangle,powerlaw_triangle_init);
CREATE_PARAMETRIC_FILE_INIT(powerlaw_tanh,powerlaw_tanh_init);
CREATE_PARAMETRIC_FILE_INIT(powerlaw_free_emissivity,powerlaw_free_emissivity_init);
CREATE_PARAMETRIC_FILE_INIT(pointsource,pointsource_init);
CREATE_PARAMETRIC_FILE_INIT(pointsource_bydet,pointsource_bydet_init);

