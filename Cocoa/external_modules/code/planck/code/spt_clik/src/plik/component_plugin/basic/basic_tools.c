#include "basic.h"

void powerlaw_triangle_norm_derivative(parametric * egl, int iv, double *Rq, double *dRq, error **err) {
  int ell,m1,m2,mell,nfreq,mv,lell;
  double l_pivot,index,v;
  double x, lnx, ln2x;
  double *A,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;
  
  A = egl->payload;
  parametric_triangle_fill_derivative(egl,iv,A,err);
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
        dRq[IDX_R(egl,ell,m1,m2)] = A[m1*egl->nfreq+m2] * v;
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
  return;
}

void powerlaw_tanh_norm_derivative(parametric * egl, int iv, double *Rq, double *dRq, error **err) {
  int ell,m1,m2,mell,nfreq,mv,lell;
  double l_pivot,index,v;
  double x, lnx, ln2x;
  double *A,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;
  
  A = egl->payload;
  parametric_tanh_fill_derivative(egl,iv,A,err);
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
        dRq[IDX_R(egl,ell,m1,m2)] = A[m1*egl->nfreq+m2] * v;
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
  return;
}

double powerlaw_tensor_step_index(parametric* egl, int m1,int m2) {
  return fabs(log(egl->freqlist[m1])-log(egl->freqlist[m2]))*2.33;
}

void pointsource_A_derivative(parametric* egl, int iv,double *Rq, double *dRq, error **err) {
  int ell,m1,m2,mell,nfreq,mv;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;
  double nrm;

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

  l_pivot = egl->l_pivot;

  nrm = l_pivot*(l_pivot+1)/2/M_PI;
  stop = 0;
  v = 1/nrm;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
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

