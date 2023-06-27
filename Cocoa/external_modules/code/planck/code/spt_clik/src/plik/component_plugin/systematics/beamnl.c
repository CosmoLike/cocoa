#include "clik_parametric.h"
#include "clik_parametric_addon.h"

void fill_offset_freq(int idreq,double *dreq, parametric *egl,int *mv,int def, error **err);

void mul0_compute(parametric *egl, double *rq, error **err);

parametric *mul0_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &mul0_compute;
  
  parametric_set_default(egl,"mul0_value",100,err); 
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}

void mul0_compute(parametric *egl, double *Rq, error **err) {
  double l_pivot,delta_l,v;
  int m1,m2,ell;
  double dip_nrm;
  int lmin,lmax;

  v = parametric_get_value(egl,"mul0_value",err);
  forwardError(*err,__LINE__,);
  lmin =  egl->lmin;
  lmax =  egl->lmax;
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=0;m2<egl->nfreq;m2++) {
      for(ell=lmin;ell<=lmax;ell++) {
        Rq[IDX_R(egl,ell,m1,m2)] = v; 
      }
    }
  }
}
CREATE_PARAMETRIC_FILE_INIT(mul0,mul0_init);

void beamnl_compute(parametric *egl, double *rq, error **err);

parametric *beamnl_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;
  double *Ip;
  int *mv;
  int nfreq_template,lmax_template;
  pfchar name;
  int nmod,ntrunc;
  double *dreq;
  int i,m1,p;

  testErrorRet(template==NULL,-11111,"template must be defined",*err,__LINE__,NULL);

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &beamnl_compute;
  egl->eg_free = &parametric_simple_payload_free;

  nfreq_template = parametric_get_value(egl,"beamnl_nfreq_template",err);
  forwardError(*err,__LINE__,NULL);

  lmax_template = parametric_get_value(egl,"beamnl_lmax_template",err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"beamnl_nmode",5 ,err); 
  forwardError(*err,__LINE__,NULL);
  nmod = parametric_get_value(egl,"beamnl_nmode",err);
  forwardError(*err,__LINE__,NULL);
  nmod+=1;

  parametric_set_default(egl,"beamnl_ptrunc",14 ,err); 
  forwardError(*err,__LINE__,NULL);
  ntrunc = parametric_get_value(egl,"beamnl_ptrunc",err);
  forwardError(*err,__LINE__,NULL);
  ntrunc+=1;

  egl->payload = malloc_err(sizeof(double)* ((lmax_template+1)*(nfreq_template)) + sizeof(double)*egl->nfreq*ntrunc + sizeof(double)*egl->nfreq+sizeof(int)*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  memcpy(egl->payload,template,sizeof(double)* ((lmax_template+1)*(nfreq_template)));
  
  mv = egl->payload + sizeof(double)* ((lmax_template+1)*(nfreq_template)) + sizeof(double)*egl->nfreq*ntrunc + sizeof(double)*egl->nfreq;
  
  dreq = malloc_err(sizeof(double)*nfreq_template,err);
  forwardError(*err,__LINE__,NULL);

  for(i=0;i<nfreq_template;i++) {
    sprintf(name,"beamnl_freq_%d",i);
    dreq[i] = parametric_get_value(egl,name,err);
    forwardError(*err,__LINE__,NULL);
  }  

  fill_offset_freq(nfreq_template,dreq, egl,mv,-1,err);
  forwardError(*err,__LINE__,NULL);

  for(m1=0;m1<egl->nfreq;m1++) {
    for(p=2;p<nmod;p++) {
      sprintf(name,"beamnl_r_%d_%d",(int)egl->freqlist[m1],p);
      parametric_set_default(egl,name,0,err); 
      forwardError(*err,__LINE__,NULL);
    }
  }

  free(dreq);

  return egl;
}

void beamnl_compute(parametric *egl, double *Rq, error **err) {
  double l_pivot,delta_l,v;
  int m1,m2,ell,p,p1,p2,t,ntrunc;
  double dip_nrm;
  int lmin,lmax;
  double *template,*It,*d;
  int *mv;
  int nfreq_template,lmax_template;
  pfchar name;
  int nmod;
  int strict;
  double t1,t2,lt1,lt2,rt,r,pt,d1,tp1;


  nfreq_template = parametric_get_value(egl,"beamnl_nfreq_template",err);
  forwardError(*err,__LINE__,);

  lmax_template = parametric_get_value(egl,"beamnl_lmax_template",err);
  forwardError(*err,__LINE__,);
  
  nmod = parametric_get_value(egl,"beamnl_nmode",err);
  forwardError(*err,__LINE__,);
  nmod+=1;
  ntrunc = parametric_get_value(egl,"beamnl_ptrunc",err);
  forwardError(*err,__LINE__,);
  ntrunc+=1;
  
  template = egl->payload;
  It = egl->payload + sizeof(double)* ((lmax_template+1)*(nfreq_template));
  d = egl->payload + sizeof(double)* ((lmax_template+1)*(nfreq_template)) + sizeof(double)*egl->nfreq*ntrunc;
  mv = egl->payload + sizeof(double)* ((lmax_template+1)*(nfreq_template)) + sizeof(double)*egl->nfreq*ntrunc + sizeof(double)*egl->nfreq;

  for(m1=0;m1<egl->nfreq;m1++) {
    rt = 0;
    for(t=1;t<ntrunc;t++) {
      It[m1*ntrunc+t] = 0;
    }
    for(p=2;p<nmod;p++) {
      sprintf(name,"beamnl_r_%d_%d",(int)egl->freqlist[m1],p);
      r =  parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      //_DEBUGHERE_("r %d %d = %g",m1,p,r);
      rt +=r;
      pt = 1;
      for(t=1;t<ntrunc;t++) {
        pt *= p;
        It[m1*ntrunc+t] += r/pt; 
        //_DEBUGHERE_("I %d(%d) %d = %g",m1,t,p,It[m1*ntrunc+t]);
      }

    }
    r = 1-rt;
    //_DEBUGHERE_("r %d %d = %g",m1,1,r);
    for(t=1;t<ntrunc;t++) {
      It[m1*ntrunc+t] += r; 
      //_DEBUGHERE_("I %d %d = %g",m1,t,It[m1*ntrunc+t]);
    }
  }


  lmin =  egl->lmin;
  lmax =  egl->lmax;
  
  for(ell=lmin;ell<=lmax;ell++) {
    for(m1=0;m1<egl->nfreq;m1++) {
      d1 = 0;
      t1 = template[ell+(lmax_template+1)*mv[m1]];
      lt1 = log(t1);
      tp1 = 1./t1;
      for(t=1;t<ntrunc;t++) {
        tp1 *= lt1/t;
        d1 += tp1 * (It[m1*ntrunc+t]-1);
        //if (ell==lmin) _DEBUGHERE_("%d(%d) %g",m1,t,d1);
      }
      d[m1]=d1;
      //if (ell==lmin) _DEBUGHERE_("%d %g",m1,d1);
    }
    for(m1=0;m1<egl->nfreq;m1++) {    
      for(m2=m1;m2<egl->nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = 1 - d[m1] - d[m2] + d[m1]*d[m2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
}
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(beamnl,beamnl_init);
