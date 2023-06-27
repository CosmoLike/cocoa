#include "clik_parametric.h"
#include "clik_parametric_addon.h"

void fill_offset_freq_TP(int idreq,double *dreq, int nfreq, double* freqlist,int *mv,int off, error **err);

void nslb_compute(parametric* egl, double *Rq, error **err);

parametric *nslb_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int *mv,m1,m2,f1;
  double dreq[4];
  pfchar name;
  char tp[2];
  
  testErrorRet(template==NULL,-11111,"template must be defined",*err,__LINE__,NULL);
  testErrorRet(ndet_T!=0 && has_TEB[0]==0,-11111,"incoherent definition, ndet_T !0 and hascl[0]=1",*err,__LINE__,NULL);
  testErrorRet(ndet_P!=0 && (has_TEB[1]==0 && has_TEB[2]==0),-11111,"incoherent definition, ndet_P !0 and hascl[1]=0 and hascl[2]=0",*err,__LINE__,NULL);


  // init
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);

  egl->payload = malloc_err(sizeof(double)* (3001*12*12 + 3001*8 + 2*8 + 12) + sizeof(int)*egl->nfreq*2,err);
  forwardError(*err,__LINE__,NULL);

  memcpy(egl->payload,template,sizeof(double)* (3001*12*12 + 3001*8));

  mv = egl->payload + sizeof(double)* (3001*12*12 + 3001*8 +2*8);

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  fill_offset_freq_TP(4,dreq, egl->nfreq_T*has_TEB[0],egl->freqlist_T,mv,0,err);
  forwardError(*err,__LINE__,NULL);
  fill_offset_freq_TP(4,dreq, egl->nfreq_P*has_TEB[1],egl->freqlist_P,mv + egl->nfreq_T*has_TEB[0],4,err);
  forwardError(*err,__LINE__,NULL);
  fill_offset_freq_TP(4,dreq, egl->nfreq_P*has_TEB[2],egl->freqlist_P,mv + egl->nfreq_P*has_TEB[1]+egl->nfreq_T*has_TEB[0],4,err);
  forwardError(*err,__LINE__,NULL);

  if (has_TEB[1]!=0 || has_TEB[2]!=0) {
    fill_offset_freq_TP(4,dreq, egl->nfreq_P,egl->freqlist_P,mv + egl->nfreq_T*has_TEB[0],4,err);
    forwardError(*err,__LINE__,NULL);
  }

  tp[0] = 'T';
  tp[1] = 'P';

  for(m1=0;m1<4;m1++) {
    for(f1=0;f1<2;f1++) {
      sprintf(name,"nslb_epsilon_%d_%c",(int)dreq[m1],tp[f1]);
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);      
      sprintf(name,"nslb_fwhm_%d_%c",(int)dreq[m1],tp[f1]);
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);      
    }
  }


  egl->eg_compute = &nslb_compute;
  egl->eg_free = &parametric_simple_payload_free;
    
 return egl;
}


void nslb_compute(parametric* egl, double *Rq, error **err) {
  double *template,*bl,*epsilon,*sigma;
  int *mv;
  double b,v;
  int m1,m2,f1,f2,ell;
  int dreq[4];
  char tp[2];
  pfchar name;
  int m1p,m2p;
  double bp;        

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  tp[0] = 'T';
  tp[1] = 'P';

  mv = egl->payload + sizeof(double)* (3001*12*12 + 3001*8 + 2*8);
  epsilon = egl->payload + sizeof(double)* (3001*12*12 + 3001*8);
  sigma = egl->payload + sizeof(double)* (3001*12*12+ + 3001*8 +8);

  template = egl->payload;
  bl = egl->payload + sizeof(double)* (3001*12*12);

  /*for(m1=0;m1<egl->nfreq;m1++) {
    _DEBUGHERE_("%d %d %d",m1,mv[m1],egl->nfreq);
  }*/
  for(m1=0;m1<4;m1++) {
    for(f1=0;f1<2;f1++) {
      sprintf(name,"nslb_epsilon_%d_%c",(int)dreq[m1],tp[f1]);
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      //_DEBUGHERE_("epsilon %d %c %g",dreq[m1],tp[f1],v)
      epsilon[m1+f1*4] = v;
      sprintf(name,"nslb_fwhm_%d_%c",(int)dreq[m1],tp[f1]);
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      //_DEBUGHERE_("fwhm %d %c %g",dreq[m1],tp[f1],v)
      sigma[m1+f1*4] = v/sqrt(8*log(2))/60/180*M_PI;
    }
  }  
  
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        /*if(ell==100) {
            _DEBUGHERE_("%d %d %d %d ----------",m1,m2,mv[m1],mv[m2]);
            _DEBUGHERE_("%g %g %g %g",epsilon[mv[m1]],sigma[mv[m1]],epsilon[mv[m2]],sigma[mv[m2]]);
        }*/
        b = epsilon[mv[m1]]*exp(-.5*ell*(ell+1)*sigma[mv[m1]]*sigma[mv[m1]])/bl[ell*8+mv[m1]] + epsilon[mv[m2]]*exp(-.5*ell*(ell+1)*sigma[mv[m2]]*sigma[mv[m2]])/bl[ell*8+mv[m2]];
        if (m1<egl->nfreq_T*egl->has_TEB[0] && m2>=egl->nfreq_T*egl->has_TEB[0]) {
          // TP case, need to symetrise;
          m1p = m1+egl->nfreq_T*egl->has_TEB[0];
          m2p = m2-egl->nfreq_T*egl->has_TEB[0];
          bp = epsilon[mv[m1p]]*exp(-.5*ell*(ell+1)*sigma[mv[m1p]]*sigma[mv[m1p]])/bl[ell*8+mv[m1p]] + epsilon[mv[m2p]]*exp(-.5*ell*(ell+1)*sigma[mv[m2p]]*sigma[mv[m2p]])/bl[ell*8+mv[m2p]];
          /*if(ell==100) {
            _DEBUGHERE_("ext %d %d %d %d",m1p,m2p,mv[m1p],mv[m2p]);
            _DEBUGHERE_("%g %g %g %g -> %g %g",epsilon[mv[m1p]],sigma[mv[m1p]],epsilon[mv[m2p]],sigma[mv[m2p]],b,bp);
          } */  
          b = .5*(b+bp);
          
        }
        /*if(ell==100) {
          _DEBUGHERE_("%g ------------",b);
          //_DEBUGHERE_("%d %d %d %g",ell,mv[m1],mv[m2],template[ell*12*12+mv[m1]*12+mv[m2]]);  
        }*/
        
        Rq[IDX_R(egl,ell,m1,m2)] = b*template[ell*12*12+mv[m1]*12+mv[m2]];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }

  return;
  
}

CREATE_PARAMETRIC_POL_TEMPLATE_FILE_INIT(nslb,nslb_init);

void trans_compute(parametric* egl, double *Rq, error **err);

parametric *trans_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int *mv,m1,m2,f1;
  double dreq[4];
  pfchar name;
  char tp[2];
  
  testErrorRet(template==NULL,-11111,"template must be defined",*err,__LINE__,NULL);
  testErrorRet(ndet_T!=0 && has_TEB[0]==0,-11111,"incoherent definition, ndet_T !0 and hascl[0]=1",*err,__LINE__,NULL);
  testErrorRet(ndet_P!=0 && (has_TEB[1]==0 && has_TEB[2]==0),-11111,"incoherent definition, ndet_P !0 and hascl[1]=0 and hascl[2]=0",*err,__LINE__,NULL);


  // init
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);

  egl->payload = malloc_err(sizeof(double)* (3001*12*12 + 3001*8 + 2*8 + 12) + sizeof(int)*egl->nfreq*2,err);
  forwardError(*err,__LINE__,NULL);

  memcpy(egl->payload,template,sizeof(double)* (3001*12*12 + 3001*8));

  mv = egl->payload + sizeof(double)* (3001*12*12 + 3001*8 +2*8);

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  fill_offset_freq_TP(4,dreq, egl->nfreq_T*has_TEB[0],egl->freqlist_T,mv,0,err);
  forwardError(*err,__LINE__,NULL);
  fill_offset_freq_TP(4,dreq, egl->nfreq_P*has_TEB[1],egl->freqlist_P,mv + egl->nfreq_T*has_TEB[0],4,err);
  forwardError(*err,__LINE__,NULL);
  fill_offset_freq_TP(4,dreq, egl->nfreq_P*has_TEB[2],egl->freqlist_P,mv + egl->nfreq_P*has_TEB[1]+egl->nfreq_T*has_TEB[0],4,err);
  forwardError(*err,__LINE__,NULL);

  if (has_TEB[1]!=0 || has_TEB[2]!=0) {
    fill_offset_freq_TP(4,dreq, egl->nfreq_P,egl->freqlist_P,mv + egl->nfreq_T*has_TEB[0],4,err);
    forwardError(*err,__LINE__,NULL);
  }

  tp[0] = 'T';
  tp[1] = 'P';

  for(m1=0;m1<4;m1++) {
    for(f1=0;f1<2;f1++) {
      sprintf(name,"trans_epsilon0_%d_%c",(int)dreq[m1],tp[f1]);
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);      
      sprintf(name,"trans_epsilon2_%d_%c",(int)dreq[m1],tp[f1]);
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);      
    }
  }


  egl->eg_compute = &trans_compute;
  egl->eg_free = &parametric_simple_payload_free;
    
 return egl;
}


void trans_compute(parametric* egl, double *Rq, error **err) {
  double *template,*bl,*epsilon0,*epsilon2;
  int *mv;
  double b,v;
  int m1,m2,f1,f2,ell;
  int dreq[4];
  char tp[2];
  pfchar name;
  int m1p,m2p;
  double bp;        

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  tp[0] = 'T';
  tp[1] = 'P';

  mv = egl->payload + sizeof(double)* (3001*12*12 + 3001*8 + 2*8);
  epsilon0 = egl->payload + sizeof(double)* (3001*12*12 + 3001*8);
  epsilon2 = egl->payload + sizeof(double)* (3001*12*12+ + 3001*8 +8);

  template = egl->payload;
  bl = egl->payload + sizeof(double)* (3001*12*12);

  /*for(m1=0;m1<egl->nfreq;m1++) {
    _DEBUGHERE_("%d %d %d",m1,mv[m1],egl->nfreq);
  }*/
  for(m1=0;m1<4;m1++) {
    for(f1=0;f1<2;f1++) {
      sprintf(name,"trans_epsilon0_%d_%c",(int)dreq[m1],tp[f1]);
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      //_DEBUGHERE_("epsilon %d %c %g",dreq[m1],tp[f1],v)
      epsilon0[m1+f1*4] = v;
      sprintf(name,"trans_epsilon2_%d_%c",(int)dreq[m1],tp[f1]);
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      //_DEBUGHERE_("fwhm %d %c %g",dreq[m1],tp[f1],v)
      epsilon2[m1+f1*4] = v;
    }
  }  
  
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        /*if(ell==100) {
            _DEBUGHERE_("%d %d %d %d ----------",m1,m2,mv[m1],mv[m2]);
            _DEBUGHERE_("%g %g %g %g",epsilon[mv[m1]],sigma[mv[m1]],epsilon[mv[m2]],sigma[mv[m2]]);
        }*/
        b = (epsilon0[mv[m1]]+epsilon2[mv[m1]]*ell*ell + epsilon0[mv[m2]]+epsilon2[mv[m2]]*ell*ell)/2.;
        if (m1<egl->nfreq_T*egl->has_TEB[0] && m2>=egl->nfreq_T*egl->has_TEB[0]) {
          // TP case, need to symetrise;
          m1p = m1+egl->nfreq_T*egl->has_TEB[0];
          m2p = m2-egl->nfreq_T*egl->has_TEB[0];
          bp = (epsilon0[mv[m1p]]+ell*ell*epsilon2[mv[m1p]] + epsilon0[mv[m2p]]+epsilon2[mv[m2p]]*ell*ell)/2.;
          /*if(ell==100) {
            _DEBUGHERE_("ext %d %d %d %d",m1p,m2p,mv[m1p],mv[m2p]);
            _DEBUGHERE_("%g %g %g %g -> %g %g",epsilon[mv[m1p]],sigma[mv[m1p]],epsilon[mv[m2p]],sigma[mv[m2p]],b,bp);
          } */  
          b = .5*(b+bp);
          
        }
        /*if(ell==100) {
          _DEBUGHERE_("%g ------------",b);
          //_DEBUGHERE_("%d %d %d %g",ell,mv[m1],mv[m2],template[ell*12*12+mv[m1]*12+mv[m2]]);  
        }*/
        //_DEBUGHERE_("%g",b);
        Rq[IDX_R(egl,ell,m1,m2)] = b*template[ell*12*12+mv[m1]*12+mv[m2]];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }

  return;
  
}

CREATE_PARAMETRIC_POL_TEMPLATE_FILE_INIT(trans,trans_init);

void trans2_compute(parametric* egl, double *Rq, error **err);

parametric *trans2_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int *mv,m1,m2,f1;
  double dreq[4];
  pfchar name;
  char tp[2];
  
  testErrorRet(template==NULL,-11111,"template must be defined",*err,__LINE__,NULL);
  testErrorRet(ndet_T!=0 && has_TEB[0]==0,-11111,"incoherent definition, ndet_T !0 and hascl[0]=1",*err,__LINE__,NULL);
  testErrorRet(ndet_P!=0 && (has_TEB[1]==0 && has_TEB[2]==0),-11111,"incoherent definition, ndet_P !0 and hascl[1]=0 and hascl[2]=0",*err,__LINE__,NULL);



  // init
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);

  egl->payload = malloc_err(sizeof(double)* (3001*12*12 + 3001*8 + 2*8 + 12) + sizeof(int)*egl->nfreq*2,err);
  forwardError(*err,__LINE__,NULL);

  memcpy(egl->payload,template,sizeof(double)* (3001*12*12 + 3001*8));

  mv = egl->payload + sizeof(double)* (3001*12*12 + 3001*8 +2*8);

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  fill_offset_freq_TP(4,dreq, egl->nfreq_T*has_TEB[0],egl->freqlist_T,mv,0,err);
  forwardError(*err,__LINE__,NULL);
  fill_offset_freq_TP(4,dreq, egl->nfreq_P*has_TEB[1],egl->freqlist_P,mv + egl->nfreq_T*has_TEB[0],4,err);
  forwardError(*err,__LINE__,NULL);
  fill_offset_freq_TP(4,dreq, egl->nfreq_P*has_TEB[2],egl->freqlist_P,mv + egl->nfreq_P*has_TEB[1]+egl->nfreq_T*has_TEB[0],4,err);
  forwardError(*err,__LINE__,NULL);

  if (has_TEB[1]!=0 || has_TEB[2]!=0) {
    fill_offset_freq_TP(4,dreq, egl->nfreq_P,egl->freqlist_P,mv + egl->nfreq_T*has_TEB[0],4,err);
    forwardError(*err,__LINE__,NULL);
  }

  tp[0] = 'T';
  tp[1] = 'P';

  for(m1=0;m1<4;m1++) {
    for(f1=0;f1<2;f1++) {
      sprintf(name,"trans_epsilon0_%d_%c",(int)dreq[m1],tp[f1]);
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);      
      sprintf(name,"trans_epsilon2_%d_%c",(int)dreq[m1],tp[f1]);
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);      
    }
  }


  egl->eg_compute = &trans2_compute;
  egl->eg_free = &parametric_simple_payload_free;
    
 return egl;
}


void trans2_compute(parametric* egl, double *Rq, error **err) {
  double *template,*bl,*epsilon0,*epsilon2;
  int *mv;
  double b,v;
  int m1,m2,f1,f2,ell;
  int dreq[4];
  char tp[2];
  pfchar name;
  int m1p,m2p;
  double bp;        

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  tp[0] = 'T';
  tp[1] = 'P';

  mv = egl->payload + sizeof(double)* (3001*12*12 + 3001*8 + 2*8);
  epsilon0 = egl->payload + sizeof(double)* (3001*12*12 + 3001*8);
  epsilon2 = egl->payload + sizeof(double)* (3001*12*12+ + 3001*8 +8);

  template = egl->payload;
  bl = egl->payload + sizeof(double)* (3001*12*12);

  /*for(m1=0;m1<egl->nfreq;m1++) {
    _DEBUGHERE_("%d %d %d",m1,mv[m1],egl->nfreq);
  }*/
  for(m1=0;m1<4;m1++) {
    for(f1=0;f1<2;f1++) {
      sprintf(name,"trans_epsilon0_%d_%c",(int)dreq[m1],tp[f1]);
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      //_DEBUGHERE_("epsilon %d %c %g",dreq[m1],tp[f1],v)
      epsilon0[m1+f1*4] = v;
      sprintf(name,"trans_epsilon2_%d_%c",(int)dreq[m1],tp[f1]);
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      //_DEBUGHERE_("fwhm %d %c %g",dreq[m1],tp[f1],v)
      epsilon2[m1+f1*4] = v;
    }
  }  
  
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        /*if(ell==100) {
            _DEBUGHERE_("%d %d %d %d ----------",m1,m2,mv[m1],mv[m2]);
            _DEBUGHERE_("%g %g %g %g",epsilon[mv[m1]],sigma[mv[m1]],epsilon[mv[m2]],sigma[mv[m2]]);
        }*/
        b = sqrt((epsilon0[mv[m1]]+epsilon2[mv[m1]]*ell*ell) * (epsilon0[mv[m2]]+epsilon2[mv[m2]]*ell*ell));

        if (m1<egl->nfreq_T*egl->has_TEB[0] && m2>=egl->nfreq_T*egl->has_TEB[0]) {
          // TP case, need to symetrise;
          m1p = m1+egl->nfreq_T*egl->has_TEB[0];
          m2p = m2-egl->nfreq_T*egl->has_TEB[0];
          bp = sqrt((epsilon0[mv[m1p]]+ell*ell*epsilon2[mv[m1p]]) * (epsilon0[mv[m2p]]+epsilon2[mv[m2p]]*ell*ell));
          /*if(ell==100) {
            _DEBUGHERE_("ext %d %d %d %d",m1p,m2p,mv[m1p],mv[m2p]);
            _DEBUGHERE_("%g %g %g %g -> %g %g",epsilon[mv[m1p]],sigma[mv[m1p]],epsilon[mv[m2p]],sigma[mv[m2p]],b,bp);
          } */  
          b = .5*(b+bp);
          
        }
        /*if(ell==100) {
          _DEBUGHERE_("%g ------------",b);
          //_DEBUGHERE_("%d %d %d %g",ell,mv[m1],mv[m2],template[ell*12*12+mv[m1]*12+mv[m2]]);  
        }*/
        //_DEBUGHERE_("%g",b);
        Rq[IDX_R(egl,ell,m1,m2)] = b*template[ell*12*12+mv[m1]*12+mv[m2]];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }

  return;
  
}

CREATE_PARAMETRIC_POL_TEMPLATE_FILE_INIT(trans2,trans2_init);

void transG_compute(parametric* egl, double *Rq, error **err);

parametric *transG_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int *mv,m1,m2,f1;
  double dreq[4];
  pfchar name;
  char tp[2];
  int deg,ideg;

  testErrorRet(template==NULL,-11111,"template must be defined",*err,__LINE__,NULL);
  testErrorRet(ndet_T!=0 && has_TEB[0]==0,-11111,"incoherent definition, ndet_T !0 and hascl[0]=1",*err,__LINE__,NULL);
  testErrorRet(ndet_P!=0 && (has_TEB[1]==0 && has_TEB[2]==0),-11111,"incoherent definition, ndet_P !0 and hascl[1]=0 and hascl[2]=0",*err,__LINE__,NULL);



  // init
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"transG_deg",2,err);
  forwardError(*err,__LINE__,NULL);   
  deg = parametric_get_value(egl,"transG_deg",err);
  forwardError(*err,__LINE__,NULL);


  egl->payload = malloc_err(sizeof(double)* (3001*12*12 + 3001*8 + deg*8 + 12) + sizeof(int)*egl->nfreq*2,err);
  forwardError(*err,__LINE__,NULL);

  memcpy(egl->payload,template,sizeof(double)* (3001*12*12 + 3001*8));

  mv = egl->payload + sizeof(double)* (3001*12*12 + 3001*8 +deg*8);

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  fill_offset_freq_TP(4,dreq, egl->nfreq_T*has_TEB[0],egl->freqlist_T,mv,0,err);
  forwardError(*err,__LINE__,NULL);
  fill_offset_freq_TP(4,dreq, egl->nfreq_P*has_TEB[1],egl->freqlist_P,mv + egl->nfreq_T*has_TEB[0],4,err);
  forwardError(*err,__LINE__,NULL);
  fill_offset_freq_TP(4,dreq, egl->nfreq_P*has_TEB[2],egl->freqlist_P,mv + egl->nfreq_P*has_TEB[1]+egl->nfreq_T*has_TEB[0],4,err);
  forwardError(*err,__LINE__,NULL);

  if (has_TEB[1]!=0 || has_TEB[2]!=0) {
    fill_offset_freq_TP(4,dreq, egl->nfreq_P,egl->freqlist_P,mv + egl->nfreq_T*has_TEB[0],4,err);
    forwardError(*err,__LINE__,NULL);
  }

  tp[0] = 'T';
  tp[1] = 'P';

      
  for(m1=0;m1<4;m1++) {
    for(f1=0;f1<2;f1++) {
      for(ideg=1;ideg<deg;ideg++) {
        sprintf(name,"transG_epsilon%d_%d_%c",ideg,(int)dreq[m1],tp[f1]);
        parametric_set_default(egl,name,0,err);
        forwardError(*err,__LINE__,NULL);      
      } 
    }
  }


  egl->eg_compute = &transG_compute;
  egl->eg_free = &parametric_simple_payload_free;
    
 return egl;
}


void transG_compute(parametric* egl, double *Rq, error **err) {
  double *template,*bl,*epsilon;
  int *mv;
  double b,v;
  int m1,m2,f1,f2,ell;
  int dreq[4];
  char tp[2];
  pfchar name;
  int m1p,m2p;
  double bp,b1,b2,b1p,b2p,ln;        
  int deg,ideg;

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  tp[0] = 'T';
  tp[1] = 'P';

  deg = parametric_get_value(egl,"transG_deg",err);
  forwardError(*err,__LINE__,);

  mv = egl->payload + sizeof(double)* (3001*12*12 + 3001*8 + deg*8);
  epsilon = egl->payload + sizeof(double)* (3001*12*12 + 3001*8);
  
  template = egl->payload;
  bl = egl->payload + sizeof(double)* (3001*12*12);


  /*for(m1=0;m1<egl->nfreq;m1++) {
    _DEBUGHERE_("%d %d %d",m1,mv[m1],egl->nfreq);
  }*/
  for(m1=0;m1<4;m1++) {
    for(f1=0;f1<2;f1++) {
      for(ideg=1;ideg<deg;ideg++) {
        sprintf(name,"transG_epsilon%d_%d_%c",ideg,(int)dreq[m1],tp[f1]);
        v = parametric_get_value(egl,name,err);
        forwardError(*err,__LINE__,);
        epsilon[ideg+m1*deg+f1*4*deg] = v;
        //_DEBUGHERE_("%s %d %d %d %g",name,m1,f1,deg,v);
      }
    }
  }  
  
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    for(m1=0;m1<egl->nfreq;m1++) {
      b1 = 1;
      ln=1;
      for (ideg=1;ideg<deg;ideg++){
        ln*=ell;
        b1 += epsilon[ideg+m1*deg]*ln;
        if(ell==100) {
        //_DEBUGHERE_("%d %d %d -> %g %d %d %g",ell,m1,m2,b1,ideg,ln,epsilon[ideg+m1*deg]*ln)
        }
      }
      for(m2=m1;m2<egl->nfreq;m2++) {
        b2 = 1;
        ln=1;
        for (ideg=1;ideg<deg;ideg++){
          ln*=ell;
          b2 += epsilon[ideg+m2*deg]*ln;
        }
        
        if(ell==100) {
            //_DEBUGHERE_("%d %d %d %d %g %g----------",m1,m2,mv[m1],mv[m2],b1,b2);
            //  =_DEBUGHERE_("%g %g %g %g",epsilon[mv[m1]],sigma[mv[m1]],epsilon[mv[m2]],sigma[mv[m2]]);
        }
        b = b1*b2-1;

        //_DEBUGHERE_("%g %g %g",b,b1,b2);
        
        if (m1<egl->nfreq_T*egl->has_TEB[0] && m2>=egl->nfreq_T*egl->has_TEB[0]) {
          // TP case, need to symetrise;
          m1p = m1+egl->nfreq_T*egl->has_TEB[0];
          m2p = m2-egl->nfreq_T*egl->has_TEB[0];
          b1p = 1;
          ln=1;
          for (ideg=1;ideg<deg;ideg++){
            ln*=ell;
            b1p += epsilon[ideg+m1p*deg]*ln;
          }   
          b2p = 1;
          ln=1;
          for (ideg=1;ideg<deg;ideg++){
            ln*=ell;
            b2p += epsilon[ideg+m2p*deg]*ln;
          }   

          bp = b1p*b2p-1;
          /*if(ell==100) {
            _DEBUGHERE_("ext %d %d %d %d",m1p,m2p,mv[m1p],mv[m2p]);
            _DEBUGHERE_("%g %g %g %g -> %g %g",epsilon[mv[m1p]],sigma[mv[m1p]],epsilon[mv[m2p]],sigma[mv[m2p]],b,bp);
          } */  
          b = .5*(b+bp);
          //_DEBUGHERE_("%g %g %g %g ",b,bp,b1p,b2p);
          
        }
        if(ell==100) {
          //_DEBUGHERE_("%g ------------",b);
          //_DEBUGHERE_("%d %d %d %g",ell,mv[m1],mv[m2],template[ell*12*12+mv[m1]*12+mv[m2]]);  
        }
        Rq[IDX_R(egl,ell,m1,m2)] = b*template[ell*12*12+mv[m1]*12+mv[m2]];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }

  return;
  
}

CREATE_PARAMETRIC_POL_TEMPLATE_FILE_INIT(transG,transG_init);

