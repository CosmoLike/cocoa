#include "clik_parametric.h"
#include "clik_parametric_addon.h"

void ffp8_gal_compute(parametric *egl, double *rq, error **err);

void ffp8_fill_offset_freq_TP(int idreq,double *dreq, int nfreq, double* freqlist,int *mv,int off, error **err) {
  int m,i;
  
  for(m=0;m<nfreq;m++) {
    double f;
    f = freqlist[m];
    mv[m]=-1;
    for(i=0;i<idreq;i++) {
      if (fabs(f-dreq[i])<1e-6) {
        mv[m] = i + off;
        break;
      }  
    }
    testErrorRetVA(mv[m]==-1,-431432,"Don't know how to compute component for freq %g",*err,__LINE__,,f);
  }    
}

void fill_offset_freq(int idreq,double *dreq, parametric *egl,int *mv,int def, error **err);


parametric *ffp8_gal_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int i,m,*mv,m1,m2;
  double dreq[4];
  pfchar name;
  double *conv,*A;
  char teb[3];
  int f1,f2;

  testErrorRet(template==NULL,-11111,"template must be defined",*err,__LINE__,NULL);
  testErrorRet(ndet_T!=0 && has_TEB[0]==0,-11111,"incoherent definition, ndet_T !0 and hascl[0]=1",*err,__LINE__,NULL);
  testErrorRet(ndet_P!=0 && (has_TEB[1]==0 && has_TEB[2]==0),-11111,"incoherent definition, ndet_P !0 and hascl[1]=0 and hascl[2]=0",*err,__LINE__,NULL);


  // init
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  
  egl->eg_compute = &ffp8_gal_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)* (3001*12*12 + 12*12) + sizeof(int)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);


  memcpy(egl->payload,template,sizeof(double)* (3001*12*12));

  mv = egl->payload + sizeof(double)* (3001*12*12 + 12*12);

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  ffp8_fill_offset_freq_TP(4,dreq, egl->nfreq_T*has_TEB[0],egl->freqlist_T,mv,0,err);
  forwardError(*err,__LINE__,NULL);
  ffp8_fill_offset_freq_TP(4,dreq, egl->nfreq_P*has_TEB[1],egl->freqlist_P,mv + egl->nfreq_T*has_TEB[0],4,err);
  forwardError(*err,__LINE__,NULL);
  ffp8_fill_offset_freq_TP(4,dreq, egl->nfreq_P*has_TEB[2],egl->freqlist_P,mv + egl->nfreq_P*has_TEB[1]+egl->nfreq_T*has_TEB[0],8,err);
  forwardError(*err,__LINE__,NULL);

  teb[0] = 'T';
  teb[1] = 'E';
  teb[2] = 'B';

  for(m1=0;m1<4;m1++) {
    for(m2=m1;m2<4;m2++) {
      for(f1=0;f1<3;f1++) {
        for(f2=f1;f2<3;f2++) {
          sprintf(name,"A_ffp8_gal_%d_%d_%c%c",(int)dreq[m1],(int)dreq[m2],teb[f1],teb[f2]);
          parametric_set_default(egl,name,0,err);
          forwardError(*err,__LINE__,NULL);      
        }
      }
    }
  }  

  parametric_set_default(egl,"ffp8_gal_l_pivot",2000,err); 
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void ffp8_gal_compute(parametric *egl, double *Rq, error **err) {
  double *template;
  double l_pivot,index,v;
  int m1,m2,ell;
  double nrm;
  int *mv;
  int rigid,irigid;
  double *conv,*A;
  pfchar name;
  double dreq[4];
  char teb[3];
  int rm1,rm2;
  int f1,f2;
  int abso;

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  teb[0] = 'T';
  teb[1] = 'E';
  teb[2] = 'B';

  mv = egl->payload + sizeof(double)* (3001*12*12 + 12*12);
  A = egl->payload + sizeof(double)* (3001*12*12);

  template = egl->payload;
  l_pivot = parametric_get_value(egl,"ffp8_gal_l_pivot",err);
  forwardError(*err,__LINE__,);;

  if (abso!=0) {
    abso = 1;
  }

  for(m1=0;m1<4;m1++) {
    for(m2=m1;m2<4;m2++) {
      for(f1=0;f1<3;f1++) {
        for(f2=f1;f2<3;f2++) {
          sprintf(name,"A_ffp8_gal_%d_%d_%c%c",(int)dreq[m1],(int)dreq[m2],teb[f1],teb[f2]);
          v = parametric_get_value(egl,name,err);
          forwardError(*err,__LINE__,);
          rm1 = m1 + f1*4;
          rm2 = m2 + f2*4;
          A[rm1*12+rm2] = (template[((int) l_pivot)*12*12+rm1*12+rm2]*l_pivot*(l_pivot+1)/2/M_PI)==0 ? 0. : v/(template[((int) l_pivot)*12*12+rm1*12+rm2]*l_pivot*(l_pivot+1)/2/M_PI);  
          A[rm2*12+rm1] = A[rm1*12+rm2];
          //_DEBUGHERE_("%d %d %g %g",rm1,rm2,v,template[((int) l_pivot)*12*12+rm1*12+rm2]);
          //_DEBUGHERE_("%s %g %g %d %d %g",name,v,(1-abso+abso*template[((int) l_pivot)*12*12+rm1*12+rm2]/l_pivot/(l_pivot+1)*2*M_PI),rm1,rm2,A[rm1*12+rm2]);
          rm1 = m1 + f2*4;
          rm2 = m2 + f1*4;
          //_DEBUGHERE_("%d %d %g %g",rm1,rm2,v,template[((int) l_pivot)*12*12+rm1*12+rm2]);
          A[rm1*12+rm2] = (template[((int) l_pivot)*12*12+rm1*12+rm2]*l_pivot*(l_pivot+1)/2/M_PI)==0 ? 0. : v/(template[((int) l_pivot)*12*12+rm1*12+rm2]*l_pivot*(l_pivot+1)/2/M_PI);  
          A[rm2*12+rm1] = A[rm1*12+rm2];
          //_DEBUGHERE_("%s %g %g %d %d %g",name,v,(1-abso+abso*template[((int) l_pivot)*12*12+rm1*12+rm2]/l_pivot/(l_pivot+1)*2*M_PI),rm1,rm2,A[rm1*12+rm2]);
          if ((f1!=f2) && (m1!=m2)) {
            rm1 = m2 + f1*4;
            rm2 = m1 + f2*4;
            A[rm1*12+rm2] = (template[((int) l_pivot)*12*12+rm1*12+rm2]*l_pivot*(l_pivot+1)/2/M_PI)==0 ? 0. : v/(template[((int) l_pivot)*12*12+rm1*12+rm2]*l_pivot*(l_pivot+1)/2/M_PI);  
            A[rm2*12+rm1] = A[rm1*12+rm2];
            //_DEBUGHERE_("%d %d %g %g",rm1,rm2,v,template[((int) l_pivot)*12*12+rm1*12+rm2]);
            //_DEBUGHERE_("%s %g %g %d %d %g",name,v,(1-abso+abso*template[((int) l_pivot)*12*12+rm1*12+rm2]/l_pivot/(l_pivot+1)*2*M_PI),rm1,rm2,A[rm1*12+rm2]);
            rm1 = m2 + f2*4;
            rm2 = m1 + f1*4;
            //_DEBUGHERE_("%d %d %g %g",rm1,rm2,v,template[((int) l_pivot)*12*12+rm1*12+rm2]);
            A[rm1*12+rm2] = (template[((int) l_pivot)*12*12+rm1*12+rm2]*l_pivot*(l_pivot+1)/2/M_PI)==0 ? 0. : v/(template[((int) l_pivot)*12*12+rm1*12+rm2]*l_pivot*(l_pivot+1)/2/M_PI);  
            A[rm2*12+rm1] = A[rm1*12+rm2];
            
          }      
        }
      }
    }
  }  
  
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        if(ell==egl->lmin) {
          //_DEBUGHERE_("%d %d %d %g",ell,mv[m1],mv[m2],template[ell*12*12+mv[m1]*12+mv[m2]]);  
        }
        
        Rq[IDX_R(egl,ell,m1,m2)] = template[ell*12*12+mv[m1]*12+mv[m2]] * A[mv[m1]*12+mv[m2]];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
}

CREATE_PARAMETRIC_POL_TEMPLATE_FILE_INIT(ffp8_gal,ffp8_gal_init);

void ffp8_gcib_compute(parametric *egl, double *rq, error **err);

parametric *ffp8_gcib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int i,m,*mv,m1,m2;
  double dreq[4];
  pfchar name;
  double *conv,*A;

  testErrorRet(template==NULL,-11111,"template must be defined",*err,__LINE__,NULL);
  

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);

  egl->eg_compute = &ffp8_gcib_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)* (10001*4*4) + sizeof(double)*(4+egl->nfreq*egl->nfreq)+sizeof(int)*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);


  memcpy(egl->payload,template,sizeof(double)* (10001*4*4));

  mv = egl->payload + sizeof(double)* (10001*4*4) + sizeof(double)*(4+egl->nfreq*egl->nfreq);
  conv = egl->payload + sizeof(double)* (10001*4*4);

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;
  fill_offset_freq(4,dreq, egl,mv,-1,err);
  forwardError(*err,__LINE__,NULL);

  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"A_gib_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"A_gib_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      parametric_set_default(egl,name,70,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  

  parametric_set_default(egl,"gib_index",-1.3,err); 
  forwardError(*err,__LINE__,NULL);

  conv[0] = 1;
  conv[1] = 1;
  conv[2] = 1;
  conv[3] = 1;

  
  parametric_set_default(egl,"gib_rigid",217,err); 
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gib_l_pivot",3000,err); 
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void ffp8_gcib_compute(parametric *egl, double *Rq, error **err) {
  double *template;
  double l_pivot,index,v;
  int m1,m2,ell;
  double nrm;
  int *mv;
  int rigid,irigid;
  double *conv,*A;
  pfchar name;
  

  mv = egl->payload + sizeof(double)* (10001*4*4) + sizeof(double)*(4+egl->nfreq*egl->nfreq);
  conv = egl->payload + sizeof(double)* (10001*4*4);
  A = egl->payload + sizeof(double)* (10001*4*4)+sizeof(double)*4;

  template = egl->payload;
  l_pivot = parametric_get_value(egl,"gib_l_pivot",err);
  forwardError(*err,__LINE__,);;

  rigid = parametric_get_value(egl,"gib_rigid",err);
  forwardError(*err,__LINE__,);
  //_DEBUGHERE_("rigid %d",rigid);

  index = parametric_get_value(egl,"gib_index",err);
  forwardError(*err,__LINE__,);

  if (rigid==0) {
    sprintf(name,"A_gib_%d",217);    
    irigid = 2;
  } else {
    double dreq[4];

    dreq[0] = 100;
    dreq[1] = 143;
    dreq[2] = 217;
    dreq[3] = 353;

    sprintf(name,"A_gib_%d",(int) rigid);  
    irigid=-1;
    for(m1=0;m1<4;m1++) {
      if (fabs(rigid-dreq[m1])<1e-6) {
        irigid=m1;
      }
    }
    testErrorRet(irigid==-1,-55214,"AAAAAAA",*err,__LINE__,);
    rigid = 1;
    
  }


  nrm = parametric_get_value(egl,name,err);
  forwardError(*err,__LINE__,);
  //_DEBUGHERE_("%s %g",name,nrm);
    
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"A_gib_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"A_gib_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*egl->nfreq+m2] = (v/template[((int) l_pivot)*16+mv[m1]*4+mv[m2]]*(1-rigid) +  nrm/template[((int)l_pivot)*16+irigid*4+irigid]*rigid*conv[mv[m1]]*conv[mv[m2]]/conv[irigid]/conv[irigid]) /l_pivot/(l_pivot+1)*2*M_PI ;
      A[m2*egl->nfreq+m1] = A[m1*egl->nfreq+m2];
    }
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index-(-1.3));
    //_DEBUGHERE_("%d %g",ell,v);
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        //_DEBUGHERE_("%d %d %d %g",ell,mv[m1],mv[m2],template[ell*16+mv[m1]*4+mv[m2]]);
        Rq[IDX_R(egl,ell,m1,m2)] = v*template[ell*16+mv[m1]*4+mv[m2]] * A[m2*egl->nfreq+m1];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
}

CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(ffp8_gcib,ffp8_gcib_init);

void grcT_compute(parametric *egl, double *rq, error **err);

parametric *grcT_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int i,m,*mv,m1,m2;
  pfchar name;
  double *conv,*A;
  int nfreq_template,lmax_template;
  double *dreq;

  testErrorRet(template==NULL,-11111,"template must be defined",*err,__LINE__,NULL);
  
  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);

  egl->eg_compute = &grcT_compute;
  egl->eg_free = &parametric_simple_payload_free;

  nfreq_template = parametric_get_value(egl,"grcT_nfreq_template",err);
  forwardError(*err,__LINE__,NULL);

  lmax_template = parametric_get_value(egl,"grcT_lmax_template",err);
  forwardError(*err,__LINE__,NULL);


  egl->payload = malloc_err(sizeof(double)* ((lmax_template+1)*(nfreq_template)*(nfreq_template)) + sizeof(double)*((nfreq_template)+egl->nfreq*egl->nfreq)+sizeof(double)*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);


  memcpy(egl->payload,template,sizeof(double)* ((lmax_template+1)*(nfreq_template)*(nfreq_template)));
  
  mv = egl->payload + sizeof(double)* ((lmax_template+1)*(nfreq_template)*(nfreq_template)) + sizeof(double)*(nfreq_template+egl->nfreq*egl->nfreq);
  conv = egl->payload + sizeof(double)* ((lmax_template+1)*(nfreq_template)*(nfreq_template));
  
  dreq = malloc_err(sizeof(double)*nfreq_template,err);
  forwardError(*err,__LINE__,NULL);

  for(i=0;i<nfreq_template;i++) {
    sprintf(name,"grcT_freq_%d",i);
    dreq[i] = parametric_get_value(egl,name,err);
    forwardError(*err,__LINE__,NULL);
    sprintf(name,"grcT_conv_%d",i);
    parametric_set_default(egl,name,1,err); 
    forwardError(*err,__LINE__,NULL);
    conv[i] = parametric_get_value(egl,name,err);
    forwardError(*err,__LINE__,NULL);
  }  

  fill_offset_freq(nfreq_template,dreq, egl,mv,-1,err);
  forwardError(*err,__LINE__,NULL);

  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"A_grcT_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"A_grcT_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  

  parametric_set_default(egl,"grcT_index",-1.3,err); 
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"grcT_index_0",-1.3,err); 
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"grcT_rigid",0,err); 
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"grcT_l_pivot",3000,err); 
  forwardError(*err,__LINE__,NULL);

  free(dreq);
  return egl;
}

void grcT_compute(parametric *egl, double *Rq, error **err) {
  double *template;
  double l_pivot,index,v,index0;
  int m1,m2,ell,i;
  double nrm;
  int *mv;
  int rigid,irigid;
  double *conv,*A;
  pfchar name;
  int nfreq_template,lmax_template;
  double *dreq;
  
  nfreq_template = parametric_get_value(egl,"grcT_nfreq_template",err);
  forwardError(*err,__LINE__,);

  lmax_template = parametric_get_value(egl,"grcT_lmax_template",err);
  forwardError(*err,__LINE__,);
  
  mv = egl->payload + sizeof(double)* ((lmax_template+1)*(nfreq_template)*(nfreq_template)) + sizeof(double)*(nfreq_template+egl->nfreq*egl->nfreq);
  conv = egl->payload + sizeof(double)* ((lmax_template+1)*(nfreq_template)*(nfreq_template));
  A = egl->payload + ((lmax_template+1)*(nfreq_template)*(nfreq_template))*sizeof(double) + sizeof(double)*(nfreq_template);

  template = egl->payload;
  l_pivot = parametric_get_value(egl,"grcT_l_pivot",err);
  forwardError(*err,__LINE__,);;
  
  rigid = parametric_get_value(egl,"grcT_rigid",err);
  forwardError(*err,__LINE__,);
  //_DEBUGHERE_("rigid %d",rigid);
  
  index = parametric_get_value(egl,"grcT_index",err);
  forwardError(*err,__LINE__,);
  index0 = parametric_get_value(egl,"grcT_index_0",err);
  forwardError(*err,__LINE__,);
  
  if (rigid==0) {
    sprintf(name,"A_grcT_%d",217);    
    irigid = 2;
  } else {
    dreq = malloc_err(sizeof(double)*nfreq_template,err);
    forwardError(*err,__LINE__,);

    for(i=0;i<nfreq_template;i++) {
      sprintf(name,"grcT_freq_%d",i);
      dreq[i] = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
    }
    
    sprintf(name,"A_grcT_%d",(int) rigid);  
    irigid=-1;
    for(m1=0;m1<nfreq_template;m1++) {
      if (fabs(rigid-dreq[m1])<1e-6) {
        irigid=m1;
      }
    }
    testErrorRet(irigid==-1,-55214,"AAAAAAA",*err,__LINE__,);
    rigid = 1;
    
  }

  
  nrm = parametric_get_value(egl,name,err);
  forwardError(*err,__LINE__,);
  //_DEBUGHERE_("%s %g",name,nrm);
     
  
  for(m1=0;m1<egl->nfreq;m1++) {
    
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"A_grcT_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"A_grcT_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*egl->nfreq+m2] = (v/template[((int) l_pivot)*(nfreq_template)*(nfreq_template)+mv[m1]*(nfreq_template)+mv[m2]]*(1-rigid) +  nrm/template[((int)l_pivot)*(nfreq_template)*(nfreq_template)+irigid*(nfreq_template)+irigid]*rigid*conv[mv[m1]]*conv[mv[m2]]/conv[irigid]/conv[irigid]) /l_pivot/(l_pivot+1)*2*M_PI ;
      A[m2*egl->nfreq+m1] = A[m1*egl->nfreq+m2];
    }
  }
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index-index0);
    //_DEBUGHERE_("%d %g",ell,v);
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        //_DEBUGHERE_("%d %d %d %d %g %g %g %d",ell,mv[m1],mv[m2],nfreq_template,template[ell*(nfreq_template)*(nfreq_template)+mv[m1]*(nfreq_template)+mv[m2]],v,A[m2*egl->nfreq+m1],ell*(nfreq_template)*(nfreq_template)+mv[m1]*(nfreq_template)+mv[m2]);
        Rq[IDX_R(egl,ell,m1,m2)] = v*template[ell*(nfreq_template)*(nfreq_template)+mv[m1]*(nfreq_template)+mv[m2]] * A[m2*egl->nfreq+m1];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
}

CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(grcT,grcT_init);

