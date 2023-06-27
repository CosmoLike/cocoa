#include "cibsz.h"

void parametric_sz_cib_payload_free(void **pp) {
  sz_cib_payload *p;
  p = *pp;
  if (p!=NULL) {
    if(p->template!=NULL) {
      free(p->template);
    }
    free(p->ind_freq);
    free(p);
  }
  *pp=NULL;
}

int CIB_FREQ[] = {100,143,217};
double CIB_FIX_COLOR[] = {1.122, 1.134, 1.33};
double SQ_SZ_FIX_COLOR[] = {0.9619, 0.95, 0.0};


void parametric_sz_cib_payload_init(parametric *egl, double *template, error **err) {
  int m1,m2;
  int l;
  double fac;
  int lnorm = 3000;
  sz_cib_payload *payload;
  int template_size;
  long fix_sz_color;

  testErrorRetVA(egl->lmax>lmax_sz_template,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,,egl->lmax,lmax_sz_template);
  testErrorRetVA(egl->lmax>lmax_corr_template,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,,egl->lmax,lmax_corr_template);

  egl->payload = malloc_err(sizeof(sz_cib_payload), err);
  forwardError(*err, __LINE__,);

  payload = egl->payload;
  
  payload ->ind_freq = malloc_err(sizeof(int)*egl->nfreq,err);
  forwardError(*err,__LINE__,);
  
  for (m1=0;m1<egl->nfreq;m1++) {
    payload->ind_freq[m1]=-1;
    for (m2=0;m2<N_FREQ_CIB;m2++) {
      if (fabs(egl->freqlist[m1]-CIB_FREQ[m2])<1e-6) {
         payload->ind_freq[m1]=m2;
      }
    }
  }

  template_size = lmax_sz_template - lmin_sz_template + 1 + lmax_corr_template - lmin_corr_template + 1;
  payload->template = malloc_err(sizeof(double)*(template_size+egl->nfreq),err);
  forwardError(*err,__LINE__,);
  if (template!=NULL) {
    memcpy(payload->template,template,sizeof(double)*template_size);
    payload->sz_template = payload->template;
    payload->corr_template = payload->template + lmax_sz_template - lmin_sz_template + 1;
    fac = payload->sz_template[lnorm-lmin_sz_template];
    for (l=lmin_sz_template;l<=lmax_sz_template;l++) {
      payload->sz_template[l-lmin_sz_template] /= fac;
    }
    fac = payload->corr_template[lnorm-lmin_corr_template];
    for (l=lmin_corr_template;l<=lmax_corr_template;l++) {
      payload->corr_template[l-lmin_corr_template] /= fac;
    }
  }
  payload->fnu = payload->template + template_size;
  fix_sz_color = 1;
  fix_sz_color = pflist_get_int_value(egl->pf,"use_sz_fix_color",&fix_sz_color,err);             
  forwardError(*err,__LINE__,);     
  for (m1=0;m1<egl->nfreq;m1++) {                                        
    payload->fnu[m1] = sz_spectrum((double)egl->freqlist[m1],PRM_NU0);  
    if (fix_sz_color !=0) {
      payload->fnu[m1] *= sqrt(SQ_SZ_FIX_COLOR[payload->ind_freq[m1]]);
    }
  }
  for(m1=0;m1<N_FREQ_CIB;m1++) {
    payload->r_cib[m1*N_FREQ_CIB+m1] = 1;
    for(m2=m1+1;m2<N_FREQ_CIB;m2++) {
      payload->r_cib[m1*N_FREQ_CIB+m2] = 0;
      payload->r_cib[m2*N_FREQ_CIB+m1] = 0;
    }
  }                                                                 
}

void sz_cib_fill_cib_struct(parametric *egl,error **err) {
  int i,j;
  char tmpl[1000];
  long fix_color;
  sz_cib_payload *payload;
  
  payload = egl->payload;

  fix_color = 1;
  fix_color = pflist_get_int_value(egl->pf,"use_cib_fix_color",&fix_color,err);             
  forwardError(*err,__LINE__,);                                       
  
  
  for(i=0;i<N_FREQ_CIB;i++) {
    sprintf(tmpl,"A_cib_%d",CIB_FREQ[i]);
    payload->a_cib[i] = parametric_get_value(egl,tmpl,err);
    forwardError(*err,__LINE__,);
    if (fix_color !=0) {
      payload->a_cib[i] *= CIB_FIX_COLOR[i];
    }
    for(j=i+1;j<N_FREQ_CIB;j++) {
      sprintf(tmpl,"r_cib_%d_%d",CIB_FREQ[i],CIB_FREQ[j]);
      payload->r_cib[i*N_FREQ_CIB+j] = parametric_get_value(egl,tmpl,err);
      forwardError(*err,__LINE__,);
      payload->r_cib[j*N_FREQ_CIB+i] =payload->r_cib[i*N_FREQ_CIB+j];
    }
  }
}


void parametric_zero_rq(parametric *egl, double* rrq) {
  int ell,m1,m2;
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {  
    for(m1=0;m1<egl->nfreq;m1++) {                
      for(m2=m1;m2<egl->nfreq;m2++) {             
        rrq[IDX_R(egl,ell,m1,m2)] = 0;       
        rrq[IDX_R(egl,ell,m2,m1)] = 0;       
      }                                      
    }                                        
  }                                          
}                                             

void cib_A_derivative(parametric* egl, int iv,double *Rq, double *dRq, error **err) {
  int ell,m1,m2,mell,nfreq,mv,nrm;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;


  l_pivot = parametric_get_value(egl,"cib_l_pivot",err);
  forwardError(*err,__LINE__,);
  nrm = l_pivot*(l_pivot+1)/2/M_PI;

  index = parametric_get_value(egl,"cib_index",err);
  forwardError(*err,__LINE__,);

  //A = egl->payload;
  nfreq = egl->nfreq;
  //for(m1=0;m1<nfreq;m1++) {
  //  for(m2=m1;m2<nfreq;m2++) {
  //    sprintf(name,"cib_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
  //    v = 1;
  //    v = parametric_get_value(egl,name,err);
  //    A[m1*nfreq+m2] = v;
  //    A[m2*nfreq+m1] = v;
  //  }
  //}

  stop = 0;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"cib_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      if (strcmp(egl->varkey[iv],name)==0) {
        stop=1;
        memset(dRq,0,sizeof(double)*(egl->lmax+1-egl->lmin)*nfreq*nfreq);
        for(ell=egl->lmin;ell<=egl->lmax;ell++) {
          v = pow(ell/l_pivot,index);
          dRq[IDX_R(egl,ell,m1,m2)] = v/nrm;
          dRq[IDX_R(egl,ell,m2,m1)] = v/nrm;
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

void sz_cib_A_cib_XXX_derivative(parametric *egl, int II,int iv, double *Rq, double *dRq, error **err) {
  SZ_CIB_DEFS;

  SZ_CIB_INITS;
  
  parametric_zero_rq(egl,dRq);
  
  // NO SZ part
  // Add the CIB part
  if (do_cib == 1) {
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {
      for(m1=0;m1<nfreq;m1++) {
        ind1 = ind_freq[m1];
        for (m2=m1;m2<nfreq;m2++){
          ind2 = ind_freq[m2];
          if ((ind1 >= 0) && (ind2 >= 0)) { // 100, 143 or 217 GHz
            if ((ind1==II) || (ind2==II)) { // freq1 or freq2 = 100 GHz
              if (ind1==ind2) { // freq1 = freq2 = 100 GHz
                dRq[IDX_R(egl,ell,m1,m2)] += r_cib[ind1*N_FREQ_CIB+ind2]/d3000 * pow(((double)ell/(double)lnorm),cib_index);
              } else { // freq1 != freq2 = 143, 217 GHz
                if (ind1==II) {
                  dRq[IDX_R(egl,ell,m1,m2)] += r_cib[ind1*N_FREQ_CIB+ind2]*sqrt(a_cib[ind2]/a_cib[ind1])/2.0/d3000 * pow(((double)ell/(double)lnorm),cib_index);  
                } else {
                  dRq[IDX_R(egl,ell,m1,m2)] += r_cib[ind1*N_FREQ_CIB+ind2]*sqrt(a_cib[ind1]/a_cib[ind2])/2.0/d3000 * pow(((double)ell/(double)lnorm),cib_index);
                }
                dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
              }
            }
          }
        }
      }
    }
  }
  
  // Add the CIB-SZ correlation part
  if (do_szxcib == 1) {
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {
      dell = (double)ell;
      mell = ell - lmin_corr_template;
      for (m1=0;m1<nfreq;m1++){
        ind1 = ind_freq[m1];
        for (m2=m1;m2<nfreq;m2++) {
          ind2 = ind_freq[m2];
          if ((ind1 >= 0) && (ind2 >= 0)) { // 100, 143 or 217 GHz
            if ((ind1==II) || (ind2==II)) { // freq1 or freq2 = 100 GHz
              if (ind1==ind2) { // freq1 = freq2 = 100 GHz
                dRq[IDX_R(egl,ell,m1,m2)] -= xi_sz_cib * sqrt(a_sz) * sqrt(fnu[m1]/a_cib[ind1]) * corr_template[mell] * 2.0*M_PI/(dell*(dell+1.0));
              } else { 
                if (ind1==II) {
                  dRq[IDX_R(egl,ell,m1,m2)] -= xi_sz_cib * sqrt(a_sz) * sqrt(fnu[m2]/a_cib[ind1])/2.0 *  corr_template[mell] * 2.0*M_PI/(dell*(dell+1.0));
                } else {
                  dRq[IDX_R(egl,ell,m1,m2)] -= xi_sz_cib * sqrt(a_sz) * sqrt(fnu[m1]/a_cib[ind2])/2.0 *  corr_template[mell] * 2.0*M_PI/(dell*(dell+1.0));
                }
                dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
              }
            }
          }
        }
      }
    }
  }
  return;
}

void sz_cib_A_cib_100_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  sz_cib_A_cib_XXX_derivative(egl,0,iv,Rq,dRq,err);
  forwardError(*err,__LINE__,);
  return;
}

void sz_cib_A_cib_143_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  sz_cib_A_cib_XXX_derivative(egl,1,iv,Rq,dRq,err);
  forwardError(*err,__LINE__,);
  return;
}

void sz_cib_A_cib_217_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  sz_cib_A_cib_XXX_derivative(egl,2,iv,Rq,dRq,err);
  forwardError(*err,__LINE__,);
  return;
}

void sz_cib_r_cib_XXX_YYY_derivative(parametric *egl, int II, int JJ, int iv, double *Rq, double *dRq, error **err) {
  SZ_CIB_DEFS;

  SZ_CIB_INITS;
  
  parametric_zero_rq(egl,dRq);

  /* NO SZ part*/                                                                                                          
  /* Add the CIB part */                                                                                                   
  if (do_cib == 1)  {                                                                                                      
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {                                                                             
      for(m1=0;m1<nfreq;m1++) {                                                                                            
        ind1 = ind_freq[m1];
        if (ind1!=II && ind1!=JJ) {
          continue;
        }                                                                                               
        for (m2=m1;m2<nfreq;m2++){                                                                                         
          ind2 = ind_freq[m2];                                                                                             
          if (((ind1 == II) && (ind2 == JJ)) || ((ind1 == JJ) && (ind2 == II))) {                                                                              
            dRq[IDX_R(egl,ell,m1,m2)] += sqrt(a_cib[ind1]*a_cib[ind2])/d3000 * pow(((double)ell/(double)lnorm),cib_index); 
            dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];                                                         
          }                                                                                                                
        }                                                                                                                  
      }                                                                                                                    
    }                                                                                                                      
  }                                                                                                                        
  return; 
}

void sz_cib_r_cib_100_143_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  sz_cib_r_cib_XXX_YYY_derivative(egl,0,1,iv,Rq,dRq,err);
  forwardError(*err,__LINE__,);
  return;
}

void sz_cib_r_cib_100_217_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  sz_cib_r_cib_XXX_YYY_derivative(egl,0,2,iv,Rq,dRq,err);
  forwardError(*err,__LINE__,);
  return;
}

void sz_cib_r_cib_143_217_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
sz_cib_r_cib_XXX_YYY_derivative(egl,1,2,iv,Rq,dRq,err);
  forwardError(*err,__LINE__,);
  return;
}

void sz_cib_xi_sz_cib_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  SZ_CIB_DEFS;

  SZ_CIB_INITS;  

  parametric_zero_rq(egl,dRq);

  // NO SZ part, no CIB part

  // Add the CIB-SZ correlation part
  if (do_szxcib ==1) {
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {
      dell = (double)ell;
      mell = ell - lmin_corr_template;
      for (m1=0;m1<nfreq;m1++){
        ind1 = ind_freq[m1];
        for (m2=m1;m2<nfreq;m2++) {
          ind2 = ind_freq[m2];
          if ((ind1 >= 0) && (ind2 >= 0)) { // 100, 143 or 217 GHz
            dRq[IDX_R(egl,ell,m1,m2)] -= sqrt(a_sz) * ( sqrt(fnu[m1]*a_cib[ind2]) + sqrt(fnu[m2]*a_cib[ind1]) ) *
              corr_template[mell] * 2.0*M_PI/(dell*(dell+1.0));
            dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
          }
        }
      }
    }
  }
  return;
}

void sz_cib_cib_index_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  SZ_CIB_DEFS;

  SZ_CIB_INITS;

  parametric_zero_rq(egl,dRq);

  // NO SZ part

  // Add the CIB part
  if (do_cib == 1) {
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {
      for(m1=0;m1<nfreq;m1++) {
        ind1 = ind_freq[m1];
        for (m2=m1;m2<nfreq;m2++){
          ind2 = ind_freq[m2];
          if ((ind1 >= 0) && (ind2 >= 0)) { // 100, 143 or 217 GHz
            dRq[IDX_R(egl,ell,m1,m2)] += r_cib[ind1*N_FREQ_CIB+ind2]*sqrt(a_cib[ind1]*a_cib[ind2])/d3000 * log((double)ell/(double)lnorm)*pow(((double)ell/(double)lnorm),cib_index);
            dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
          }
        }
      }
    }
  }

  // No CIB-SZ correlation part

  return;
}

void sz_cib_A_sz_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  SZ_CIB_DEFS;

  SZ_CIB_INITS;
    
  parametric_zero_rq(egl,dRq);
  
  // Compute the SZ part first
  if (do_sz ==1) {
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {
      dell = (double)ell;
      mell = (ell-lmin_sz_template);
      for (m1=0;m1<nfreq;m1++) {
        for (m2=m1;m2<nfreq;m2++) {
          dRq[IDX_R(egl,ell,m1,m2)] += 2.0*M_PI/(dell*(dell+1.0)) * sz_template[mell] * fnu[m1] * fnu[m2];
          dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
        }
      }
    }
  }

  // NO CIB part
  // Add the CIB-SZ correlation part
  if (do_szxcib ==1) {
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {
      dell = (double)ell;
      mell = ell - lmin_corr_template;
      for (m1=0;m1<nfreq;m1++){
        ind1 = ind_freq[m1];
        for (m2=m1;m2<nfreq;m2++) {
          ind2 = ind_freq[m2];
          if ((ind1 >= 0) && (ind2 >= 0)) { // 100, 143 or 217 GHz for now
            dRq[IDX_R(egl,ell,m1,m2)] -= xi_sz_cib * 1./(2.0*sqrt(a_sz)) * ( sqrt(fnu[m1]*a_cib[ind2]) + sqrt(fnu[m2]*a_cib[ind1]) ) *
               corr_template[mell] * 2.0*M_PI/(dell*(dell+1.0));
            dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
          }
        }
      }
    }
  }
  return;
}  


void fill_offset_freq(int idreq,double *dreq, parametric *egl,int *mv,int def, error **err) {
  int m,i;

  for(m=0;m<egl->nfreq;m++) {
    double f;
    f = egl->freqlist[m];
    mv[m]=def;
    for(i=0;i<idreq;i++) {
      //_DEBUGHERE_("%g %d",f,dreq[i]);
      if (fabs(f-dreq[i])<1e-6) {
        mv[m]=i;
        break;
      }  
    }
    testErrorRetVA(mv[m]==-1,-431432,"Don't know how to compute component for freq %g",*err,__LINE__,,f);
  }    
}
