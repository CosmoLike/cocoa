/*
 *  smica.c
 *  lowly_project
 *
 *  Created by Karim Benabed on 30/10/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#include "smica_ext.h"


// modified gaussian criterion
typedef struct {
  int vec_size;
  int *vec_select;
  double *vec_buf,*vec;
  double *sigma_inverse;
} gauss_lkl_data;

double smica_crit_mgauss(void *vsmic, error **err) {
  int iq,im1,im2,iv,m2,m;
  double res;
  Smica *smic;
  char uplo;
  double done,dzero;
  int one,i,j;
  gauss_lkl_data *pld;

  smic = vsmic;
  pld = smic->lkl_data;
  //_DEBUGHERE_("%d %d",smic->nq,smic->m);
  //write_bin_vector(smic->rq, "rq.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);   
  //forwardError(*err,__LINE__,0);
  //write_bin_vector(smic->rq_hat, "rq_hat.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);   
  //forwardError(*err,__LINE__,0);
  // reorganize data
  //write_bin_vector(smic->quad_mask, "quad_mask.dat", sizeof(int)*(smic->quad_sn), err);   
  m = smic->m;
  m2 = m*m;
  iv = 0;
  for (iv=0;iv<pld->vec_size;iv++) {
    //int civ;
    pld->vec[iv] = smic->rq[pld->vec_select[iv]] - smic->rq_hat[pld->vec_select[iv]];
    //civ = pld->vec_select[iv];
    //_DEBUGHERE_("%d (%d %d %d) %d %g %g %g",iv,civ/m2,(civ-(civ/m2)*m2)/m,(civ-(civ/m2)*m2)%m, pld->vec_select[iv],smic->rq[pld->vec_select[iv]], smic->rq_hat[pld->vec_select[iv]],pld->vec[iv]);
  }
  //_DEBUGHERE_("%d",iv);
  
  one = 1;
  done = 1;
  dzero = 0;
  uplo = 'L';
  //printMat(smic->crit_cor,pld->vec_size,pld->vec_size);
  //_DEBUGHER_("%d",pld->vec_size);
  //write_bin_vector(pld->vec, "gvec.dat", sizeof(double)*(pld->vec_size), err);   
  //write_bin_vector(smic->crit_cor, "crit_cor.dat", sizeof(double)*(pld->vec_size)*(pld->vec_size), err);   
  //_DEBUGHERE_("","");
  //for(i=0;i<pld->vec_size;i++) {
  //  pld->vec_buf[i] = smic->crit_cor[i*pld->vec_size]*pld->vec[0];
  //  for(j=1;j<pld->vec_size;j++) {
  //    pld->vec_buf[i] += smic->crit_cor[i*pld->vec_size+j]*pld->vec[j];
  //  }
  //}
  dsymv(&uplo, &pld->vec_size, &done, pld->sigma_inverse, &pld->vec_size, pld->vec, &one, &dzero, pld->vec_buf, &one);

  //write_bin_vector(pld->vec, "gvecCC.dat", sizeof(double)*(pld->vec_size), err);   
  
  res = 0;
  for(iq=0;iq<pld->vec_size;iq++) {
    //_DEBUGHERE_("%g %g %g",res,pld->vec[iq],pld->vec_buf[iq])
    res += pld->vec[iq]*pld->vec_buf[iq];
  }
  return -.5*res;  
}

void free_smica_crit_mgauss(void **ppld) {
  gauss_lkl_data *pld;
  pld = *ppld;
  free(pld->vec);
  free(pld->sigma_inverse);
  free(pld);
  
  *ppld = NULL;
}

void smica_set_crit_mgauss(Smica *smic, int select_size, double *sigma, int *ordering,int* qmin, int* qmax, error **err) {
  int iv,jv,iq,ic;
  int nv,info;
  gauss_lkl_data *pld;
  char uplo;

  pld = malloc_err(sizeof(gauss_lkl_data),err);
  forwardError(*err,__LINE__,);
  
  pld->vec_size = select_size;

  pld->vec_select = malloc_err(sizeof(int)*select_size,err);
  forwardError(*err,__LINE__,);
  
  nv = 0;
  for (ic=0;ic<(smic->m*(smic->m+1))/2;ic++) {
    iv = ordering[ic*2];
    jv = ordering[ic*2+1];
    for (iq=qmin[iv*smic->m+jv];iq<qmax[iv*smic->m+jv];iq++) {
      pld->vec_select[nv] = iq*smic->m*smic->m+iv*smic->m+jv;
      nv++;
      testErrorRet(nv>pld->vec_size,-2422,"number of line in covariance matrix smaller than number of selected Cl !",*err,__LINE__,);
    }
  }

  pld->vec = malloc_err(sizeof(double)*pld->vec_size*2,err);
  forwardError(*err,__LINE__,);
  pld->vec_buf = pld->vec + pld->vec_size;
  
  pld->sigma_inverse = malloc_err(sizeof(double)*pld->vec_size*pld->vec_size,err);
  forwardError(*err,__LINE__,);
  memcpy(pld->sigma_inverse,sigma,sizeof(double)*pld->vec_size*pld->vec_size);
  uplo = 'L';
  dpotri(&uplo,&select_size,pld->sigma_inverse,&select_size,&info);
  testErrorRetVA(info!=0,-432432,"cannot inverse covariance with dpotri (info = %d)",*err,__LINE__,,info);

  smic->crit = &smica_crit_mgauss;
  smic->lkl_data = pld;
  smic->lkl_data_free = free_smica_crit_mgauss;
}

//void smica_set_crit_mgauss(Smica *smic, double *crit_cor, int *mask,int *ordering,error **err) {
//  int iv,jv,iq;
//  int nv;
//  gauss_lkl_data *pld;
//
//  pld = malloc_err(sizeof(gauss_lkl_data),err);
//  forwardError(*err,__LINE__,);
//
//  _DEBUGHERE_("","");
//  nv = (smic->nq * smic->m * (smic->m+1))/2;
//  pld->vec_select = malloc_err(sizeof(int)*smic->nq*smic->m*smic->m,err);
//  forwardError(*err,__LINE__,);
//
//  // mask is both in spectra an ell !!!
//  nv = 0;
//  if (ordering == NULL) {
//    for (iv=0;iv<smic->m;iv++) {
//      for (jv=iv;jv<smic->m;jv++) {
//        for (iq=0;iq<smic->nq;iq++) {
//          //_DEBUGHERE_("%d %d %d %d (%d)",iq,iv,jv,mask[iq*smic->m*smic->m+iv*smic->m+jv],nv);
//          if (mask == NULL || mask[iq*smic->m*smic->m+iv*smic->m+jv]!=0) {
//            pld->vec_select[nv] = iq*smic->m*smic->m+iv*smic->m+jv;
//            nv++;
//          }
//        }
//      }
//    }
//  } else {
//    int ic;
//    for (ic=0;ic<(smic->m*(smic->m+1))/2;ic++) {
//      iv = ordering[ic*2];
//      jv = ordering[ic*2+1];
//      for (iq=0;iq<smic->nq;iq++) {
//        if (mask == NULL || mask[iq*smic->m*smic->m+iv*smic->m+jv]!=0) {
//          //_DEBUGHERE_("%d %d %d %d (%d)",iq,iv,jv,mask[iq*smic->m*smic->m+iv*smic->m+jv],nv);
//          pld->vec_select[nv] = iq*smic->m*smic->m+iv*smic->m+jv;
//          nv++;
//        }
//      }
//    }
//  }
//
//  pld->vec_size = nv;
//
//  pld->vec = malloc_err(sizeof(double)*pld->vec_size*2,err);
//  forwardError(*err,__LINE__,);
//  pld->vec_buf = pld->vec + pld->vec_size;
//
//  pld->sigma_inverse = malloc_err(sizeof(double)*pld->vec_size*pld->vec_size,err);
//  forwardError(*err,__LINE__,);
//  
//  memcpy(pld->sigma_inverse,crit_cor,sizeof(double)*pld->vec_size*pld->vec_size);
//
//  smic->crit = &smica_crit_mgauss;
//  smic->lkl_data = pld;
//  smic->lkl_data_free = free_smica_crit_mgauss;
//}

// classic criterion


void smica_set_crit_eig(Smica *smic, double *nrm, error **err) {
  char uplo,jobz;
  int info, m;
  double rq, w,wrk;


  uplo = 'L';
  jobz = 'V';
  info = 0;
  smic->eig_lwork = -1;
  m = smic->m;
  
  dsyev(&jobz, &uplo, &m, &rq, &m, &w, &wrk, &smic->eig_lwork, &info);
  smic->eig_lwork = wrk;

  smic->eig_buf = malloc_err(sizeof(double)*(m+smic->eig_lwork),err);
  forwardError(*err,__LINE__,);

  smic->eig_nrm = malloc_err(sizeof(double)*smic->nq,err);
  forwardError(*err,__LINE__,);

  if (nrm!=NULL) {
    memcpy(smic->eig_nrm,nrm,sizeof(double)*smic->nq);
  } else {
    memset(smic->eig_nrm,0,sizeof(double)*smic->nq);
  }
  
  smic->crit = smica_crit_eig;
  //_DEBUGHERE_("","");
  
}


double smica_crit_eig(void* vsmic, error **err) {
  char jobz,uplo;
  double *wrk, *w;
  int info;
  double res,pes;
  int i,j,k,m,lwork,iq;
  Smica *smic;
  double jes,les,dg;
  double *rq,*rq_hat;
  smic = vsmic;

  m = smic->m;
  // assumes rq_hat is rq_hat !
  // nrm is -.5 (log(det(rq_hat)) + ndim) or any arbitrary number
  // buf is big enough to hold wrk and w
  // compute the eigenvalues and modes of rq

  uplo = 'L';
  jobz = 'V';
  info = 0;
  w = smic->eig_buf;
  wrk = w + m;
  lwork = smic->eig_lwork;
  res = 0;
  
  for (iq=0;iq<smic->nq;iq++) {
    rq = &(smic->rq[iq*m*m]);
    rq_hat = &(smic->rq_hat[iq*m*m]);

    dsyev(&jobz, &uplo, &m, rq, &m, w, wrk, &lwork, &info);
    testErrorRetVA(info !=0,-1234321,"dsyev failed (%d)",*err,__LINE__,0,info);
  
    les = 0;

    // ici peut etre un peu de regul du resultat, un jour
    for(i=0;i<m;i++) {
      les += log(w[i]); 
      pes = 0;
      dg=0;
      for(j=0;j<m;j++) {
        jes = 0;
        for(k=0;k<m;k++) {
          jes += rq_hat[j*m+k] * rq[k+i*m];
        }
        dg += rq[i+j*m] * rq[i+j*m];
        pes += jes * rq[j+i*m];
      }
      les += pes/w[i];
    }

    //_DEBUGHERE_("%g",les*.5 - smic->eig_nrm[iq])
    
    res += (les*.5 - smic->eig_nrm[iq]) * smic->wq[iq];
}
  
  return -res;
}


// quad criterion

void smica_set_crit_quad(Smica *smic, double *fid,int *mask,error **err) {
  int iq,m,nq;
  char uplo,trans,diag,side;
  int info;
  double *fdq;
  int sz;
  int i,j,sn;

  //_DEBUGHERE_("","");
  smic->crit = smica_crit_quad;
  smic->quad_sn = smic->m;
  smic->quad_mask = NULL;
  if (mask!=NULL) {
    smic->quad_mask = malloc_err(sizeof(int)*smic->m*smic->m,err);
    forwardError(*err,_LINE__,);
    memcpy(smic->quad_mask,mask,sizeof(int)*smic->m*smic->m);
    sn = 0;
    
    for(i=0;i<smic->m;i++) {
      for(j=i;j<smic->m;j++) {
        sn += mask[i*smic->m+j];
      }
    }
    smic->quad_sn = sn;
    smic->eig_buf = malloc_err(sizeof(double)*sn*sn*smic->nq,err);
    forwardError(*err,_LINE__,);
    memcpy(smic->eig_buf,fid,sizeof(double)*sn*sn*smic->nq);
    smic->crit = smica_crit_quad_mask;
    
    return;
  }
  if (fid!=NULL) {
    //_DEBUGHERE_("","");
    smic->crit = smica_crit_quadfid;
    smic->eig_buf = malloc_err(sizeof(double)*smic->m*smic->m*smic->nq,err);
    forwardError(*err,_LINE__,);
    memcpy(smic->eig_buf,fid,sizeof(double)*smic->m*smic->m*smic->nq);
    uplo = 'L';
    m = smic->m;
    sz = m*m;

    for(iq = 0;iq<smic->nq;iq++) {
      fdq = smic->eig_buf + iq*sz;
      dpotrf(&uplo,&m,fdq,&m,&info);
      testErrorRetVA(info!=0,lowly_chol,"Could not cholesky decompose rq_hat using dpotrf (%d)",*err,__LINE__,,info);   
    }
  }
}

double smica_crit_quad_mask(void *vsmic,error **err) {
  double res;
  int iq,m,nq;
  Smica *smic;
  int one,sz;
  double dminusone;
  double *c_rq_hat,*c_rq,*c_rq_fid;
  char uplo,trans,diag,side;
  int info;
  int i,j,p;
  double *vecq;
  double done,dzero;
  int sn;
  double kdd;

  smic = vsmic;
  m = smic->m;
  nq = smic->nq;

  one = 1;
  sz = m*m;
  dminusone = -1;
  uplo = 'L';
  done = 1;
  dzero = 1;
  sn = smic->quad_sn;

  res = 0;
      
  for(iq=0;iq<nq;iq++) {
    double kdd;
    c_rq_hat = smic->rq_hat + sz*iq;
    c_rq = smic->rq + sz*iq;
    c_rq_fid = smic->eig_buf + sn*sn*iq;
    vecq = smic->z_buf;

    //printMat(c_rq,m,m);
    // c_rq = c_rq - c_rq_hat  
    daxpy(&sz, &dminusone,c_rq_hat,&one, c_rq, &one);
    
    p = 0;
    for(i=0;i<m;i++) {
      for(j=i;j<m;j++) {
        
        if (smic->quad_mask[i*m+j]==0) {
          continue;
        }
        if (i==j) {
          vecq[p] = c_rq[i*m+j];
        } else {
          vecq[p] = c_rq[i*m+j];
        }
        
        p++;
      }
    }
    
    
    dsymv(&uplo, &sn, &done, c_rq_fid, &sn, vecq, &one, &dzero, c_rq, &one);
    
        
    kdd = 0;
    for(i=0;i<sn;i++) {
      kdd += c_rq[i]*vecq[i];
    }
    
    kdd=kdd/4.;
    

    res += smic->wq[iq] * kdd;
  }

  return -res;  
}

double smica_crit_quadfid(void *vsmic,error **err) {
  double res;
  int iq,m,nq;
  Smica *smic;
  int one,sz;
  double dminusone;
  double *c_rq_hat,*c_rq,*c_rq_fid;
  char uplo,trans,diag,side;
  int info;
  int i,j;

  
  smic = vsmic;
  m = smic->m;
  nq = smic->nq;

  one = 1;
  sz = m*m;
  dminusone = -1;
  uplo = 'L';
    
  res = 0;

  for(iq=0;iq<nq;iq++) {
    double kdd;
    c_rq_hat = smic->rq_hat + sz*iq;
    c_rq = smic->rq + sz*iq;
    c_rq_fid = smic->eig_buf + sz*iq;
    
    //printMat(c_rq,m,m);
    // c_rq = c_rq - c_rq_hat  
    daxpy(&sz, &dminusone,c_rq_hat,&one, c_rq, &one);
    

    dpotrs(&uplo, &m, &m, c_rq_fid, &m, c_rq, &m, &info);
    testErrorRetVA(info!=0,lowly_chol,"Could not solve rq using dpotrs (%d)",*err,__LINE__,0,info);

    kdd = 0;
    for(i=0;i<m;i++) {
      for(j=0;j<m;j++) {
        kdd += c_rq[i*m+j]*c_rq[j*m+i];
      }
    }
    kdd=kdd/4.;
    
    res += smic->wq[iq] * kdd;
  }

  return -res;  
}


double smica_crit_quad(void *vsmic,error **err) {
  double res;
  int iq,m,nq;
  Smica *smic;
  int one,sz;
  double dminusone;
  double *c_rq_hat,*c_rq;
  char uplo,trans,diag,side;
  int info;
  int i,j;

  smic = vsmic;
  m = smic->m;
  nq = smic->nq;

  one = 1;
  sz = m*m;
  dminusone = -1;
  uplo = 'L';
    
  res = 0;

  for(iq=0;iq<nq;iq++) {
    double kdd;
    c_rq_hat = smic->z_buf;
    c_rq = smic->rq + sz*iq;
    memcpy(c_rq_hat,smic->rq_hat+sz*iq,m*m*sizeof(double));
    
    // c_rq_hat = c_rq_hat - c_rq  
    daxpy(&sz, &dminusone,c_rq,&one, c_rq_hat, &one);

    dpotrf(&uplo,&m,c_rq,&m,&info);
    testErrorRetVA(info!=0,lowly_chol,"Could not cholesky decompose rq using dpotrf (%d)",*err,__LINE__,0,info);
  
    dpotrs(&uplo, &m, &m, c_rq, &m, c_rq_hat, &m, &info);
    testErrorRetVA(info!=0,lowly_chol,"Could not solve rq using dpotrs (%d)",*err,__LINE__,0,info);

    kdd = 0;
    for(i=0;i<m;i++) {
      for(j=0;j<m;j++) {
        kdd += c_rq_hat[i*m+j]*c_rq_hat[j*m+i];
      }
    }
    kdd=kdd/4.;

    res += smic->wq[iq] * kdd;
  }

  return -res;  
}




// Simple components



//diagonal

SmicaComp* comp_diag_init(int nq, int m, error **err) {
  SmicaComp* SC;
  
  SC = alloc_SC(nq*m,nq,m,NULL,&comp_diag_update,NULL,err);
  forwardError(*err,__LINE__,NULL);
  return SC;
}




void comp_diag_update(void* data,double* locpars, double* rq, error **err) {
  int iq,im,m2,ii,m;
  SmicaComp* SC;
  
  SC = data;
  m = SC->m;
  
  m2=m*m;
  ii=0;
  //_DEBUGHERE_("pars","");
  //printMat(locpars, SC->nq, m);    
  for(iq=0;iq<SC->nq;iq++) {
    //_DEBUGHERE_("rq a","");
    //printMat(rq+m2*iq, m, m);    
    for(im=0;im<m;im++) {
      rq[iq*m2+im*m+im] += locpars[ii];
      ii++;
    }
    //_DEBUGHERE_("rq b","");
    //printMat(rq+m2*iq, m, m);    
  }
}

SmicaComp* amp_diag_init(int nq, int m, double* tmpl, error **err) {
  SmicaComp *SC;
  double *itmpl;

  itmpl = malloc_err(sizeof(double)*nq*m,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(itmpl,tmpl,sizeof(double)*nq*m);

  SC = alloc_SC(m,nq,m,itmpl,&amp_diag_update,&amp_diag_free,err);
  forwardError(*err,__LINE__,NULL);
  
  return SC;
}

void amp_diag_update(void* data,double* locpars, double* rq, error **err) {
  int iq,im,m2,ii,m,q2;
  SmicaComp* SC;
  double *tmpl;

  SC = data;
  m = SC->m;
  tmpl = SC->data;
  m2 = m*m;
  q2 = SC->nq*SC->nq;

  ii=0;
  //_DEBUGHERE_("pars","");
  //printMat(locpars, SC->nq, m);
  for(im=0;im<SC->m;im++) {
    double pp;
    pp = locpars[im];
    for(iq=0;iq<SC->nq;iq++) {
      rq[iq*m2+im*m+im] += pp*tmpl[im*q2+iq];
    }
  }     
}

void amp_diag_free(void** data) {
  SmicaComp *SC;
  SC = *data;
  free(SC->data);
  free(SC);
  *data = NULL;
}



SC_gcal* comp_gcal_gen_init(int q, int m, int* ngcal, double *gcaltpl, int nell, double *bins, error **err) {
  SC_gcal *gc;
  int ntpl,im;
  
  gc = malloc_err(sizeof(SC_gcal),err);
  forwardError(*err,__LINE__,NULL);
  
  gc->ngcal = malloc_err(sizeof(int)*m,err);
  forwardError(*err,__LINE__,NULL);
  
  ntpl=0;
  for(im=0;im<m;im++) {
    ntpl+=ngcal[im];
    gc->ngcal[im]=ngcal[im];
  }
  gc->ntpl = ntpl;
  
  gc->gcaltpl = malloc_err(sizeof(double)*q*ntpl,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->gcaltpl,gcaltpl,sizeof(double)*q*ntpl);

  if (nell>0) {
    int iq,il;
    gc->gcaltpl = malloc_err(sizeof(double)*nell*ntpl,err);
    forwardError(*err,__LINE__,NULL);
    memcpy(gc->gcaltpl,gcaltpl,sizeof(double)*nell*ntpl);
    gc->nell = nell;
    gc->tpl = malloc_err(sizeof(double)*nell*m,err);
    forwardError(*err,__LINE__,NULL); 
    gc->bins = malloc_err(sizeof(double)*(nell*q),err);
    forwardError(*err,__LINE__,NULL);
    memcpy(gc->bins,bins,sizeof(double)*q*nell);
    gc->tpll = malloc_err(sizeof(double)*q,err);
    forwardError(*err,__LINE__,NULL);
    for (iq=0;iq<q;iq++) {
      gc->tpll[iq] = 0;
      for(il=0;il<nell;il++) {
        gc->tpll[iq] += bins[iq*nell+il]*bins[iq*nell+il];
      }
    }    
  } else {
    gc->nell =0;
    gc->gcaltpl = malloc_err(sizeof(double)*q*ntpl,err);
    forwardError(*err,__LINE__,NULL);
    memcpy(gc->gcaltpl,gcaltpl,sizeof(double)*q*ntpl);
    gc->tpl = malloc_err(sizeof(double)*q*m,err);
    forwardError(*err,__LINE__,NULL);
    gc->tpll=NULL;
    gc->bins = NULL;
  }
  return gc;
}

void comp_gcal_apply_update(int m,int nq,double* rq,double *tpl,int nell,double *tpll, double *bins) {
  int iq,im1,im2,il;

  if (nell>0) {
    for (iq=0;iq<nq;iq++) {
      for(im1=0;im1<m;im1++) {  
        for(im2=0;im2<m;im2++) {
          double cr;
          cr=0;
          for(il=0;il<nell;il++) {
            cr += bins[iq*nell+il] * bins[iq*nell+il] * tpl[im1*nq + il] * tpl[im2*nq + il];
          }
          rq[iq*m*m+im1*m+im2] *= cr/tpll[iq];
        }
      }
    }  
  } else {
    for (iq=0;iq<nq;iq++) {
      for(im1=0;im1<m;im1++) {  
        for(im2=0;im2<m;im2++) {
          rq[iq*m*m+im1*m+im2] *= tpl[im1*nq + iq] * tpl[im2*nq + iq];
        }
      }
    }
  }
}


SmicaComp* comp_gcal_lin_init(int q,int m, int *ngcal, double* gcaltpl, int nell, double *bins,error **err) {
  SC_gcal *gc;
  SmicaComp *SC;

  gc = comp_gcal_gen_init(q,m,ngcal,gcaltpl,nell,bins,err);
  forwardError(*err,__LINE__,NULL);

  SC = alloc_SC(gc->ntpl,q,m,gc, &comp_gcal_lin_update, &comp_gcal_free,err);
  forwardError(*err,__LINE__,NULL);

  return SC;
}

void comp_gcal_lin_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  SC_gcal *gc;
  int im1,im2, iq,ig;
  int m2,m,nq,ipars;

  SC = data;
  gc = SC->data;
  m2 = SC->m*SC->m;
  m = SC->m;
  nq = SC->nq;
  ipars = 0;

  if (gc->nell>0) {
    nq = gc->nell;
  }
  ipars = 0;
  // first build the template 
  for(im1=0;im1<m;im1++) {
    if (gc->ngcal[im1]==0) {
      for(iq=0;iq<nq;iq++) {
          gc->tpl[im1*nq + iq] = 1;
      }
      continue;
    }
    for(iq=0;iq<nq;iq++) {
      gc->tpl[im1*nq + iq] = locpars[ipars] * gc->gcaltpl[ipars*nq+iq];
    }
    ipars++;
    for(ig=1;ig<gc->ngcal[im1];ig++) {
      for(iq=0;iq<nq;iq++) {
        gc->tpl[im1*nq + iq] += locpars[ipars] * gc->gcaltpl[ipars*nq+iq];
       
      }
      ipars++;
    }
  }

  comp_gcal_apply_update(m,SC->nq,rq,gc->tpl,gc->nell,gc->tpll,gc->bins);
}

SmicaComp* comp_gcal_log_init(int q,int m, int* ngcal,  double* gcaltpl,int nell, double *bins,error **err) {
  SC_gcal *gc;
  SmicaComp *SC;

  gc = comp_gcal_gen_init(q,m,ngcal,gcaltpl,nell,bins,err);
  forwardError(*err,__LINE__,NULL);

  SC = alloc_SC(gc->ntpl,q,m,gc, &comp_gcal_log_update, &comp_gcal_free,err);
  forwardError(*err,__LINE__,NULL);

  return SC;
}

void comp_gcal_log_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  SC_gcal *gc;
  int im1,im2, iq,ig;
  int m2,m,nq,ipars;

  SC = data;
  gc = SC->data;
  m2 = SC->m*SC->m;
  m = SC->m;
  nq = SC->nq;
  ipars = 0;

  if (gc->nell>0) {
    nq = gc->nell;
  }
  ipars = 0;
  // first build the template 
  for(im1=0;im1<m;im1++) {
    if (gc->ngcal[im1]==0) {
      for(iq=0;iq<nq;iq++) {
          gc->tpl[im1*nq + iq] = 1;
      }
      continue;
    }
    for(iq=0;iq<nq;iq++) {
      gc->tpl[im1*nq + iq] = exp(locpars[ipars] * gc->gcaltpl[ipars*nq+iq]);
    }
    ipars++;
    for(ig=1;ig<gc->ngcal[im1];ig++) {
      for(iq=0;iq<nq;iq++) {
        gc->tpl[im1*nq + iq] *= exp(locpars[ipars] * gc->gcaltpl[ipars*nq+iq]);
       
      }
      ipars++;
    }
  }

  comp_gcal_apply_update(m,SC->nq,rq,gc->tpl,gc->nell,gc->tpll,gc->bins);
}

void comp_gcal_free(void** data) {
  SmicaComp *SC;
  SC_gcal *gc;
  
  SC = *data;
  gc = SC->data;
  free(gc->ngcal);
  free(gc->gcaltpl);
  free(gc->tpl);
  if (gc->nell>0) {
    free(gc->tpll);
    free(gc->bins);
  }
  free(gc);

  free(SC);
  *data = NULL;
}

SmicaComp* comp_gcal2_init(int q,int m, int npar, int *im, int *jm, double* tpl, error **err) {
 SC_gcal2 *gc;
  SmicaComp *SC;

  gc = malloc_err(sizeof(SC_gcal2),err);
  forwardError(*err,__LINE__,NULL);

  gc->npar = npar;

  gc->rtpl = malloc_err(sizeof(double)*npar*q,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->rtpl,tpl,sizeof(double)*npar*q);

  gc->rqbuf = malloc_err(sizeof(double)*m*m*q,err);
  forwardError(*err,__LINE__,NULL);
  
  gc->im = malloc_err(sizeof(int)*npar,err);
  forwardError(*err,__LINE__,NULL);
  gc->jm = malloc_err(sizeof(int)*npar,err);
  forwardError(*err,__LINE__,NULL);
  
  SC = alloc_SC(npar,q,m,gc, &comp_gcal2_update, &comp_gcal2_free,err);
  forwardError(*err,__LINE__,NULL);

  return SC;

}


void comp_gcal2_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  SC_gcal2 *gc;
  int ipar;
  int m2, q,m;
  int im,jm,iq,z;
  double *tpl;

  SC = data;
  gc = SC->data;

  m = SC->m;
  m2 = SC->m*SC->m;
  m = SC->m;
  q = SC->nq;
  
  memset(gc->rqbuf,0,sizeof(double)*m2*q);

  for(ipar=0;ipar<gc->npar;ipar++) {
    im = gc->im[ipar];
    jm = gc->jm[ipar];
    tpl = &(gc->rtpl[ipar*q]);
    for(iq=0;iq<q;iq++) {
      gc->rqbuf[iq*m2+im*m+jm] += 2*locpars[ipar] * tpl[iq];  
      gc->rqbuf[iq*m2+jm*m+im] = gc->rqbuf[iq*m2+im*m+jm];
    }    
  }
  for(z=0;z<m2*q;z++) {
    rq[z] *= exp(gc->rqbuf[z]);
  }
}

void comp_gcal2_free(void** data) {
  SmicaComp *SC;
  SC_gcal2 *gc;
  
  SC = *data;
  gc = SC->data;
  free(gc->im);
  free(gc->jm);
  free(gc->rtpl);
  free(gc->rqbuf);
  free(gc);

  free(SC);
  *data = NULL;
}

