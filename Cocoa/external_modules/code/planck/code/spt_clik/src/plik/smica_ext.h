/*
 *  smica.h
 *  lowly_project
 *
 *  Created by Karim Benabed on 30/10/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */


// some includes
#include "smica.h"

#ifndef _SMICA_EXT_
#define _SMICA_EXT_


double smica_crit_classic(void *vsmic,error **err);
double smica_crit_mgauss(void *vsmic, error **err);
double smica_crit_eig(void *vsmic, error **err);
double smica_crit_quad(void *vsmic,error **err);
double smica_crit_quadfid(void *vsmic,error **err);
double smica_crit_quad_mask(void *vsmic,error **err);

void smica_set_crit_mgauss(Smica *smic, int select_size, double *sigma, int *ordering,int* qmin, int* qmax, error **err);
void smica_set_crit_eig(Smica *smic, double *nrm, error **err);
void smica_set_crit_quad(Smica *smic,double *fid, int *mask, error **err);

void quad_mask_me(Smica *smic,double *delta);

double kld(int n, double* rq_hat, double* rq, error **err);


SmicaComp* comp_diag_init(int q, int m, error **err);
void comp_diag_update(void* data,double* locpars, double* rq, error **err);


typedef struct {
  int ntpl;
  int *ngcal;
  double *gcaltpl;
  double *tpl;
  double *tpll;
  double *bins;
  int nell;
} SC_gcal;


SmicaComp* comp_gcal_log_init(int q,int m, int *ngcal, double* galtpl, int nell, double *bins,error **err);
SmicaComp* comp_gcal_lin_init(int q,int m, int *ngcal, double* galtpl, int nell, double *bins,error **err);
void comp_gcal_log_update(void* data,double* locpars, double* rq, error **err);
void comp_gcal_lin_update(void* data,double* locpars, double* rq, error **err);
void comp_gcal_free(void** data);

SmicaComp* comp_cst_init(int nq, int m, double *rq_0, error **err);
void comp_cst_update(void* data,double* locpars, double* rq, error **err);
void free_comp_cst(void** data);

void amp_diag_free(void** data);
void amp_diag_update(void* data,double* locpars, double* rq, error **err);
SmicaComp* amp_diag_init(int nq, int m, double* tmpl, error **err);



typedef struct {
  double *rtpl;
  double *rqbuf;
  int *im,*jm;
  int npar;

} SC_gcal2;

SmicaComp* comp_gcal2_init(int q,int m, int npar, int *im, int *jm, double* tpl, error **err);
void comp_gcal2_update(void* data,double* locpars, double* rq, error **err);
void comp_gcal2_free(void** data);

#endif
