#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(NAGLIB)
  #include <nag.h>
  #include <nag_stdlib.h>
  #include <nage01.h>
  #include <nage02.h>
#endif
#include "log.c/src/log.h"
#include "basics.h"
#include "cosmo3D.h"
#include "recompute.h"
#include "structs.h"

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Background
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void setup_chi(int* io_nz, double** io_z, double** io_chi, int io) {
  static int nz;
  static double* z = NULL;
  static double* chi = NULL;
  // IO == 1 IMPLES THAT IO_CHI(Z) WILL COPIED TO LOCAL H(Z)
  if (io == 1) {
    if (z != NULL) {
      #if defined(NAGLIB)
      NAG_FREE(z);
      #else
      free(z);
      #endif
    }
    if (chi != NULL) {
      #if defined(NAGLIB)
      NAG_FREE(chi);
      #else
      free(chi);
      #endif
    }
    nz = (*io_nz);
    if (!(nz > 5)) {
      log_fatal("array to small");
      exit(1);
    }
    #if defined(NAGLIB)
    z = (double*) NAG_ALLOC(nz, double);
    chi = (double *) NAG_ALLOC(nz, double);
    #else
    z = (double*) malloc(nz*sizeof(double));
    chi = (double *) malloc(nz*sizeof(double));
    #endif
    if (z == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    if (chi == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    for (int i=0; i<nz; i++) {
      z[i] = (*io_z)[i];
      chi[i] = (*io_chi)[i];
    }
  }
  else {
    // IO != 1 IMPLES THAT LOCAL H(Z) WILL BE COPIED TO IO_chi(Z)
    if (z == NULL || chi == NULL) {
      printf("cosmo3D::setup_hubble: z or H not allocated\n");
      exit(1);
    }
    (*io_nz) = nz;
    if (io_z != NULL) {
      if((*io_z) != NULL) {
        #if defined(NAGLIB)
        NAG_FREE((*io_z));
        #else
        free((*io_z));
        #endif
        (*io_z) = NULL;
      }
    }
    else {
      printf("cosmo3D::setup_hubble: input pointer not allocated\n");
      exit(1);
    }
    if (io_chi != NULL) {
      if((*io_chi) != NULL) {
        #if defined(NAGLIB)
        NAG_FREE((*io_chi));
        #else
        free((*io_chi));
        #endif
        (*io_chi) = NULL;
      }
    }
    else {
      log_fatal("input pointer not allocated\n");
      exit(1);
    }
    #if defined(NAGLIB)
    (*io_z) = (double*) NAG_ALLOC(nz, double);
    (*io_chi) = (double*) NAG_ALLOC(nz, double);
    #else
    (*io_z) = (double*) malloc(nz*sizeof(double));
    (*io_chi) = (double*) malloc(nz*sizeof(double));
    #endif
    if ((*io_z) == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    if ((*io_chi) == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    for (int i = 0; i < nz; i++) {
      (*io_z)[i] = z[i];
      (*io_chi)[i] = chi[i];
    }
  }
}

double chi(double io_a) {
  static cosmopara C;
  static int first = 0;
  double redshift, out_chi;
  #if defined(NAGLIB)
  static Nag_Spline spline;
  if (first == 0) {
    spline.lamda = NULL;
    spline.c = NULL;
    first = 1;
  }
  #else
  static int nz;
  static double** z;
  static double** chi;
  if (first == 0) {
    z = NULL;
    chi = NULL;
    first = 1;
  }
  #endif
  if (recompute_cosmo3D(C)) {
    update_cosmopara(&C);
    #if defined(NAGLIB)
    int nz;
    double** z = (double**) NAG_ALLOC(1, double*);
    double** chi = (double**) NAG_ALLOC(1, double*);
    #else
    if(z != NULL) {
      free(z);
      z = NULL;
    }
    z = (double**) malloc(1*sizeof(double*));
    if(chi != NULL) {
      free(chi);
      chi = NULL;
    }
    chi = (double**) malloc(1*sizeof(double*));
    #endif
    if(z == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*z) = NULL;
    if(chi == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*chi) = NULL;
    #if defined(NAGLIB)
    if(spline.lamda != NULL) {
      NAG_FREE(spline.lamda);
      spline.lamda = NULL;
    }
    if(spline.c != NULL) {
      NAG_FREE(spline.c);
      spline.c = NULL;
    }
    #endif
    setup_chi(&nz, z, chi, 0);
    #if defined(NAGLIB)
    NagError fail;
    INIT_FAIL(fail);
    nag_1d_spline_interpolant(nz, (*z), (*chi), &spline, &fail);
    if(fail.code != NE_NOERROR) {
      log_fatal("NAG interpolation failed %s", fail.message);
      exit(1);
    }
    NAG_FREE((*z));
    NAG_FREE(z);
    NAG_FREE((*chi));
    NAG_FREE(chi);
    #endif
  }
  redshift = 1.0/io_a - 1.0;
  #if defined(NAGLIB)
  NagError fail;
  INIT_FAIL(fail);
  nag_1d_spline_evaluate(redshift, &out_chi, &spline, &fail);
  if(fail.code != NE_NOERROR) {
    log_fatal("NAG interpolation failed %s", fail.message);
    exit(1);
  }
  #else
  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = nz-1;
    while (ihi>ilo+1) {
      size_t ll = (ihi+ilo)/2;
      if((*z)[ll]>redshift)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }
  double dy = (redshift-(*z)[j])/((*z)[j+1]-(*z)[j]);
  out_chi = (*chi)[j] + dy*((*chi)[j+1]-(*chi)[j]);
  #endif
  // convert from (Mpc/h) to (Mpc/h)/(c/H0=100)^3 (dimensioneless)
  return (out_chi/(cosmology.coverH0));
}

double dchi_da(double io_a) {
  static cosmopara C;
  static int first = 0;
  double redshift;
  double out_dchi[4];
  #if defined(NAGLIB)
  static Nag_Spline spline;
  if (first == 0) {
    spline.lamda = NULL;
    spline.c = NULL;
    first = 1;
  }
  #else
  static int nz;
  static double** z;
  static double** chi;
  if (first == 0) {
    z = NULL;
    chi = NULL;
    first = 1;
  }
  #endif
  if (recompute_cosmo3D(C)) {
    update_cosmopara(&C);
    #if defined(NAGLIB)
    int nz;
    double** z = (double**) NAG_ALLOC(1, double*);
    double** chi = (double**) NAG_ALLOC(1, double*);
    #else
    if(z != NULL) {
      free(z);
      z = NULL;
    }
    z = (double**) malloc(1*sizeof(double*));
    if(chi != NULL) {
      free(chi);
      chi = NULL;
    }
    chi = (double**) malloc(1*sizeof(double*));
    #endif
    if(z == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*z) = NULL;
    if(chi == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*chi) = NULL;
    #if defined(NAGLIB)
    if(spline.lamda != NULL) {
      NAG_FREE(spline.lamda);
      spline.lamda = NULL;
    }
    if(spline.c != NULL) {
      NAG_FREE(spline.c);
      spline.c = NULL;
    }
    #endif
    setup_chi(&nz, z, chi, 0);
    #if defined(NAGLIB)
    NagError fail;
    INIT_FAIL(fail);
    nag_1d_spline_interpolant(nz, (*z), (*chi), &spline, &fail);
    if(fail.code != NE_NOERROR) {
      log_fatal("NAG interpolation failed %s", fail.message);
      exit(1);
    }
    NAG_FREE((*z));
    NAG_FREE(z);
    NAG_FREE((*chi));
    NAG_FREE(chi);
    #endif
  }
  redshift = 1.0/io_a - 1.0;
  #if defined(NAGLIB)
  NagError fail;
  Nag_DerivType derivs;
  INIT_FAIL(fail);
  derivs = Nag_LeftDerivs;
  nag_1d_spline_deriv(derivs,redshift,out_dchi,&spline,&fail);
  if(fail.code != NE_NOERROR) {
    log_fatal("NAG interpolation failed %s", fail.message);
    exit(1);
  }
  #else
  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = nz-1;
    while (ihi>ilo+1) {
      size_t ll = (ihi+ilo)/2;
      if((*z)[ll]>redshift)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }
  if (j>0) {
    double dy = (redshift-(*z)[j])/((*z)[j+1]-(*z)[j]);
    double tmp_up = ((*chi)[j+2]-(*chi)[j])/((*z)[j+2]-(*z)[j]);
    double tmp_down = ((*chi)[j+1]-(*chi)[j-1])/((*z)[j+1]-(*z)[j-1]);
    out_dchi[1] = tmp_down + dy*(tmp_up-tmp_down);
  } else {
    double dy = (redshift-(*z)[j])/((*z)[j+1]-(*z)[j]);
    double tmp_up = ((*chi)[j+2]-(*chi)[j])/((*z)[j+2]-(*z)[j]);
    double tmp_down = ((*chi)[j+1]-(*chi)[j])/((*z)[j+1]-(*z)[j]);
    out_dchi[1] = tmp_down + dy*(tmp_up-tmp_down);
  }
  #endif
  // convert from (Mpc/h) to (Mpc/h)/(c/H0=100)^3 (dimensioneless)
  out_dchi[1] = (out_dchi[1]/cosmology.coverH0);
  // convert from d\chi/dz to d\chi/da
  return out_dchi[1]/(io_a*io_a);
}

struct chis chi_all(double io_a) {
  static cosmopara C;
  static int first = 0;
  double redshift;
  double out_dchi[4];
  #if defined(NAGLIB)
  static Nag_Spline spline;
  if (first == 0) {
    spline.lamda = NULL;
    spline.c = NULL;
    first = 1;
  }
  #else
  static int nz;
  static double** z;
  static double** chi;
  if (first == 0) {
    z = NULL;
    chi = NULL;
    first = 1;
  }
  #endif
  if (recompute_cosmo3D(C)) {
    update_cosmopara(&C);
    #if defined(NAGLIB)
    int nz;
    double** z = (double**) NAG_ALLOC(1, double*);
    double** chi = (double**) NAG_ALLOC(1, double*);
    #else
    if(z != NULL) {
      free(z);
      z = NULL;
    }
    z = (double**) malloc(1*sizeof(double*));
    if(chi != NULL) {
      free(chi);
      chi = NULL;
    }
    chi = (double**) malloc(1*sizeof(double*));
    #endif
    if(z == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*z) = NULL;
    if(chi == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*chi) = NULL;
    #if defined(NAGLIB)
    if(spline.lamda != NULL) {
      NAG_FREE(spline.lamda);
      spline.lamda = NULL;
    }
    if(spline.c != NULL) {
      NAG_FREE(spline.c);
      spline.c = NULL;
    }
    #endif
    setup_chi(&nz, z, chi, 0);
    #if defined(NAGLIB)
    NagError fail;
    INIT_FAIL(fail);
    nag_1d_spline_interpolant(nz, (*z), (*chi), &spline, &fail);
    if(fail.code != NE_NOERROR) {
      log_fatal("NAG interpolation failed %s", fail.message);
      exit(1);
    }
    NAG_FREE((*z));
    NAG_FREE(z);
    NAG_FREE((*chi));
    NAG_FREE(chi);
    #endif
  }
  redshift = 1.0/io_a - 1.0;
  #if defined(NAGLIB)
  NagError fail;
  Nag_DerivType derivs;
  INIT_FAIL(fail);
  derivs = Nag_LeftDerivs;
  nag_1d_spline_deriv(derivs,redshift,out_dchi,&spline,&fail);
  if(fail.code != NE_NOERROR) {
    log_fatal("NAG interpolation failed %s", fail.message);
    exit(1);
  }
  #else
  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = nz-1;
    while (ihi>ilo+1) {
      size_t ll = (ihi+ilo)/2;
      if((*z)[ll]>redshift)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }
  double dy = (redshift-(*z)[j])/((*z)[j+1]-(*z)[j]);
  out_dchi[0] = (*chi)[j] + dy*((*chi)[j+1]-(*chi)[j]);
  if (j>0) {
    double tmp_up = ((*chi)[j+2]-(*chi)[j])/((*z)[j+2]-(*z)[j]);
    double tmp_down = ((*chi)[j+1]-(*chi)[j-1])/((*z)[j+1]-(*z)[j-1]);
    out_dchi[1] = tmp_down + dy*(tmp_up-tmp_down);
  } else {
    double tmp_up = ((*chi)[j+2]-(*chi)[j])/((*z)[j+2]-(*z)[j]);
    double tmp_down = ((*chi)[j+1]-(*chi)[j])/((*z)[j+1]-(*z)[j]);
    out_dchi[1] = tmp_down + dy*(tmp_up-tmp_down);
  }
  #endif
  // convert from (Mpc/h) to (Mpc/h)/(c/H0=100)^3 (dimensioneless)
  out_dchi[1] = (out_dchi[1]/cosmology.coverH0);
  // convert from d\chi/dz to d\chi/da
  out_dchi[1] = out_dchi[1]/(io_a*io_a);
  /* output result */
  struct chis result;
  result.chi = out_dchi[0]/cosmology.coverH0;
  result.dchida = out_dchi[1];
  return result;
}

double dchi_dz(double io_a) {
  return (io_a*io_a)*dchi_da(io_a);
}

double hoverh0(double io_a) {
  return 1.0/dchi_dz(io_a);
}

double hoverh0v2(double io_a, double dchida) {
  return 1.0/((io_a*io_a)*dchida);
}

double a_chi(double io_chi) {
  static cosmopara C;
  static int first = 0;
  double out_redshift;
  #if defined(NAGLIB)
  static Nag_Spline spline;
  if (first == 0) {
    spline.lamda = NULL;
    spline.c = NULL;
    first = 1;
  }
  #else
  static int nz;
  static double** z;
  static double** chi;
  if (first == 0) {
    z = NULL;
    chi = NULL;
    first = 1;
  }
  #endif

  // convert from (Mpc/h)/(c/H0=100)^3 (dimensioneless) to (Mpc/h)
  io_chi *= cosmology.coverH0;

  if (recompute_cosmo3D(C)) {
    update_cosmopara(&C);
    #if defined(NAGLIB)
    int nz;
    double** z = (double**) NAG_ALLOC(1, double*);
    double** chi = (double**) NAG_ALLOC(1, double*);
    #else
    if(z != NULL) {
      free(z);
      z = NULL;
    }
    z = (double**) malloc(1*sizeof(double*));
    if(chi != NULL) {
      free(chi);
      chi = NULL;
    }
    chi = (double**) malloc(1*sizeof(double*));
    #endif
    if(z == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*z) = NULL;
    if(chi == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*chi) = NULL;
    #if defined(NAGLIB)
    if(spline.lamda != NULL) {
      NAG_FREE(spline.lamda);
      spline.lamda = NULL;
    }
    if(spline.c != NULL) {
      NAG_FREE(spline.c);
      spline.c = NULL;
    }
    #endif
    setup_chi(&nz, z, chi, 0);
    #if defined(NAGLIB)
    NagError fail;
    INIT_FAIL(fail);
    nag_1d_spline_interpolant(nz, (*chi), (*z), &spline, &fail);
    if(fail.code != NE_NOERROR) {
      log_fatal("NAG interpolation failed %s", fail.message);
      exit(1);
    }
    NAG_FREE((*z));
    NAG_FREE(z);
    NAG_FREE((*chi));
    NAG_FREE(chi);
    #endif
  }

  #if defined(NAGLIB)
  NagError fail;
  INIT_FAIL(fail);
  nag_1d_spline_evaluate(io_chi, &out_redshift, &spline, &fail);
  if(fail.code != NE_NOERROR) {
    log_fatal("NAG interpolation failed %s", fail.message);
    exit(1);
  }
  #else
  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = nz-1;
    while (ihi>ilo+1) {
      size_t ll = (ihi+ilo)/2;
      if((*chi)[ll]>io_chi)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }
  double dy = (io_chi-(*chi)[j])/((*chi)[j+1]-(*chi)[j]);
  out_redshift = (*z)[j] + dy*((*z)[j+1]-(*z)[j]);
  #endif
  return 1.0/(1.0+out_redshift);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Growth Factor
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void setup_growth(int* io_nz, double** io_z, double** io_G, int io) {
  static int nz;
  static double* z = NULL;
  static double* G = NULL;
  if (io == 1) {
    // IO == 1 IMPLES THAT IO_G(Z) WILL COPIED TO LOCAL D(Z)
    if (z != NULL) {
      #if defined(NAGLIB)
      NAG_FREE(z);
      #else
      free(z);
      #endif
    }
    if (G != NULL) {
      #if defined(NAGLIB)
      NAG_FREE(G);
      #else
      free(G);
      #endif
    }
    nz = (*io_nz);
    if (!(nz > 5)) {
      log_fatal("array to small");
      exit(1);
    }
    #if defined(NAGLIB)
    z = (double*) NAG_ALLOC(nz, double);
    G = (double*) NAG_ALLOC(nz, double);
    #else
    z = (double*) malloc(nz*sizeof(double));
    G = (double*) malloc(nz*sizeof(double));
    #endif
    if (z == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    if (G == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    // Set Values
    for (int i=0; i<nz; i++) {
      z[i] = (*io_z)[i];
      G[i] = (*io_G)[i];
    }
  }
  else {
    // IO != 1 IMPLES THAT LOCAL D(Z) WILL BE COPIED TO IO_G(Z)
    if (z == NULL || G == NULL) {
      log_fatal("growth not setup from Cobaya\n");
      exit(1);
    }
    (*io_nz) = nz;
    if (io_z != NULL) {
      if((*io_z) != NULL) {
        #if defined(NAGLIB)
        NAG_FREE((*io_z));
        #else
        free((*io_z));
        #endif
        (*io_z) = NULL;
      }
    }
    else {
      log_fatal("input pointer not allocated\n");
      exit(1);
    }
    if (io_G != NULL) {
      if((*io_G) != NULL) {
        #if defined(NAGLIB)
        NAG_FREE((*io_G));
        #else
        free((*io_G));
        #endif
        (*io_G) = NULL;
      }
    }
    else {
      log_fatal("input pointer not allocated\n");
      exit(1);
    }
    #if defined(NAGLIB)
    (*io_z) = (double*) NAG_ALLOC(nz, double);
    (*io_G) = (double*) NAG_ALLOC(nz, double);
    #else
    (*io_z) = (double*) malloc(nz*sizeof(double));
    (*io_G) = (double*) malloc(nz*sizeof(double));
    #endif
    if ((*io_z) == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    if ((*io_G) == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    for (int i=0; i<nz; i++) {
      (*io_z)[i] = z[i];
      (*io_G)[i] = G[i];
    }
  }
}

double growfac(double io_a) {
  return norm_growfac(io_a, true);
}

double norm_growfac(double io_a, bool normalize_z0) {
  static cosmopara C;
  static double growfact1;
  static int first = 0;
  double redshift, out_G;
  #if defined(NAGLIB)
  static Nag_Spline spline;
  if (first == 0) {
    spline.lamda = NULL;
    spline.c = NULL;
    first = 1;
  }
  #else
  static int nz;
  static double** z;
  static double** G;
  if (first == 0) {
    z = NULL;
    G = NULL;
    first = 1;
  }
  #endif
  if (recompute_cosmo3D(C)) {
    update_cosmopara(&C);
    #if defined(NAGLIB)
    int nz;
    double** z = (double**) NAG_ALLOC(1, double*);
    double** G = (double**) NAG_ALLOC(1,double*);
    #else
    if(z != NULL) {
      free(z);
      z = NULL;
    }
    z = (double**) malloc(1*sizeof(double*));
    if(G != NULL) {
      free(G);
      G = NULL;
    }
    G = (double**) malloc(1*sizeof(double*));
    #endif
    if (z == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*z) = NULL;
    if (G == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*G) = NULL;
    #if defined(NAGLIB)
    if (spline.lamda != NULL) {
      NAG_FREE(spline.lamda);
      spline.lamda = NULL;
    }
    if (spline.c != NULL) {
      NAG_FREE(spline.c);
      spline.c = NULL;
    }
    #endif
    setup_growth(&nz, z, G, 0);
    #if defined(NAGLIB)
    NagError fail;
    INIT_FAIL(fail);
    nag_1d_spline_interpolant(nz, (*z), (*G), &spline, &fail);
    if (fail.code != NE_NOERROR) {
      log_fatal("NAG interpolation failed %s", fail.message);
      exit(1);
    }
    NAG_FREE((*z));
    NAG_FREE(z);
    NAG_FREE((*G));
    NAG_FREE(G);
    {
      NagError fail;
      INIT_FAIL(fail);
      redshift = 0.0;
      nag_1d_spline_evaluate(redshift, &growfact1, &spline, &fail);
      if (fail.code != NE_NOERROR) {
        log_fatal("NAG interpolation failed %s", fail.message);
        exit(1);
      }
    }
    #else
    {
      redshift = 0.0;
      int j = 0;
      {
        size_t ilo = 0;
        size_t ihi = nz-1;
        while (ihi>ilo+1) {
          size_t ll = (ihi+ilo)/2;
          if((*z)[ll]>redshift)
            ihi = ll;
          else
            ilo = ll;
        }
        j = ilo;
      }
      double dy = (redshift-(*z)[j])/((*z)[j+1]-(*z)[j]);
      growfact1 = (*G)[j] + dy*((*G)[j+1]-(*G)[j]);
    }
    #endif
  }
  redshift = 1.0/io_a-1.0;
  #if defined(NAGLIB)
  NagError fail;
  INIT_FAIL(fail);
  nag_1d_spline_evaluate(redshift, &out_G, &spline, &fail);
  if (fail.code != NE_NOERROR) {
    log_fatal("NAG interpolation failed %s", fail.message);
    exit(1);
  }
  #else
  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = nz-1;
    while (ihi>ilo+1) {
      size_t ll = (ihi+ilo)/2;
      if((*z)[ll]>redshift)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }
  double dy = (redshift-(*z)[j])/((*z)[j+1]-(*z)[j]);
  out_G = (*G)[j] + dy*((*G)[j+1]-(*G)[j]);
  #endif
  if(normalize_z0) {
    return (out_G*io_a)/growfact1; // Growth D = G * a
  } else {
    return out_G*io_a; // Growth D = G * a
  }
}

double f_growth(double redshift) {
  static cosmopara C;
  static int first = 0;
  double out_G, out_dlnGdlna;
  #if defined(NAGLIB)
  static Nag_Spline spline;
  static Nag_DerivType derivs;
  if (first == 0) {
    spline.lamda = NULL;
    spline.c = NULL;
    first = 1;
  }
  #else
  static int nz;
  static double** z;
  static double** G;
  if (first == 0) {
    z = NULL;
    G = NULL;
    first = 1;
  }
  #endif
  if (recompute_cosmo3D(C)) {
    update_cosmopara(&C);
    #if defined(NAGLIB)
    int nz;
    double** z = (double**) NAG_ALLOC(1, double*);
    double** G = (double**) NAG_ALLOC(1,double*);
    #else
    if(z != NULL) {
      free(z);
      z = NULL;
    }
    z = (double**) malloc(1*sizeof(double*));
    if(G != NULL) {
      free(G);
      G = NULL;
    }
    G = (double**) malloc(1*sizeof(double*));
    #endif
    if (z == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*z) = NULL;
    if (G == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*G) = NULL;
    #if defined(NAGLIB)
    if (spline.lamda != NULL) {
      NAG_FREE(spline.lamda);
      spline.lamda = NULL;
    }
    if (spline.c != NULL) {
      NAG_FREE(spline.c);
      spline.c = NULL;
    }
    #endif
    setup_growth(&nz, z, G, 0);
  }
  #if defined(NAGLIB)
  NagError fail;
  INIT_FAIL(fail);
  derivs = Nag_LeftDerivs;
  nag_1d_spline_deriv(derivs, redshift, s, &spline, &fail);
  if (fail.code != NE_NOERROR) {
    log_fatal("NAG interpolation failed %s", fail.message);
    exit(1);
  }
  out_G = s[0]
  out_dGdlna = -s[1]*(1.0+redshift)
  #else
  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = nz-1;
    while (ihi>ilo+1) {
      size_t ll = (ihi+ilo)/2;
      if((*z)[ll]>redshift)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }
  double dy = (redshift-(*z)[j])/((*z)[j+1]-(*z)[j]);
  out_G = (*G)[j] + dy*((*G)[j+1]-(*G)[j]);
  const double dlnGdlnz =
    (((*G)[j+1]-(*G)[j])/((*z)[j+1]-(*z)[j]))*redshift/out_G;
  out_dlnGdlna = -dlnGdlnz*(1+redshift)/redshift;
  #endif
  return 1 + out_dlnGdlna; // Growth D = G * a
}

struct growths norm_growfac_all(double io_a, bool normalize_z0) {
  static cosmopara C;
  static double growfact1;
  static int first = 0;
  double redshift;
  #if defined(NAGLIB)
  static Nag_Spline spline;
  if (first == 0) {
    spline.lamda = NULL;
    spline.c = NULL;
    first = 1;
  }
  #else
  static int nz;
  static double** z;
  static double** G;
  if (first == 0) {
    z = NULL;
    G = NULL;
    first = 1;
  }
  #endif
  if (recompute_cosmo3D(C)) {
    update_cosmopara(&C);
    #if defined(NAGLIB)
    int nz;
    double** z = (double**) NAG_ALLOC(1, double*);
    double** G = (double**) NAG_ALLOC(1,double*);
    #else
    if(z != NULL) {
      free(z);
      z = NULL;
    }
    z = (double**) malloc(1*sizeof(double*));
    if(G != NULL) {
      free(G);
      G = NULL;
    }
    G = (double**) malloc(1*sizeof(double*));
    #endif
    if (z == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*z) = NULL;
    if (G == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*G) = NULL;
    #if defined(NAGLIB)
    if (spline.lamda != NULL) {
      NAG_FREE(spline.lamda);
      spline.lamda = NULL;
    }
    if (spline.c != NULL) {
      NAG_FREE(spline.c);
      spline.c = NULL;
    }
    #endif
    setup_growth(&nz, z, G, 0);
    #if defined(NAGLIB)
    NagError fail;
    INIT_FAIL(fail);
    nag_1d_spline_interpolant(nz, (*z), (*G), &spline, &fail);
    if (fail.code != NE_NOERROR) {
      log_fatal("NAG interpolation failed %s", fail.message);
      exit(1);
    }
    NAG_FREE((*z));
    NAG_FREE(z);
    NAG_FREE((*G));
    NAG_FREE(G);
    {
      NagError fail;
      INIT_FAIL(fail);
      redshift = 0.0;
      nag_1d_spline_evaluate(redshift, &growfact1, &spline, &fail);
      if (fail.code != NE_NOERROR) {
        log_fatal("NAG interpolation failed %s", fail.message);
        exit(1);
      }
    }
    #else
    {
      redshift = 0.0;
      int j = 0;
      {
        size_t ilo = 0;
        size_t ihi = nz-1;
        while (ihi>ilo+1) {
          size_t ll = (ihi+ilo)/2;
          if((*z)[ll]>redshift)
            ihi = ll;
          else
            ilo = ll;
        }
        j = ilo;
      }
      double dy = (redshift-(*z)[j])/((*z)[j+1]-(*z)[j]);
      growfact1 = (*G)[j] + dy*((*G)[j+1]-(*G)[j]);
    }
    #endif
  }
  redshift = 1.0/io_a-1.0;
  #if defined(NAGLIB)
  NagError fail;
  INIT_FAIL(fail);
  nag_1d_spline_evaluate(redshift, &out_G, &spline, &fail);
  if (fail.code != NE_NOERROR) {
    log_fatal("NAG interpolation failed %s", fail.message);
    exit(1);
  }
  #else
  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = nz-1;
    while (ihi>ilo+1) {
      size_t ll = (ihi+ilo)/2;
      if((*z)[ll]>redshift)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }
  double dy = (redshift-(*z)[j])/((*z)[j+1]-(*z)[j]);
  const double out_G = (*G)[j] + dy*((*G)[j+1]-(*G)[j]);
  const double dlnGdlnz =
    (((*G)[j+1]-(*G)[j])/((*z)[j+1]-(*z)[j]))*redshift/out_G;
  const double out_dlnGdlna = -dlnGdlnz*(1+redshift)/redshift;
  #endif
  struct growths out_Gf;
  out_Gf.f = 1.0 + out_dlnGdlna;
  if(normalize_z0) {
    out_Gf.D = (out_G*io_a)/growfact1;
  } else {
    out_Gf.D = (out_G*io_a);
  }
  return out_Gf;
}

struct growths growfac_all(double io_a) {
  return norm_growfac_all(io_a, true);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Power Spectrum
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void setup_p_lin(int* io_nlog10k, int* io_nz, double** io_log10k,
                                  double** io_z, double** io_lnP, int io) {
  static int nlog10k;
  static int nz;
  static double* log10k = NULL;
  static double* z = NULL;
  static double* lnP = NULL;
  if(io == 1) {
    // IO == 1 IMPLES THAT IO_PK(Z,LOGK) WILL COPIED TO LOCAL PK(Z,LOGK)
    if (log10k != NULL) {
      #if defined(NAGLIB)
      NAG_FREE(log10k);
      #else
      free(log10k);
      #endif
    }
    if (z != NULL) {
      #if defined(NAGLIB)
      NAG_FREE(z);
      #else
      free(z);
      #endif
    }
    if (lnP != NULL) {
      #if defined(NAGLIB)
      NAG_FREE(lnP);
      #else
      free(lnP);
      #endif
    }
    nlog10k = (*io_nlog10k);
    nz = (*io_nz);
    if (!(nlog10k > 5)) {
      log_fatal("array to small");
      exit(1);
    }
    if (!(nz > 5)) {
      log_fatal("array to small");
      exit(1);
    }
    #if defined(NAGLIB)
    log10k = (double*) NAG_ALLOC(nlog10k, double);
    z = (double*) NAG_ALLOC(nz, double);
    lnP = (double*) NAG_ALLOC(nz*nlog10k, double);
    #else
    log10k = (double*) malloc(nlog10k*sizeof(double));
    z = (double*) malloc(nz*sizeof(double));
    lnP = (double*) malloc(nz*nlog10k*sizeof(double));
    #endif
    if (log10k == NULL || z == NULL || lnP == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    // Set Values
    for (int j = 0; j < nz; j++) {
      z[j] = (*io_z)[j];
    }
    for (int i=0; i<nlog10k; i++) {
      log10k[i] = (*io_log10k)[i];
      for(int j=0; j<nz; j++) {
        lnP[i*nz+j] = (*io_lnP)[i*nz+j];
      }
    }
  }
  else {
    // IO != 1 IMPLES THAT LOCAL PK(Z,LOGK) WILL COPIED TO IO_PK(Z,LOGK)
    if (log10k == NULL || z == NULL || lnP == NULL) {
      log_fatal("pk linear not setup from cobaya\n");
      exit(1);
    }
    (*io_nlog10k) = nlog10k;
    if (io_log10k != NULL) {
      if((*io_log10k) != NULL) {
        #if defined(NAGLIB)
        NAG_FREE((*io_log10k));
        #else
        free((*io_log10k));
        #endif
        (*io_log10k) = NULL;
      }
    }
    else {
      log_fatal("input pointer not allocated\n");
      exit(1);
    }
    (*io_nz) = nz;
    if (io_z != NULL) {
      if((*io_z) != NULL) {
        #if defined(NAGLIB)
        NAG_FREE((*io_z));
        #else
        free((*io_z));
        #endif
        (*io_z) = NULL;
      }
    }
    else {
      log_fatal("input pointer not allocated\n");
      exit(1);
    }
    if (io_lnP != NULL) {
      if ((*io_lnP) != NULL) {
        #if defined(NAGLIB)
        NAG_FREE((*io_lnP));
        #else
        free((*io_lnP));
        #endif
        (*io_lnP) = NULL;
      }
    }
    else {
      log_fatal("input pointer not allocated\n");
      exit(1);
    }
    #if defined(NAGLIB)
    (*io_log10k) = (double*) NAG_ALLOC(nlog10k, double);
    (*io_z) = (double*) NAG_ALLOC(nz, double);
    (*io_lnP) = (double*) NAG_ALLOC(nz*nlog10k, double);
    #else
    (*io_log10k) = (double*) malloc(nlog10k*sizeof(double));
    (*io_z) = (double*) malloc(nz*sizeof(double));
    (*io_lnP) = (double*) malloc(nz*nlog10k*sizeof(double));
    #endif
    if ((*io_log10k) == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    if ((*io_z) == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    if ((*io_lnP) == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    for (int j=0; j<nz; j++) {
      (*io_z)[j] = z[j];
    }
    for (int i=0; i<nlog10k; i++) {
      (*io_log10k)[i] = log10k[i];
      for (int j=0; j<nz; j++) {
        (*io_lnP)[i*nz+j] = lnP[i*nz+j];
      }
    }
  }
}

void setup_p_nonlin(int* io_nlog10k, int* io_nz, double** io_log10k,
                                  double** io_z, double** io_lnP, int io) {
  static int nlog10k;
  static int nz;
  static double* log10k = NULL;
  static double* z = NULL;
  static double* lnP = NULL;
  if (io == 1) {
    // IO == 1 IMPLES THAT IO_PK(Z,K) WILL COPIED TO LOCAL PK(Z,K)
    if (log10k != NULL) {
      #if defined(NAGLIB)
      NAG_FREE(log10k);
      #else
      free(log10k);
      #endif
    }
    if (z != NULL) {
      #if defined(NAGLIB)
      NAG_FREE(z);
      #else
      free(z);
      #endif
    }
    if (lnP != NULL) {
      #if defined(NAGLIB)
      NAG_FREE(lnP);
      #else
      free(lnP);
      #endif
    }
    nlog10k = (*io_nlog10k);
    nz = (*io_nz);
    if (!(nlog10k > 5)) {
      log_fatal("array to small");
      exit(1);
    }
    if (!(nz > 5)) {
      log_fatal("array to small");
      exit(1);
    }
    #if defined(NAGLIB)
    log10k = (double*) NAG_ALLOC(nlog10k, double);
    z = (double*) NAG_ALLOC(nz, double);
    lnP = (double*) NAG_ALLOC(nz*nlog10k, double);
    #else
    log10k = (double*) malloc(nlog10k*sizeof(double));
    z = (double*) malloc(nz*sizeof(double));
    lnP = (double*) malloc(nz*nlog10k*sizeof(double));
    #endif
    if (log10k == NULL || z == NULL || lnP == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    for (int j=0; j<nz; j++) {
      z[j] = (*io_z)[j];
    }
    for (int i=0; i<nlog10k; i++) {
      log10k[i] = (*io_log10k)[i];
      for (int j=0; j<nz; j++) {
        lnP[i*nz+j] = (*io_lnP)[i*nz+j];
      }
    }
  }
  else {
    // IO != 1 IMPLES THAT LOCAL PK(Z,K) WILL COPIED TO IO_PK(Z,K)
    if (log10k == NULL || z == NULL || lnP == NULL) {
      log_fatal("input pointer not setup\n");
      exit(1);
    }
    (*io_nlog10k) = nlog10k;
    if (io_log10k != NULL) {
      if((*io_log10k) != NULL) {
        #if defined(NAGLIB)
        NAG_FREE((*io_log10k));
        #else
        free((*io_log10k));
        #endif
        (*io_log10k) = NULL;
      }
    }
    else {
      log_fatal("input pointer not allocated\n");
      exit(1);
    }
    (*io_nz) = nz;
    if (io_z != NULL) {
      if((*io_z) != NULL) {
        #if defined(NAGLIB)
        NAG_FREE((*io_z));
        #else
        free((*io_z));
        #endif
        (*io_z) = NULL;
      }
    }
    else {
      log_fatal("input pointer not allocated\n");
      exit(1);
    }
    if (io_lnP != NULL) {
      if((*io_lnP) != NULL) {
        #if defined(NAGLIB)
        NAG_FREE((*io_lnP));
        #else
        free((*io_lnP));
        #endif
        (*io_lnP) = NULL;
      }
    }
    else {
      log_fatal("input pointer not allocated\n");
      exit(1);
    }
    #if defined(NAGLIB)
    (*io_log10k) = (double*) NAG_ALLOC(nlog10k, double);
    (*io_z) = (double*) NAG_ALLOC(nz, double);
    (*io_lnP) = (double*) NAG_ALLOC(nz*nlog10k, double);
    #else
    (*io_log10k) = (double*) malloc(nlog10k*sizeof(double));
    (*io_z) = (double*) malloc(nz*sizeof(double));
    (*io_lnP) = (double*) malloc(nz*nlog10k*sizeof(double));
    #endif
    if ((*io_log10k) == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    if ((*io_z) == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    if ((*io_lnP) == NULL ) {
      log_fatal("fail allocation");
      exit(1);
    }
    for (int j = 0; j<nz; j++) {
      (*io_z)[j] = z[j];
    }
    for (int i=0; i<nlog10k; i++) {
      (*io_log10k)[i] = log10k[i];
      for (int j=0; j<nz; j++) {
        (*io_lnP)[i*nz+j] = lnP[i*nz+j];
      }
    }
  }
}

double p_lin(double io_k, double io_a) {
  // k in units H_0/c (dimensioneless)
  static cosmopara C;
  static int first = 0;
  const int size = 1;
  double tmp_log10k[size];
  double tmp_z[size];
  double out_lnP[size];
  static double limits_z[2];
  #if defined(NAGLIB)
  static Nag_2dSpline spline;
  if (first == 0) {
    spline.lamda = NULL;
    spline.mu = NULL;
    spline.c = NULL;
    first = 1;
  }
  #else
  static double** log10k;
  static double** z;
  static double** lnP;
  static int nlog10k;
  static int nz;
  if (first == 0) {
    z = NULL;
    lnP = NULL;
    log10k = NULL;
    first = 1;
  }
  #endif
  if (recompute_cosmo3D(C)) {
    update_cosmopara(&C);
    #if defined(NAGLIB)
    int nlog10k;
    int nz;
    double** log10k = (double**) NAG_ALLOC(1, double*);
    double** z = (double**) NAG_ALLOC(1, double*);
    double** lnP = (double**) NAG_ALLOC(1, double*);
    #else
    if (log10k != NULL) {
      free((*log10k));
      free(log10k);
    }
    log10k = (double**) malloc(1*sizeof(double*));
    if (z != NULL) {
      free((*z));
      free(z);
    }
    z = (double**) malloc(1*sizeof(double*));
    if (lnP != NULL) {
      free((*lnP));
      free(lnP);
    }
    lnP = (double**) malloc(1*sizeof(double*));
    #endif
    if (log10k == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*log10k) = NULL;
    if (z == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*z) = NULL;
    if (lnP == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*lnP) = NULL;
    #if defined(NAGLIB)
    if (spline.lamda != NULL) {
      NAG_FREE(spline.lamda);
      spline.lamda = NULL;
    }
    if (spline.mu != NULL) {
      NAG_FREE(spline.mu);
      spline.mu = NULL;
    }
    if (spline.c != NULL) {
      NAG_FREE(spline.c);
      spline.c = NULL;
    }
    #endif
    setup_p_lin(&nlog10k, &nz, log10k, z, lnP, 0);
    limits_z[0] = z[0][0];
    limits_z[1] = z[0][nz-1];
    #if defined(NAGLIB)
    NagError fail;
    INIT_FAIL(fail);
    nag_2d_spline_interpolant((Integer) nlog10k, (Integer) nz, (*log10k),
                                             (*z), (*lnP), &spline, &fail);
    if (fail.code != NE_NOERROR) {
      log_fatal("NAG interpolation failed %s", fail.message);
      exit(1);
    }
    // Deallocation of 2D array
    NAG_FREE((*log10k));
    NAG_FREE(log10k);
    NAG_FREE((*z));
    NAG_FREE(z);
    NAG_FREE((*lnP));
    NAG_FREE(lnP);
    #endif
  }
  // convert from (x/Mpc/h - dimensioneless) to h/Mpc with x = c/H0 (Mpc)
  tmp_log10k[0] = log10(io_k/cosmology.coverH0);
  tmp_z[0] = 1.0/io_a-1.0;

  if ((tmp_z[0]<limits_z[0]) || (tmp_z[0]>limits_z[1]))  {
    log_fatal("redshift out of bounds z_min = %4.3e, z_max = %4.3e and z = %4.3e",
      limits_z[0], limits_z[1], tmp_z[0]);
    exit(1);
  }

  #if defined(NAGLIB)
  NagError fail;
  INIT_FAIL(fail);
  nag_2d_spline_eval((Integer) size, tmp_log10k, tmp_z, out_lnP,
                                                           &spline, &fail);
  if (fail.code != NE_NOERROR) {
    log_fatal("NAG interpolation failed %s", fail.message);
    exit(1);
  }
  #else

  int i = 0;
  {
    size_t ilo = 0;
    size_t ihi = nlog10k-2;
    while (ihi>ilo+1) {
      size_t ll = (ihi+ilo)/2;
      if((*log10k)[ll]>tmp_log10k[0])
        ihi = ll;
      else
        ilo = ll;
    }
    i = ilo;
  }

  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = nz-2;
    while (ihi>ilo+1) {
      size_t ll = (ihi+ilo)/2;
      if((*z)[ll]>tmp_z[0])
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }

  double dx = (tmp_log10k[0]-(*log10k)[i])/((*log10k)[i+1]-(*log10k)[i]);

  double dy = (tmp_z[0]-(*z)[j])/((*z)[j+1]-(*z)[j]);

  out_lnP[0]=(1-dx)*(1-dy)*(*lnP)[i*nz+j]+(1-dx)*dy*(*lnP)[i*nz+(j+1)]
                   +dx*(1-dy)*(*lnP)[(i+1)*nz+j]+dx*dy*(*lnP)[(i+1)*nz+(j+1)];
  #endif

  // convert from (Mpc/h)^3 to (Mpc/h)^3/(c/H0=100)^3 (dimensioneless)
  return exp(out_lnP[0])/pow(cosmology.coverH0,3.0);
}

double p_nonlin(double io_k, double io_a) {
  // k in units H_0/c
  static cosmopara C;
  const int size = 1;
  double tmp_log10k[size];
  double tmp_z[size];
  double out_lnP[size];
  static double limits_z[2];
  static int first = 0;
  #if defined(NAGLIB)
  static Nag_2dSpline spline;
  if (first == 0) {
    spline.lamda = NULL;
    spline.mu = NULL;
    spline.c = NULL;
    first = 1;
  }
  #else
  static double** log10k;
  static double** z;
  static double** lnP;
  static int nlog10k;
  static int nz;
  if (first == 0) {
    z = NULL;
    lnP = NULL;
    log10k = NULL;
    first = 1;
  }
  #endif
  if (recompute_cosmo3D(C)) {
    update_cosmopara(&C);
    #if defined(NAGLIB)
    int nlog10k;
    int nz;
    double** log10k = (double**) NAG_ALLOC(1, double*);
    double** z = (double**) NAG_ALLOC(1, double*);
    double** lnP = (double**) NAG_ALLOC(1, double*);
    #else
    if (log10k != NULL) {
      free((*log10k));
      free(log10k);
    }
    log10k = (double**) malloc(1*sizeof(double*));
    if (z != NULL) {
      free((*z));
      free(z);
    }
    z = (double**) malloc(1*sizeof(double*));
    if (lnP != NULL) {
      free((*lnP));
      free(lnP);
    }
    lnP = (double**) malloc(1*sizeof(double*));
    #endif
    if (log10k == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*log10k) = NULL;
    if (z == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*z) = NULL;
    if (lnP == NULL) {
      log_fatal("fail allocation");
      exit(1);
    }
    (*lnP) = NULL;
    #if defined(NAGLIB)
    if (spline.lamda != NULL) {
      NAG_FREE(spline.lamda);
      spline.lamda = NULL;
    }
    if (spline.mu != NULL) {
      NAG_FREE(spline.mu);
      spline.mu = NULL;
    }
    if (spline.c != NULL) {
      NAG_FREE(spline.c);
      spline.c = NULL;
    }
    #endif
    setup_p_nonlin(&nlog10k, &nz, log10k, z, lnP, 0);
    limits_z[0] = z[0][0];
    limits_z[1] = z[0][nz-1];
    #if defined(NAGLIB)
    NagError fail;
    INIT_FAIL(fail);
    nag_2d_spline_interpolant((Integer) nlog10k, (Integer) nz, (*log10k),
                                             (*z), (*lnP), &spline, &fail);
    if (fail.code != NE_NOERROR) {
      log_fatal("NAG interpolation failed %s", fail.message);
      exit(1);
    }
    // Deallocation of the 2D arrays
    NAG_FREE((*log10k));
    NAG_FREE(log10k);
    NAG_FREE((*z));
    NAG_FREE(z);
    NAG_FREE((*lnP));
    NAG_FREE(lnP);
    #endif
  }
  // convert to h/Mpc
  tmp_log10k[0] = log10(io_k/(cosmology.coverH0));
  tmp_z[0] = 1.0/io_a-1.0;
  if ((tmp_z[0]<limits_z[0]) || (tmp_z[0]>limits_z[1]))  {
    log_fatal("redshift out of bounds z_min = %4.3e, z_max = %4.3e and z = %4.3e",
      limits_z[0], limits_z[1], tmp_z[0]);
    exit(1);
  }
  #if defined(NAGLIB)
  NagError fail;
  INIT_FAIL(fail);
  nag_2d_spline_eval((Integer)size,tmp_log10k,tmp_z,out_lnP,&spline,&fail);
  if (fail.code != NE_NOERROR) {
    log_fatal("NAG interpolation failed %s", fail.message);
    exit(1);
  }
  #else
  int i = 0;
  {
    size_t ilo = 0;
    size_t ihi = nlog10k-2;
    while (ihi>ilo+1) {
      size_t ll = (ihi+ilo)/2;
      if((*log10k)[ll]>tmp_log10k[0])
        ihi = ll;
      else
        ilo = ll;
    }
    i = ilo;
  }
  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = nz-2;
    while (ihi>ilo+1) {
      size_t ll = (ihi+ilo)/2;
      if((*z)[ll]>tmp_z[0])
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }
  double dx = (tmp_log10k[0]-(*log10k)[i])/((*log10k)[i+1]-(*log10k)[i]);
  double dy = (tmp_z[0]-(*z)[j])/((*z)[j+1]-(*z)[j]);
  out_lnP[0]=(1-dx)*(1-dy)*(*lnP)[i*nz+j]+(1-dx)*dy*(*lnP)[i*nz+(j+1)]
                   +dx*(1-dy)*(*lnP)[(i+1)*nz+j]+dx*dy*(*lnP)[(i+1)*nz+(j+1)];
  #endif
  // convert from (Mpc/h)^3 to (Mpc/h)^3/(c/H0=100)^3 (dimensioneless)
  return exp(out_lnP[0])/pow(cosmology.coverH0,3);
}

double Pdelta(double io_kNL, double io_a) {
  double out_PK;
  static int P_type = -1;
  if (P_type == -1) {
    if (strcmp(pdeltaparams.runmode,"linear") == 0) {
      P_type = 3;
    }
  }
  switch (P_type) {
    case 3:
      out_PK = p_lin(io_kNL, io_a);
      break;
    default:
      out_PK = p_nonlin(io_kNL, io_a);
      break;
  }
  return out_PK;
}

// calculating the angular diameter distance f_K
// BS01 2.4, 2.30: f_K is a radial function that, depending on the curvature of
// the Universe, is a trigonometric, linear, or hyperbolic function of chi
double f_K(double chi) {
  double K, K_h, f;
  K = (cosmology.Omega_m + cosmology.Omega_v - 1.);
  if (K > precision.medium) { // open
    K_h = sqrt(K); // K in units H0/c see BS eq. 2.30
    f = 1. / K_h * sin(K_h * chi);
  } else if (K < -precision.medium) { // closed
    K_h = sqrt(-K);
    f = 1. / K_h * sinh(K_h * chi);
  } else { // flat
    f = chi;
  }
  return f;
}

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// Baryons
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

// return P(k)_bary/P(k)_DMO from hydro sims ; kintern in unit [h/Mpc]
double PkRatio_baryons(double kintern, double a) {
  FILE *infile;
  static barypara B;
  double logkin = log10(kintern);
  double res;
  static double *logk_bins = 0;
  static double *a_bins = 0;
  static double **TblogPkR = 0;

  if (bary.isPkbary == 0) {
    return 1.;
  }

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  gsl_interp2d *interp2d = gsl_interp2d_alloc(T, bary.Nkbins, bary.Nabins);
  double *GSLPKR = malloc(bary.Nkbins * bary.Nabins * sizeof(double));

  if (recompute_PkRatio(B)) {
    update_PkRatio(&B);

    printf("in recompute PkRatio \n");

    if (TblogPkR != 0) {
      free_double_matrix(TblogPkR, 0, bary.Nkbins - 1, 0, bary.Nabins - 1);
    }
    TblogPkR = create_double_matrix(0, bary.Nkbins - 1, 0, bary.Nabins - 1);

    if (logk_bins != 0) {
      free_double_vector(logk_bins, 0, bary.Nkbins - 1);
    }
    logk_bins = create_double_vector(0, bary.Nkbins - 1);

    if (a_bins != 0) {
      free_double_vector(logk_bins, 0, bary.Nabins - 1);
    }
    a_bins = create_double_vector(0, bary.Nabins - 1);

    for (int i = 0; i < bary.Nabins; i++) {
      a_bins[i] = 1. / (1 + bary.z_bins[i]);
      // printf("a: %le,z: %le\n",a[i],z[i]);
    }

    infile = fopen(bary.FILE_logPkR, "r");
    if (infile == NULL) {
      printf("Error opening logPkRatio file\n");
      exit(1);
    }

    fscanf(infile, "%*[^\n]"); // Read and discard the 1st line

    for (int i = 0; i < bary.Nkbins; i++) {
      fscanf(infile, "%le ", &logk_bins[i]);
      for (int j = 0; j < bary.Nabins; j++) {
        fscanf(infile, "%le ", &TblogPkR[i][j]);
      }
    }
    fclose(infile);
  }

  for (int i = 0; i < bary.Nkbins; i++) {
    for (int j = 0; j < bary.Nabins; j++) {
      gsl_interp2d_set(interp2d, GSLPKR, i, j, TblogPkR[i][j]);
    }
  }

  gsl_interp2d_init(interp2d, logk_bins, a_bins, GSLPKR, bary.Nkbins,
                    bary.Nabins);

  // allow extrapolation beyond k>1500 Pk_ratio
  res = gsl_interp2d_eval_extrap(interp2d, logk_bins, a_bins, GSLPKR, logkin, a,
                                 NULL, NULL);
  res = pow(10, res);

  gsl_interp2d_free(interp2d);
  free(GSLPKR);
  return res;
}

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// MODIFIED GRAVITY
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

double MG_Sigma(double a) {
  return 0.0;
}