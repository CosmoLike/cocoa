#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../cfftlog/cfftlog.h"
#include <fftw3.h>

#include "bias.h"
#include "basics.h"
#include "cosmo3D.h"
#include "cosmo2D.h"
#include "IA.h"
#include "pt_cfastpt.h"
#include "radial_weights.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"

static int GSL_WORKSPACE_SIZE = 250;
static int use_linear_ps_limber = 0; /* 0 or 1 */
static int include_RSD_GS = 0; /* 0 or 1 */
static int include_RSD_GG = 1; /* 0 or 1 */
static int include_RSD_GK = 0; /* 0 or 1 */
static double w_l_min = 0.0001;
static double w_l_max = 5.0e6;

double beam_cmb(const double l)
{
  const double sigma = cmb.fwhm/sqrt(8.0*constants.ln2);
  return exp(-0.5*l*l*sigma*sigma);
}

double C_gk_tomo_limber_nointerp_wrapper(double l, int ni, int use_linear_ps)
{
  return C_gk_tomo_limber_nointerp(l, ni, use_linear_ps)*beam_cmb(l);
}

double C_gk_tomo_limber_wrapper(double l, int ni)
{
  return C_gk_tomo_limber(l, ni)*beam_cmb(l);
}

double C_ks_tomo_limber_nointerp_wrapper(double l, int ni, int use_linear_ps)
{
  return C_ks_tomo_limber_nointerp(l, ni, use_linear_ps)*beam_cmb(l);
}

double C_ks_tomo_limber_wrapper(double l, int ni)
{
  return C_ks_tomo_limber(l, ni)*beam_cmb(l);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real Space) - Full Sky - bin average
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

typedef struct 
{
  double xmin;
  double xmax;
  double Pmin;
  double Pmax;
  double dPmin;
  double dPmax;
} bin_avg;

bin_avg set_bin_average(int i_theta, int j_L)
{
  static double** Pmin  = 0;
  static double** Pmax  = 0;
  static double** dPmin = 0;
  static double** dPmax = 0;
  static double* xmin = 0;
  static double* xmax = 0;
  
  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }

  if (Pmin == 0)
  {
    Pmin = (double**) malloc(like.Ntheta*sizeof(double*));
    Pmax = (double**) malloc(like.Ntheta*sizeof(double*));
    dPmin = (double**) malloc(like.Ntheta*sizeof(double*));
    dPmax = (double**) malloc(like.Ntheta*sizeof(double*));
    for(int i=0; i<like.Ntheta ; i++)
    {
      Pmin[i] = (double*) calloc(limits.LMAX, sizeof(double));
      Pmax[i] = (double*) calloc(limits.LMAX, sizeof(double));
      dPmin[i] = (double*) calloc(limits.LMAX, sizeof(double));
      dPmax[i] = (double*) calloc(limits.LMAX, sizeof(double));
    }
    xmin = (double*) calloc(like.Ntheta, sizeof(double));
    xmax = (double*) calloc(like.Ntheta, sizeof(double));
    const double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    for(int i=0; i<like.Ntheta ; i++)
    {
      xmin[i] = cos(exp(log(like.vtmin) + (i + 0.0)*logdt));
      xmax[i] = cos(exp(log(like.vtmin) + (i + 1.0)*logdt));
    }
    #pragma omp parallel for
    for (int i=0; i<like.Ntheta; i++)
    {
      gsl_sf_legendre_Pl_deriv_array(limits.LMAX, xmin[i], Pmin[i], dPmin[i]);
      gsl_sf_legendre_Pl_deriv_array(limits.LMAX, xmax[i], Pmax[i], dPmax[i]);
    }
  }
  bin_avg r;
  r.xmin = xmin[i_theta];
  r.xmax = xmax[i_theta];
  r.Pmin = Pmin[i_theta][j_L];
  r.Pmax = Pmax[i_theta][j_L];
  r.dPmin = dPmin[i_theta][j_L];
  r.dPmax = dPmax[i_theta][j_L];
  return r;
}

double xi_pm_tomo(int pm, int nt, int ni, int nj, int limber)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }

  static double** Glplus = 0;
  static double** Glminus = 0;
  static double* xi_vec_plus = 0;
  static double* xi_vec_minus = 0;
  static cosmopara C;
  static nuisancepara N;
  
  const int NSIZE = tomo.shear_Npowerspectra;
  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;

  if (Glplus == 0)
  {
    Glplus = create_double_matrix(0, ntheta-1, 0, nell-1);
    Glminus = create_double_matrix(0, ntheta-1, 0, nell-1);
    xi_vec_plus = create_double_vector(0, NSIZE*ntheta-1);
    xi_vec_minus = create_double_vector(0, NSIZE*ntheta-1);

    double xmin[ntheta];
    double xmax[ntheta];

    for (int i=0; i<ntheta; i++)
    { // Cocoa: dont thread (init of static variables inside set_bin_average)
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }
    #pragma omp parallel for
    for (int i=0; i<ntheta; i ++)
    {
      double* Pmin = create_double_vector(0, nell + 1);
      double* Pmax = create_double_vector(0, nell + 1);
      double* dPmin = create_double_vector(0, nell + 1);
      double* dPmax = create_double_vector(0, nell + 1);
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[l] = r.Pmin;
        Pmax[l] = r.Pmax;
        dPmin[l] = r.dPmin;
        dPmax[l] = r.dPmax;
      }
      for (int l=3; l<nell; l++)
      {
        Glplus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(
        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1.)/(2.*l+1)           * (Pmin[l+1]-Pmax[l+1])
        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])
        +2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        -2*(l+2) * (dPmin[l-1]-dPmax[l-1])
        )/(xmin[i]-xmax[i]);

        Glminus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(
        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1.)/(2.*l+1)           * (Pmin[l+1]-Pmax[l+1])
        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])
        -2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        +2*(l+2) * (dPmin[l-1]-dPmax[l-1])
        )/(xmin[i]-xmax[i]);
      }
      free_double_vector(Pmin, 0, nell + 1);
      free_double_vector(Pmax, 0, nell + 1);
      free_double_vector(dPmin, 0, nell + 1);
      free_double_vector(dPmax, 0, nell + 1);
    }
  }

  if (recompute_shear(C, N))
  {
    if(limber == 1)
    {
      if (like.IA == 5 || like.IA == 6)
      { // NEW TATT MODELING
        double** Cl_EE = (double**) malloc(NSIZE*sizeof(double*));
        double** Cl_BB = (double**) malloc(NSIZE*sizeof(double*));
        for (int nz = 0; nz<NSIZE; nz++)
        {
          Cl_EE[nz] = calloc(nell, sizeof(double));
          Cl_BB[nz] = calloc(nell, sizeof(double));
        }

        for (int nz=0; nz<NSIZE; nz++) // no threading because B mode condition
        {
          const int Z1NZ = Z1(nz);
          const int Z2NZ = Z2(nz);
          
          // only compute BB if the TATT parameters allow for B-mode terms
          const int BM = (nuisance.b_ta_z[0] || nuisance.b_ta_z[Z1NZ] ||
              nuisance.b_ta_z[Z2NZ] || nuisance.A2_ia ||
              nuisance.A2_z[Z1NZ] || nuisance.A2_z[Z2NZ]) ? 1 : 0;

          {
            int l = 2;
            Cl_EE[nz][l] = C_ss_tomo_TATT_EE_limber_nointerp((double) l, Z1NZ, Z2NZ);
            Cl_BB[nz][l] = (BM == 1) ? C_ss_tomo_TATT_BB_limber_nointerp((double) l, Z1NZ, Z2NZ) : 0;
            l = 3;
            Cl_EE[nz][l] = C_ss_tomo_TATT_EE_limber_nointerp((double) l, Z1NZ, Z2NZ);
            Cl_BB[nz][l] = (BM == 1) ? C_ss_tomo_TATT_BB_limber_nointerp((double) l, Z1NZ, Z2NZ) : 0;
          }
          #pragma omp parallel for
          for (int l=4; l<limits.LMIN_tab; l++)
          {
            Cl_EE[nz][l] = C_ss_tomo_TATT_EE_limber_nointerp((double) l, Z1NZ, Z2NZ);
            Cl_BB[nz][l] = (BM == 1) ? C_ss_tomo_TATT_BB_limber_nointerp((double) l, Z1NZ, Z2NZ) : 0;
          }
          {
            const int l = limits.LMIN_tab;
            Cl_EE[nz][l] = C_ss_tomo_TATT_EE_limber((double) l,Z1NZ, Z2NZ);
            Cl_BB[nz][l] = (BM == 1) ? C_ss_tomo_TATT_BB_limber((double) l, Z1NZ, Z2NZ) : 0.0;
          }
          #pragma omp parallel for
          for (int l=limits.LMIN_tab+1; l<nell; l++)
          {
            Cl_EE[nz][l] = C_ss_tomo_TATT_EE_limber((double) l,Z1NZ, Z2NZ);
            Cl_BB[nz][l] = (BM == 1) ?  C_ss_tomo_TATT_BB_limber((double) l, Z1NZ, Z2NZ) : 0.0;
          }
        }

        #pragma omp parallel for
        for (int nz=0; nz<NSIZE; nz++)
        {
          for (int i=0; i<ntheta; i++)
          {
            xi_vec_plus[nz * ntheta + i] = 0;
            xi_vec_minus[nz * ntheta + i] = 0;
            for (int l=2; l<nell; l++)
            {
              xi_vec_plus[nz*ntheta+i] +=  Glplus[i][l] * (Cl_EE[nz][l] + Cl_BB[nz][l]);
              xi_vec_minus[nz*ntheta+i] += Glminus[i][l] * (Cl_EE[nz][l] - Cl_BB[nz][l]);
            }
          }
        }

        for (int nz = 0; nz<NSIZE; nz++)
        {
          free(Cl_EE[nz]);
          free(Cl_BB[nz]);
        }
        free(Cl_EE);
        free(Cl_BB);
      }
      else
      {
        double** Cl = malloc(NSIZE*sizeof(double*));
        for (int nz = 0; nz<NSIZE; nz++)
        {
          Cl[nz] = calloc(nell, sizeof(double));
        }

        {
          const int nz=0;
          const int Z1NZ = Z1(nz);
          const int Z2NZ = Z2(nz);
          {
            int l = 2;
            Cl[nz][l] = C_ss_tomo_limber_nointerp((double) l, Z1NZ, Z2NZ, use_linear_ps_limber);  
            l = 3;
            Cl[nz][l] = C_ss_tomo_limber_nointerp((double) l, Z1NZ, Z2NZ, use_linear_ps_limber);
          }
          #pragma omp parallel for
          for (int l=4; l<limits.LMIN_tab; l++) 
          {
            Cl[nz][l] = C_ss_tomo_limber_nointerp((double) l, Z1NZ, Z2NZ, use_linear_ps_limber);
          }
          {
            const int l=limits.LMIN_tab;
            Cl[nz][l] = C_ss_tomo_limber((double) l, Z1NZ, Z2NZ);
          }
          #pragma omp parallel for
          for (int l=limits.LMIN_tab+1; l<nell; l++)
          {
            Cl[nz][l] = C_ss_tomo_limber((double) l, Z1NZ, Z2NZ);
          }
        }
        #pragma omp parallel for
        for (int nz=1; nz<NSIZE; nz++)
        {
          const int Z1NZ = Z1(nz);
          const int Z2NZ = Z2(nz);
          for (int l=2; l<limits.LMIN_tab; l++)
          {
            Cl[nz][l] = C_ss_tomo_limber_nointerp((double) l, Z1NZ, Z2NZ, use_linear_ps_limber);
          }
          for (int l=limits.LMIN_tab; l<nell; l++)
          {
            Cl[nz][l] = C_ss_tomo_limber((double) l, Z1NZ, Z2NZ);
          }
        }

        #pragma omp parallel for
        for (int nz=0; nz<NSIZE; nz++)
        {
          for (int i=0; i<ntheta; i++)
          {
            xi_vec_plus[nz*ntheta+i] = 0;
            xi_vec_minus[nz*ntheta+i] = 0;
            for (int l=2; l<nell; l++)
            {
              xi_vec_plus[nz*ntheta+i] += Glplus[i][l]*Cl[nz][l];
              xi_vec_minus[nz*ntheta+i] += Glminus[i][l]*Cl[nz][l];
            }
          }
        }

        for (int nz=0; nz<NSIZE; nz++)
        {
          free(Cl[nz]);
        }
        free(Cl);
      }
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  const int q = N_shear(ni,nj);
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  if (pm > 0)
  {
    return xi_vec_plus[q*ntheta + nt];
  }
  else
  {
    return xi_vec_minus[q*ntheta + nt];
  }
}

// ---------------------------------------------------------------------------

double w_gammat_tomo(int nt, int ni, int nj, int limber)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }

  static double** Pl = 0;
  static double* w_vec = 0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
    
  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;
  const int NSIZE = tomo.ggl_Npowerspectra;

  if (Pl == 0)
  {
    Pl = create_double_matrix(0, ntheta - 1, 0, nell - 1);
    w_vec = create_double_vector(0,NSIZE*ntheta-1);
    double xmin[ntheta];
    double xmax[ntheta];
    // Cocoa: dont thread (init of static variables inside set_bin_average)
    for (int i=0; i<ntheta; i++)
    {
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }
    #pragma omp parallel for
    for (int i=0; i<like.Ntheta; i++)
    {
      double* Pmin = create_double_vector(0, nell + 1);
      double* Pmax = create_double_vector(0, nell + 1);
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[l] = r.Pmin;
        Pmax[l] = r.Pmax;
      }
      for (int l=1; l<nell; l++)
      {
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
          *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
          +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
          -2./(2*l+1.)*(Pmin[l+1]-Pmax[l+1]));
      }
      free_double_vector(Pmin, 0, nell+1);
      free_double_vector(Pmax, 0, nell+1);
    }
  }

  if (recompute_gs(C, G, N))
  {
    double** Cl = malloc(NSIZE*sizeof(double*));
    for (int nz = 0; nz<NSIZE; nz++)
    {
      Cl[nz] = calloc(nell, sizeof(double));
    }
    if(limber == 1)
    {
      {
        const int nz = 0;
        const int ZLNZ = ZL(nz);
        const int ZSNZ = ZS(nz);
        {
          int l = 2;
          Cl[nz][l] = C_gs_tomo_limber_nointerp((double) l, ZLNZ, ZSNZ, use_linear_ps_limber);
          l = 3;
          Cl[nz][l] = C_gs_tomo_limber_nointerp((double) l, ZLNZ, ZSNZ, use_linear_ps_limber);
        }
        #pragma omp parallel for
        for (int l=4; l<limits.LMIN_tab; l++)
        {
          Cl[nz][l] = C_gs_tomo_limber_nointerp((double) l, ZLNZ, ZSNZ, use_linear_ps_limber);
        }
        {
          const int l=limits.LMIN_tab;
          Cl[nz][l] = C_gs_tomo_limber((double) l, ZLNZ, ZSNZ);
        }
        #pragma omp parallel for
        for (int l=limits.LMIN_tab+1; l<nell; l++)
        {
          Cl[nz][l] = C_gs_tomo_limber((double) l, ZLNZ, ZSNZ);
        }
      }
      #pragma omp parallel for
      for (int nz=1; nz<NSIZE; nz++)
      {
        const int ZLNZ = ZL(nz);
        const int ZSNZ = ZS(nz);
        for (int l=2; l<limits.LMIN_tab; l++)
        {
          Cl[nz][l] = C_gs_tomo_limber_nointerp((double) l, ZLNZ, ZSNZ, use_linear_ps_limber);
        }
        #pragma omp parallel for
        for (int l=limits.LMIN_tab; l<nell; l++)
        {
          Cl[nz][l] = C_gs_tomo_limber((double) l, ZLNZ, ZSNZ);
        }
      }
    }
    else
    {
      const int L = 1;
      const double tolerance = 0.01; //required fractional accuracy in C(l)
      const double dev = 10. * tolerance; // will be diff  exact vs Limber
                                          // init to large value in order to
                                          // start while loop
      // Cocoa: no threading allowed here - (fftw allocation @C_gl_mixed)
      for (int nz=0; nz<NSIZE; nz++)
      {
        C_gl_mixed(L, nell, ZL(nz), ZS(nz), Cl[nz], dev, tolerance);
      }
    }
    #pragma omp parallel for
    for (int nz=0; nz<NSIZE; nz++)
    {
      for (int i=0; i<like.Ntheta; i++)
      {
        w_vec[nz*like.Ntheta+i] = 0;
        for (int l = 1; l<nell; l++)
        {
          w_vec[nz*like.Ntheta+i] += Pl[i][l]*Cl[nz][l];
        }
      }
    }

    for (int nz=0; nz<NSIZE; nz++)
    {
      free(Cl[nz]);
    }
    free(Cl);
    
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }

  const int q = N_ggl(ni, nj);
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  return w_vec[q*like.Ntheta + nt];
}

// ---------------------------------------------------------------------------

double w_gg_tomo(int nt, int ni, int nj, int limber)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static double **Pl = 0;
  static double *w_vec = 0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;

  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  
  const int nell = limits.LMAX;

  if (Pl == 0)
  {
    Pl = create_double_matrix(0, like.Ntheta - 1, 0, nell - 1);
    w_vec = create_double_vector(0, tomo.clustering_Nbin * like.Ntheta - 1);
    double xmin[like.Ntheta];
    double xmax[like.Ntheta];
    // Cocoa: dont thread (init of static variables inside set_bin_average)
    for (int i=0; i<like.Ntheta; i ++)
    {
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }
    #pragma omp parallel for
    for (int i=0; i<like.Ntheta; i ++)
    {
      double* Pmin = create_double_vector(0, nell + 1);
      double* Pmax = create_double_vector(0, nell + 1);
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[l] = r.Pmin;
        Pmax[l] = r.Pmax;
      }
      Pl[i][0] = 1.0;
      const double a = (1.0/(xmin[i] - xmax[i]))*(1. / (4. * M_PI));
      for (int l=1; l<nell; l++)
      {
        Pl[i][l] = a*(Pmin[l + 1] - Pmax[l + 1] - Pmin[l - 1] + Pmax[l - 1]);
      }
      free_double_vector(Pmin, 0, nell + 1);
      free_double_vector(Pmax, 0, nell + 1);
    }
  }

  if (recompute_gg(C, G, N))
  {
    double** Cl = malloc(tomo.clustering_Nbin*sizeof(double*));
    for (int nz=0; nz<tomo.clustering_Nbin; nz++)
    {
      Cl[nz] = calloc(nell, sizeof(double));
    }
    if(limber == 1)
    {
      {
        const int nz = 0;
        {
          int l = 1;
          Cl[nz][l] = C_gg_tomo_limber_nointerp((double) l, nz, nz, use_linear_ps_limber);
          l = 2;
          Cl[nz][l] = C_gg_tomo_limber_nointerp((double) l, nz, nz, use_linear_ps_limber);
        }
        #pragma omp parallel for
        for (int l=3; l<limits.LMIN_tab; l++)
        {
          Cl[nz][l] = C_gg_tomo_limber_nointerp((double) l, nz, nz, use_linear_ps_limber);
        }
        {
          const int l = limits.LMIN_tab;
          Cl[nz][l] = C_gg_tomo_limber(1.0*l, nz, nz);
        }
        #pragma omp parallel for
        for (int l=limits.LMIN_tab+1; l<nell; l++)
        {
          Cl[nz][l] = C_gg_tomo_limber((double) l, nz, nz);
        }
      }
      #pragma omp parallel for
      for (int nz=1; nz<tomo.clustering_Nbin; nz++)
      {
        for (int l=1; l<limits.LMIN_tab; l++)
        {
          Cl[nz][l] = C_gg_tomo_limber_nointerp((double) l, nz, nz, use_linear_ps_limber);
        }
        for (int l=limits.LMIN_tab; l<nell; l++)
        {
          Cl[nz][l] = C_gg_tomo_limber((double) l, nz, nz);
        }
      }
    }
    else
    {
      const int L = 0;
      const double tolerance = 0.01; // required fractional accuracy in C(l)
      const double dev = 10. * tolerance; // will be diff  exact vs Limber
                                          // init to large value in order to
                                          // start while loop
      
      // Cocoa: no threading allowed here - (fftw allocation @C_cl_mixed)
      for (int nz=0; nz<tomo.clustering_Nbin; nz++)
      {
        C_cl_mixed(L, nell, nz, nz, Cl[nz], dev, tolerance);
      }
    }
    #pragma omp parallel for
    for (int nz=0; nz<tomo.clustering_Nbin; nz++)
    {
      for (int i=0; i<like.Ntheta; i++)
      {
        w_vec[nz*like.Ntheta + i] = 0;
        for (int l=1; l<nell; l++)
        {
          w_vec[nz*like.Ntheta + i] += Pl[i][l]*Cl[nz][l];
        }
      }
    }
    for (int nz=0; nz<tomo.clustering_Nbin; nz++)
    {
      free(Cl[nz]);
    }
    free(Cl);
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  
  if (ni != nj)
  {
    log_fatal("ni != nj tomography not supported");
    exit(1);
  }
  
  return w_vec[ni * like.Ntheta + nt];
}

// ---------------------------------------------------------------------------

double w_gk_tomo(int nt, int ni, int limber)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  
  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }

  static double **Pl =0;
  static double *w_vec =0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  const int nell = limits.LMAX;

  if (Pl == 0)
  {
    Pl = create_double_matrix(0, like.Ntheta-1, 0, nell-1);
    w_vec = create_double_vector(0,tomo.clustering_Nbin*like.Ntheta-1);
    double xmin[like.Ntheta];
    double xmax[like.Ntheta];
    // Cocoa: dont thread (init of static variables inside set_bin_average)
    for (int i=0; i<like.Ntheta; i++)
    {
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }
    #pragma omp parallel for
    for (int i=0; i<like.Ntheta; i++)
    {
      double* Pmin = create_double_vector(0, nell + 1);
      double* Pmax = create_double_vector(0, nell + 1);
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[l] = r.Pmin;
        Pmax[l] = r.Pmax;
      }
      for (int l=1; l<nell; l++)
      {
        Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/
          (xmin[i]-xmax[i]);
      }
      free_double_vector(Pmin, 0, nell+1);
      free_double_vector(Pmax, 0, nell+1);
    }
  }
  
  if (recompute_gk(C, G, N))
  {
    if(limber == 1)
    {
      double** Cl = malloc(tomo.clustering_Nbin*sizeof(double*));
      for (int nz = 0; nz<tomo.clustering_Nbin; nz++)
      {
        Cl[nz] = calloc(nell, sizeof(double));
      }

      {
        const int nz = 0;
        {
          int l = 1;
          Cl[nz][l] = C_gk_tomo_limber_nointerp_wrapper((double) l, nz, use_linear_ps_limber);
          l = 2;
          Cl[nz][l] = C_gk_tomo_limber_nointerp_wrapper((double) l, nz, use_linear_ps_limber);
        }
        #pragma omp parallel for
        for (int l=3; l<limits.LMIN_tab; l++)
        {
          Cl[nz][l] = C_gk_tomo_limber_nointerp_wrapper((double) l, nz, use_linear_ps_limber);
        }
        {
          int l = limits.LMIN_tab;
          Cl[nz][l] = C_gk_tomo_limber_wrapper((double) l, nz);
        }
        #pragma omp parallel for
        for (int l=limits.LMIN_tab+1; l<nell; l++)
        {
          Cl[nz][l] = C_gk_tomo_limber_wrapper((double) l, nz);
        }
      }
      #pragma omp parallel for
      for (int nz=1; nz<tomo.clustering_Nbin; nz++)
      {
        for (int l=1; l<limits.LMIN_tab; l++)
        {
          Cl[nz][l] = C_gk_tomo_limber_nointerp_wrapper((double) l, nz, use_linear_ps_limber);
        }
        for (int l=limits.LMIN_tab; l<nell; l++)
        {
          Cl[nz][l] = C_gk_tomo_limber_wrapper((double) l, nz);
        }
      }
      #pragma omp parallel for
      for (int nz=0; nz<tomo.clustering_Nbin; nz++)
      {
        for (int i=0; i<like.Ntheta; i++)
        {
          w_vec[nz*like.Ntheta+i] = 0;
          for (int l=1; l<nell; l++)
          {
            w_vec[nz*like.Ntheta+i] += Pl[i][l]*Cl[nz][l];
          }
        }
      }
      for (int nz=0; nz<tomo.clustering_Nbin; nz++)
      {
        free(Cl[nz]);
      }
      free(Cl);
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  
  return w_vec[ni*like.Ntheta+nt];
}

// ---------------------------------------------------------------------------

double w_ks_tomo(int nt, int ni, int limber)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  if (like.Ntheta ==0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }

  static double** Pl =0;
  static double* w_vec =0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  const int nell = limits.LMAX;

  if (Pl == 0)
  {
    Pl = create_double_matrix(0, like.Ntheta - 1, 0, nell - 1);
    w_vec = create_double_vector(0, tomo.shear_Nbin*like.Ntheta - 1);
    double xmin[like.Ntheta];
    double xmax[like.Ntheta];
    // Cocoa: dont thread (init of static variables inside set_bin_average)
    for (int i=0; i<like.Ntheta; i++)
    {
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }
    #pragma omp parallel for
    for (int i=0; i<like.Ntheta; i++)
    {
      double* Pmin = create_double_vector(0, nell + 1);
      double* Pmax = create_double_vector(0, nell + 1);
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[l] = r.Pmin;
        Pmax[l] = r.Pmax;
      }
      for (int l=1; l<nell; l++)
      {
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
          *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
          +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
          -2./(2*l+1.)*(Pmin[l+1]-Pmax[l+1]));
      }
      free_double_vector(Pmin, 0, nell+1);
      free_double_vector(Pmax, 0, nell+1);
    }
  }
  
  if (recompute_ks(C, G, N))
  {
    if(limber == 1)
    {
      double** Cl = malloc(tomo.shear_Nbin*sizeof(double*));
      for (int nz=0; nz<tomo.shear_Nbin; nz++)
      {
        Cl[nz] = calloc(nell, sizeof(double));
      }

      {
        const int nz = 0;
        {
          int l = 1;
          Cl[nz][l] = C_ks_tomo_limber_nointerp_wrapper(l, nz, use_linear_ps_limber);
          l = 2;
          Cl[nz][l] = C_ks_tomo_limber_nointerp_wrapper(l, nz, use_linear_ps_limber);
        }
        #pragma omp parallel for
        for (int l=3; l<limits.LMIN_tab; l++)
        {
          Cl[nz][l] = C_ks_tomo_limber_nointerp_wrapper(l, nz, use_linear_ps_limber);
        }
        {
          int l = limits.LMIN_tab;
          Cl[nz][l] = C_ks_tomo_limber_wrapper(l, nz);
        }
        #pragma omp parallel for
        for (int l=limits.LMIN_tab+1; l<nell; l++)
        {
          Cl[nz][l] = C_ks_tomo_limber_wrapper(l, nz);
        }
      }
      #pragma omp parallel for
      for (int nz=1; nz<tomo.shear_Nbin; nz++)
      {
        for (int l=1; l<limits.LMIN_tab; l++)
        {
          Cl[nz][l] = C_ks_tomo_limber_nointerp_wrapper(l, nz, use_linear_ps_limber);
        }
        for (int l=limits.LMIN_tab; l<nell; l++)
        {
          Cl[nz][l] = C_ks_tomo_limber_wrapper(l, nz);
        }
      }
      #pragma omp parallel for
      for (int nz=0; nz<tomo.shear_Nbin; nz++)
      {
        for (int i=0; i<like.Ntheta; i++)
        {
          w_vec[nz*like.Ntheta+i] = 0;
          for (int l=2; l<nell; l++)
          {
            w_vec[nz*like.Ntheta+i] += Pl[i][l]*Cl[nz][l];
          }
        }
      }

      for (int nz=0; nz<tomo.shear_Nbin; nz++)
      {
        free(Cl[nz]);
      }
      free(Cl);
    } else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return w_vec[ni*like.Ntheta + nt];
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real space) - flat sky
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double xi_pm_tomo_flatsky(int pm, double theta, int ni, int nj, int limber)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int NSIZE = tomo.shear_Npowerspectra;
  const int ntheta = Ntable.N_thetaH;
  const double l_min = w_l_min;
  const double l_max = w_l_max;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax - lnlmin)/(1.0*ntheta);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2 + 1;
  
  const double lnthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double) ntheta);
  const double lntheta = log(theta);

  if (table == 0)
  {
    table = create_double_matrix(0, 2*NSIZE - 1, 0, ntheta - 1); // 2 NSIZE = {xi+, xi-}
  }
  if (recompute_shear(C, N))
  {
    if (limber == 1)
    {
      if (like.IA == 5 || like.IA == 6)
      { // NEW TATT MODELING
        log_fatal("Limber && TATT not implemented");
        exit(1);
      }
      else
      {
        typedef fftw_complex fftwZ;

        fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
        for (int j=0; j<NSIZE; j++)
        {
          flP[j] = (fftwZ*) fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftwZ));
        }
        {
          double** lP = (double**) malloc(sizeof(double*)*NSIZE);
          fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
          for (int j=0; j<NSIZE; j++)
          {
            lP[j] = (double*) malloc(ntheta*sizeof(double));
            plan[j] = fftw_plan_dft_r2c_1d(ntheta,lP[j],flP[j],FFTW_ESTIMATE);
          }
          
          // Power spectrum on logarithmic bins          
          {
            const int k = 0;
            const int Z1NZ = Z1(k);
            const int Z2NZ = Z2(k);
            {
              const int p = 0;
              const double l = exp(lnrc + (p - nc)*dlnl);
              lP[k][p] = l*C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, use_linear_ps_limber);
              lP[k][limits.LMIN_tab] = C_ss_tomo_limber(limits.LMIN_tab, Z1NZ, Z2NZ); 
            }
            #pragma omp parallel for
            for(int p=1; p<ntheta; p++)
            {
              const double l = exp(lnrc + (p - nc)*dlnl);
              if(l > limits.LMIN_tab)
              {
                lP[k][p] = l*C_ss_tomo_limber(l, Z1NZ, Z2NZ);
              } 
              else if (l == limits.LMIN_tab)
              {
                // do nothing
              }
              else
              {
                lP[k][p] = l*C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, use_linear_ps_limber);
              }
            }
          }
          #pragma omp parallel for
          for (int k=1; k<NSIZE; k++)
          {
            const int Z1NZ = Z1(k);
            const int Z2NZ = Z2(k);
            for(int p=0; p<ntheta; p++)
            {
              const double l = exp(lnrc+(p - nc)*dlnl);
              if(l > limits.LMIN_tab)
              {
                lP[k][p] = l*C_ss_tomo_limber(l, Z1NZ, Z2NZ);
              } 
              else
              {
                lP[k][p] = l*C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, use_linear_ps_limber);
              }
            }
          }

          #pragma omp parallel for
          for (int j=0; j<NSIZE; j++)
          {
            fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)
          }
          
          for (int j = 0; j < NSIZE; j++)
          {
            fftw_free(lP[j]);
            fftw_destroy_plan(plan[j]);
          }
          free(lP);
          free(plan);
        }
        
        double*** lP = (double***) malloc(sizeof(double**)*NSIZE);
        fftwZ*** kernel = (fftwZ***) malloc(sizeof(fftwZ**)*NSIZE);
        fftwZ*** conv = (fftwZ***) malloc(sizeof(fftwZ**)*NSIZE);
        fftw_plan** plan = (fftw_plan**) malloc(sizeof(fftw_plan*)*NSIZE);
        double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
        for (int j=0; j<NSIZE; j++)
        {
          lP[j] = (double**) malloc(sizeof(double*)*2);
          kernel[j] = (fftwZ**) malloc(sizeof(fftwZ*)*2);
          conv[j] = (fftwZ**) malloc(sizeof(fftwZ*)*2);
          plan[j] = (fftw_plan*) malloc(sizeof(fftw_plan)*2);
          tab[j] = create_double_matrix(0, 1, 0, ntheta-1);
          for (int m=0; m<2; m++)
          {
            const int COVSZ = (ntheta/2+1);
            lP[j][m] = (double*) malloc(ntheta*sizeof(double));
            kernel[j][m] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
            conv[j][m] = (fftwZ*) fftw_malloc(COVSZ*sizeof(fftwZ));
            plan[j][m] = fftw_plan_dft_c2r_1d(ntheta,conv[j][m],lP[j][m],FFTW_ESTIMATE);
          }
        }
        
        #pragma omp parallel for
        for (int j=0; j<NSIZE; j++)
        {
          for (int m=0; m<2; m++)
          {
            double arg[2];
            arg[0] = 0; // bias
            arg[1] = (m == 0 ? 0 : 4); // order of Bessel function

            // perform the convolution, negative sign for kernel (complex conj.)
            for(int i=0; i<(ntheta/2+1); i++)
            {
              const double k = 2*constants.pi*i/(dlnl*ntheta);
              hankel_kernel_FT(k, kernel[j][m], arg, 2);
              conv[j][m][i][0] = flP[j][i][0]*(kernel[j][m][0][0]) - flP[j][i][1]*(kernel[j][m][0][1]);
              conv[j][m][i][1] = flP[j][i][1]*(kernel[j][m][0][0]) + flP[j][i][0]*(kernel[j][m][0][1]);
            }

            // force Nyquist- and 0-frequency-components to be double
            conv[j][m][0][1] = 0;
            conv[j][m][ntheta/2][1] = 0;

            fftw_execute(plan[j][m]);

            for(int k=0; k<ntheta; k++)
            {
              const double t = exp((nc-k)*dlnl-lnrc); // theta=1/l
              tab[j][m][ntheta-k-1] = lP[j][m][k]/(t*2*constants.pi*ntheta);
            }
          }
          for (int k=0; k<ntheta; k++)
          {
            table[2*j][k] = tab[j][0][k];
            table[2*j+1][k] = tab[j][1][k];
          }
        }
        
        for (int j=0; j<NSIZE; j++)
        {
          for (int m=0; m<2; m++)
          {
            fftw_free(lP[j][m]);
            fftw_free(kernel[j][m]);
            fftw_free(conv[j][m]);
            fftw_destroy_plan(plan[j][m]);
          }
          fftw_free(flP[j]);
          free(lP[j]);
          free(kernel[j]);
          free(conv[j]);
          free(plan[j]);
          free_double_matrix(tab[j], 0, 1, 0, ntheta);
        }
        free(flP);
        free(lP);
        free(kernel);
        free(conv);
        free(plan);
        free(tab);
        // --------------------------------------------------------------------
        // Cocoa: code extracted (& adapted) from xipm_via_hankel (ends)
        // --------------------------------------------------------------------
      }
      update_cosmopara(&C);
      update_nuisance(&N);
    } else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
  }
  const int q = 2*N_shear(ni, nj) + (1 - pm)/2;
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }

  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 0, 0);
}

// ---------------------------------------------------------------------------

double w_gammat_tomo_flatsky(double theta, int ni, int nj, int limber)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int NSIZE = tomo.ggl_Npowerspectra;
  const int ntheta = Ntable.N_thetaH;

  const double l_min = w_l_min;
  const double l_max = w_l_max;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax-lnlmin)/(ntheta - 1.);
  const double lnrc = 0.5*(lnlmax+lnlmin);
  const double nc = ntheta/2 + 1;
  
  const double lnthetamin = (nc-ntheta + 1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double)ntheta);
  const double lntheta = log(theta);

  if (table == 0)
  {
    table = create_double_matrix(0, NSIZE - 1, 0, ntheta - 1);
  }
  if (recompute_gs(C, G, N)) 
  {
    if (limber == 1)
    {
      typedef fftw_complex fftwZ;

      fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE); // go to log-Fourier-space
      for (int j=0; j<NSIZE; j++)
      {
        flP[j] = fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
      }
      {
        double** lP = (double**) malloc(sizeof(double*)*NSIZE);
        fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
        for (int j=0; j<NSIZE; j++)
        {
          lP[j] = (double*) malloc(ntheta*sizeof(double));
          plan[j] = fftw_plan_dft_r2c_1d(ntheta, lP[j], flP[j], FFTW_ESTIMATE);
        }

        { // Power spectrum on logarithmic bins (begins)
          {
            const int k = 0;
            const int ZLNZ = ZL(k);
            const int ZSNZ = ZS(k);
            {
              const int p = 0;
              const double l = exp(lnrc + (p - nc)*dlnl);
              lP[k][p] = l*C_gs_tomo_limber_nointerp(l, ZLNZ, ZSNZ, use_linear_ps_limber);
              lP[k][limits.LMIN_tab] = C_gs_tomo_limber(limits.LMIN_tab, ZLNZ, ZSNZ);
            }
            #pragma omp parallel for
            for(int p=1; p<ntheta; p++)
            {
              const double l = exp(lnrc + (p - nc)*dlnl);
              if(l > limits.LMIN_tab)
              {
                lP[k][p] = l*C_gs_tomo_limber(l, ZLNZ, ZSNZ);
              }
              else if (l == limits.LMIN_tab)
              {
                // do nothing
              }
              else
              {
                lP[k][p] = l*C_gs_tomo_limber_nointerp(l, ZLNZ, ZSNZ, use_linear_ps_limber);
              }
            }
          }
          #pragma omp parallel for
          for (int k=1; k<NSIZE; k++)
          {
            const int ZLNZ = ZL(k);
            const int ZSNZ = ZS(k);
            for(int p=0; p<ntheta; p++)
            {
              const double l = exp(lnrc + (p - nc)*dlnl);
              if(l > limits.LMIN_tab)
              {
                lP[k][p] = l*C_gs_tomo_limber(l, ZLNZ, ZSNZ);
              }
              else
              {
                lP[k][p] = l*C_gs_tomo_limber_nointerp(l, ZLNZ, ZSNZ, use_linear_ps_limber);
              }
            }
          }
        } // Power spectrum on logarithmic bins (ends)
        
        #pragma omp parallel for
        for (int j=0; j<NSIZE; j++)
        {
          fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)
        }
        
        for (int j=0; j<NSIZE; j++)
        {
          fftw_free(lP[j]);
          fftw_destroy_plan(plan[j]);
        }
        
        free(lP);
        free(plan);
      }

      double** lP = (double**) malloc(sizeof(double*)*NSIZE);
      fftwZ** kernel = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
      fftwZ** conv = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
      for (int j=0; j<NSIZE; j++)
      {
        const int COVSZ = (ntheta/2+1);
        lP[j] = (double*) malloc(ntheta*sizeof(double));
        kernel[j] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
        conv[j] = (fftwZ*) fftw_malloc(COVSZ*sizeof(fftwZ));
        plan[j] = fftw_plan_dft_c2r_1d(ntheta, conv[j], lP[j], FFTW_ESTIMATE);
        tab[j] = create_double_matrix(0, 1, 0, ntheta-1);
      }

      #pragma omp parallel for
      for (int j=0; j<NSIZE; j++)
      {
        double arg[2];
        arg[0] = 0; // bias
        arg[1] = 2; // order of Bessel function

        for(int i=0; i<(ntheta/2+1); i++)
        {
          const double kk = 2*constants.pi*i/(dlnl*ntheta);
          hankel_kernel_FT(kk, kernel[j], arg, 2);
          conv[j][i][0] = flP[j][i][0]*kernel[j][0][0] - flP[j][i][1]*kernel[j][0][1];
          conv[j][i][1] = flP[j][i][1]*kernel[j][0][0] + flP[j][i][0]*kernel[j][0][1];
        }

        // force Nyquist- and 0-frequency-components to be double
        conv[j][0][1] = 0;
        conv[j][ntheta/2][1] = 0;

        fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)

        for(int k=0; k<ntheta; k++)
        {
          const double t = exp((nc-k)*dlnl-lnrc); // theta=1/l
          tab[j][0][ntheta-k-1] = lP[j][k]/(t*2*constants.pi*ntheta);
        }
        
        for (int k=0; k<ntheta; k++)
        {
          table[j][k] = tab[j][0][k];
        }
      }

      for (int j=0; j<NSIZE; j++)
      {
        fftw_free(flP[j]);
        fftw_free(lP[j]);
        fftw_free(conv[j]);
        fftw_free(kernel[j]);
        fftw_destroy_plan(plan[j]);
        free_double_matrix(tab[j], 0, 1, 0, ntheta);
      }
      free(flP);
      free(lP);
      free(conv);
      free(kernel);
      free(plan);
      free(tab);
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }

  if(test_zoverlap(ni, nj))
  {
    const int q = N_ggl(ni, nj);
    if(q > NSIZE - 1)
    {
      log_fatal("error in selecting bin number");
      exit(1);
    }

    if (lntheta < lnthetamin || lntheta > lnthetamax)
    {
      const double theta = exp(lntheta);
      const double theta_min = exp(lnthetamin);
      const double theta_max = exp(lnthetamax);
      log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
      exit(1);
    }

    return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 0, 0);
  } 
  else 
  {
    return 0.0;
  }
}

// ---------------------------------------------------------------------------

double w_gg_tomo_flatsky(double theta, int ni, int nj, int limber)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int NSIZE = tomo.clustering_Nbin; // cross redshift bin not supported so not using
                                          // tomo.cc_clustering_Npowerspectra
  const int ntheta = Ntable.N_thetaH;

  const double lmin = w_l_min;
  const double lmax = w_l_max;
  const double lnlmax = log(lmax);
  const double lnlmin = log(lmin);
  const double dlnl = (lnlmax - lnlmin)/(1.0*ntheta - 1.);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2 + 1;

  const double lntheta = log(theta);
  const double lnthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax-lnthetamin)/((double) ntheta);

  if (table == 0)
  {
    table = create_double_matrix(0, NSIZE - 1, 0, ntheta - 1);
  }  
  if (recompute_gg(C, G, N)) 
  {
    if (limber == 1)
    {
      typedef fftw_complex fftwZ;

      fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE); // go to log-Fourier-space
      for (int j=0; j<NSIZE; j++)
      {
        flP[j] = fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
      }
      {
        double** lP = (double**) malloc(sizeof(double*)*NSIZE);
        fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
        for (int j=0; j<NSIZE; j++)
        {
          int ARRAYSZ = ntheta;
          lP[j] = (double*) malloc(ARRAYSZ*sizeof(double));
          plan[j] = fftw_plan_dft_r2c_1d(ARRAYSZ, lP[j], flP[j], FFTW_ESTIMATE);
        }

        { // Power spectrum on logarithmic bins (begins)
          {
            const int k = 0;
            const int ZCL1 = k; // cross redshift bin not supported so not using ZCL1(k)
            const int ZCL2 = k; // cross redshift bin not supported so not using ZCL2(k) 
            {
              const int p = 0;
              const double l = exp(lnrc + (p - nc)*dlnl);
              if(l > limits.LMIN_tab)
              {
                lP[k][p] = l*C_gg_tomo_limber(l, ZCL1, ZCL2);
              } 
              else
              {
                lP[k][p] = l*C_gg_tomo_limber_nointerp(l, ZCL1, ZCL2, use_linear_ps_limber);
                lP[k][limits.LMIN_tab] = C_gg_tomo_limber(limits.LMIN_tab, ZCL1, ZCL2);
              }
            }
            #pragma omp parallel for
            for(int p=1; p<ntheta; p++)
            {
              const double l = exp(lnrc + (p - nc)*dlnl);
              if(l > limits.LMIN_tab)
              {
                lP[k][p] = l*C_gg_tomo_limber(l, ZCL1, ZCL2);
              } 
              else if (l == limits.LMIN_tab)
              {
                // do nothing
              }
              else
              {
                lP[k][p] = l*C_gg_tomo_limber_nointerp(l, ZCL1, ZCL2, use_linear_ps_limber);
              }
            }
          }
          #pragma omp parallel for
          for (int k=1; k<NSIZE; k++) 
          {
            const int ZCL1 = k; // cross redshift bin not supported so not using ZCL1(k)
            const int ZCL2 = k; // cross redshift bin not supported so not using ZCL2(k) 
            for(int p=0; p<ntheta; p++)
            {
              const double l = exp(lnrc + (p - nc)*dlnl);
              if(l > limits.LMIN_tab)
              {
                lP[k][p] = l*C_gg_tomo_limber(l, ZCL1, ZCL2);
              } 
              else
              {
                lP[k][p] = l*C_gg_tomo_limber_nointerp(l, ZCL1, ZCL2, use_linear_ps_limber);
              }
            }
          }
        } // Power spectrum on logarithmic bins (ends)

        #pragma omp parallel for
        for (int j=0; j<NSIZE; j++)
        {
          fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)
        }

        for (int j=0; j<NSIZE; j++)
        {
          fftw_free(lP[j]);
          fftw_destroy_plan(plan[j]);
        }
        free(lP);
        free(plan);
      }

      double** lP = (double**) malloc(sizeof(double*)*NSIZE);
      fftwZ** kernel = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
      fftwZ** conv = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
      for (int j=0; j<NSIZE; j++)
      {
        const int ARRAYSZ = ntheta;
        const int COVSZ = (ntheta/2+1);
        lP[j] = (double*) malloc(ARRAYSZ*sizeof(double));
        kernel[j] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
        conv[j] = (fftwZ*) fftw_malloc(COVSZ*sizeof(fftwZ));
        plan[j] = fftw_plan_dft_c2r_1d(ARRAYSZ, conv[j], lP[j], FFTW_ESTIMATE);
        tab[j] = create_double_matrix(0, 1, 0, ntheta-1);
      }

      #pragma omp parallel for
      for (int j=0; j<NSIZE; j++)
      {
        double arg[2];
        arg[0] = 0; // bias
        arg[1] = 0; // order of Bessel function

        for(int i=0; i<(ntheta/2+1); i++)
        {
          const double kk = 2*constants.pi*i/(dlnl*ntheta);
          hankel_kernel_FT(kk, kernel[j], arg, 2);
          conv[j][i][0] = flP[j][i][0]*kernel[j][0][0] - flP[j][i][1]*kernel[j][0][1];
          conv[j][i][1] = flP[j][i][1]*kernel[j][0][0] + flP[j][i][0]*kernel[j][0][1];
        }

        // force Nyquist- and 0-frequency-components to be double
        conv[j][0][1] = 0;
        conv[j][ntheta/2][1] = 0;

        fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)

        for(int k=0; k<ntheta; k++)
        {
          const double t = exp((nc-k)*dlnl-lnrc); // theta=1/l
          tab[j][0][ntheta-k-1] = lP[j][k]/(t*2*constants.pi*ntheta);
        }
        for (int k=0; k<ntheta; k++)
        {
          table[j][k] = tab[j][0][k];
        }
      }
      
      for (int j=0; j<NSIZE; j++)
      {
        fftw_free(flP[j]);
        fftw_free(lP[j]);
        fftw_free(conv[j]);
        fftw_free(kernel[j]);
        fftw_destroy_plan(plan[j]);
        free_double_matrix(tab[j], 0, 1, 0, ntheta);
      }
      free(flP);
      free(lP);
      free(conv);
      free(kernel);
      free(plan);
      free(tab);
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  
  if (ni != nj) 
  {
    log_fatal("cross-tomography not supported");
    exit(1);
  }
    
  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }

  const int q = ni; // cross redshift bin not supported so not using N_CL(ni, nj) instead of ni
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 1, 1);
}

// ---------------------------------------------------------------------------

double w_gk_tomo_flatsky(double theta, int ni, int limber)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int nisize = tomo.clustering_Nbin;
  const int ntheta = Ntable.N_thetaH;
  const int NSIZE = nisize;

  const double lmin = w_l_min;
  const double lmax = w_l_max;
  const double lnlmax = log(lmax);
  const double lnlmin = log(lmin);
  const double dlnl = (lnlmax - lnlmin)/(1.0*ntheta - 1.);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2 + 1;
  
  const double lntheta = log(theta);
  const double lnthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax-lnthetamin)/((double)ntheta);

  if (table==0)
  {
    table = create_double_matrix(0, NSIZE - 1, 0, ntheta - 1);
  }
  if (recompute_gk(C, G, N))
  {
    if (limber == 1)
    {
      typedef fftw_complex fftwZ;

      // go to log-Fourier-space
      fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      for (int i=0; i<NSIZE; i++)
      {
        flP[i] = fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
      }
      {
        double** lP = (double**) malloc(sizeof(double*)*NSIZE);
        fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
        for (int i=0; i<NSIZE; i++)
        {
          lP[i] = (double*) malloc(ntheta*sizeof(double));
          plan[i] = fftw_plan_dft_r2c_1d(ntheta, lP[i], flP[i], FFTW_ESTIMATE);
        }

        { // Power spectrum on logarithmic bins (begins)
          {
            const int i = 0;
            {
              const int k = 0;
              const double l = exp(lnrc + (k - nc)*dlnl);
              lP[i][k] = l*C_gk_tomo_limber_nointerp_wrapper(l, i, use_linear_ps_limber);
              lP[i][limits.LMIN_tab] = C_gk_tomo_limber_wrapper(limits.LMIN_tab, i);
            }
            #pragma omp parallel for
            for(int k=1; k<ntheta; k++)
            {
              const double l = exp(lnrc + (k - nc)*dlnl);
              if(l > limits.LMIN_tab) 
              {
                lP[i][k] = l*C_gk_tomo_limber_wrapper(l, i);
              } 
              else if (l == limits.LMIN_tab)
              {
                // do nothing
              }
              else 
              {
                lP[i][k] = l*C_gk_tomo_limber_nointerp_wrapper(l, i, use_linear_ps_limber);
              }
            }
          }
          #pragma omp parallel for
          for (int i=1; i<NSIZE; i++)
          {
            for(int k=0; k<ntheta; k++)
            {
              const double l = exp(lnrc + (k - nc)*dlnl);
              if(l > limits.LMIN_tab) 
              {
                lP[i][k] = l*C_gk_tomo_limber_wrapper(l, i);
              } 
              else 
              {
                lP[i][k] = l*C_gk_tomo_limber_nointerp_wrapper(l, i, use_linear_ps_limber);
              }
            }
          }
        } // Power spectrum on logarithmic bins (ends)
        
        #pragma omp parallel for
        for (int i=0; i<NSIZE; i++)
        {
          fftw_execute(plan[i]); // Execute FFTW in parallel (thread-safe)
        }
        
        for (int i=0; i<NSIZE; i++)
        {
          fftw_free(lP[i]);
          fftw_destroy_plan(plan[i]);
        } 
        free(lP);
        free(plan);
      }

      double** lP = (double**) fftw_malloc(sizeof(double**)*NSIZE);
      fftwZ** kernel = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      fftwZ** conv = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
      for (int i=0; i<NSIZE; i++)
      {
        lP[i] = (double*) fftw_malloc(ntheta*sizeof(double));
        kernel[i] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
        conv[i] = (fftwZ*) fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
        plan[i] = fftw_plan_dft_c2r_1d(ntheta, conv[i], lP[i], FFTW_ESTIMATE);
        tab[i] = create_double_matrix(0, 1, 0, ntheta-1);
      }
      
      #pragma omp parallel for
      for (int i=0; i<NSIZE; i++)
      {
        double arg[2];
        arg[0] = 0; // bias
        arg[1] = 0; // order of Bessel function

        for(int k=0; k<(ntheta/2+1); k++)
        {
          const double kk = 2*constants.pi*i/(dlnl*ntheta);
          hankel_kernel_FT(kk, kernel[i], arg, 2);
          
          conv[i][k][0] = flP[i][k][0]*kernel[i][0][0] - flP[i][k][1]*kernel[i][0][1];
          conv[i][k][1] = flP[i][k][1]*kernel[i][0][0] + flP[i][k][0]*kernel[i][0][1];
        }

        // force Nyquist- and 0-frequency-components to be double
        conv[i][0][1] = 0;
        conv[i][ntheta/2][1] = 0;

        fftw_execute(plan[i]); // Execute FFTW in parallel (thread-safe)

        for(int k=0; k<ntheta; k++)
        {
          const double t = exp((nc-k)*dlnl-lnrc); // theta=1/l
          tab[i][0][ntheta-k-1] = lP[i][k]/(t*2*constants.pi*ntheta);
        }
        for (int k=0; k<ntheta; k++)
        {
          table[i][k] = tab[i][0][k];
        }
      }
      
      for (int i=0; i<NSIZE; i++)
      {
        fftw_free(flP[i]);
        fftw_free(lP[i]);
        fftw_free(kernel[i]);
        fftw_free(conv[i]);
        fftw_destroy_plan(plan[i]);
        free_double_matrix(tab[i], 0, 1, 0, ntheta);
      }
      free(flP);
      free(lP);
      free(kernel);
      free(conv);
      free(plan);
      free(tab);
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }

  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }

  return interpol(table[ni], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 1, 1);
}

// ---------------------------------------------------------------------------

double w_ks_tomo_flatsky(double theta, int ni, int limber)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int nisize = tomo.shear_Nbin;
  const int ntheta = Ntable.N_thetaH;
  const int NSIZE = nisize;

  const double l_min = w_l_min;
  const double l_max = w_l_max;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax-lnlmin)/(ntheta - 1.);
  const double lnrc = 0.5*(lnlmax+lnlmin);
  const double nc = ntheta/2 + 1;
  
  const double lnthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax-lnthetamin)/((double) ntheta);
  const double lntheta = log(theta);

  if (table == 0) 
  {
    table = create_double_matrix(0, NSIZE - 1, 0, ntheta - 1);
  }
  if (recompute_shear(C, N))
  {
    if (limber == 1)
    {
      typedef fftw_complex fftwZ;

      // go to log-Fourier-space
      fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      for (int i=0; i<NSIZE; i++)
      {
        flP[i] = fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
      }
      {
        double** lP = (double**) malloc(sizeof(double*)*NSIZE);
        fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
        for (int i=0; i<NSIZE; i++)
        {
          lP[i] = (double*) malloc(ntheta*sizeof(double));
          plan[i] = fftw_plan_dft_r2c_1d(ntheta, lP[i], flP[i], FFTW_ESTIMATE);
        }
        { // Power spectrum on logarithmic bins (begins)
          {
            const int i = 0;
            {
              const int k = 0;
              const double l = exp(lnrc + (k - nc)*dlnl);
              lP[i][k] = l*C_ks_tomo_limber_nointerp_wrapper(l, i, use_linear_ps_limber);
              lP[i][limits.LMIN_tab] = C_ks_tomo_limber_wrapper(limits.LMIN_tab, i);
            }
            #pragma omp parallel for
            for(int k=1; k<ntheta; k++)
            {
              const double l = exp(lnrc + (k - nc)*dlnl);
              if(l > limits.LMIN_tab)
              {
                lP[i][k] = l*C_ks_tomo_limber_wrapper(l, i);
              }
              else if(l == limits.LMIN_tab)
              {
                // do nothing
              } 
              else
              {
                lP[i][k] = l*C_ks_tomo_limber_nointerp_wrapper(l, i,
                  use_linear_ps_limber);
              }
            }
          }
          #pragma omp parallel for
          for (int i=1; i<NSIZE; i++)
          {
            for(int k=0; k<ntheta; k++)
            {
              const double l = exp(lnrc + (k - nc)*dlnl);
              if(l > limits.LMIN_tab)
              {
                lP[i][k] = l*C_ks_tomo_limber_wrapper(l, i);
              }
              else
              {
                lP[i][k] = l*C_ks_tomo_limber_nointerp_wrapper(l, i, use_linear_ps_limber);
              }
            }
          }
        } // Power spectrum on logarithmic bins (ends)
        
        #pragma omp parallel for
        for (int i=0; i<NSIZE; i++)
        {
          fftw_execute(plan[i]); // Execute FFTW in parallel (thread-safe)
        }
        
        for (int i=0; i<NSIZE; i++)
        {
          fftw_free(lP[i]);
          fftw_destroy_plan(plan[i]);
        }
        free(lP);
        free(plan);
      }

      double** lP = (double**) fftw_malloc(sizeof(double**)*NSIZE);
      fftwZ** kernel = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      fftwZ** conv = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
      for (int i=0; i<NSIZE; i++)
      {
        lP[i] = (double*) fftw_malloc(ntheta*sizeof(double));
        kernel[i] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
        conv[i] = (fftwZ*) fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
        plan[i] = fftw_plan_dft_c2r_1d(ntheta, conv[i], lP[i], FFTW_ESTIMATE);
        tab[i] = create_double_matrix(0, 1, 0, ntheta-1);
      }
      
      #pragma omp parallel for
      for (int i=0; i<NSIZE; i++)
      {
        double arg[2];
        arg[0] = 0; // bias
        arg[1] = 0; // order of Bessel function

        for(int k=0; k<(ntheta/2+1); k++)
        {
          const double kk = 2*constants.pi*i/(dlnl*ntheta);
          hankel_kernel_FT(kk, kernel[i], arg, 2);
          conv[i][k][0] = flP[i][k][0]*kernel[i][0][0] - flP[i][k][1]*kernel[i][0][1];
          conv[i][k][1] = flP[i][k][1]*kernel[i][0][0] + flP[i][k][0]*kernel[i][0][1];
        }

        // force Nyquist- and 0-frequency-components to be double
        conv[i][0][1] = 0;
        conv[i][ntheta/2][1] = 0;

        fftw_execute(plan[i]); // Execute FFTW in parallel (thread-safe)

        for(int k=0; k<ntheta; k++)
        {
          const double t = exp((nc-k)*dlnl-lnrc); // theta=1/l
          tab[i][0][ntheta-k-1] = lP[i][k]/(t*2*constants.pi*ntheta);
        }
        for (int k=0; k<ntheta; k++)
        {
          table[i][k] = tab[i][0][k];
        }
      }
      
      for (int i=0; i<NSIZE; i++)
      {
        fftw_free(flP[i]);
        fftw_free(lP[i]);
        fftw_free(kernel[i]);
        fftw_free(conv[i]);
        fftw_destroy_plan(plan[i]);
        free_double_matrix(tab[i], 0, 1, 0, ntheta);
      }  
      free(flP);
      free(lP);
      free(kernel);
      free(conv);
      free(plan);
      free(tab);
      }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  
  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }

  return interpol(table[ni], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 1, 1);
}


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Limber Approximation (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// SS ANGULAR CORRELATION FUNCTION - TATT
// -----------------------------------------------------------------------------

// NLA/TA amplitude C1, nz argument only need if per-bin amplitude
double C1_TA(double a, double nz, double growfac_a)
{
  // per-bin IA parameters
  if (like.IA == 3 || like.IA == 5)
  {
    return -nuisance.A_z[(int)nz]*cosmology.Omega_m*nuisance.c1rhocrit_ia/ growfac_a;
  }
  // power law evolution
  return -cosmology.Omega_m * nuisance.c1rhocrit_ia /
         growfac_a * nuisance.A_ia *
         pow(1. / (a * nuisance.oneplusz0_ia), nuisance.eta_ia);
}

// TA source bias parameter, nz argument only need if per-bin amplitude
double b_TA(double a __attribute__((unused)), double nz)
{
  // per-bin IA parameters
  if (like.IA == 5) {
    return nuisance.b_ta_z[(int)nz];
  }
  // power law evolution
  return nuisance.b_ta_z[0];
}

// TT amplitude C2, nz argument only need if per-bin amplitude
double C2_TT(double a, double nz, double growfac_a)
{
  // per-bin IA parameters
  if (like.IA == 5)
  {
    return 5. * nuisance.A2_z[(int)nz] * cosmology.Omega_m *
           nuisance.c1rhocrit_ia * pow(1.0/growfac_a, 2.0);
  }
  // power law evolution
  return 5. * nuisance.A2_ia * cosmology.Omega_m * nuisance.c1rhocrit_ia *
         (1.0 /(growfac_a*growfac_a)) *
         pow(1. / (a * nuisance.oneplusz0_ia), nuisance.eta_ia_tt);
}

double int_for_C_ss_tomo_TATT_EE_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double* ) params;

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  // radial n_z weight for first source bin (for use with IA term)
  const double ws1 = W_source(a, ar[0], hoverh0);

  // radial n_z weight for second source bin (for use with IA term)
  const double ws2 = W_source(a, ar[1], hoverh0);
  // radial lens efficiency for first source bin
  const double wk1 = W_kappa(a, fK, ar[0]);
  // radial lens efficiency for second source bin
  const double wk2 = W_kappa(a, fK, ar[1]);

  // IA parameters for first source bin
  const double C1 = C1_TA(a, ar[0], growfac_a);
  const double b_ta = b_TA(a, ar[0]);
  const double C2 = C2_TT(a, ar[0], growfac_a);

  // IA parameters for second source bin
  const double C1_2 = C1_TA(a, ar[1], growfac_a);
  const double b_ta_2 = b_TA(a, ar[1]);
  const double C2_2 = C2_TT(a, ar[1], growfac_a);

  // GG cosmic shear
  const double pdelta_ak = Pdelta(k, a);
  double res = wk1 * wk2 * pdelta_ak;
  //COCOA: Took these evaluations of the parenthesis - to force them to update
  //COCOA: the static variables in the first call that is done outside OpenMP loop
  const double tmp1 =
    TATT_II_EE(k, a, C1, C2, b_ta, C1_2, C2_2, b_ta_2, growfac_a, pdelta_ak);
  const double tmp2 =
    TATT_GI_E(k, a, C1, C2, b_ta, growfac_a, pdelta_ak);
  const double tmp3 =
    TATT_GI_E(k, a, C1_2, C2_2, b_ta_2, growfac_a, pdelta_ak);
  if (C1 || C1_2 || C2 || C2_2) {
    // II contribution
    res += ws1 * ws2 * tmp1;
    // GI contribution
    res += ws1 * wk2 * tmp2 + ws2 * wk1 * tmp3;
  }
  return res * chidchi.dchida / (fK * fK);
}

double int_for_C_ss_tomo_TATT_BB_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double* ) params;

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  // radial n_z weight for first source bin (for use with IA term)
  const double ws1 = W_source(a, ar[0], hoverh0);

  // radial n_z weight for second source bin (for use with IA term)
  const double ws2 = W_source(a, ar[1], hoverh0);

  // IA parameters for first source bin
  const double C1 = C1_TA(a, ar[0], growfac_a);
  const double b_ta = b_TA(a, ar[0]);
  const double C2 = C2_TT(a, ar[0], growfac_a);

  // IA parameters for second source bin
  const double C1_2 = C1_TA(a, ar[1], growfac_a);
  const double b_ta_2 = b_TA(a, ar[1]);
  const double C2_2 = C2_TT(a, ar[1], growfac_a);

  //COCOA: Took these evaluations of the parenthesis - to force them to update
  //COCOA: the static variables in the first call that is done outside OpenMP loop
  double res = 0.;
  const double tmp1 = TATT_II_BB(k, a, C1, C2, b_ta, C1_2, C2_2, b_ta_2, growfac_a);
  if ((b_ta || C2) && (b_ta_2 || C2_2))
  {
    res = ws1 * ws2 * tmp1;
  }
  return res * chidchi.dchida / (fK * fK);
}


double C_ss_tomo_TATT_EE_limber_nointerp(double l, int ni, int nj)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  double array[3] = {(double) ni, (double) nj, l};

  return int_gsl_integrate_low_precision(int_for_C_ss_tomo_TATT_EE_limber, (void*) array, 
    fmax(amin_source(ni), amin_source(nj)), amax_source(ni), NULL, GSL_WORKSPACE_SIZE);
}

double C_ss_tomo_TATT_BB_limber_nointerp(double l, int ni, int nj)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  double array[3] = {(double) ni, (double) nj, l};

  return int_gsl_integrate_low_precision(int_for_C_ss_tomo_TATT_BB_limber, (void*) array,
    fmax(amin_source(ni), amin_source(nj)), fmin(amax_source_IA(ni), amax_source_IA(nj)),
    NULL, GSL_WORKSPACE_SIZE);
}

double C_ss_tomo_TATT_EE_limber(double l, int ni, int nj)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static double** table;
  static double* sig;
  static int osc[100];

  const int NSIZE = tomo.shear_Npowerspectra;
  const int nell = Ntable.N_ell_TATT;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin) / (nell - 1.);

  if (table == 0)
  {
    table = create_double_matrix(0, NSIZE - 1, 0, nell - 1);
    sig = create_double_vector(0, NSIZE - 1);
  }
  if (recompute_shear(C, N))
  {
    {
      const int k = 0;
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ss_tomo_TATT_EE_limber_nointerp(500., Z1NZ, Z2NZ) < 0)
      {
        sig[k] = -1.;
      }
      #pragma omp parallel for
      for (int i=0; i<nell; i++)
      {
        const double llog = lnlmin + i*dlnl;
        table[k][i] = C_ss_tomo_TATT_EE_limber_nointerp(exp(llog), Z1NZ, Z2NZ);
        if (table[k][i] * sig[k] < 0.)
        {
          #pragma omp critical
          {
            osc[k] = 1;
          }
        }
      }
      if (osc[k] == 0)
      {
        for (int i = 0; i<nell; i++)
        {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    #pragma omp parallel for
    for (int k=1; k<NSIZE; k++)
    {
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ss_tomo_TATT_EE_limber_nointerp(500., Z1NZ, Z2NZ) < 0)
      {
        sig[k] = -1.;
      }
      for (int i=0; i<nell; i++)
      {
        const double llog = lnlmin + i*dlnl;
        table[k][i] = C_ss_tomo_TATT_EE_limber_nointerp(exp(llog), Z1NZ, Z2NZ);
        if (table[k][i] * sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        for (int i = 0; i<nell; i++)
        {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }

  int k = N_shear(ni, nj);
  if(k > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  double f1;
  if (osc[k] == 0)
  {
    f1 = sig[k] * exp(interpol(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 0., 0.));
  }
  else
  {
    f1 = interpol(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 0., 0.);
  }
  if (isnan(f1))
  {
    f1 = 0.;
  }
  return f1;
}

double C_ss_tomo_TATT_BB_limber(double l, int ni, int nj)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static double** table;
  static double* sig;
  static int osc[100];

  const int nell = Ntable.N_ell_TATT;
  const int NSIZE = tomo.shear_Npowerspectra;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin) / (nell - 1.);

  if (table == 0)
  {
    table = create_double_matrix( 0, NSIZE - 1, 0, nell - 1);
    sig = create_double_vector(0, NSIZE - 1);
  }
  if (recompute_shear(C, N))
  {
    {
      const int k = 0;
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ss_tomo_TATT_BB_limber_nointerp(500., Z1NZ, Z2NZ) < 0)
      {
        sig[k] = -1.;
      }
      #pragma omp parallel for
      for (int i=0; i<nell; i++)
      {
        const double llog = lnlmin + i*dlnl;
        table[k][i] = C_ss_tomo_TATT_BB_limber_nointerp(exp(llog), Z1NZ, Z2NZ);
        if (table[k][i] * sig[k] < 0.)
        {
          #pragma omp critical
          {
            osc[k] = 1;
          }
        }
      }
      if (osc[k] == 0)
      {
        for (int i=0; i<nell; i++)
        {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    #pragma omp parallel for
    for (int k=1; k<NSIZE; k++)
    {
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ss_tomo_TATT_BB_limber_nointerp(500., Z1NZ, Z2NZ) < 0)
      {
        sig[k] = -1.;
      }
      for (int i=0; i<nell; i++)
      {
        const double llog = lnlmin + i*dlnl;
        table[k][i] = C_ss_tomo_TATT_BB_limber_nointerp(exp(llog), Z1NZ, Z2NZ);
        if (table[k][i] * sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        for (int i=0; i<nell; i++)
        {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }

  const int k = N_shear(ni, nj);
  if(k > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  double f1;
  if (osc[k] == 0)
  {
    f1 = sig[k] * exp(interpol(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 0, 0));
  }
  else
  {
    f1 = interpol(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 0, 0);
  }
  if (isnan(f1)) 
  {
    f1 = 0.;
  }
  return f1;
}

// -----------------------------------------------------------------------------
// SS ANGULAR CORRELATION FUNCTION - NON TATT
// -----------------------------------------------------------------------------

double int_for_C_ss_tomo_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double* ) params;

  // Cocoa: added extra options to reduce code duplication
  const int use_linear_ps = ar[3];

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor = (ar[2] - 1.)*(ar[2])*(ar[2] + 1.)*(ar[2] + 2.);

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double ell4 = ell*ell*ell*ell;

  const double ws1 = W_source(a, ar[0], hoverh0);
  const double ws2 = W_source(a, ar[1], hoverh0);
  const double wk1 = W_kappa(a, fK, ar[0]);
  const double wk2 = W_kappa(a, fK, ar[1]);

  double res = 0.0;

  switch(like.IA)
  {
    case 0:
    {
      res = wk1*wk2;

      break;
    }
    case 1:
    {
      const double norm =
        A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      res = ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm + wk1*wk2;

      break;
    }
    case 3:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      res = ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm + wk1*wk2;

      break;
    }
    case 4:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a*
        nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia), nuisance.eta_ia);

      res = ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm + wk1*wk2;

      break;
    }
    default:
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
  }

  if(use_linear_ps == 1)
  {
    res *= p_lin(k,a);
  }
  else
  {
    res *= Pdelta(k,a);
  }

  return res*(chidchi.dchida/(fK*fK))*ell_prefactor/ell4;
}

double C_ss_tomo_limber_nointerp(double l, int ni, int nj, int use_linear_ps)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  double array[4] = {(double) ni, (double) nj, l, (double) use_linear_ps};

  int j,k;
  if (ni <= nj)
  {
    j = nj;
    k = ni;
  }
  else
  {
    j = ni;
    k = nj;
  }

  switch(like.IA)
  { // different IA might require different integrator precision
    case 0:
    {
      return int_gsl_integrate_low_precision(
        int_for_C_ss_tomo_limber,
        (void*) array,
        amin_source(j),
        0.99999,
        NULL,
        GSL_WORKSPACE_SIZE
      );

      break;
    }
    case 1:
    {
      return int_gsl_integrate_low_precision(
        int_for_C_ss_tomo_limber,
        (void*) array,
        amin_source(j),
        amax_source(k),
        NULL,
        GSL_WORKSPACE_SIZE
      );

      break;
    }
    case 3:
    {
      return int_gsl_integrate_low_precision(
        int_for_C_ss_tomo_limber,
        (void*) array,
        amin_source(j),
        amax_source(k),
        NULL,
        GSL_WORKSPACE_SIZE
      );

      break;
    }
    case 4:
    {
      return int_gsl_integrate_low_precision(
        int_for_C_ss_tomo_limber,
        (void*) array,
        amin_source(j),
        amax_source(k),
        NULL,
        GSL_WORKSPACE_SIZE
      );

      break;
    }
    default:
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
  }
}

double C_ss_tomo_limber(double l, int ni, int nj)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static double **table;  

  const int nell = Ntable.N_ell;
  const int NSIZE = tomo.shear_Npowerspectra;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin)/(nell - 1.);

  if (table == 0)
  {
    table = create_double_matrix(0, NSIZE - 1, 0, nell - 1);
  }

  if (recompute_shear(C, N))
  {
    {
      const int k = 0;
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      {
        const int i = 0;
        const double lnl = lnlmin + i*dlnl;
        const double l = exp(lnl);
        table[k][i]= log(C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, use_linear_ps_limber));
      }
      #pragma omp parallel for
      for (int i=1; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        const double l = exp(lnl);
        table[k][i]= log(C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, use_linear_ps_limber));
      }
    }
    #pragma omp parallel for
    for (int k=1; k<NSIZE; k++)
    {
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      for (int i=0; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        const double l = exp(lnl);
        table[k][i]= log(C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, use_linear_ps_limber));
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }

  const int q = N_shear(ni, nj);
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  
  const double f1 = exp(interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1.));
  if (isnan(f1)) 
  {
    return 0.0;
  }
  return f1;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// GS ANGULAR CORRELATION FUNCTION - TATT
// -----------------------------------------------------------------------------

double int_for_C_gs_tomo_limber_TATT(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }

  double* ar = (double*) params;
  const double l = ar[2];
  const int ni = (int) ar[0];
  const int nj = (int) ar[1];

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = l + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  // radial n_z weight for source bin (for use with IA term)
  const double ws = W_source(a, nj, hoverh0);
  // radial lens efficiency for source bin
  const double wk = W_kappa(a, fK, nj);
  // IA parameters for first source bin
  const double C1 = C1_TA(a, nj, growfac_a);
  const double b_ta = b_TA(a, nj);
  const double C2 = C2_TT(a, nj, growfac_a);
  // radial n_z weight for lens bin (for use with clustering term)
  const double w_density = W_HOD(a, ni, hoverh0);
  // lens efficiency *b_mag for lens bin (for lens magnification)
  const double w_mag = W_mag(a, fK, ni) * gbias.b_mag[ni];
  // galaxy bias parameters for lens bin
  const double b1 = gbias.b1_function(1. / a - 1., ni);
  const double b2 = gbias.b2[ni];
  const double bs2 = gbias.bs2[ni];

  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;
  const double Pnl = Pdelta(k, a);
  double P_1loop = b1 * Pnl;
  if (w_density * b2 != 0)
  {
    P_1loop += g4 * (0.5 * b2 * PT_d1d2(k) + 0.5 * bs2 * PT_d1s2(k) +
                     0.5 * b3nl_from_b1(b1) * PT_d1d3(k));
  }

  // 1-loop P_gm ggl terms
  double res = w_density * wk * P_1loop;
  // lens magnification x G term
  res += w_mag * wk * Pnl;
  // (linear bias lens density + lens magnification) with TATT_GI terms

  //COCOA: Took these evaluations of the parenthesis - to force them to update
  //COCOA: the static variables in the first call that is done outside OpenMP loop
  const double tmp1 = TATT_GI_E(k, a, C1, C2, b_ta, growfac_a, Pnl);
  if (C1 || C2)
  {
    res += (b1 * w_density + w_mag) * ws * tmp1;
  }
  return res * chidchi.dchida / fK / fK;
}

// -----------------------------------------------------------------------------
// GS ANGULAR CORRELATION FUNCTION - NON TATT
// -----------------------------------------------------------------------------

double int_for_C_gs_tomo_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int use_linear_ps = (int) ar[3];

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor1 = (ar[2])*(ar[2]+1.);
  double ell_prefactor2 = (ar[2]-1.)*ell_prefactor1*(ar[2]+2.);
  if (ell_prefactor2 <= 0.)
  {
    ell_prefactor2 = 0.;
  }
  else
  {
    ell_prefactor2 = sqrt(ell_prefactor2);
  }

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2]+0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;

  const double wgal = W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[0])*(ell_prefactor1/ell/ell -1.)*gbias.b_mag[(int) ar[0]];

  double res = 0.0;

  if(include_RSD_GS == 1)
  {
    static double chi_a_min = 0;
    if (chi_a_min == 0)
    {
      chi_a_min = chi(limits.a_min);
    }
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell+1.)/k);
    if (chi_1 > chi_a_min)
    {
      return 0;
    }
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);

    res = (wgal + W_RSD(ell, a_0, a_1, ar[0]));
  }
  else
  {
    res = wgal;
  }

  switch(like.IA)
  {
    case 0:
    {
      res *= W_kappa(a, fK, ar[1]);

      break;
    }
    case 1:
    {
      const double norm =
        A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      res *= (W_kappa(a, fK, ar[1]) - W_source(a, ar[1], hoverh0)*norm);

      break;
    }
    case 3:
    {
      const double norm = nuisance.A_z[(int)ar[1]]*cosmology.Omega_m*
        nuisance.c1rhocrit_ia/growfac_a;

      res *= (W_kappa(a, fK, ar[1]) - W_source(a, ar[1], hoverh0)*norm);

      break;
    }
    case 4:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a*
        nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);

      res *= (W_kappa(a, fK, ar[1]) - W_source(a, ar[1], hoverh0)*norm);

      break;
    }
    default:
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
  }

  if (use_linear_ps == 1)
  {
    res *= p_lin(k,a);
  }
  else
  {
    res *= Pdelta(k,a);
  }

  return res*(chidchi.dchida/(fK*fK))*(ell_prefactor2/(ell*ell));
}

double int_for_C_gs_tomo_limber_withb2(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor1 = (ar[2])*(ar[2]+1.);
  double ell_prefactor2 = (ar[2]-1.)*ell_prefactor1*(ar[2]+2.);
  if (ell_prefactor2 <= 0.)
  {
    ell_prefactor2 = 0.;
  }
  else
  {
    ell_prefactor2 = sqrt(ell_prefactor2);
  }

  const double b1 = gbias.b1_function(1./a - 1., (int)ar[0]);
  const double b2 = gbias.b2[(int)ar[0]];
  const double bs2 = gbias.bs2[(int)ar[0]];

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;

  const double wgal = W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[0])*(ell_prefactor1/ell/ell -1.)*gbias.b_mag[(int) ar[0]];

  double linear_part;
  double non_linear_part = W_HOD(a, ar[0], hoverh0);

  if(include_RSD_GS == 1)
  {
    static double chi_a_min = 0;
    if (chi_a_min == 0)
    {
      chi_a_min = chi(limits.a_min);
    }
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell+1.)/k);
    if (chi_1 > chi_a_min)
    {
      return 0;
    }
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);

    linear_part = (wgal + W_RSD(ell, a_0, a_1, ar[0]))*Pdelta(k,a);
  }
  else
  {
    linear_part = wgal*Pdelta(k,a);
  }

  const double ws = W_source(a, ar[1], hoverh0);
  const double wk = W_kappa(a, fK, ar[1]);

  switch(like.IA)
  {
    case 0:
    {
      linear_part *= wk;
      non_linear_part *= wk;

      break;
    }
    case 1:
    {
      const double norm =
        A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      linear_part *= (wk - ws*norm);
      non_linear_part *= (wk - ws*norm);

      break;
    }
    case 3:
    {
      const double norm = nuisance.A_z[(int)ar[1]]*cosmology.Omega_m*
        nuisance.c1rhocrit_ia/growfac_a;

      linear_part *= (wk - ws*norm);
      non_linear_part *= (wk - ws*norm);

      break;
    }
    case 4:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a*
        nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);

      linear_part *= (wk - ws*norm);
      non_linear_part *= (wk - ws*norm);

      break;
    }
    default:
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
  }

  non_linear_part *= g4*(0.5*b2*PT_d1d2(k) +
      0.5*bs2*PT_d1s2(k) + 0.5*b3nl_from_b1(b1)*PT_d1d3(k));

  return (linear_part + non_linear_part)*
    (chidchi.dchida/(fK*fK))*(ell_prefactor2/(ell*ell));
}

double C_gs_tomo_limber_nointerp(double l, int nl, int ns, int use_linear_ps)
{
  if(nl < -1 || nl > tomo.clustering_Nbin -1 || ns < -1 || ns > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", nl, ns);
    exit(1);
  }

  double array[4] = {(double) nl, (double) ns, l, use_linear_ps};

  switch(like.IA)
  { // different IA might require different integrator precision
    case 0:
    {
      if (gbias.b2[nl] && use_linear_ps == 0)
      {
        return int_gsl_integrate_low_precision(int_for_C_gs_tomo_limber_withb2, (void*) array,
          amin_lens(nl), 0.99999, NULL, GSL_WORKSPACE_SIZE);
      }
      return int_gsl_integrate_medium_precision(int_for_C_gs_tomo_limber, (void*) array,
        amin_lens(nl), 0.99999, NULL, GSL_WORKSPACE_SIZE);

      break;
    }
    case 1:
    {
     if (gbias.b2[nl] && use_linear_ps == 0)
      {
        return int_gsl_integrate_medium_precision(int_for_C_gs_tomo_limber_withb2, (void*) array,
          amin_lens(nl), 0.99999, NULL, GSL_WORKSPACE_SIZE);
      }
      return int_gsl_integrate_low_precision(int_for_C_gs_tomo_limber, (void*) array,
        amin_lens(nl), 0.99999, NULL, GSL_WORKSPACE_SIZE);

      break;
    }
    case 3:
    {
      if (gbias.b2[nl] && use_linear_ps == 0)
      {
        return int_gsl_integrate_medium_precision(int_for_C_gs_tomo_limber_withb2, (void*) array,
          amin_lens(nl), 0.99999, NULL, GSL_WORKSPACE_SIZE);
      }
      return int_gsl_integrate_low_precision(int_for_C_gs_tomo_limber, (void*) array,
        amin_lens(nl), 0.99999, NULL, GSL_WORKSPACE_SIZE);

      break;
    }
    case 4:
    {
      if (gbias.b2[nl] && use_linear_ps == 0)
      {
        return int_gsl_integrate_medium_precision(int_for_C_gs_tomo_limber_withb2, (void*) array,
          amin_lens(nl), 0.99999, NULL, GSL_WORKSPACE_SIZE);
      }
      return int_gsl_integrate_low_precision(int_for_C_gs_tomo_limber, (void*) array, 
        amin_lens(nl), 0.9999, NULL,GSL_WORKSPACE_SIZE);

      break;
    }
    case 5:
    {
      if(use_linear_ps == 1)
      {
        log_fatal("use linear power spectrum option not implemented with TATT");
        exit(1);
      }
      return int_gsl_integrate_low_precision(int_for_C_gs_tomo_limber_TATT,
        (void*) array, amin_lens(nl), 0.9999, NULL, GSL_WORKSPACE_SIZE
      );

      break;
    }
    case 6:
    {
      if(use_linear_ps == 1)
      {
        log_fatal("use linear power spectrum option not implemented with TATT");
        exit(1);
      }

      return int_gsl_integrate_low_precision(int_for_C_gs_tomo_limber_TATT, (void*) array,
        amin_lens(nl), 0.9999, NULL, GSL_WORKSPACE_SIZE);

      break;
    }
    default:
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
  }
}

double C_gs_tomo_limber(double l, int ni, int nj)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;
  static double* sig;
  static int osc[100];

  const int NSIZE = tomo.ggl_Npowerspectra;
  const int nell = (like.IA == 5 || like.IA == 6) ? Ntable.N_ell_TATT :  Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin) / (nell - 1.);

  if (table == 0) 
  {
    table = create_double_matrix(0, NSIZE - 1, 0, nell - 1);
    sig = create_double_vector(0, NSIZE - 1);
  }
  if (recompute_gs(C, G, N))
  {
    if (like.IA == 5 || like.IA == 6) // TATT MODELING
    { 
      {
        const int k = 0;
        const int ZLNZ = ZL(k);
        const int ZSNZ = ZS(k);
        sig[k] = 1.;
        osc[k] = 0;
        if (C_gs_tomo_limber_nointerp(500., ZLNZ, ZSNZ, use_linear_ps_limber) < 0)
        {
          sig[k] = -1.;
        }
        #pragma omp parallel for
        for (int i=0; i<nell; i++)
        {
          const double lnl = lnlmin + i*dlnl;
          const double l = exp(lnl);
          table[k][i] = C_gs_tomo_limber_nointerp(l, ZLNZ, ZSNZ, use_linear_ps_limber);
          if (table[k][i] * sig[k] < 0.)
          {
            #pragma omp critical
            {
              osc[k] = 1;
            }
          }
        }
        if (osc[k] == 0)
        {
          for (int i=0; i<nell; i++)
          {
            table[k][i] = log(sig[k] * table[k][i]);
          }
        }
      }
      #pragma omp parallel for
      for (int k=1; k<NSIZE; k++)
      {
        const int ZLNZ = ZL(k);
        const int ZSNZ = ZS(k);
        sig[k] = 1.;
        osc[k] = 0;
        if (C_gs_tomo_limber_nointerp(500., ZLNZ, ZSNZ, use_linear_ps_limber) < 0)
        {
          sig[k] = -1.;
        }
        for (int i=0; i<nell; i++)
        {
          const double lnl = lnlmin + i*dlnl;
          const double l = exp(lnl);
          table[k][i] = C_gs_tomo_limber_nointerp(l, ZLNZ, ZSNZ, use_linear_ps_limber);
          if (table[k][i] * sig[k] < 0.)
          {
            #pragma omp critical
            {
              osc[k] = 1;
            }
          }
        }
        if (osc[k] == 0)
        {
          for (int i=0; i<nell; i++)
          {
            table[k][i] = log(sig[k] * table[k][i]);
          }
        }
      }
    }
    else
    {
      {
        const int k = 0;
        const int ZLNZ = ZL(k);
        const int ZSNZ = ZS(k);
        {
          int i = 0;
          const double lnl = lnlmin + i*dlnl;
          const double l = exp(lnl);
          table[k][i] = log(C_gs_tomo_limber_nointerp(l, ZLNZ, ZSNZ, use_linear_ps_limber));
        }
        #pragma omp parallel for
        for (int i=1; i<nell; i++)
        {
          const double lnl = lnlmin + i*dlnl;
          const double l = exp(lnl);
          table[k][i] = log(C_gs_tomo_limber_nointerp(l, ZLNZ, ZSNZ, use_linear_ps_limber));
        }
      }
      #pragma omp parallel for
      for (int k=1; k<NSIZE; k++)
      {
        const int ZLNZ = ZL(k);
        const int ZSNZ = ZS(k);
        for (int i=0; i<nell; i++)
        {
          const double lnl = lnlmin + i*dlnl;
          const double l = exp(lnl);
          table[k][i] = log(C_gs_tomo_limber_nointerp(l, ZLNZ, ZSNZ, use_linear_ps_limber));
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
    update_galpara(&G);
  }
  
  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }

  const int k = N_ggl(ni, nj);
  if(k > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  if (like.IA == 5 || like.IA == 6) // TATT MODELING
  {     
    double f1 = 0.;
    if (test_zoverlap(ni, nj) && osc[k] == 0)
    {
      f1 = sig[k] * exp(interpol(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 0, 0));
    }
    else if (test_zoverlap(ni, nj) && osc[k] == 1)
    {
      f1 = interpol(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 0, 0);
    }
    if (isnan(f1))
    {
      f1 = 0;
    }
    return f1;
  }
  else
  {
    if(test_zoverlap(ni, nj))
    {
      const double f1 = exp(interpol_fitslope(table[k], nell, lnlmin, lnlmax, dlnl, lnl, 1));
      if (isnan(f1))
      {
       return 0.0;
      }
      return f1;
    }
    else
    {
      return 0.0;
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gg_tomo_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int use_linear_ps = (int) ar[3];

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor = (ar[1])*(ar[1] + 1.);

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  const double wgal1 = W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[0])*(ell_prefactor/(ell*ell) - 1.0)*gbias.b_mag[(int) ar[0]];

  const double wgal2 = W_gal(a, ar[1], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[1])*(ell_prefactor/(ell*ell) - 1.0)*gbias.b_mag[(int) ar[1]];

  double res = 0.0;
  if(include_RSD_GG == 1)
  {
    static double chi_a_min = 0;
    if (chi_a_min == 0)
    {
      chi_a_min = chi(limits.a_min);
    }
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell + 1.)/k);
    if (chi_1 > chi_a_min)
    {
      return 0;
    }
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);

    res = (wgal1 + W_RSD(ell, a_0, a_1, ar[0]))*
      (wgal2 + W_RSD(ell, a_0, a_1, ar[1]));
  }
  else
  {
    res = wgal1*wgal2;
  }

  if (use_linear_ps == 1)
  {
    res *= p_lin(k,a);
  }
  else
  {
    res *= Pdelta(k,a);
  }

  return (res*chidchi.dchida/(fK*fK));
}

double int_for_C_gg_tomo_limber_withb2(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor = (ar[1])*(ar[1] + 1.);

  const double b1 = gbias.b1_function(1./a - 1., (int) ar[0]);
  const double b2 = gbias.b2[(int) ar[0]];
  const double bs2 = gbias.bs2[(int) ar[0]];

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double s4 = 0.; // PT_sigma4(k);

  const double wgal1 = W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[0])*(ell_prefactor/(ell*ell) - 1.0)*gbias.b_mag[(int) ar[0]];

  const double wgal2 = W_gal(a, ar[1], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[1])*(ell_prefactor/(ell*ell) - 1.0)*gbias.b_mag[(int) ar[1]];

  double linear_part = 0.0;
  if(include_RSD_GG == 1)
  {
    static double chi_a_min = 0;
    if (chi_a_min == 0)
    {
      chi_a_min = chi(limits.a_min);
    }
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell+1.)/k);
    if (chi_1 > chi_a_min)
    {
      return 0;
    }
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);

    linear_part = (wgal1 + W_RSD(ell, a_0, a_1, ar[0]))*
      (wgal2 + W_RSD(ell, a_0, a_1, ar[1]))*Pdelta(k, a);
  }
  else
  {
    linear_part = wgal1*wgal2*Pdelta(k, a);
  }

  const double non_linear_part = g4*W_HOD(a, ar[0], hoverh0)*
    W_HOD(a, ar[1], hoverh0)*
    (b1 * b2 * PT_d1d2(k) + 0.25 * b2 * b2 * (PT_d2d2(k) - 2. * s4) +
    b1 * bs2 * PT_d1s2(k) + 0.5 * b2 * bs2 * (PT_d2s2(k) - 4. / 3. * s4) +
    0.25 * bs2 * bs2 * (PT_s2s2(k) - 8. / 9. * s4) +
    b1 * b3nl_from_b1(b1) * PT_d1d3(k));

  return (linear_part + non_linear_part)*(chidchi.dchida/(fK*fK));
}

// WARNING: C_gg beyond linear bias for cross-tomography bins not yet supported
double C_gg_tomo_limber_nointerp(double l, int ni, int nj, int use_linear_ps)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  double array[4] = {(double) ni, (double) nj, l, (double) use_linear_ps};

  if ((gbias.b2[ni] || gbias.b2[nj]) && use_linear_ps == 0)
  {
    if (ni != nj)
    {
      // WARNING: C_gg beyond linear bias for cross-tomography bins not yet supported
      return int_gsl_integrate_low_precision(int_for_C_gg_tomo_limber, (void*) array,
        fmax(amin_lens(ni), amin_lens(nj)), fmin(amax_lens(ni), amax_lens(nj)), NULL,
        GSL_WORKSPACE_SIZE);
    }
    return int_gsl_integrate_low_precision(int_for_C_gg_tomo_limber_withb2, (void*) array,
      amin_lens(ni), amax_lens(ni), NULL, GSL_WORKSPACE_SIZE);
  }
  else
  {
    if (ni == nj)
    {
      return int_gsl_integrate_medium_precision(int_for_C_gg_tomo_limber, (void*) array,
        amin_lens(ni), 0.999999, NULL, GSL_WORKSPACE_SIZE);
    }
    return int_gsl_integrate_low_precision(int_for_C_gg_tomo_limber, (void*) array,
      amin_lens(nj), 0.99999, NULL, GSL_WORKSPACE_SIZE); // zi<=zj
  }
}

double C_gg_tomo_limber(double l, int ni, int nj) // cross redshift bin not supported 
{ 
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX+1);
  const double dlnl = (lnlmax - lnlmin) / (nell);
  const int NSIZE = tomo.clustering_Nbin;

  if (table == 0)
  {
    table = create_double_matrix(0, NSIZE - 1, 0, nell - 1);
  }
  
  if (recompute_gg(C, G, N))
  {
    {
      const int k = 0;
      const int ZCL1 = k; // cross redshift bin not supported so not using ZCL1(k)
      const int ZCL2 = k; // cross redshift bin not supported so not using ZCL2(k)
      {
        const int p = 0;
        const double lnl = lnlmin + p*dlnl;
        const double l = exp(lnl);
        const double result = C_gg_tomo_limber_nointerp(l, ZCL1, ZCL2, use_linear_ps_limber);
        if (result <= 0)
        {
          table[k][p] = -100;
        }
        else
        {
          table[k][p] = log(result);
        }
      }
      #pragma omp parallel for
      for(int p=1; p<nell; p++)
      {
        const double lnl = lnlmin + p*dlnl;
        const double l = exp(lnl);
        const double result = C_gg_tomo_limber_nointerp(l, ZCL1, ZCL2, use_linear_ps_limber);
        if (result <= 0)
        {
          table[k][p] = -100;
        }
        else
        {
          table[k][p] = log(result);
        }
      }
    }
    #pragma omp parallel for
    for (int k=1; k<NSIZE; k++)  
    {
      const int ZCL1 = k; // cross redshift bin not supported so not using ZCL1(k)
      const int ZCL2 = k; // cross redshift bin not supported so not using ZCL2(k)
      for(int p=0; p<nell; p++)
      {
        const double lnl = lnlmin + p*dlnl;
        const double l = exp(lnl);
        const double result = C_gg_tomo_limber_nointerp(l, ZCL1, ZCL2, use_linear_ps_limber);
        if (result <= 0)
        {
          table[k][p] = -100;
        }
        else
        {
          table[k][p] = log(result);
        }
      }
    }
    update_galpara(&G);
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  if (ni != nj) 
  {
    log_fatal("cross-tomography not supported");
    exit(1);
  }

  const int q = ni; // cross redshift bin not supported so not using N_CL(ni, nj) instead of ni
  if (q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }

  const double f1 = exp(interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1));
  if (isnan(f1))
  {
    return 0.0;
  }
  return f1;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gk_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int use_linear_ps = (int) ar[2];

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor = (ar[1])*(ar[1] + 1.);

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[1] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;

  const double wgal = W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[0])*(ell_prefactor/(ell*ell) - 1.0)*gbias.b_mag[(int) ar[0]];

  double res;
  if(include_RSD_GK == 1)
  {
    static double chi_a_min = 0;
    if (chi_a_min == 0)
    {
      chi_a_min = chi(limits.a_min);
    }
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell+1.)/k);
    if (chi_1 > chi_a_min)
    {
      return 0;
    }
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);

    res = (wgal + W_RSD(ell, a_0, a_1, ar[0]))*W_k(a, fK);
  }
  else
  {
    res = wgal*W_k(a, fK);
  }

  if (use_linear_ps == 1)
  {
    res *= p_lin(k,a);
  }
  else
  {
    res *= Pdelta(k,a);
  }

  return (res*chidchi.dchida/(fK*fK))*(ell_prefactor/(ell*ell));
}

double int_for_C_gk_limber_withb2(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor = (ar[1])*(ar[1] + 1.);

  //const double b1 = gbias.b1_function(1./a - 1.,(int) ar[0]);
  const double b2 = gbias.b2[(int) ar[0]];
  const double bs2 = gbias.bs2[(int) ar[0]];

  const double growfac_a = growfac(a);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double ell = ar[1] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double WK = W_k(a, fK);

  const double wgal = W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[0])*(ell_prefactor/(ell*ell) - 1.0)*gbias.b_mag[(int) ar[0]];

  double linear_part = 0.0;
  if(include_RSD_GK == 1)
  {
    static double chi_a_min = 0;
    if (chi_a_min == 0)
    {
      chi_a_min = chi(limits.a_min);
    }
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell+1.)/k);
    if (chi_1 > chi_a_min)
    {
      return 0;
    }
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);

    linear_part = (wgal + W_RSD(ell, a_0, a_1, ar[0]))*WK*Pdelta(k, a);
  }
  {
    linear_part = wgal*WK*Pdelta(k, a);
  }

  const double non_linear_part = g4*W_HOD(a, ar[0], hoverh0)*WK*(
    0.5*b2*PT_d1d2(k) + 0.5*bs2*PT_d1s2(k));

  return (linear_part + non_linear_part)*(chidchi.dchida/(fK*fK))*
    (ell_prefactor/(ell*ell));
}

double C_gk_tomo_limber_nointerp(double l, int nl, int use_linear_ps)
{
  if(nl < -1 || nl > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", nl);
    exit(1);
  }

  double array[3] = {(double) nl, l, (double) use_linear_ps};

  if ((gbias.b2[nl] || gbias.b2[nl]) && use_linear_ps == 0)
  {
    return int_gsl_integrate_medium_precision(
      int_for_C_gk_limber_withb2,
      (void*) array,
      amin_lens(nl),
      amax_lens(nl),
      NULL,
      GSL_WORKSPACE_SIZE
    );
  }
  return int_gsl_integrate_medium_precision(
    int_for_C_gk_limber,
    (void*) array,
    amin_lens(nl),
    0.99999,
    NULL,
    GSL_WORKSPACE_SIZE
  );
}

double C_gk_tomo_limber(double l, int ni)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int NSIZE = tomo.clustering_Nbin;
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin)/(nell);

  if (table == 0)
  {
    table = create_double_matrix(0, NSIZE - 1, 0, nell - 1);
  }
  if (recompute_gk(C, G, N))
  {
    {
      const int k = 0;
      {
        const int i = 0;
        const double lnl = lnlmin + i*dlnl;
        const double l = exp(lnl);
        table[k][i] = log(C_gk_tomo_limber_nointerp(l, k, use_linear_ps_limber));
      }
      #pragma omp parallel for
      for (int i=1; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        const double l = exp(lnl);
        table[k][i]= log(C_gk_tomo_limber_nointerp(l, k, use_linear_ps_limber));
      }
    }
    #pragma omp parallel for
    for (int k=1; k<NSIZE; k++)
    {
      for (int i=0; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        const double l = exp(lnl);
        table[k][i]= log(C_gk_tomo_limber_nointerp(l, k, use_linear_ps_limber));
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
    update_galpara(&G);
  }

  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }
 
  const int q =  ni; 
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  } 

  const double f1 = exp(interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
  if (isnan(f1)) 
  {
    return 0.0;
  }
  return f1;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_ks_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int use_linear_ps = (int) ar[2];

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor1 = (ar[1])*(ar[1] + 1.);
  double ell_prefactor2 = (ar[1] - 1.)*ell_prefactor1*(ar[1] + 2.);
  if(ell_prefactor2 <= 0.)
  {
    ell_prefactor2 = 0.;
  }
  else
  {
    ell_prefactor2 = sqrt(ell_prefactor2);
  }

  const double ell = ar[1] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double ell4 = ell*ell*ell*ell;

  const double ws1 = W_source(a, ar[0], hoverh0);
  const double wk1 = W_kappa(a, fK, ar[0]);
  const double wk2 = W_k(a, fK);

  double res = 0;

  switch(like.IA)
  {
    case 0:
    {
      res = wk1*wk2;

      break;
    }
    case 1:
    {
      const double norm =
        A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      res = (-ws1*wk2*norm + wk1*wk2);

      break;
    }
    case 3:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      res = (-ws1*wk2*norm + wk1*wk2);

      break;
    }
    case 4:
    {
      const double norm = (cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a)*
        nuisance.A_ia*pow(1.0/(a*nuisance.oneplusz0_ia), nuisance.eta_ia);

      res = (-ws1*wk2*norm + wk1*wk2);

      break;
    }
    default:
    {
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
    }
  }

  if (use_linear_ps == 1)
  {
    res *= p_lin(k,a);
  }
  else
  {
    res *= Pdelta(k,a);
  }

  return (res*chidchi.dchida/(fK*fK))*ell_prefactor1*ell_prefactor2/ell4;
}

double C_ks_tomo_limber_nointerp(double l, int nj, int use_linear_ps)
{
  if(nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input nj = %d", nj);
    exit(1);
  }

  double array[3] = {(double) nj, l, (double) use_linear_ps};

  switch(like.IA)
  { // different IA might require different integrator precision
    case 0:
    {
      return int_gsl_integrate_medium_precision(
        int_for_C_ks_limber,
        (void*) array,
        amin_source(nj),
        0.99999,
        NULL,
        GSL_WORKSPACE_SIZE
      );

      break;
    }
    case 1:
    {
      return int_gsl_integrate_medium_precision(
        int_for_C_ks_limber,
        (void*) array,
        amin_source(nj),
        amax_source(nj),
        NULL,
        GSL_WORKSPACE_SIZE
      );

      break;
    }
    case 3:
    { 
      return int_gsl_integrate_medium_precision(
        int_for_C_ks_limber,
        (void*) array,
        amin_source(nj),
        amax_source(nj),
        NULL,
        GSL_WORKSPACE_SIZE
      );

      break;
    }
    case 4:
    {
      return int_gsl_integrate_medium_precision(
        int_for_C_ks_limber,
        (void*) array,
        amin_source(nj),
        amax_source(nj),
        NULL,
        GSL_WORKSPACE_SIZE
      );

      break;
    }
    default:
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
  }
}

double C_ks_tomo_limber(double l, int ni)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static double** table;
  static double* sig;
  static int osc[100];

  const int NSIZE = tomo.shear_Nbin;
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin)/(nell);

  if (table == 0)
  {
    table = create_double_matrix(0, NSIZE - 1, 0, nell - 1);
    sig = create_double_vector(0,tomo.shear_Nbin-1);
  }
  if (recompute_shear(C, N))
  {
    {
      const int k = 0;
      sig[k] = 1.;
      osc[k] = 0;
      const double res = C_ks_tomo_limber_nointerp(500., k, use_linear_ps_limber);
      if (res < 0)
      {
        sig[k] = -1.;
      }
      #pragma omp parallel for
      for (int i=0; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        table[k][i]= C_ks_tomo_limber_nointerp(exp(lnl), k, use_linear_ps_limber);
        if (res*sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        for(int i=0; i<nell; i++)
        {
          table[k][i] = log(sig[k]*table[k][i]);
        }
      }
    }
    #pragma omp parallel for
    for (int k=1; k<NSIZE; k++)
    {
      sig[k] = 1.;
      osc[k] = 0;
      const double res = C_ks_tomo_limber_nointerp(500., k, use_linear_ps_limber);
      if (res < 0)
      {
        sig[k] = -1.;
      }
      for (int i=0; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        table[k][i]= C_ks_tomo_limber_nointerp(exp(lnl), k, use_linear_ps_limber);
        if (res*sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        for(int i=0; i<nell; i++)
        {
          table[k][i] = log(sig[k]*table[k][i]);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  
  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }

  const int q =  ni; 
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  } 

  double f1 = 0.;
  if (osc[ni] == 0)
  {
    f1 = sig[ni]*exp(interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1));
  }
  if (osc[ni] == 1)
  {
    f1 = interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1);
  }
  if (isnan(f1))
  {
    f1 = 0.;
  }
  return f1;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_kk_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int use_linear_ps = (int) ar[1];

  struct chis chidchi = chi_all(a);

  // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor = (ar[1])*(ar[1] + 1.);

  const double ell = ar[0] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double ell4 = ell*ell*ell*ell;
  const double WK = W_k(a, fK);

  double res = WK*WK;

  if(use_linear_ps == 1)
  {
    res *= p_lin(k,a);
  }
  else
  {
    res *= Pdelta(k,a);
  }

  return res*(chidchi.dchida/(fK*fK))*ell_prefactor*ell_prefactor/ell4;
}

double C_kk_limber_nointerp(double l, int use_linear_ps)
{
  double array[2] = {l, (double) use_linear_ps};

  return int_gsl_integrate_medium_precision(int_for_C_kk_limber, (void*) array,
    limits.a_min*(1.+1.e-5), 1.-1.e-5, NULL, GSL_WORKSPACE_SIZE);
}

double C_kk_limber(double l)
{
  static cosmopara C;
  static double* table;
  
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(limits.LMAX + 1);
  const double dlnl = (lnlmax - lnlmin)/(nell);

  if (table == 0)
  {
    table   = create_double_vector(0, nell - 1);
  }
  if (recompute_cosmo3D(C))
  {
    {
      const int i = 0;
      const double lnl = lnlmin + i*dlnl;
       table[i] = log(C_kk_limber_nointerp(exp(lnl), use_linear_ps_limber));
    }
    #pragma omp parallel for
    for (int i=1; i<nell; i++)
    {
      const double lnl = lnlmin + i*dlnl;
       table[i] = log(C_kk_limber_nointerp(exp(lnl), use_linear_ps_limber));
    }
    update_cosmopara(&C);
  }

  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }
  
  double f1 = exp(interpol(table, nell, lnlmin, lnlmax, dlnl, lnl, 1., 1.));
  if (isnan(f1))
  {
    f1 = 0.;
  }
  return f1;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Non Limber (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void f_chi_for_Psi_sh(double* chi_ar, int Nchi, double* f_chi_ar, int nj)
{
  if(nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input nj = %d", nj);
    exit(1);
  }

  const double real_coverH0 = cosmology.coverH0/cosmology.h0;

  for (int i=0; i<Nchi; i++)
  {
    // first convert unit of chi from Mpc to c/H0
    const double a = a_chi(chi_ar[i] / real_coverH0) ;
    const double z = 1./a - 1.;
    const double fK = f_K(chi_ar[i]/real_coverH0);

    if (z > tomo.shear_zmax[nj])
    {
      f_chi_ar[i] = 0.;
    }
    else
    {
      const double wkappa = W_kappa(a, fK, (double) nj);
      const double window_L = wkappa/fK/(real_coverH0*real_coverH0);
      f_chi_ar[i] = window_L*growfac(a); // unit [Mpc^-2]
    }
  }
}

void f_chi_for_Psi_sh_IA(double* chi_ar, int Nchi, double* f_chi_IA_ar, int nj)
{
  if(nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input nj = %d", nj);
    exit(1);
  }

  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;

  for(int i=0; i<Nchi; i++)
  {
    // first convert unit of chi from Mpc to c/H0
    const double a = a_chi(chi_ar[i] / real_coverH0) ;
    const double z = 1./a - 1.;
    struct chis chidchi = chi_all(a);
    const double hoverh0 = hoverh0v2(a, chidchi.dchida);
    const double fK = f_K(chi_ar[i]/real_coverH0);

    const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac(a)*
      nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);

    if( (z<tomo.shear_zmin[nj]) || (z>tomo.shear_zmax[nj]) )
    {
      f_chi_IA_ar[i] = 0.;
    }
    else
    {
      double wsource = W_source(a, (double) nj, hoverh0);
      wsource = (wsource > 0.) ? wsource : 0.;
      const double window_ia = -wsource*norm/fK/(real_coverH0*real_coverH0);
      f_chi_IA_ar[i] = window_ia*growfac(a); // unit [Mpc^-2]
    }
  }
}

// Integrand for galaxy density
void f_chi_for_Psi_cl(double* chi_ar, int Nchi, double* f_chi_ar, int ni)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  const double real_coverH0 = cosmology.coverH0/cosmology.h0; // unit Mpc

  for (int i=0; i<Nchi; i++)
  {
    const double a = a_chi(chi_ar[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1. / a - 1.;
    if ((z < tomo.clustering_zmin[ni]) || (z > tomo.clustering_zmax[ni]))
    {
      f_chi_ar[i] = 0.;
    }
    else
    {
      // get rid of unphysical negatives
      const double tmp1 = pf_photoz(z, ni);
      const double pf = (tmp1 < 0.) ? 0 : tmp1;
      f_chi_ar[i] = chi_ar[i] * pf * growfac(a) *
                    gbias.b1_function(z, ni) * hoverh0(a)/real_coverH0;
    }
  }
}

// Integrand for galaxy density RSD
void f_chi_for_Psi_cl_RSD(double* chi_ar, int Nchi, double* f_chi_RSD_ar, int ni)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;

  for (int i = 0; i < Nchi; i++) {
    const double a = a_chi(chi_ar[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1. / a - 1.;
    if ((z < tomo.clustering_zmin[ni]) || (z > tomo.clustering_zmax[ni]))
    {
      f_chi_RSD_ar[i] = 0.;
    } else
    {
      // get rid of unphysical negatives
      const double tmp1 = pf_photoz(z, ni);
      const double pf = (tmp1 < 0.) ? 0 : tmp1;
      struct growths tmp2 = growfac_all(a);
      f_chi_RSD_ar[i] = -chi_ar[i] * pf * tmp2.D * tmp2.f *
                        hoverh0(a) / real_coverH0;
    }
  }
}

// Integrand for lensing magnification of galaxy density
void f_chi_for_Psi_cl_Mag(double* chi_ar, int Nchi, double* f_chi_Mag_ar, int ni)
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;

  for (int i = 0; i < Nchi; i++)
  {
    const double a = a_chi(chi_ar[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1. / a - 1.;
    const double fK = f_K(chi_ar[i]/real_coverH0 /* convert unit to c/H0 */);
    if (z > tomo.clustering_zmax[ni])
    {
      f_chi_Mag_ar[i] = 0.;
    } else
    {
      const double wmag = W_mag(a, fK, (double) ni);
      const double window_M = wmag / fK / (real_coverH0 * real_coverH0);
      f_chi_Mag_ar[i] = window_M * growfac(a); // unit [Mpc^-2]
    }
  }
}

// Mixture of non-Limber and Limber of C_cl (galaxy clustering)
void C_cl_mixed(int L, int max, int ni, int nj, double* Cl, double dev, double tolerance) 
{
  if(ni < -1 || ni > tomo.clustering_Nbin -1 || nj < -1 || nj > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static double** k1_ar;
  static double** k2_ar;
  static double** Fk1_ar;
  static double** Fk2_ar;
  static double** Fk1_Mag_ar;
  static double** Fk2_Mag_ar;
  static double* chi_ar;

  const int Nell_block = 20; //COCOA: original value = 100
  const int Nchi = 500; //COCOA: original value = 1000
  int ell_ar[Nell_block];
  double f1_chi_ar[Nchi];
  double f2_chi_ar[Nchi];
  double f1_chi_RSD_ar[Nchi];
  double f2_chi_RSD_ar[Nchi];
  double f1_chi_Mag_ar[Nchi];
  double f2_chi_Mag_ar[Nchi];

  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  const double chi_min = chi(1./(1.+0.002))*real_coverH0; // DIMENSIONELESS
  const double chi_max = chi(1./(1.+4.))*real_coverH0; // DIMENSIONELESS
  const double dlnchi = log(chi_max/chi_min) / (Nchi - 1.);
  const double dlnk = dlnchi;

  // COCOA: no need to create/destroy arrays with same size at every call
  if(k1_ar == 0)
  {
    k1_ar = (double**) malloc(Nell_block * sizeof(double*));
    k2_ar = (double**) malloc(Nell_block * sizeof(double*));
    Fk1_ar = (double**) malloc(Nell_block * sizeof(double*));
    Fk2_ar = (double**) malloc(Nell_block * sizeof(double*));
    Fk1_Mag_ar = (double**) malloc(Nell_block * sizeof(double*));
    Fk2_Mag_ar = (double**) malloc(Nell_block * sizeof(double*));
    for (int i = 0; i < Nell_block; i++) {
      k1_ar[i] = (double*) malloc(Nchi * sizeof(double));
      k2_ar[i] = (double*) malloc(Nchi * sizeof(double));
      Fk1_ar[i] = (double*) malloc(Nchi * sizeof(double));
      Fk2_ar[i] = (double*) malloc(Nchi * sizeof(double));
      Fk1_Mag_ar[i] = (double*) malloc(Nchi * sizeof(double));
      Fk2_Mag_ar[i] = (double*) malloc(Nchi * sizeof(double));
    }
    chi_ar = (double*) malloc(Nchi * sizeof(double));
  }
  for (int i = 0; i < Nchi; i++)
  {
    // chi_min and chi_max are cosmology dependent
    chi_ar[i] = chi_min * exp(dlnchi * i);
  }
  #pragma omp parallel for
  for (int i = 0; i < Nell_block; i++)
  {
    for (int j = 0; j < Nchi; j++)
    {
      k1_ar[i][j] = 0.0;
      k2_ar[i][j] = 0.0;
      Fk1_ar[i][j] = 0.0;
      Fk2_ar[i][j] = 0.0;
      Fk1_Mag_ar[i][j] = 0.0;
      Fk2_Mag_ar[i][j] = 0.0;
    }
  }
  f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi_ar, ni);
  
  if (ni != nj)
  {
    f_chi_for_Psi_cl(chi_ar, Nchi, f2_chi_ar, nj);
  }

  f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD_ar, ni);
  if (ni != nj)
  {
    f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f2_chi_RSD_ar, nj);
  }

  f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag_ar, ni);
  if (ni != nj)
  {
    f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f2_chi_Mag_ar, nj);
  }

  int i_block = 0;

  config my_config, my_config_RSD, my_config_Mag;
  my_config.nu = 1.;
  my_config.c_window_width = 0.25;
  my_config.derivative = 0;
  my_config.N_pad = 200;
  my_config.N_extrap_low = 0;
  my_config.N_extrap_high = 0;

  my_config_RSD.nu = 1.01;
  my_config_RSD.c_window_width = 0.25;
  my_config_RSD.derivative = 2;
  my_config_RSD.N_pad = 500;
  my_config_RSD.N_extrap_low = 0;
  my_config_RSD.N_extrap_high = 0;

  my_config_Mag.nu = 1.;
  my_config_Mag.c_window_width = 0.25;
  my_config_Mag.derivative = 0;
  my_config_Mag.N_pad = 500;
  my_config_Mag.N_extrap_low = 0;
  my_config_Mag.N_extrap_high = 0;

  // COCOA: LMAX_NOLIMBER might avoid infinite loop in weird models
  while ((fabs(dev) > tolerance) && (L < limits.LMAX_NOLIMBER))
  {
    for (int i = 0; i < Nell_block; i++)
    {
      ell_ar[i] = i + i_block * Nell_block;
    }

    cfftlog_ells(chi_ar, f1_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k1_ar, Fk1_ar);
    if (ni != nj)
    {
      cfftlog_ells(chi_ar, f2_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k2_ar, Fk2_ar);
    }

    cfftlog_ells_increment(chi_ar, f1_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, 
      k1_ar, Fk1_ar);
    if (ni != nj)
    {
      cfftlog_ells_increment(chi_ar, f2_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, 
        k2_ar, Fk2_ar);
    }

    // Add in lensing magnification contribution
    cfftlog_ells(chi_ar, f1_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar,  Nell_block, 
      k1_ar, Fk1_Mag_ar);
    if (ni != nj)
    {
      cfftlog_ells(chi_ar, f2_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, 
        k2_ar, Fk2_Mag_ar);
    }

    #pragma omp parallel for
    for (int i=0; i<Nell_block; i++)
    {
      const double ell_prefactor = ell_ar[i] * (ell_ar[i] + 1.);
      for (int j = 0; j < Nchi; j++)
      {
        Fk1_ar[i][j] += (ell_prefactor / (k1_ar[i][j] * k1_ar[i][j]) *
                         (gbias.b_mag[ni]) * Fk1_Mag_ar[i][j]);
        if (ni != nj)
        {
          Fk2_ar[i][j] += (ell_prefactor / (k2_ar[i][j] * k2_ar[i][j]) *
                           (gbias.b_mag[nj]) * Fk2_Mag_ar[i][j]);
        }
      }
    }
    {
      const int i = 0;
      double tmp[Nchi];
      #pragma omp parallel for
      for (int j=0; j<Nchi; j++)
      {
        const double k1_cH0 = k1_ar[i][j] * real_coverH0;
        if (ni == nj)
        {
          tmp[j] = Fk1_ar[i][j]*Fk1_ar[i][j]*k1_cH0*k1_cH0*k1_cH0 * p_lin(k1_cH0,1.0);
        } else
        {
          tmp[j] = Fk1_ar[i][j]*Fk2_ar[i][j]*k1_cH0*k1_cH0*k1_cH0*p_lin(k1_cH0,1.0);
        }
      }
      double cl_temp = 0.;
      for (int j=0; j<Nchi; j++) 
      {
        cl_temp += tmp[j];
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2. / M_PI +
         C_gg_tomo_limber_nointerp((double) ell_ar[i], ni, nj, use_linear_ps_limber)
        -C_gg_tomo_limber_nointerp((double) ell_ar[i], ni, nj, 1 /*linear PS*/);
    }
    #pragma omp parallel for
    for (int i=1; i<Nell_block; i++)
    {
      double cl_temp = 0.;
      for (int j = 0; j < Nchi; j++)
      {
        const double k1_cH0 = k1_ar[i][j] * real_coverH0;
        if (ni == nj)
        {
          cl_temp += (Fk1_ar[i][j]) * (Fk1_ar[i][j]) * k1_cH0 * k1_cH0 *
                     k1_cH0 * p_lin(k1_cH0, 1.0);
        } else
        {
          cl_temp += (Fk1_ar[i][j]) * (Fk2_ar[i][j]) * k1_cH0 * k1_cH0 *
                     k1_cH0 * p_lin(k1_cH0, 1.0);
        }
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2. / M_PI +
        C_gg_tomo_limber_nointerp((double) ell_ar[i], ni, nj, use_linear_ps_limber)
        -C_gg_tomo_limber_nointerp((double) ell_ar[i], ni, nj, 1 /*linear PS*/);
    }

    i_block++;

    if(L >= limits.LMAX_NOLIMBER - Nell_block)
    { //Xiao: break before memory leak in next iteration
      break;
    }

    L = i_block * Nell_block - 1;
    dev =
    Cl[L]/C_gg_tomo_limber_nointerp((double) L, ni, nj, use_linear_ps_limber)-1;
  }
  L++;

  Cl[L] = C_gg_tomo_limber((double) L, ni, nj);
  #pragma omp parallel for
  for (int l=L+1; l<max; l++)
  {
    Cl[l] = C_gg_tomo_limber((double) l, ni, nj);
  }
}

// Mixture of non-Limber and Limber of C_cl (G-G lensing)
void C_gl_mixed(int L, int max, int nl, int ns, double* Cl, double dev, double tolerance)
{
  if(nl < -1 || nl > tomo.clustering_Nbin -1 || ns < -1 || ns > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", nl, ns);
    exit(1);
  }

  static double** k1_ar;
  static double** k2_ar;
  static double** Fk1_ar;
  static double** Fk2_ar;
  static double** Fk1_Mag_ar;
  static double** Fk2_Mag_ar;
  static double* chi_ar;

  const int Nell_block = 20; //COCOA: original value = 100
  const int Nchi = 500; //COCOA: original value = 1000
  int ell_ar[Nell_block];
  double f1_chi_ar[Nchi];
  double f1_chi_RSD_ar[Nchi];
  double f1_chi_Mag_ar[Nchi];
  double f2_chi_ar[Nchi];
  double f2_chi_IA_ar[Nchi];

  k1_ar = malloc(Nell_block * sizeof(double* ));
  k2_ar = malloc(Nell_block * sizeof(double* ));
  Fk1_ar = malloc(Nell_block * sizeof(double* ));
  Fk2_ar = malloc(Nell_block * sizeof(double* ));

  Fk1_Mag_ar = malloc(Nell_block * sizeof(double* ));
  for(int i=0; i<Nell_block; i++)
  {
    k1_ar[i] = malloc(Nchi * sizeof(double));
    k2_ar[i] = malloc(Nchi * sizeof(double));
    Fk1_ar[i] = malloc(Nchi * sizeof(double));
    Fk2_ar[i] = malloc(Nchi * sizeof(double));
    Fk1_Mag_ar[i] = malloc(Nchi * sizeof(double));
    for(int j=0; j<Nchi;j ++)
    {
      Fk1_ar[i][j] = 0.;
      Fk2_ar[i][j] = 0.;
      Fk1_Mag_ar[i][j] = 0.;
    }
  }

  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  const double chi_min = 10.;
  const double chi_max = 7000.;
  const double dlnchi = log(chi_max/chi_min) / (Nchi - 1.);
  const double dlnk = dlnchi;

  // COCOA: no need to create/destroy arrays with same size at every call
  if(k1_ar == 0)
  {
    k1_ar = (double**) malloc(Nell_block * sizeof(double*));
    k2_ar = (double**) malloc(Nell_block * sizeof(double*));
    Fk1_ar = (double**) malloc(Nell_block * sizeof(double*));
    Fk2_ar = (double**) malloc(Nell_block * sizeof(double*));
    Fk1_Mag_ar = (double**) malloc(Nell_block * sizeof(double*));
    Fk2_Mag_ar = (double**) malloc(Nell_block * sizeof(double*));
    for (int i=0; i<Nell_block; i++)
    {
      k1_ar[i] = (double*) malloc(Nchi * sizeof(double));
      k2_ar[i] = (double*) malloc(Nchi * sizeof(double));
      Fk1_ar[i] = (double*) malloc(Nchi * sizeof(double));
      Fk2_ar[i] = (double*) malloc(Nchi * sizeof(double));
      Fk1_Mag_ar[i] = (double*) malloc(Nchi * sizeof(double));
      Fk2_Mag_ar[i] = (double*) malloc(Nchi * sizeof(double));
    }
    chi_ar = (double*) malloc(Nchi * sizeof(double));
  }
  for (int i=0; i<Nchi; i++)
  {
    // chi_min and chi_max may be cosmology dependent
    chi_ar[i] = chi_min * exp(dlnchi * i);
  }
  #pragma omp parallel for
  for (int i=0; i<Nell_block; i++)
  {
    for (int j=0; j<Nchi; j++)
    {
      k1_ar[i][j] = 0.0;
      k2_ar[i][j] = 0.0;
      Fk1_ar[i][j] = 0.0;
      Fk2_ar[i][j] = 0.0;
      Fk1_Mag_ar[i][j] = 0.0;
      Fk2_Mag_ar[i][j] = 0.0;
    }
  }

  f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi_ar, nl);

  f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD_ar, nl);

  f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag_ar, nl);

  f_chi_for_Psi_sh(chi_ar, Nchi, f2_chi_ar, ns);

  f_chi_for_Psi_sh_IA(chi_ar, Nchi, f2_chi_IA_ar, ns);

  for(int j=0; j<Nchi; j++)
  {
    f2_chi_ar[j] += f2_chi_IA_ar[j];
  }

  int i_block = 0;

  config my_config, my_config_RSD, my_config_Mag, my_config_L;
  my_config.nu = 1.;
  my_config.c_window_width = 0.25;
  my_config.derivative = 0;
  my_config.N_pad = 200;
  my_config.N_extrap_low = 0;
  my_config.N_extrap_high = 0;

  my_config_RSD.nu = 1.01;
  my_config_RSD.c_window_width = 0.25;
  my_config_RSD.derivative = 2;
  my_config_RSD.N_pad = 200;
  my_config_RSD.N_extrap_low = 0;
  my_config_RSD.N_extrap_high = 0;

  my_config_Mag.nu = 1.;
  my_config_Mag.c_window_width = 0.25;
  my_config_Mag.derivative = 0;
  my_config_Mag.N_pad = 1000;
  my_config_Mag.N_extrap_low = 0;
  my_config_Mag.N_extrap_high = 0;

  my_config_L.nu = 1.;
  my_config_L.c_window_width = 0.25;
  my_config_L.derivative = 0;
  my_config_L.N_pad = 1000.;
  my_config_L.N_extrap_low = 0;
  my_config_L.N_extrap_high = 0;

  // COCOA: LMAX_NOLIMBER might avoid infinite loop in weird models
  while ((fabs(dev) > tolerance) && (L < limits.LMAX_NOLIMBER))
  {
    for(int i=0; i<Nell_block; i++)
    {
      ell_ar[i]=i+i_block*Nell_block;
    }
    cfftlog_ells(chi_ar, f1_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k1_ar,
                Fk1_ar);

    cfftlog_ells_increment(chi_ar, f1_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar,
                Nell_block, k1_ar, Fk1_ar);

    // Add in lensing magnification contribution
    cfftlog_ells(chi_ar, f1_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar,
                Nell_block, k1_ar, Fk1_Mag_ar);

    #pragma omp parallel for
    for(int i=0; i<Nell_block; i++)
    {
      const double ell_prefactor = ell_ar[i]*(ell_ar[i]+1);
      for(int j=0; j<Nchi; j++)
      {
        Fk1_ar[i][j] += (ell_prefactor/(k1_ar[i][j]*k1_ar[i][j])*
                        (gbias.b_mag[nl])*Fk1_Mag_ar[i][j]);
      }
    }

    // shear part
    cfftlog_ells(chi_ar, f2_chi_ar, Nchi, &my_config_L, ell_ar, Nell_block,
                 k2_ar, Fk2_ar);

    for(int i=0; i<Nell_block; i++)
    {
      double ell_prefactor2 =
        (ell_ar[i]-1.)*ell_ar[i]*(ell_ar[i]+1.)*(ell_ar[i]+2.);
      if(ell_prefactor2 <= 0.)
      {
        ell_prefactor2 = 0.;
      }
      else
      {
        ell_prefactor2 = sqrt(ell_prefactor2);
      }
      for(int j=0; j<Nchi; j++)
      {
        Fk2_ar[i][j] *= (ell_prefactor2 / (k1_ar[i][j]*k1_ar[i][j]));
      }
    }

    {
      const int i=0;
      double tmp[Nchi];
      {
        const int j=0;
        const double k1_cH0 = k1_ar[i][j] * real_coverH0;
        tmp[j] = Fk1_ar[i][j]*Fk2_ar[i][j]*k1_cH0*k1_cH0*k1_cH0*p_lin(k1_cH0, 1.0);
      }
      #pragma omp parallel for
      for(int j=1; j<Nchi; j++)
      {
        const double k1_cH0 = k1_ar[i][j] * real_coverH0;
        tmp[j] = Fk1_ar[i][j]*Fk2_ar[i][j]*k1_cH0*k1_cH0*k1_cH0*p_lin(k1_cH0,1.0);
      }
      double cl_temp = 0.;
      for(int j=0; j<Nchi; j++)
      {
        cl_temp += tmp[j];
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI +
       C_gs_tomo_limber_nointerp((double) ell_ar[i], nl, ns, use_linear_ps_limber)
       -C_gs_tomo_limber_nointerp((double) ell_ar[i], nl, ns, 1 /*linear ps*/);
    }
    #pragma omp parallel for
    for(int i=1; i<Nell_block; i++)
    {
      double cl_temp = 0.;
      for(int j=0; j<Nchi; j++)
      {
        const double k1_cH0 = k1_ar[i][j] * real_coverH0;
        cl_temp +=
          (Fk1_ar[i][j])*(Fk2_ar[i][j])*k1_cH0*k1_cH0*k1_cH0*p_lin(k1_cH0,1.0);
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI +
        C_gs_tomo_limber_nointerp((double) ell_ar[i], nl, ns, use_linear_ps_limber)
        -C_gs_tomo_limber_nointerp((double) ell_ar[i], nl, ns, 1 /*linear ps*/);
    }

    i_block++;

    if(L >= limits.LMAX_NOLIMBER - Nell_block)
    { //Xiao: break before memory leak in next iteration
      break;
    }

    L = i_block*Nell_block -1 ;
    dev =
    Cl[L]/C_gs_tomo_limber_nointerp((double) L, nl, ns, use_linear_ps_limber)-1;
  }
  L++;

  Cl[L] = C_gg_tomo_limber((double) L, nl, ns);
  #pragma omp parallel for
  for (int l=L+1; l<max; l++)
  {
    Cl[l] = C_gg_tomo_limber((double) l, nl, ns);
  }
}