#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basics.h"
#include "bias.h"
#include "cosmo3D.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"

double gb1(const double z, const int ni)
{ // g = galaxy, b1 = linear galaxy bias
  if (ni < -1 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  double ans;

  switch(like.galaxy_bias_model[0])
  {
    case B1_PER_BIN:
    {
      ans = nuisance.gb[0][ni];
      break;
    }
    case B1_PER_BIN_EVOLV:
    {
      ans = nuisance.gb[0][ni] * pow((1.0 + z) / (1. + zmean(ni)), 1.0);
      break;
    }
    case B1_PER_BIN_PASS_EVOLV:
    {
      const double z_evolv_passiv = growfac(1. / (z + 1.)) /
                growfac(1. / (1. + 0.5 * (redshift.clustering_zdist_zmin[ni] +
                                          redshift.clustering_zdist_zmax[ni])));
      
      ans = (nuisance.gb[0][ni] - 1.) / z_evolv_passiv + 1.;
      break;
    }
    case B1_GROWTH_SCALING:
    {
      ans = nuisance.gb[0][0] / growfac(1. / (z + 1.));
      break;
    }
    case B1_POWER_LAW:
    {
      ans = nuisance.gb[0][0] * pow(1 + z, nuisance.gb[0][1]);
      break;
    }
    default:
    {
      log_fatal("like.galaxy_bias_model[0] = %d not supported", 
        like.galaxy_bias_model[0]);
      exit(1);  
    }
  }
  return ans;
}

double gb2(const double z, const int ni)
{ // g = galaxy, b2 = non-linear galaxy bias
  if (ni < -1 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  double ans;

  switch(like.galaxy_bias_model[1])
  {
    case B2_PER_BIN:
    {
      ans = nuisance.gb[1][ni];
      break;
    }
    case B2_FROM_B1:
    {
      // fitting formula from Lazeyras et al. 2016 (Eq. 5.2) - 1511.01096
      const double b1 = gb1(z, ni);
      ans = 0.412 - 2.143 * b1 + 0.929 * b1 * b1 + 0.008 * b1 * b1 * b1;
      break;
    }
    default:
    {
      log_fatal("like.galaxy_bias_model[1] = %d not supported", 
        like.galaxy_bias_model[1]);
      exit(1);  
    }
  }
  return ans;
}

double gbs2(const double z, const int ni)
{
  if (ni < -1 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  double ans;

  switch(like.galaxy_bias_model[2])
  {
    case BS2_PER_BIN:
    {
      ans = nuisance.gb[2][ni];
      break;
    }
    case BS2_FROM_B1:
    {
      ans = -4. / 7. * (gb1(z, ni) - 1.0); // 1405.1447v4.pdf
      break;
    }
    default:
    {
      log_fatal("like.galaxy_bias_model[2] = %d not supported", 
        like.galaxy_bias_model[2]);
      exit(1);  
    }
  }
  return ans;
}

double gb3(const double z, const int ni)
{
  if (ni < -1 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  double ans;

  switch(like.galaxy_bias_model[3])
  {
    case B3_PER_BIN:
    {
      ans = nuisance.gb[3][ni];
      break;
    }
    case B3_FROM_B1:
    {
      ans = (gb1(z, ni) - 1.0);
      break;
    }
    default:
    {
      log_fatal("like.galaxy_bias_model[3] = %d not supported", 
        like.galaxy_bias_model[3]);
      exit(1);  
    }
  }
  return ans;
}

double gbmag(const double z, const int ni)
{
  if (ni < -1 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  double ans;

  switch(like.galaxy_bias_model[4])
  {
    case BMAG_PER_BIN:
    {
      ans = nuisance.gb[4][ni];
      break;
    }
    default:
    {
      log_fatal("like.galaxy_bias_model[4] = %d not supported", 
        like.galaxy_bias_model[4]);
      exit(1);  
    }
  }

  return ans;
}


