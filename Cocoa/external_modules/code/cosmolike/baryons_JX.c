#include <assert.h>
#include <gsl/gsl_interp2d.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>

#include "basics.h"
#include "baryons.h"
#include "structs.h"

#include "log.c/src/log.h"
#include "H5Cpp.h"

void set_bary_parameters_to_scenario(const char* scenario, const char* lib_file, int sim_id)
{
  // Initialize baryon feedback suppression & binning
  if(bary.a_bins != NULL || bary.logk_bins != NULL || bary.log_PkR != NULL || bary.interp2d != NULL)
  {
    reset_bary_struct();
  }
  bary.is_Pk_bary = 1;

  // Read data set from HDF5 library file
  char sim_zBins[100]; sprintf(sim_zBins,"/%s/zBins",scenario);
  char sim_logkBins[100]; sprintf(sim_logkBins,"/%s/logkBins",scenario);
  char sim_logPkR[100]; sprintf(sim_logPkR,"/%s/logPkR/sim%d",scenario,sim_id);
  const H5std_string FILE_NAME(lib_file);
  const H5std_string DST_ZBINS(sim_zBins);
  const H5std_string DST_LOGKBINS(sim_logkBins);
  const H5std_string DST_LOGPKR(sim_logPkR);
  H5::H5File file(FILE_NAME, H5F_ACC_RDONLY);
  H5::DataSet dst_zBins  = file.openDataSet(DST_ZBINS);
  H5::DataSet dst_logkBins  = file.openDataSet(DST_LOGKBINS);
  H5::DataSet dst_logPkR  = file.openDataSet(DST_LOGPKR);
  H5::DataSpace sp_zBins = dst_zBins.getSpace();
  H5::DataSpace sp_logkBins = dst_logkBins.getSpace();
  H5::DataSpace sp_logPkR = dst_logPkR.getSpace();

  // Read metadata from HDF5 library 
  const H5std_string NA_BINS("Na_bins");
  const H5std_string NK_BINS("Nk_bins");
  H5::Attribute Na_bins, Nk_bins;
  try {
  	Na_bins = dst_zBins.openAttribute(NA_BINS);
  	Nk_bins = dst_logkBins.openAttribute(NK_BINS);
  } catch(const H5::Exception& e)
  {
  	log_fatal("\x1b[90m{}\x1b[0m: Attribute {} or {} not found!",
        "set_bary_parameters_to_scenario",
        NA_BINS, NK_BINS);
  }
  H5T_class_t Na_bins_type = Na_bins.getTypeClass();
  H5T_class_t Nk_bins_type = Nk_bins.getTypeClass();
  // Check the data type of meta info
  if (Na_bins_type == H5T_INTEGER) {
  	H5::IntType intType = Na_bins.getIntType();
  	size_t size = intType.getSize();
  	if (size == sizeof(int)){
  		log_info("\x1b[90m{}\x1b[0m: Attribute {} is int.",
        "set_bary_parameters_to_scenario",
        NA_BINS);
  		Na_bins.read(intType, &bary.Na_bins);
  	} else if (size == sizeof(long)){
  		log_info("\x1b[90m{}\x1b[0m: Attribute {} is long.",
        "set_bary_parameters_to_scenario",
        NA_BINS);
  		long longValue;
  		Na_bins.read(intType, &longValue);
  		bary.Na_bins = static_cast<int>(longValue);
  	} else if (size == sizeof(long long)) {
  		log_info("\x1b[90m{}\x1b[0m: Attribute {} is long long.",
        "set_bary_parameters_to_scenario",
        NA_BINS);
  		long long longlongValue;
  		Na_bins.read(intType, &longlongValue);
  		bary.Na_bins = static_cast<int>(longlongValue);
  	}
  } else{
  	log_fatal("\x1b[90m{}\x1b[0m: Attribute {} must be integer!",
        "set_bary_parameters_to_scenario",
        NA_BINS);
  }
  if (Nk_bins_type == H5T_INTEGER) {
  	H5::IntType intType = Nk_bins.getIntType();
  	size_t size = intType.getSize();
  	if (size == sizeof(int)){
  		log_info("\x1b[90m{}\x1b[0m: Attribute {} is int.",
        "set_bary_parameters_to_scenario",
        NK_BINS);
  		Nk_bins.read(intType, &bary.Nk_bins);
  	} else if (size == sizeof(long)){
  		log_info("\x1b[90m{}\x1b[0m: Attribute {} is long.",
        "set_bary_parameters_to_scenario",
        NK_BINS);
  		long longValue;
  		Nk_bins.read(intType, &longValue);
  		bary.Nk_bins = static_cast<int>(longValue);
  	} else if (size == sizeof(long long)) {
  		log_info("\x1b[90m{}\x1b[0m: Attribute {} is long long.",
        "set_bary_parameters_to_scenario",
        NK_BINS);
  		long long longlongValue;
  		Nk_bins.read(intType, &longlongValue);
  		bary.Nk_bins = static_cast<int>(longlongValue);
  	}
  } else{
  	log_fatal("\x1b[90m{}\x1b[0m: Attribute {} must be integer!",
        "set_bary_parameters_to_scenario",
        NK_BINS);
  }
  // End data type check
  // Read data into zBins, logkBins, and logPkR
  float* zBins = (float*) calloc(bary.Na_bins * sizeof(float));
  float* logkBins = (float *) calloc(bary.Nk_bins * sizeof(float));
  float* logPkR = (float *) calloc(bary.Nk_bins*bary.Na_bins * sizeof(float));
  sp_zBins.selectAll();
  dst_zBins.read(zBins, H5::PredType::NATIVE_FLOAT, sp_zBins);
  sp_logkBins.selectAll();
  dst_logkBins.read(logkBins, H5::PredType::NATIVE_FLOAT, sp_logkBins);
  sp_logPkR.selectAll();
  dst_logPkR.read(logPkR, H5::PredType::NATIVE_FLOAT, sp_logPkR);
  // Close data set
  Na_bins.close();Nk_bins.close();
  dst_zBins.close();dst_logkBins.close();dst_logPkR.close();
  file.close();

  // Allocate memory
  if(bary.a_bins == NULL)
  {
    bary.a_bins = (double*) malloc(sizeof(double)*bary.Na_bins);
    if (bary.a_bins == NULL)
    {
      log_fatal("\x1b[90m{}\x1b[0m: Failed Allocation of {}",
        "set_bary_parameters_to_TNG100",
        "a_bins array");
      exit(1);
    }
  }
  if(bary.logk_bins == NULL)
  {
    bary.logk_bins = (double*) malloc(sizeof(double)*bary.Nk_bins);
    if (bary.logk_bins == NULL)
    {
      log_fatal("\x1b[90m{}\x1b[0m: Failed Allocation of {}",
        "set_bary_parameters_to_TNG100",
        "logk_bins array");
      exit(1);
    }
  }
  if(bary.log_PkR == NULL)
  {
    bary.log_PkR = (double*) malloc(sizeof(double)*bary.Nk_bins*bary.Na_bins);
    if (bary.log_PkR == NULL)
    {
      log_fatal("\x1b[90m{}\x1b[0m: Failed Allocation of {}",
        "set_bary_parameters_to_TNG100",
        "log_PkR array");
      exit(1);
    }
  }

  bary.T = (gsl_interp2d_type*) gsl_interp2d_bilinear;
  if(bary.interp2d == NULL)
  {
    bary.interp2d = gsl_interp2d_alloc((const gsl_interp2d_type*) bary.T,
      bary.Nk_bins, bary.Na_bins);
    if (bary.interp2d == NULL)
    {
      log_fatal("\x1b[90m{}\x1b[0m: Failed Allocation of {}",
        "set_bary_parameters_to_TNG100", "interp2d struct");
      exit(1);
    }
  }

  // Filling a, k, and PkR arrays, and initialize interpolation
  #pragma omp parallel for
  for (int i=0; i<bary.Na_bins; i++)
  {
  	bary.a_bins[i] =  1./(1 + zBins[i]);
  	for (int j=0; j<bary.Nk_bins; j++)
    {
  	  if (i == 0) bary.logk_bins[j] = logkBins[j];
      int status = gsl_interp2d_set(bary.interp2d, bary.log_PkR, j, i, logPkR[j*bary.Na_bins+i]);
      if (status)
      {
        log_fatal("\x1b[90m{}\x1b[0m: gsl error {}",
          "set_bary_parameters_to_scenario", gsl_strerror(status));
        exit(1);
      }
  	}
  }

  int status = gsl_interp2d_init(bary.interp2d, bary.logk_bins, bary.a_bins,
 		bary.log_PkR, bary.Nk_bins, bary.Na_bins);
  if (status)
  {
    log_fatal("\x1b[90m{}\x1b[0m: gsl error {}",
      "set_bary_parameters_to_TNG100", gsl_strerror(status));
    exit(1);
  }
}

void init_baryons(const char* scenario, int sim_id)
{
  char lib_file[500] = "baryons_logPkR.h5";
  set_bary_parameters_to_scenario(scenario, lib_file, sim_id);
}

int main(){
	init_baryons("BAHAMAS", 1);
}