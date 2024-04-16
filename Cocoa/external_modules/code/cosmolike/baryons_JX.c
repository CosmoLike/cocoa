#include <assert.h>
#include <gsl/gsl_interp2d.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basics.h"
#include "baryons_JX.h"
#include "structs.h"

#include "log.c/src/log.h"
#include <hdf5.h>

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
  hid_t file_id;
  hid_t dst_id_zBins, dsp_id_zBins;
  hid_t dst_id_logkBins, dsp_id_logkBins;
  hid_t dst_id_logPkR, dsp_id_logPkR;
  hid_t attr_id_Na_bins, attr_id_Nk_bins;
  herr_t status_h5;
  file_id = H5Fopen(lib_file, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0){
    log_fatal("Failed to open file %s!", lib_file);
    exit(1);
  }
  dst_id_zBins = H5Dopen2(file_id, sim_zBins, H5P_DEFAULT);
  dst_id_logkBins = H5Dopen2(file_id, sim_logkBins, H5P_DEFAULT);
  dst_id_logPkR = H5Dopen2(file_id, sim_logPkR, H5P_DEFAULT);
  if ((dst_id_zBins < 0) || (dst_id_logkBins < 0) || (dst_id_logPkR < 0)){
    log_fatal("Failed to read dataset %s, %s, or %s!", sim_zBins, sim_logkBins, 
      sim_logPkR);
    H5Fclose(file_id);
    exit(1);
  }
  dsp_id_zBins = H5Dget_space(dst_id_zBins);
  dsp_id_logkBins = H5Dget_space(dst_id_logkBins);
  dsp_id_logPkR = H5Dget_space(dst_id_logPkR);
  if ((dst_id_zBins < 0) || (dsp_id_logkBins < 0) || (dsp_id_logPkR < 0)){
    log_fatal("Failed to get dataspace %s, %s, or %s!", sim_zBins, 
      sim_logkBins, sim_logPkR);
    H5Dclose(dst_id_zBins);H5Dclose(dst_id_logkBins);H5Dclose(dst_id_logPkR);
    H5Fclose(file_id);
    exit(1);
  }

  // Read metadata from HDF5 dataset attributes
  attr_id_Na_bins = H5Aopen(dst_id_zBins, "Na_bins", H5P_DEFAULT);
  attr_id_Nk_bins = H5Aopen(dst_id_logkBins, "Nk_bins", H5P_DEFAULT);
  if ((attr_id_Na_bins < 0) || (attr_id_Nk_bins < 0)){
    log_fatal("Failed to open attribute Na_bins/Nk_bins!");
    H5Sclose(dsp_id_zBins);H5Sclose(dsp_id_logkBins);H5Sclose(dsp_id_logPkR);
    H5Dclose(dst_id_zBins);H5Dclose(dst_id_logkBins);H5Dclose(dst_id_logPkR);
    H5Fclose(file_id);
    exit(1);
  }
  status_h5 = H5Aread(attr_id_Na_bins, H5T_NATIVE_INT32, &bary.Na_bins);
  if (status_h5 < 0){
    log_fatal("Failed to read attribute Na_bins!");
    H5Aclose(attr_id_Na_bins);H5Aclose(attr_id_Nk_bins);
    H5Sclose(dsp_id_zBins);H5Sclose(dsp_id_logkBins);H5Sclose(dsp_id_logPkR);
    H5Dclose(dst_id_zBins);H5Dclose(dst_id_logkBins);H5Dclose(dst_id_logPkR);
    H5Fclose(file_id);
    exit(1);
  }
  status_h5 = H5Aread(attr_id_Nk_bins, H5T_NATIVE_INT32, &bary.Nk_bins);
  if (status_h5 < 0){
    log_fatal("Failed to read attribute Na_bins!");
    H5Aclose(attr_id_Na_bins);H5Aclose(attr_id_Nk_bins);
    H5Sclose(dsp_id_zBins);H5Sclose(dsp_id_logkBins);H5Sclose(dsp_id_logPkR);
    H5Dclose(dst_id_zBins);H5Dclose(dst_id_logkBins);H5Dclose(dst_id_logPkR);
    H5Fclose(file_id);
    exit(1);
  }
  H5Aclose(attr_id_Na_bins);H5Aclose(attr_id_Nk_bins);
  log_info("Na_bins = %d; Nk_bins = %d", bary.Na_bins, bary.Nk_bins);

  // Read data into zBins, logkBins, and logPkR
  float* zBins = (float*) calloc(bary.Na_bins, sizeof(float));
  float* logkBins = (float *) calloc(bary.Nk_bins, sizeof(float));
  float* logPkR = (float *) calloc(bary.Nk_bins*bary.Na_bins, sizeof(float));
  log_info("Reading HDF5 dataset");
  
  status_h5 = H5Dread(dst_id_zBins, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, 
    H5P_DEFAULT, zBins);
    if (status_h5 < 0) {
        log_fatal("Failed to read zBins!");
        H5Sclose(dsp_id_zBins);H5Sclose(dsp_id_logkBins);H5Sclose(dsp_id_logPkR);
        H5Dclose(dst_id_zBins);H5Dclose(dst_id_logkBins);H5Dclose(dst_id_logPkR);
        H5Fclose(file_id);
        exit(1);
    }
  
  status_h5 = H5Dread(dst_id_logkBins, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, 
    H5P_DEFAULT, logkBins);
    if (status_h5 < 0) {
        log_fatal("Failed to read logkBins!");
        H5Sclose(dsp_id_zBins);H5Sclose(dsp_id_logkBins);H5Sclose(dsp_id_logPkR);
        H5Dclose(dst_id_zBins);H5Dclose(dst_id_logkBins);H5Dclose(dst_id_logPkR);
        H5Fclose(file_id);
        exit(1);
    }
  
  status_h5 = H5Dread(dst_id_logPkR, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, 
    H5P_DEFAULT, logPkR);
    if (status_h5 < 0) {
        log_fatal("Failed to read logPkR!");
        H5Sclose(dsp_id_zBins);H5Sclose(dsp_id_logkBins);H5Sclose(dsp_id_logPkR);
        H5Dclose(dst_id_zBins);H5Dclose(dst_id_logkBins);H5Dclose(dst_id_logPkR);
        H5Fclose(file_id);
        exit(1);
    }
  // Close data set
  H5Sclose(dsp_id_zBins);H5Sclose(dsp_id_logkBins);H5Sclose(dsp_id_logPkR);
  H5Dclose(dst_id_zBins);H5Dclose(dst_id_logkBins);H5Dclose(dst_id_logPkR);
  H5Fclose(file_id);

  // Allocate memory
  if(bary.a_bins == NULL)
  {
    bary.a_bins = (double*) malloc(sizeof(double)*bary.Na_bins);
    if (bary.a_bins == NULL)
    {
      log_fatal("Failed Allocation of a_bins array");
      exit(1);
    }
  }
  if(bary.logk_bins == NULL)
  {
    bary.logk_bins = (double*) malloc(sizeof(double)*bary.Nk_bins);
    if (bary.logk_bins == NULL)
    {
      log_fatal("Failed Allocation of logk_bins array");
      exit(1);
    }
  }
  if(bary.log_PkR == NULL)
  {
    bary.log_PkR = (double*) malloc(sizeof(double)*bary.Nk_bins*bary.Na_bins);
    if (bary.log_PkR == NULL)
    {
      log_fatal("Failed Allocation of log_PkR array");
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
      log_fatal("Failed Allocation of interp2d struct");
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
      // NOTE: here the axes has changed from (Nk, Na) to (Na, Nk)
      int status = gsl_interp2d_set(bary.interp2d, bary.log_PkR, j, i, logPkR[j*bary.Na_bins+i]);
      if (status)
      {
        log_fatal("gsl error %s", gsl_strerror(status));
        exit(1);
      }
  	}
  }

  int status = gsl_interp2d_init(bary.interp2d, bary.logk_bins, bary.a_bins,
 		bary.log_PkR, bary.Nk_bins, bary.Na_bins);
  if (status)
  {
    log_fatal("gsl error %s", gsl_strerror(status));
    exit(1);
  }
  free(zBins);free(logkBins);free(logPkR);
}

void init_baryons(const char* scenario, int sim_id)
{
  // Find the directory where the source code is saved & specify lib file
  char* filename;
  const char* libName = "baryons_logPkR.h5";
  const char* libPath = __FILE__;
  size_t pathLen, filenameLen;
  const char* lastSeparator = strrchr(libPath, '/');
  if (lastSeparator) {
    pathLen = lastSeparator - libPath + 1;// +1 to include the last /
    filenameLen = pathLen + strlen(libName);  
  } else{
    log_fatal("Error in reading the path of %s!", libPath);
    exit(1);
  }
  filename = (char*) malloc(filenameLen + 1);
  if (filename) {
    strncpy(filename, libPath, pathLen);
    strncpy(filename+pathLen, libName, strlen(libName));
    filename[filenameLen] = '\0';
  }
  set_bary_parameters_to_scenario(scenario, filename, sim_id);
}
