#include <string>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdexcept>
#include <array>
#include <random>
#include <variant>
#include <cmath> 
#include <tuple>

// SPDLOG
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/cfg/env.h>

// ARMADILLO LIB AND PYBIND WRAPPER (CARMA)
#include <carma.h>
#include <armadillo>

// Python Binding
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
namespace py = pybind11;

// boost library
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>

// COSMOLIKE
#include "cosmolike/basics.h"
#include "cosmolike/baryons.h"
#include "cosmolike/interface_aux.hpp"
#include "cosmolike/interface_redshift.hpp"
#include "cosmolike/interface_redshift.hpp"
#include "cosmolike/redshift_spline.h"
#include "cosmolike/structs.h"

namespace cosmolike_interface
{

void init_lens_sample(std::string multihisto_file, const int Ntomo)
{
  spdlog::debug("{}: Begins", "init_lens_sample");

  if (!(multihisto_file.size() > 0))
  {
    spdlog::critical("{}: empty {} string not supported", 
      "init_lens_sample", "multihisto_file");
    exit(1);
  }
  if (!(Ntomo > 0) || Ntomo > MAX_SIZE_ARRAYS)
  {
    spdlog::critical("{}: {} = {} not supported (max = {})", 
      "init_lens_sample", "Ntomo", Ntomo, MAX_SIZE_ARRAYS);
    exit(1);
  }

  redshift.clustering_photoz = 4;
  redshift.clustering_nbin = Ntomo;
  
  spdlog::debug("{}: {} = {} selected.", "init_lens_sample",
    "clustering_REDSHIFT_FILE", multihisto_file);

  spdlog::debug("{}: {} = {} selected.", "init_lens_sample",
    "clustering_Nbin", Ntomo);

  // READ THE N(Z) FILE BEGINS ------------
  matrix input_table = read_table(multihisto_file);

  if (!input_table.col(0).eval().is_sorted("ascend"))
  {
      spdlog::critical("bad n(z) file (z vector not monotonic)");
      exit(1);
  }

  int cache_update = 0;
  if (redshift.clustering_nzbins != input_table.n_rows ||
      redshift.clustering_zdist_table == NULL)
  {
    cache_update = 1;
  }
  else
  {
    for (int i=0; i<redshift.clustering_nzbins; i++) 
    {
      double** tab = redshift.clustering_zdist_table;   // alias
      double* z_v = redshift.clustering_zdist_table[Ntomo];  // alias

      if (fdiff(z_v[i], input_table(i,0))) 
      {
        cache_update = 1;
        break;
      }
      for (int k=0; k<Ntomo; k++) 
      {  
        if (fdiff(tab[k][i], input_table(i,k+1))) 
        {
          cache_update = 1;
          goto jump;
        }
      }
    }
  }

  jump:

  if (1 == cache_update)
  {
    redshift.clustering_nzbins = input_table.n_rows;
    const int nzbins = redshift.clustering_nzbins;    // alias

    if (redshift.clustering_zdist_table != NULL) 
        free(redshift.clustering_zdist_table);
    redshift.clustering_zdist_table = (double**) malloc2d(Ntomo + 1, nzbins);
    
    double** tab = redshift.clustering_zdist_table;   // alias
    double* z_v = redshift.clustering_zdist_table[Ntomo];  // alias
    
    for (int i=0; i<nzbins; i++) 
    {
      z_v[i] = input_table(i,0);
      for (int k=0; k<Ntomo; k++) 
      {
        tab[k][i] = input_table(i,k+1);
      }
    }
    
    redshift.clustering_zdist_zmin_all = fmax(z_v[0], 1.e-5);
    
    redshift.clustering_zdist_zmax_all = z_v[nzbins-1] + 
      (z_v[nzbins-1] - z_v[0]) / ((double) nzbins - 1.);

    for (int k=0; k<Ntomo; k++) 
    { // Set tomography bin boundaries
      auto nofz = input_table.col(k+1).eval();
      
      arma::uvec idx = arma::find(nofz > 0.999e-8*nofz.max());
      
      redshift.clustering_zdist_zmin[k] = z_v[idx(0)];
      
      redshift.clustering_zdist_zmax[k] = z_v[idx(idx.n_elem-1)];
    }
    // READ THE N(Z) FILE ENDS ------------

    redshift.random_clustering = RandomNumber::get_instance().get();

    pf_photoz(0.1, 0); // init static variables

    for (int k=0; k<Ntomo; k++)
    {
      spdlog::debug("{}: bin {} - {} = {}.", "init_lens_sample", 
        k, "<z_s>", zmean(k));
    } 
  }

  spdlog::debug("{}: Ends", "init_lens_sample");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_source_sample(std::string multihisto_file, const int Ntomo)
{
  spdlog::debug("{}: Begins", "init_source_sample");

  if (!(multihisto_file.size() > 0))
  {
    spdlog::critical("{}: empty {} string not supported",
        "init_source_sample", "multihisto_file");
    exit(1);
  }

  if (!(Ntomo > 0) || Ntomo > MAX_SIZE_ARRAYS)
  {
    spdlog::critical("{}: {} = {} not supported (max = {})",
        "init_source_sample", "Ntomo", Ntomo, MAX_SIZE_ARRAYS);
    exit(1);
  }

  redshift.shear_photoz = 4;
  redshift.shear_nbin = Ntomo;
  
  spdlog::debug("{}: tomo.shear_Npowerspectra = {}", 
      "init_source_sample", tomo.shear_Npowerspectra);

  spdlog::debug("{}: {} = {} selected.", 
      "init_source_sample", "shear_REDSHIFT_FILE", multihisto_file);

  spdlog::debug("{}: {} = {} selected.", 
      "init_source_sample", "shear_Nbin", Ntomo);

  // READ THE N(Z) FILE BEGINS ------------
  matrix input_table = read_table(multihisto_file);
  
  if (!input_table.col(0).eval().is_sorted("ascend"))
  {
    spdlog::critical("bad n(z) file (z vector not monotonic)");
    exit(1);
  }

  int cache_update = 0;
  if (redshift.shear_nzbins != input_table.n_rows ||
      redshift.shear_zdist_table == NULL)
  {
    cache_update = 1;
  }
  else
  {
    double** tab = redshift.shear_zdist_table;        // alias  
    double* z_v = redshift.shear_zdist_table[Ntomo];  // alias
    for (int i=0; i<redshift.shear_nzbins; i++) 
    {
      if (fdiff(z_v[i], input_table(i,0))) 
      {
        cache_update = 1;
        break;
      }
      for (int k=0; k<Ntomo; k++) 
      {
        if (fdiff(tab[k][i], input_table(i,k+1))) 
        {
          cache_update = 1;
          goto jump;
        }
      }
    }
  }

  jump:

  if (1 == cache_update)
  {
    redshift.shear_nzbins = input_table.n_rows;
    const int nzbins = redshift.shear_nzbins; // alias

    if (redshift.shear_zdist_table == NULL) free(redshift.shear_zdist_table);
    redshift.shear_zdist_table = (double**) malloc2d(Ntomo + 1, nzbins);

    double** tab = redshift.shear_zdist_table;        // alias  
    double* z_v = redshift.shear_zdist_table[Ntomo];  // alias
    for (int i=0; i<nzbins; i++) 
    {
      z_v[i] = input_table(i,0);
      for (int k=0; k<Ntomo; k++) 
        tab[k][i] = input_table(i,k+1);
    }
  
    redshift.shear_zdist_zmin_all = fmax(z_v[0], 1.e-5);
    
    redshift.shear_zdist_zmax_all = z_v[nzbins-1] + 
      (z_v[nzbins-1] - z_v[0]) / ((double) nzbins - 1.);

    for (int k=0; k<Ntomo; k++) 
    { // Set tomography bin boundaries
      auto nofz = input_table.col(k+1).eval();
      
      arma::uvec idx = arma::find(nofz > 0.999e-8*nofz.max());
      redshift.shear_zdist_zmin[k] = z_v[idx(0)];
      redshift.shear_zdist_zmax[k] = z_v[idx(idx.n_elem-1)];
    }
  
    // READ THE N(Z) FILE ENDS ------------
    if (redshift.shear_zdist_zmax_all < redshift.shear_zdist_zmax[Ntomo-1] || 
        redshift.shear_zdist_zmin_all > redshift.shear_zdist_zmin[0]) 
    {
      spdlog::critical("zhisto_min = {},zhisto_max = {}", 
                       redshift.shear_zdist_zmin_all, 
                       redshift.shear_zdist_zmax_all);
      
      spdlog::critical("shear_zdist_zmin[0] = {},"
                       " shear_zdist_zmax[redshift.shear_nbin-1] = {}", 
                       redshift.shear_zdist_zmin[0], 
                       redshift.shear_zdist_zmax[Ntomo-1]);
      
      spdlog::critical("%s n(z) file incompatible with tomo.shear bin choice", 
                       multihisto_file);
      exit(1);
    } 

    zdistr_photoz(0.1, 0); // init static variables

    for (int k=0; k<Ntomo; k++)
    {
      spdlog::debug("{}: bin {} - {} = {}.", "init_source_sample", 
        k, "<z_s>", zmean_source(k));
    }

    redshift.random_shear = RandomNumber::get_instance().get();
  }

  spdlog::debug("{}: Ends", "init_source_sample");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_ntomo_powerspectra()
{
  if (redshift.shear_nbin == 0)
  {
    spdlog::critical(
      "{}: {} not set prior to this function call", 
      "init_ntomo_powerspectra", "redshift.shear_nbin");
    exit(1);
  }
  if (redshift.clustering_nbin == 0)
  {
    spdlog::critical(
      "{}: {} not set prior to this function call", 
      "init_ntomo_powerspectra", "redshift.clustering_nbin");
    exit(1);
  }

  tomo.shear_Npowerspectra = redshift.shear_nbin * (redshift.shear_nbin + 1) / 2;

  int n = 0;
  for (int i=0; i<redshift.clustering_nbin; i++)
    for (int j=0; j<redshift.shear_nbin; j++)
      n += test_zoverlap(i, j);
  tomo.ggl_Npowerspectra = n;

  tomo.clustering_Npowerspectra = redshift.clustering_nbin;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_redshift_distributions_from_files(
  std::string lens_multihisto_file, const int lens_ntomo,
  std::string source_multihisto_file, const int source_ntomo)
{
  init_lens_sample(lens_multihisto_file, lens_ntomo);
  
  init_source_sample(source_multihisto_file, source_ntomo);

  init_ntomo_powerspectra();
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_nuisance_shear_photoz(arma::Col<double> SP)
{
  spdlog::debug("{}: Begins", "set_nuisance_shear_photoz");

  if (redshift.shear_nbin == 0)
  {
    spdlog::critical("{}: {} = 0 is invalid",
      "set_nuisance_shear_photoz", "shear_Nbin");
    exit(1);
  }
  if (redshift.shear_nbin != static_cast<int>(SP.n_elem))
  {
    spdlog::critical(
      "{}: incompatible input w/ size = {} (!= {})",
      "set_nuisance_shear_photoz", SP.n_elem, redshift.shear_nbin);
    exit(1);
  }

  int cache_update = 0;
  for (int i=0; i<redshift.shear_nbin; i++)
  {
    if (fdiff(nuisance.photoz[0][0][i], SP(i)))
    {
      cache_update = 1;
      nuisance.photoz[0][0][i] = SP(i);
    } 
  }

  if (1 == cache_update)
    nuisance.random_photoz_shear = RandomNumber::get_instance().get();

  spdlog::debug("{}: Ends", "set_nuisance_shear_photoz");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_nuisance_clustering_photoz(arma::Col<double> CP)
{
  spdlog::debug("{}: Begins", "set_nuisance_clustering_photoz");

  if (redshift.clustering_nbin == 0)
  {
    spdlog::critical("{}: {} = 0 is invalid",
      "set_nuisance_clustering_photoz", "clustering_Nbin");
    exit(1);
  }
  if (redshift.clustering_nbin != static_cast<int>(CP.n_elem))
  {
    spdlog::critical(
      "{}: incompatible input w/ size = {} (!= {})",
      "set_nuisance_clustering_photoz", CP.n_elem, redshift.clustering_nbin);
    exit(1);
  }

  int cache_update = 0;
  for (int i=0; i<redshift.clustering_nbin; i++)
  {
    if (fdiff(nuisance.photoz[1][0][i], CP(i)))
    { 
      cache_update = 1;
      nuisance.photoz[1][0][i] = CP(i);
    }
  }

  if(1 == cache_update)
    nuisance.random_photoz_clustering = RandomNumber::get_instance().get();
  
  spdlog::debug("{}: Ends", "set_nuisance_clustering_photoz");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*
void set_nuisance_clustering_photoz_stretch(vector CP)
{
  spdlog::debug("{}: Begins", "set_nuisance_clustering_photoz_stretch");

  if (redshift.clustering_nbin == 0)
  {
    spdlog::critical(
      "{}: {} = 0 is invalid",
      "set_nuisance_clustering_photoz_stretch",
      "clustering_Nbin"
    );
    exit(1);
  }

  if (redshift.clustering_nbin != static_cast<int>(CP.n_elem))
  {
    spdlog::critical(
      "{}: incompatible input w/ size = {} (!= {})",
      "set_nuisance_clustering_photoz_stretch",
      CP.n_elem,
      redshift.clustering_nbin
    );
    exit(1);
  }

  for (int i=0; i<redshift.clustering_nbin; i++)
  {
    nuisance.stretch_zphot_clustering[i] = CP(i);
  }

  spdlog::debug("{}: Ends", "set_nuisance_clustering_photoz_stretch");
}
*/

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
} // namespace cosmolike_interface