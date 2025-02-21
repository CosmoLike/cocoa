#include <string>
#include <vector>
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
#include "cosmolike/structs.h"
#include "cosmolike/interface_aux.hpp"

namespace cosmolike_interface
{

arma::Mat<double> read_table(const std::string file_name)
{
  std::ifstream input_file(file_name);

  if (!input_file.is_open())
  {
    const std::string etmp = "read_table";
    spdlog::critical("{}: file {} cannot be opened", etmp, file_name);
    exit(1);
  }

  // --------------------------------------------------------
  // Read the entire file into memory
  // --------------------------------------------------------
  std::string tmp;
  input_file.seekg(0,std::ios::end);
  tmp.resize(static_cast<size_t>(input_file.tellg()));
  input_file.seekg(0,std::ios::beg);
  input_file.read(&tmp[0],tmp.size());
  input_file.close();
  
  if(tmp.empty())
  {
    spdlog::critical("{}: file {} is empty", "read_table", file_name);
    exit(1);
  }
  
  // --------------------------------------------------------
  // Second: Split file into lines
  // --------------------------------------------------------
  std::vector<std::string> lines;
  lines.reserve(50000);

  boost::trim_if(tmp, boost::is_any_of("\t "));
  boost::trim_if(tmp, boost::is_any_of("\n"));
  boost::split(lines, tmp,boost::is_any_of("\n"), boost::token_compress_on);
  
  // Erase comment/blank lines
  auto check = [](std::string mystr) -> bool
  {
    return boost::starts_with(mystr, "#");
  };
  lines.erase(std::remove_if(lines.begin(), lines.end(), check), lines.end());
  
  // --------------------------------------------------------
  // Third: Split line into words
  // --------------------------------------------------------
  arma::Mat<double> result;
  size_t ncols = 0;
  { // first line
    std::vector<std::string> words;
    words.reserve(100);
    
    boost::split(words,lines[0],boost::is_any_of(" \t"),boost::token_compress_on);
    
    ncols = words.size();
    
    result.set_size(lines.size(), ncols);
    for (size_t j=0; j<ncols; j++)
      result(0,j) = std::stod(words[j]);
  }

  #pragma omp parallel for
  for (size_t i=1; i<lines.size(); i++)
  {
    std::vector<std::string> words;
    
    boost::split(
      words, 
      lines[i], 
      boost::is_any_of(" \t"),
      boost::token_compress_on
    );
    
    if (words.size() != ncols)
    {
      const std::string etmp = "read_table";
      spdlog::critical("{}: file {} is not well formatted", etmp, file_name);
      exit(1);
    }
    
    for (size_t j=0; j<ncols; j++)
      result(i,j) = std::stod(words[j]);
  };
  
  return result;
}

} // namespace cosmolike_interface