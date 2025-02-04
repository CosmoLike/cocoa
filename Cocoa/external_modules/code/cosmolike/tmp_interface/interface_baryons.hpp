#include <carma.h>
#include <armadillo>
#include <map>
#include <string>
#include <tuple>

#ifndef __COSMOLIKE_GENERIC_INTERFACE_BARYONS_HPP
#define __COSMOLIKE_GENERIC_INTERFACE_BARYONS_HPP

namespace cosmolike_interface
{

class BaryonScenario
{
public:
  static BaryonScenario& get_instance()
  {
    static BaryonScenario instance;
    return instance;
  }
  ~BaryonScenario() = default;

  int nscenarios() const;

  bool is_pcs_set() const;

  bool is_scenarios_set() const;

  void set_scenarios(std::string data_sims, std::string scenarios);

  void set_scenarios(std::string scenarios);

  void set_pcs(arma::Mat<double> eigenvectors);

  std::string get_scenario(const int i) const;

  std::tuple<std::string,int> select_baryons_sim(const std::string scenario);

  arma::Mat<double> get_pcs() const;

  double get_pcs(const int ci, const int cj) const;

private:
  bool is_pcs_set_;

  bool is_scenarios_set_;

  int nscenarios_;

  std::map<int, std::string> scenarios_;

  arma::Mat<double> eigenvectors_;

  BaryonScenario() = default;
  BaryonScenario(BaryonScenario const&) = delete;
};

void init_baryons_contamination(std::string sim); // old API

#ifdef HDF5LIB
void init_baryons_contamination(std::string sim, std::string all_sims_hdf5_file);
#endif

double get_baryon_power_spectrum_ratio(const double log10k, const double a);

arma::Mat<double> compute_baryon_pcas_3x2pt_real(arma::Col<int>::fixed<3> order);

arma::Mat<double> compute_baryon_pcas_6x2pt(arma::Col<int>::fixed<6> order);

} // namespace cosmolike_interface