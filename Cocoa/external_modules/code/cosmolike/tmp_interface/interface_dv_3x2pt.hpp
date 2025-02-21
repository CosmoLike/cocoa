#include <carma.h>
#include <armadillo>

#ifndef __COSMOLIKE_GENERIC_INTERFACE_DV_3x2PT_REAL_HPP
#define __COSMOLIKE_GENERIC_INTERFACE_DV_3x2PT_REAL_HPP

namespace cosmolike_interface
{

class PointMass
{
public:
  static PointMass& get_instance()
  {
    static PointMass instance;
    return instance;
  }
  ~PointMass() = default;

  void set_pm_vector(arma::Col<double> pm);

  arma::Col<double> get_pm_vector() const;

  double get_pm(const int zl, const int zs, const double theta) const;

private:
  arma::Col<double> pm_;
  PointMass() = default;
  PointMass(PointMass const&) = delete;
};

void init_data_vector_size_3x2pt_real_space();

void init_data_vector_size_3x2pt_fourier_space();

arma::Col<double> compute_data_vector_3x2pt_real_masked_any_order(
    arma::Col<int>::fixed<3> order
  );

arma::Col<double> compute_data_vector_3x2pt_real_masked_any_order(
    arma::Col<double> Q,                // PC amplitudes
    arma::Col<int>::fixed<3> order
  );

} // namespace cosmolike_interface