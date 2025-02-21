#include <carma.h>
#include <armadillo>

#ifndef __COSMOLIKE_GENERIC_INTERFACE_MASKED_DATA_HPP
#define __COSMOLIKE_GENERIC_INTERFACE_MASKED_DATA_HPP

class IP
{
  public:
    static IP& get_instance()
    {
      static IP instance;
      return instance;
    }

    ~IP() = default;

    // ----------------------------------------------

    bool is_mask_set() const;
    bool is_data_set() const;
    bool is_inv_cov_set() const;

    // ----------------------------------------------

    void set_data(std::string datavector_filename);

    // 3x2pt
    void set_mask(std::string mask_filename, arma::Col<int>::fixed<3> order);

    // 6x2pt
    void set_mask(std::string mask_filename, arma::Col<int>::fixed<6> order);

    void set_inv_cov(std::string covariance_filename);

    // ----------------------------------------------

    int get_mask(const int ci) const;

    double get_dv_masked(const int ci) const;

    double get_inv_cov_masked(const int ci, const int cj) const;

    int get_index_sqzd(const int ci) const;

    double get_dv_masked_sqzd(const int ci) const;

    double get_inv_cov_masked_sqzd(const int ci, const int cj) const;

    arma::Col<double> expand_theory_data_vector_from_sqzd(arma::Col<double>) const;

    arma::Col<double> sqzd_theory_data_vector(arma::Col<double>) const;

    double get_chi2(arma::Col<double> datavector) const;

    // ----------------------------------------------

    int get_ndata() const;

    arma::Col<int> get_mask() const;

    arma::Col<double> get_dv_masked() const;

    arma::Mat<double> get_cov_masked() const;

    arma::Mat<double> get_inv_cov_masked() const;

    int get_ndata_sqzd() const;

    arma::Col<double> get_dv_masked_sqzd() const;

    arma::Mat<double> get_cov_masked_sqzd() const;

    arma::Mat<double> get_inv_cov_masked_sqzd() const;

  private:

    bool is_mask_set_ = false;
    
    bool is_data_set_ = false;
    
    bool is_inv_cov_set_ = false;
    
    int ndata_ = 0;
    
    int ndata_sqzd_ = 0;
    
    std::string mask_filename_;
    
    std::string cov_filename_;
    
    std::string data_filename_;
    
    arma::Col<int> mask_;

    arma::Col<double> data_masked_;
    
    arma::Mat<double> cov_masked_;

    arma::Col<int> index_sqzd_;
    
    arma::Mat<double> inv_cov_masked_;
    
    arma::Col<double> data_masked_sqzd_;

    arma::Mat<double> cov_masked_sqzd_; 

    arma::Mat<double> inv_cov_masked_sqzd_;
    
    IP() = default;
    IP(IP const&) = delete;
};

void init_binning_fourier(const int Ncl, const int lmin, const int lmax);

void init_binning_real_space(
    const int Ntheta, 
    const double theta_min_arcmin, 
    const double theta_max_arcmin
  );

void init_data_3x2pt_real_space(
    std::string cov, 
    std::string mask, 
    std::string data,
    arma::Col<int>::fixed<3> order
  );

void init_data_6x2pt_real_space(
    std::string cov, 
    std::string mask, 
    std::string data,
    arma::Col<int>::fixed<6> order
  );

void init_probes(std::string possible_probes);

} // namespace cosmolike_interface