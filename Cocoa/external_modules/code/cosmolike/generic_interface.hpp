#include <carma.h>
#include <armadillo>
#include <map>

#ifndef __COSMOLIKE_GENERIC_INTERFACE_HPP
#define __COSMOLIKE_GENERIC_INTERFACE_HPP

namespace cosmolike_interface
{

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Class RandomNumber
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

class RandomNumber
{ // Singleton Class that holds a random number generator
  public:
    static RandomNumber& get_instance()
    {
      static RandomNumber instance;
      return instance;
    }
    ~RandomNumber() = default;

    double get()
    {
      return dist_(mt_);
    }

  protected:
    std::random_device rd_;
    std::mt19937 mt_;
    std::uniform_real_distribution<double> dist_;
  
  private:
    RandomNumber() :
      rd_(),
      mt_(rd_()),
      dist_(0.0, 1.0) {
      };
  
    RandomNumber(RandomNumber const&) = delete;
};

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Class IP (InterfaceProducts)
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

class IP
{ // InterfaceProducts: Singleton Class that holds data vector, covariance...
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

    //void set_PMmarg(std::string U_PMmarg_file);

    // ----------------------------------------------

    int get_mask(const int ci) const;

    double get_dv_masked(const int ci) const;

    double get_inv_cov_masked(const int ci, const int cj) const;

    int get_index_sqzd(const int ci) const;

    double get_dv_masked_sqzd(const int ci) const;

    double get_inv_cov_masked_sqzd(const int ci, const int cj) const;

    arma::Col<double> expand_theory_data_vector_from_sqzd(arma::Col<double>) const;

    arma::Col<double> sqzd_theory_data_vector(arma::Col<double>) const;

    double get_chi2(std::vector<double> datavector) const;

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
    
    int ndata_sqzd_ = 0;                      // for baryon PC, reduced dim
    
    std::string mask_filename_;
    
    std::string cov_filename_;
    
    std::string data_filename_;
    
    arma::Col<int> mask_;

    arma::Col<double> data_masked_;
    
    arma::Mat<double> cov_masked_;

    arma::Col<int> index_sqzd_;
    
    arma::Mat<double> inv_cov_masked_;
    
    arma::Col<double> data_masked_sqzd_;     // for baryon project, reduced dim

    arma::Mat<double> cov_masked_sqzd_;      // for baryon project, reduced dim

    arma::Mat<double> inv_cov_masked_sqzd_;  // for baryon project, reduced dim
    
    IP() = default;
    IP(IP const&) = delete;
};

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Class IPCMB (InterfaceProducts - 6x2pt (includes CMB lensing))
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

class IPCMB
{
  public:
    static IPCMB& get_instance()
    {
      static IPCMB instance;
      return instance;
    }

    ~IPCMB() = default;

    // ----------------------------------------------

    bool is_cmb_binmat_set() const;

    bool is_cmb_offset_set() const;

    // ----------------------------------------------

    void set_cmb_binning_mat(std::string cmb_binned_matrix_filename);

    void set_cmb_theory_offset(std::string cmb_theory_offset_filename);

    // ----------------------------------------------

    double get_binning_matrix_with_correction(const int ci, const int cj) const;

    double get_cmb_theory_offset(const int ci) const;

    // ----------------------------------------------

    arma::Mat<double> get_binning_matrix_with_correction() const;

    arma::Mat<double> get_cmb_theory_offset() const;

  private:   
    bool is_cmb_binmat_set_ = false;
    
    bool is_cmb_offset_set_ = false;

    int nbp_ = 0;
    
    int ncl_ = 0;

    std::string binmat_filename_;                // LSS x CMB
    
    std::string offset_filename_;                // LSS x CMB
    
    arma::Col<double> cmb_theory_offset_;                   // LSS x CMB, see Eqn 35 Planck 2018 VIII

    arma::Mat<double> cmb_binning_matrix_with_correction_;  // LSS x CMB, see Eqn 35 Planck 2018 VIII
    
    IPCMB() = default;
    IPCMB(IP const&) = delete;
};

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Class PointMass
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

class PointMass
{// Singleton Class that Evaluate Point Mass Marginalization
  public:
    static PointMass& get_instance()
    {
      static PointMass instance;
      return instance;
    }
    ~PointMass() = default;

    void set_pm_vector(std::vector<double> pm);

    std::vector<double> get_pm_vector() const;

    double get_pm(const int zl, const int zs, const double theta) const;

  private:
    std::vector<double> pm_;

    PointMass() = default;
    
    PointMass(PointMass const&) = delete;
};

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Class BaryonScenario
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

class BaryonScenario
{ // Singleton Class that map Baryon Scenario (integer to name)
  public:
    static BaryonScenario& get_instance()
    {
      static BaryonScenario instance;
      return instance;
    }
    ~BaryonScenario() = default;

    int nscenarios() const;

    void set_scenarios(std::string data_sims, std::string scenarios);

    void set_scenarios(std::string scenarios);

    std::string get_scenario(const int i) const;

    std::tuple<std::string,int> select_baryons_sim(const std::string scenario);

  private:
    int nscenarios_;
    std::map<int, std::string> scenarios_;

    BaryonScenario() = default;
    BaryonScenario(BaryonScenario const&) = delete;
};

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// AUX FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Mat<double> read_table(const std::string file_name);

// https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type 
almost_equal(T x, T y, int ulp = 100)
{
  // the machine epsilon has to be scaled to the magnitude of the values used
  // and multiplied by the desired precision in ULPs (units in the last place)
  return std::fabs(x-y) <= std::numeric_limits<T>::epsilon() * std::fabs(x+y) * ulp
      // unless the result is subnormal
      || std::fabs(x-y) < std::numeric_limits<T>::min();
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GLOBAL INIT FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_accuracy_boost(
    const double accuracy_boost, 
    const double sampling_boost,
    const int integration_accuracy
  );

void init_binning_fourier(
    const int Ncl, 
    const int lmin, 
    const int lmax
  );

void init_binning_real_space(
    const int Ntheta, 
    const double theta_min_arcmin, 
    const double theta_max_arcmin
  );

void init_binning_cmb_bandpower(
    const int Nbandpower, 
    const int lmin, 
    const int lmax
  );

void init_cosmo_runmode(
    const bool is_linear
  );

void init_cmb(
    const double lmin_kappa_cmb, 
    const double lmax_kappa_cmb, 
    const double fwhm, 
    std::string pixwin_file
  );

void init_cmb_bandpower(
    const int is_cmb_bandpower, 
    const int is_cmb_kkkk_covariance_from_simulation, 
    const double alpha
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

void init_data_vector_size(
    arma::Col<int>::fixed<6> exclude
  );

void init_data_vector_size_real_space(
    arma::Col<int>::fixed<6> exclude
  );

void init_data_vector_size_3x2pt_real_space(
  );

void init_data_vector_size_6x2pt_real_space(
  );

void init_distances(
    std::vector<double> io_z, 
    std::vector<double> io_chi
  );

void init_growth(
    std::vector<double> io_z, 
    std::vector<double> io_G
  );

void init_IA(
    const int IA_MODEL, 
    const int IA_REDSHIFT_EVOL
  );

void init_lens_sample(
    std::string multihisto_file, 
    const int Ntomo
  );

void init_linear_power_spectrum(
    std::vector<double> io_log10k,
    std::vector<double> io_z, 
    std::vector<double> io_lnP
  );

void init_non_linear_power_spectrum(
    std::vector<double> io_log10k,
    std::vector<double> io_z, 
    std::vector<double> io_lnP
  );

void init_probes(
    std::string possible_probes
  );

void initial_setup(
  );

void init_source_sample(
    std::string multihisto_file, 
    const int Ntomo
  );

void init_survey(
    std::string surveyname, 
    double area, 
    double sigma_e
  );

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GLOBAL SET FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_pm(
    std::vector<double> PM
  );

void set_nuisance_bias(
    std::vector<double> B1, 
    std::vector<double> B2, 
    std::vector<double> B_MAG
  );

void set_nuisance_clustering_photoz(
    std::vector<double> CP
  );

void set_nuisance_IA(
    std::vector<double> A1, 
    std::vector<double> A2,
    std::vector<double> BTA
  );

void set_nuisance_magnification_bias(
  std::vector<double> B_MAG
  );

void set_nuisance_nonlinear_bias(
    std::vector<double> B1,
    std::vector<double> B2
  );

void set_nuisance_shear_calib(
    std::vector<double> M
  );

void set_nuisance_shear_photoz(
    std::vector<double> SP
  );

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GLOBAL COMPUTE FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Mat<double> compute_baryon_pcas_3x2pt(
    arma::Col<int>::fixed<3> order
  );

arma::Mat<double> compute_baryon_pcas_6x2pt(
    arma::Col<int>::fixed<6> order
  );

double compute_chi2(
    std::vector<double> datavector
  );

std::vector<double> compute_data_vector_6x2pt_masked_any_order(
    arma::Col<int>::fixed<6> order
  );

std::vector<double> compute_data_vector_3x2pt_masked_any_order(
    arma::Col<int>::fixed<3> order
  );

}  // namespace cosmolike_interface
#endif // HEADER GUARD