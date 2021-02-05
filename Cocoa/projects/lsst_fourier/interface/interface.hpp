#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <carma/carma.h>

#ifndef __COSMOLIKE_INTERFACE_HPP
#define __COSMOLIKE_INTERFACE_HPP

// --- Auxiliary Code ---  
namespace interface_mpp_aux {

class RandomNumber {  
// Singleton Class that holds a random number generator  
public:    
  static RandomNumber& get_instance() {      
    static RandomNumber instance;
		return instance; 
  }	
  double get() {
    return dist_(mt_);
  }
	~RandomNumber() = default;
protected:
  std::random_device rd_;
  std::mt19937 mt_;
  std::uniform_real_distribution<double> dist_;
private:   
  RandomNumber() :
    rd_(),
    mt_(rd_()),
    dist_(0.0,1.0) {  
	};
  RandomNumber(RandomNumber const&) = delete; 
};

class RealData {    
public:    
  static RealData& get_instance() {      
    static RealData instance;
    return instance; 
  } 
  ~RealData() = default;

  void set_data(std::string DATA); 

  void set_mask(std::string MASK);

  void set_inv_cov(std::string COV);

  arma::Col<double> get_mask() const;

  arma::Col<double> get_data() const;

  int get_mask(const int ci) const;

  double get_data(const int ci) const; 
  
  double get_inv_cov(const int ci, const int cj) const; 

  double get_chi2(std::vector<double> datavector) const;

  bool is_mask_set() const;

  bool is_data_set() const;

  bool is_inv_cov_set() const;

private:
  bool is_mask_set_ = false;
  bool is_data_set_ = false; 
  bool is_inv_cov_set_ = false;
  std::string mask_filename_;
  std::string cov_filename_;
  std::string data_filename_;
  arma::Col<double> data_;
  arma::Col<double> mask_;
  arma::Mat<double> inv_cov_mask_;
  RealData() = default;
  RealData(RealData const&) = delete; 
};

arma::Mat<double> read_table(const std::string file_name);

// https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
    almost_equal(T x, T y, int ulp = 100) {
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::fabs(x-y) <= std::numeric_limits<T>::epsilon() * std::fabs(x+y) * ulp
        // unless the result is subnormal
        || std::fabs(x-y) < std::numeric_limits<T>::min();
}

}  // namespace interface_mpp_aux
#endif // HEADER GUARD