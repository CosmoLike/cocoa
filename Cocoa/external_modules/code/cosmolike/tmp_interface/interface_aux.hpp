#include <carma.h>
#include <armadillo>
#include <map>

#ifndef __COSMOLIKE_GENERIC_INTERFACE_AUX_HPP
#define __COSMOLIKE_GENERIC_INTERFACE_AUX_HPP

namespace cosmolike_interface
{

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

} // namespace cosmolike_interface