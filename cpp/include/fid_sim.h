#ifndef FID_ANALYSIS_INCLUDE_FID_SIM_H_
#define FID_ANSLYSIS_INCLUDE_FID_SIM_H_

/*===========================================================================*\

Author: Matthias W. Smith
Email : mwmsith2@uw.edu
Date  : 14/2/14

Notes :
  This simulates FID using the Bloch equations and numerical integration.
  There is a default parameter set which can be overridden by setting up a 
  config file and passing it as the command line argument.  If no argument is
  given, it loads a default config file called "sim_config.json".  If no config
  file is desired, then pass "none" to use the default params.

Example Usages:
  ./generate_fids none
  ./generate_fids
  ./generate_fids my_params.json

Dependencies:
  C++11 support
  Boost Libraries
  Armadillo Linalg Library

Compiling:
  I use clang++ on my mac, g++ also works though the -std=c++11 may need to be
  changed to -std=c++0x

  clang++ -std=c++11 -O4 generate_fids.cxx -o generate_fids

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include <algorithm>
#include <random>
#include <cmath>

//--- other includes --------------------------------------------------------//

#include <boost/numeric/odeint.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <armadillo>

//--- project includes ------------------------------------------------------//

/*--- Namespaces ------------------------------------------------------------*/

using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::vector;
using namespace std::chrono;
using namespace boost::numeric::odeint;
using namespace boost::property_tree;


namespace fid {

// A template function to handle vector addition.
template <typename T>
inline vector<T>& operator+(vector<T>& a, vector<T>& b)
{
  assert(a.size() == b.size());

  transform(a.begin(), a.end(), b.begin(), a.begin(), std::plus<T>());
  return a;
}

// A template function to handle vector subtraction.
template <typename T>
inline vector<T>& operator-(vector<T>& a, vector<T>& b)
{
  assert(a.size() == b.size());

  transform(a.begin(), a.end(), b.begin(), a.begin(), std::minus<T>());
  return a;
}

// A template function to handle vector multiplying with a scalar.
template <typename T>
inline vector<T>& operator*(T c, vector<T>& a)
{
  for (auto it = a.begin(); it != a.end(); ++it){
    *it = c * (*it);
  }

  return a;
}

class FidFactory
{
 public:

  // ctors
  FidFactory();

  // accessors

  // mutators

  // member methods
  void SimulateFid(vec& tm, vec& wf);
  void IdealFid(vec& tm, vec& wf);

 private:

  // simulation parameters
  int num_fids_;
  int num_points_;  // samples in each FID
  int reduction_; // reduction factor after mixing down the signal
  int num_steps_;   // number of steps in a gradient
  double dbfield_;
  double ti_;
  double tf_;
  double dt_;
  double t_total_;

  double freq_lar_;
  double freq_ref_;
  double snr_;

  //strength and duration pulse; omega_r * t_pulse = 1/4 for pi-pulse
  double omega_r_;
  double t_pulse_;
  double gamma_g_;
  double gamma_1_;
  double gamma_2_;

  double x_;

  vec spin_;
  vec time_vec_;
  vec spin_sum_;
  vec cos_cache_;
  vector<double> gradient_;
  
  // private member functions
  // Define the cross product for 3-vectors.
  inline void Cross(const vec& u, const vec& v, vec& res)
  {
    res[0] = u[1] * v[2] - u[2] * v[1];
    res[1] = u[2] * v[0] - u[0] * v[2];
    res[2] = u[0] * v[1] - u[1] * v[0];
  } 

  vector<double> LowPassFilter(vector<double>& s);
  void AddNoise(vector<double>& s);

  // Function which returns time dependent Bfield.
  vec Bfield(const double& x, const double& t);

  // The time evolution equation for the fields.
  void Bloch(vec const &s, vec &dsdt, double t);

}; // FidFactory

} // ::fid

#endif