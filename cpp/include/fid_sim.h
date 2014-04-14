#ifndef FID_ANALYSIS_INCLUDE_FID_SIM_H_
#define FID_ANALYSIS_INCLUDE_FID_SIM_H_

/*===========================================================================*\

author: Matthias W. Smith
email: mwmsith2@uw.edu

notes:
  The class, FidFactory simulates FIDs using the Bloch equations and 
  numerical integration.

initalization

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
#include "TFile.h"
#include "TTree.h"

//--- project includes ------------------------------------------------------//
#include "fid_params.h"

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

  // member methods
  void SimulateFid(vec& wf, vec& tm);
  void IdealFid(vec& wf, vec& tm);

 private:

  static double ti_;
  static double tf_;
  static double dt_;
  static double x_;

  static vec s_;
  static vec spin_;
  static vec time_vec_;
  static vec spin_sum_;
  static vec cos_cache_;
  static vec gradient_;
  
  // private member functions
  // Define the cross product for 3-vectors.
  static inline void Cross(const vec& u, const vec& v, vec& res)
  {
    res[0] = u[1] * v[2] - u[2] * v[1];
    res[1] = u[2] * v[0] - u[0] * v[2];
    res[2] = u[0] * v[1] - u[1] * v[0];
  } 

  // Low pass filter to extract mixed down signal.
  static vector<double> LowPassFilter(vector<double>& s);

  // Function which returns time dependent Bfield.
  static vec Bfield(const double& x, const double& t);

  // The time evolution equation for the fields.
  static void Bloch(vec const &s, vec &dsdt, double t);

  // The integration monitor function
  static void Printer(vec const &s , double t);

}; // FidFactory

class GradientFidFactory
{
 public:

  // ctor
  GradientFidFactory();

  // dtor
  ~GradientFidFactory();

  // member methods
  void ConstructFid(const vec& gradient, vec& wf);

 private:

  int num_sim_fids_;
  int zero_idx_;
  double d_grad_;

  TFile *pf_fid_;
  TTree *pt_fid_;
  vector<Double_t> wf_;

  int GetTreeIndex(double grad_strength);
};

} // ::fid

#endif