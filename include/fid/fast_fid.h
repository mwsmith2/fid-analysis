#ifndef LIBFID_FID_FAST_FID_H_
#define LIBFID_FID_FAST_FID_H_

/*============================================================================*\

author: Matthias W. Smith
email: mwmsith2@uw.edu

notes:

  Defines several frequency extraction and analysis methods for FIDs 
  as a class to encapsulate all the ideas.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <string>
#include <vector>
#include <numeric>
#include <random>
#include <complex>
#include <cmath>

//--- other includes --------------------------------------------------------//
#include <boost/algorithm/string.hpp>
#include <armadillo>
#include "TGraph.h"
#include "TF1.h"

//--- project includes ------------------------------------------------------//
#include "fid/base_fid.h"
#include "fid/params.h"
#include "fid/math.h"

namespace fid {

class FastFid : public BaseFid {

 public:
  
  // ctors
  FastFid(const std::string& fid_file);
  FastFid(const char* fid_file);
  FastFid(const std::vector<double>& wf);
  FastFid(const std::vector<double>& wf, const std::vector<double>& tm);

  // Simplified frequency extraction
  double CalcFreq();
  
 private:
  
}; // FastFid
 
} // fid

#endif
