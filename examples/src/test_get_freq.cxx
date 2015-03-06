/*===========================================================================*\

author: Matthias W. Smith
email:  mwsmith2@uw.edu
date:   2014/12/02

notes: Tests the frequency extraction ability at different sampling 
      frequencies.

usage:

./sampling_rate_test <output>

The parameters in brackets are optional.  The default output is 
sampling_rate_data.csv.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <fstream>

//--- other includes --------------------------------------------------------//
#include <boost/filesystem.hpp>

//--- project includes ------------------------------------------------------//
#include "fid.h"

using namespace fid;
using std::cout;
using std::endl;

int main(int argc, char **argv)
{
  int nfids = 10;
  cout.precision(12);

  // some necessary parameters
  vec wf;
  vec tm;

  // Create random number engine/distribution.
  std::default_random_engine gen;
  std::uniform_real_distribution<double> rand_flat_dist(35.0, 40.0);

  double final_time = sim::start_time + sim::num_samples * sim::delta_time;
  tm = construct_range(sim::start_time, final_time, sim::delta_time);

  for (int i = 0; i < nfids; ++i) {
      
    // Make ideal FID waveform
    double freq = rand_flat_dist(gen);
    ideal_fid(wf, tm, freq, 0.0, 100 * 100);

    FID myfid(wf, tm);

    cout << "zc: " << myfid.GetFreq("zc") << " ";
    cout << myfid.CalcZeroCountFreq() << endl;

    cout << "cn: " << myfid.GetFreq("cn") << " ";
    cout << myfid.CalcCentroidFreq() << endl;

    cout << "lz: " << myfid.GetFreq("lz") << " ";
    cout << myfid.CalcLorentzianFreq() << endl;

    cout << "ph: " << myfid.GetFreq("ph") << " ";
    cout << myfid.CalcPhaseFreq() << endl;

    cout << "sn: " << myfid.GetFreq("sn") << " ";
    cout << myfid.CalcSinusoidFreq() << endl;
  }

  return 0;
}