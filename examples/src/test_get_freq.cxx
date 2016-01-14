/*===========================================================================*\

author: Matthias W. Smith
email:  mwsmith2@uw.edu
date:   2014/12/02

notes: Tests the functionality of selecting methods through GetFreq

usage:

./test_get_freq

\*===========================================================================*/


//--- std includes ----------------------------------------------------------//

//--- other includes --------------------------------------------------------//

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
  std::vector<double> wf;
  std::vector<double> tm;

  // Create random number engine/distribution.
  std::default_random_engine gen;
  std::uniform_real_distribution<double> rand_flat_dist(35.0, 40.0);

  double final_time = sim::start_time + sim::num_samples * sim::sample_time;
  tm = construct_range(sim::start_time, final_time, sim::sample_time);

  FidFactory ff;
  ff.SetSignalToNoise(100 * 100);

  for (int i = 0; i < nfids; ++i) {
      
    // Make ideal Fid waveform
    ff.SetFidFreq(rand_flat_dist(gen));
    ff.IdealFid(wf, tm);

    Fid myfid(wf, tm);

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
