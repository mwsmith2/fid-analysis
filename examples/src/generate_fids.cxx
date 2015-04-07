/*===========================================================================*\

author: Matthias W. Smith
email:  mwsmith2@uw.edu
file:   generate_fids.cxx

about: This is a new test program for my FID libraries 

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <fstream>
#include <vector>
#include <map>

//--- other includes --------------------------------------------------------//
#include "TFile.h"
#include "TTree.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

//--- project includes ------------------------------------------------------//
#include "fid.h"

using namespace fid;

int main(int argc, char **argv)
{
  // initialize the configurable parameters
  if (argc > 1) load_params(argv[1]); 

  // some necessary parameters
  vec wf;
  double final_time = sim::start_time + sim::num_samples*sim::delta_time;
  vec tm(construct_range(sim::start_time, sim::delta_time, final_time));

  vec freqs;
  vec phases;
  vec snrs;
  vec vals;

  boost::property_tree::ptree conf;
  boost::property_tree::ptree pt;
  boost::property_tree::read_json("sim_params.json", conf);

  // Get the range to sweep over.
  pt = conf.get_child("sweep.freq");
  if (pt.get<bool>("in_use")) {

    vals.resize(0);
    for (auto &val : pt.get_child("values")) {
      vals.push_back(val.second.get_value<double>());
    }

    freqs = construct_range(vals);

  } else {

    freqs.push_back(sim::freq_larmor - sim::freq_ref);
  }

  pt = conf.get_child("sweep.phase");
  if (pt.get<bool>("in_use")) {

    vals.resize(0);
    for (auto &val : pt.get_child("values")) {
      vals.push_back(val.second.get_value<double>());
    }

    phases = construct_linspace(vals);

  } else {

    phases.push_back(sim::mixdown_phi);
  }

  pt = conf.get_child("sweep.snr");
  if (pt.get<bool>("in_use")) {
   
    vals.resize(0);
    for (auto &val : pt.get_child("values")) {
      vals.push_back(val.second.get_value<double>());
    }

    snrs = construct_range(vals);

  } else {

    snrs.push_back(sim::snr);
  }

  // Make FidFactory
  FidFactory ff;

  // begin sweeps
  for (auto f: freqs){

    for (auto p: phases){

      for (auto s: snrs){

        if (freqs.size() > 1) cout << "Running for frequency " << f;
        if (phases.size() > 1) cout << ", phase " << p;
        if (snrs.size() > 1) cout << ", signal-to-noise " << s;
        cout << endl;

        ff.SimulateFid(wf, tm);

        FID my_fid(wf, tm);
        draw_fid(my_fid, string("data/fig/test.pdf"), string("Test FID"));

      } // snr

    } // phi

  } // freq

  return 0;
}
