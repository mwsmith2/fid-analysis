/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   11/02/14

Detail: This is a new test program for my Fid libraries 

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using std::cout;
using std::endl;

//--- project includes ------------------------------------------------------//
#include "fid.h"

using namespace fid;

int main(int argc, char **argv)
{
  // run parameters 
  std::string data_file = "data/ideal_fid_sweep_data.csv";
  std::string sim_conf_file = "runtime/sim_params.json";
  int num_fids = 100;

  // initialize the configurable parameters
  if (argc > 1) load_params(argv[1]);

  // some necessary parameters
  std::vector<double> wf;
  std::vector<double> tm;

  FidFactory ff;

  double final_time = sim::start_time + sim::num_samples*sim::delta_time;
  tm = construct_range(sim::start_time, sim::delta_time, final_time);

  std::ofstream out;
  out.precision(10);
  out.open(data_file);

  std::vector<double> freqs;
  std::vector<double> phases;
  std::vector<double> snrs;
  std::vector<double> vals;

  // Get the range to sweep over.
  boost::property_tree::ptree conf;
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(sim_conf_file, conf);

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

  // begin sweeps
  for (auto f: freqs){

    for (auto p: phases){

      for (auto s: snrs){

        if (freqs.size() > 1) cout << "Running for frequency " << f;
        if (phases.size() > 1) cout << ", phase " << p;
        if (snrs.size() > 1) cout << ", signal-to-noise " << s;
        cout << endl;

        for (int i = 0; i < num_fids; i++){

          if (freqs.size() > 1) out << f << ", ";
          if (phases.size() > 1) out << p << ", ";
          if (snrs.size() > 1) out << s << ", ";

          sim::freq_larmor = f;
          sim::mixdown_phi = p;
          sim::snr = s;

          ff.IdealFid(wf, tm, true);
          Fid my_fid(wf, tm);

          my_fid.WriteFreqCsv(out);

        } // num_fids

      } // snr

    } // phi

  } // freq

  out.close();
  return 0;
}
