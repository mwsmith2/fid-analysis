#include "fid_params.h"


namespace fid {

  // fid parameter namespace

  namespace sweep {

    // per configuration statistics
    int n_fids = 1000;

    // default parameters (s_ for single)
    double s_freq = 23.456789;
    double s_phi  = 0.0;
    double s_grad = 0.0;
    double s_s2n  = 100.0;

    // time vector variables
    int fid_length = 10000;
    double i_time = -1.0;
    double d_time = 0.001;

    // freqeuency sweep
    double i_freq = 23.0;
    double f_freq = 23.2;
    double d_freq = 0.001;

    // phase sweep
    int n_phi = 10;
    double i_phi = 0.0;
    double f_phi = 2 * M_PI;
    double d_phi = (f_phi - i_phi) / (n_phi - 1);

    // gradient sweep
    double i_grad = 0.0;
    double f_grad = 250.0;
    double d_grad = 10.0;

    // signal-to-noise sweep
    double i_s2n = 50.0;
    double f_s2n = 300.0;
    double d_s2n = 10.0;

    // data vectors
    vector<double> wf;
    vector<double> tm;
    
    // output stream
    std::ofstream out;

    // necessary initalizations
    inline void init_params(){

      wf.reserve(fid_length);
      tm.reserve(fid_length);

      // set precision for output streams
      cout.precision(10);
      cout.setf(std::ios::fixed, std::ios::floatfield);

      std::ofstream out;
      out.precision(10);
      out.setf(std::ios::fixed, std::ios::floatfield);
    }
  } // ::params

  void load_params(int argc, char **argv);

} // ::fid
