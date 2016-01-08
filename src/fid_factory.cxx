#include "fid/fid_factory.h"

namespace fid {

//---------------------------------------------------------------------------//
//--- FID Factory -----------------------------------------------------------//
//---------------------------------------------------------------------------//


FidFactory::FidFactory()
{
  using std::cout;
  using std::endl;

  // Set the start/stop times.
  ti_ = sim::start_time;
  tf_ = ti_ + sim::delta_time * sim::num_samples;

  // Calculate the decimation ratio.
  sim_to_fid_ = (tf_ - ti_) / (sim::dt_integration * sim::num_samples) + 0.5; 

  // If zero, more less integration steps than sampling times were requested.
  if (sim_to_fid_ == 0) {

    cout << "WARNING: The given integration step was larger than the ";
    cout << "sampling time, so the sampling time, " << sim::delta_time;
    cout << ", will be used instead." << endl;

    sim_to_fid_ = 1;
    dt_ = sim::delta_time;

  } else {

    dt_ = sim::delta_time / sim_to_fid_;
  
    if (dt_ != sim::dt_integration) {
      cout << "WARNING: The given integration time step was not an even";
      cout << " divisor of the sampling rate, so it has been rounded to ";
      cout << dt_ << endl;
    }
  }

  // Set number of simulation samples.
  sim_length_ = sim_to_fid_ * sim::num_samples;
  printer_idx_ = 0;

  // open the default root file
  pf_fid_ = new TFile(grad::root_file.c_str());

  if (pf_fid_->IsOpen()) {
    pt_fid_ = (TTree *)pf_fid_->Get("t");
  
    wf_.resize(sim::num_samples);
    pt_fid_->SetBranchAddress(grad::fid_branch.c_str(), &wf_[0]);
    pt_fid_->GetEntry(0);

    num_sim_fids_ = pt_fid_->GetEntries();
    d_grad_ = (grad::max - grad::min) / num_sim_fids_;
    zero_idx_ = -1 * grad::min / d_grad_ + 0.5;
  }
}


FidFactory::~FidFactory()
{
  // closing the TFile will take care of the TTree pointer
  pf_fid_->Close();
  delete pf_fid_;
}


// Create an idealized FID with current Simulation parameters
void FidFactory::IdealFid(std::vector<double>& wf, 
                          std::vector<double>& tm, 
                          bool withnoise,
                          bool discretize)
{
  wf.reserve(tm.size());
  wf.resize(0);

  // Define the waveform
  double temp;
  double w = kTau * (sim::freq_larmor - sim::freq_ref);
  double phi = sim::mixdown_phi;
  double tau = 1.0 / sim::gamma_1;
  double amp = sim::amplitude;
  double base = sim::baseline;

  for (auto it = tm.begin(); it != tm.end(); ++it){

    if (*it >= sim::t_pulse){

      temp = amp * std::exp(-(*it - sim::t_pulse) / tau);
      temp *= std::sin((*it) * w + phi);
      wf.push_back(temp + base);

    } else {

      wf.push_back(base);

    }
  } 

  if (withnoise) addnoise(wf, sim::snr);

  if (discretize) floor(wf);
}


// Simulate an FID with the current class settings.
void FidFactory::SimulateFid(std::vector<double>& wf, 
                             std::vector<double>& tm, 
                             bool withnoise,
                             bool discretize)
{
  using namespace boost::numeric::odeint;
  using std::bind;
  using std::ref;
  using namespace std::placeholders; // _1, _2, _3

  // make sure memory is allocated for the final FIDs
  tm.reserve(sim::num_samples);
  wf.reserve(sim::num_samples);
  s_ = sim::spin_0; // The starting spin vector

  // Bind the member function and make a reference so it isn't copied.
  integrate_const(runge_kutta4<std::vector<double>>(), 
                  bind(&FidFactory::Bloch, ref(*this), _1, _2, _3), 
                  s_, 
                  ti_, 
                  tf_, 
                  dt_, 
                  bind(&FidFactory::Printer, ref(*this), _1, _2));

  // Set the results
  tm.resize(0);
  wf.resize(0);

  // Fill the FID vectors with simulation results.
  for (int i = 0; i < sim::num_samples; ++i) {
    tm.push_back(time_vec_[i * sim_to_fid_]);
    wf.push_back(spin_vec_[i * sim_to_fid_] + sim::baseline);
  }

  // Take care of optional effects.
  if (withnoise) addnoise(wf, sim::snr);
  if (discretize) floor(wf);
}


// Simulate a gradient FID by combining multiple simulated FIDs.
void FidFactory::GradientFid(const std::vector<double>& gradient, 
                             std::vector<double>& wf, 
                             bool withnoise,
                             bool discretize)
{
  // Make sure we have an example file draw from.
  if (!pf_fid_->IsOpen()) {
    std::cout << "No ROOT file loaded.  Cannot make gradient FIDs.";
    std::cout << std::endl;
    return;
  }

  // Find the appropriate FIDs and sum them
  wf.assign(sim::num_samples, 0.0);


  // Add each FID in the gradient.
  for (auto& g : gradient){

    pt_fid_->GetEntry(GetTreeIndex(g));

    for (uint i = 0; i < wf.size(); i++){
      wf[i] += wf_[i];
    }
  }

  // Calculate the proper amplitude and scale + offset the sim FID.
  double amp = sim::amplitude / gradient.size();

  for (auto& val : wf) {
    val = val * amp + sim::baseline;
  }

  // Take care of optional effects.
  if (withnoise) addnoise(wf, sim::snr);
  if (discretize) floor(wf);
}


void FidFactory::PrintDiagnosticInfo()
{
  using std::cout;
  using std::endl;

  cout << endl;
  cout << "Printing Diagnostic Info for FidFactory @" << this << endl;
  cout << "The time step, fid length: " << dt_ << ", ";
  cout << sim::num_samples << endl;
  cout << "The sim length, sim-to-fid: " << sim_length_ << ", ";
  cout << sim_to_fid_ << endl;
}


// Define the Bloch equation governing FID signals.
void FidFactory::Bloch(std::vector<double> const &s, 
                       std::vector<double> &dsdt, 
                       double t)
{
  // Again static to save time on memory allocations.
  static std::vector<double> bf = {{0., 0., 0.}};   // Bfield
  static std::vector<double> s1 = {{0., 0., 0.}};  // Cross product spin
  static std::vector<double> s2 = {{0., 0., 0.}};  // Relaxation spin

  // Update the Bfield.
  bf = Bfield(t);

  // Set the relaxtion bits of the differential.
  s2[0] = sim::gamma_2 * s[0];
  s2[1] = sim::gamma_2 * s[1];
  s2[2] = sim::gamma_1 * (s[2] - 1.0);

  // Calculate the cross product.
  cross(bf, s, s1);

  // Set the differential to be integrated.
  dsdt = s1 - s2;
}


// A function that defines the static + RF B-fields.
std::vector<double> FidFactory::Bfield(const double& t)
{
  // Made static to save on memory calls.
  static std::vector<double> a = {0., 0., 0.}; // constant external field
  static std::vector<double> b = {0., 0., 0.}; // time dependent B field

  // Return static external field if after the pulsed field.
  if (t >= sim::t_pulse){
    return a;

  // Set the fields if the simulation is just starting.
  } else if (t <= ti_ + dt_) {

    a[2] = kTau * sim::freq_larmor;
    b[2] = kTau * sim::freq_larmor;

  }

  // Return static field if the time is before the pulsed field.
  if (t < 0.0) return a;

  // If none of the above, return the time-dependent, pulsed field.
  b[0] = sim::omega_r * cos(kTau * sim::freq_ref * t);
  b[1] = sim::omega_r * sin(kTau * sim::freq_ref * t);
  return b;
}

// The function called each step of integration.
void FidFactory::Printer(std::vector<double> const &s , double t)
{
  // Cache the cosine function for mixing later.
  if (cos_cache_.size() == 0){
    cos_cache_.reserve(sim_length_);

    double temp = t;
    double val;

    for (int i = 0; i < sim_length_; i++){
      val = cos(kTau * sim::freq_ref * temp + sim::mixdown_phi);
      cos_cache_.push_back(val);
      temp += dt_;
    }

    spin_vec_.assign(sim_length_, 0.0);
    time_vec_.assign(sim_length_, 0.0);
  }

  // Make sure the index is reset
  if (t < ti_ + dt_) printer_idx_ = 0;

  // Record spin in the y-direction and mix down
  spin_vec_[printer_idx_] = sim::amplitude * s[1] * cos_cache_[printer_idx_]; 

  // Record the time and increment printer_idx_
  time_vec_[printer_idx_++] = t;

  // If the FID is done, do more stuff.
  if (t > tf_ - dt_) {

    // Final spin sum from all different gradients
    spin_vec_ = LowPassFilter(spin_vec_);

    // Reset the counters
    printer_idx_ = 0; // They'll get incremented after.
  }
}

// Low pass filter suppresses higher frequencies from mixing down.
std::vector<double> FidFactory::LowPassFilter(std::vector<double>& s)
{
  // Allocate the filter and set the central frequency.
  std::vector<double> filter;
  double freq_cut = sim::freq_cut_ratio * sim::freq_larmor;

  // Define the filter if not defined.  Using 3rd order Butterworth filter.
  if (filter.size() == 0) {

    filter.resize(sim_length_);
    int i = 0;
    int j = sim_length_ - 1;
    double temp;

    // The filter is symmetric, so we can fill both sides in tandem.
    while (i < sim_length_ / 2){
      // scaling term
      temp = pow(i / (dt_ * sim_length_ * freq_cut), 6);
      filter[i] = pow(1.0 / (1.0 + temp), 0.5);
      filter[j--] = filter[i++];
    }
  }

  // Copy the vectors to Armadillo formats first.
  arma::vec v1(&s[0], s.size(), false); 
  arma::vec v2(&filter[0], filter.size(), false);

  // Now get the FFT
  arma::cx_vec fft = arma::fft(v1);

  // Multiply the FFT and the filter element wise
  fft = fft % v2;

  // This confusing oneliner, get the inverse FFT and converts to a std::vector
  return arma::conv_to<std::vector<double>>::from(arma::real(arma::ifft(fft)));
}

int FidFactory::GetTreeIndex(double grad_strength)
{
  return (int)(std::nearbyint(grad_strength / d_grad_ + zero_idx_) + 0.5);
}

} // ::fid