#include "fid_sim.h"

namespace fid {

// static variables get initialized here
double FidFactory::ti_ = i_time;
double FidFactory::tf_ = ti_ + d_time * len_fids;
double FidFactory::dt_ = sim::dt_integration;
double FidFactory::x_ = 0.0;

vec FidFactory::s_;
vec FidFactory::spin_;
vec FidFactory::time_vec_;
vec FidFactory::spin_sum_;
vec FidFactory::cos_cache_;
vec FidFactory::gradient_(1, 0.0);


//---------------------------------------------------------------------------//
//--- FID Factory -----------------------------------------------------------//
//---------------------------------------------------------------------------//

FidFactory::FidFactory(){
  
  ti_ = i_time;
  tf_ = ti_ + d_time * len_fids;
  dt_ = sim::dt_integration;
  x_ = 0.0;
}

void FidFactory::SimulateFid(vec& wf, vec& tm)
{
  // make sure memory is allocated for the final FIDs
  tm.reserve(len_fids);
  wf.reserve(len_fids);

  s_ = sim::spin_0; // The starting spin vector
  integrate_const(runge_kutta4<vec>(), Bloch, s_, ti_, tf_, dt_, Printer);

  // Set the results
  tm.resize(0);
  wf.resize(0);
  for (int i = 0; i < sim::num_points; i += sim::reduction){
    tm.push_back(time_vec_[i]);
    wf.push_back(spin_sum_[i]);
  }
}

// Create an idealized FID with current Simulation parameters
void FidFactory::IdealFid(vec& wf, vec& tm)
{
  wf.reserve(tm.size());
  wf.resize(0);

  // Define the waveform
  double temp;

  for (auto it = tm.begin(); it != tm.end(); ++it){

    if (*it >= sim::t_pulse){
      temp = std::exp(-(*it - sim::t_pulse) * s_tau);
      temp *= std::sin((*it) * kTau * s_freq + s_phase);
      wf.push_back(temp);

    } else {
      wf.push_back(0.0);

    }
  } 

  // Add some noise
  static std::default_random_engine gen(sim::seed);
  std::normal_distribution<double> norm(0.0, 1.0 / s_snr);
  for (auto it = wf.begin(); it != wf.end(); ++it){
    *it += norm(gen);
  }
}

void FidFactory::Bloch(vec const &s, vec &dsdt, double t)
{
  // Again static to save time on memory allocations.
  static vec b = {{0., 0., 0.}};   // Bfield
  static vec s1 = {{0., 0., 0.}};  // Cross product piece of spin
  static vec s2 = {{0., 0., 0.}};  // Relaxation piece of spin
  static double a1 = kTau * sim::gamma_1; // Relaxation time #1
  static double a2 = kTau * sim::gamma_2; // Relaxation time #2

  // Only need to update the Bfield before the pulsed kick has finished.
  if (t <= sim::t_pulse + dt_) b = Bfield(x_, t);

  // Set the relaxtion bits of the field.
  s2[0] = a2 * s[0];
  s2[1] = a2 * s[1];
  s2[2] = a1 * s[2] - a1;

  // Calculate the cross product.
  Cross(b, s, s1);

  // Set the differential to be integrated.
  dsdt = s1 - s2;
}

vec FidFactory::Bfield(const double& x, const double& t)
{
  // Made static to save on memory calls.
  static vec a = {0., 0., 0.}; // holds constant external field
  static vec b = {0., 0., 0.}; // for time dependent B field
  static double w = kTau * sim::freq_ref;

  // Return static external field if after the pulsed field.
  if (t >= sim::t_pulse){
    return a;

  // Set the fields if the simulation is just starting.
  } else if (t <= ti_ + dt_) {
    a[2] = kTau * (1. + x_) * sim::freq_larmor;
    b[2] = kTau * (1. + x_) * sim::freq_larmor;
  }

  // Return static field if the time is before the pulsed field.
  if (t < 0.0) return a;

  // If none of the above, return the time-dependent, pulsed field.
  b[0] = sim::omega_r * cos(w * t);
  b[1] = sim::omega_r * sin(w * t);
  return b;
}

// Printer function is called to do stuff each step of integration.
void FidFactory::Printer(vec const &s , double t){

  // Initialized static memory cache.
  static int count = 0;
  static int index = 0;
  static int fid_num = 0;
  static int step_num = 0;
  static int ref_count = (tf_ - ti_) / (dt_ * sim::num_points) + 0.5; 

  // Cache the cosine function for mixing later.
  if (cos_cache_.size() == 0){
    cos_cache_.reserve(sim::num_points);
    double temp = t;
    for (int i = 0; i < sim::num_points; i++){
      cos_cache_.push_back(cos(kTau * sim::freq_ref * temp));
      temp += dt_ * ref_count;
    }
    spin_.assign(sim::num_points, 0.0);
    spin_sum_.assign(sim::num_points, 0.0);
    time_vec_.assign(sim::num_points, 0.0);
  }

  // Time is equal to one of our time steps, do stuff.
  if (count++ % ref_count == 0){ 

    // Record spin in the y-direction and mix down
    spin_[index] = s[1] * cos_cache_[index]; 

    // Record the time and increment index
    time_vec_[index++] = t;

    // If the FID is done, do more stuff.
    if (t == tf_){

      // Reset the spin_sum_ if we are starting a new FID
      if (step_num == 0) spin_sum_.assign(sim::num_points, 0.0);

      // If the FID has no more gradient steps left
      if (step_num == gradient_.size() - 1){

        // Final spin sum from all different gradients
        spin_sum_ = spin_sum_ + spin_;
        spin_sum_ = LowPassFilter(spin_sum_);

        // Reset the counters
        index = 0; // They'll get incremented after.
        step_num = 0;
        count = 0;

        // Progress report
        if (fid_num % (num_fids/10) == 0){
          cout << fid_num << " FIDs have been simulated." << endl;
        }

      } else {

        // Add the spin to the set and do the next gradient.
        spin_sum_ = spin_sum_ + spin_;
        step_num++;
        index = 0;
        count = 0;

      }
    }
  }
}

// Low pass filter to suppress the higher frequency introducing in mixing down.
vec FidFactory::LowPassFilter(vec& s){

  // Store the filter statically though this might be a minimal speed boost.
  static vec filter;
  static double freq_cut = kTau * sim::freq_ref;

  // Define the filter if not defined.  Using 3rd order Butterworth filter.
  if (filter.size() == 0){

    filter.resize(sim::num_points);
    int i = 0;
    int j = sim::num_points-1;
    double temp;

    // The filter is symmetric, so we can fill both sides in tandem.
    while (i < sim::num_points / 2){
      // scaling term
      temp = pow((kTau * i) / (dt_ * sim::num_points * freq_cut), 6);
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
  return arma::conv_to<vec>::from(arma::real(arma::ifft(fft)));
}


//---------------------------------------------------------------------------//
//--- Gradient FID Factory --------------------------------------------------//
//---------------------------------------------------------------------------//
GradientFidFactory::GradientFidFactory()
{
  // open the default root file
  pf_fid_ = new TFile(grad::root_file.c_str());
  pt_fid_ = (TTree *)pf_fid_->Get("t");
  
  pt_fid_->GetEntry(0);
  pt_fid_->SetBranchAddress(grad::fid_branch.c_str(), &wf_[0]);

  num_sim_fids_ = pt_fid_->GetEntries();
  d_grad_ = (grad::max - grad::min) / num_sim_fids_;
  zero_idx_ = -1 * grad::min / d_grad_ + 0.5;
}

GradientFidFactory::~GradientFidFactory()
{
  // closing the TFile will take care of the TTree pointer
  pf_fid_->Close();
  delete pf_fid_;
}

void GradientFidFactory::ConstructFid(const vec& gradient, vec& wf)
{
  // Find the appropriate FIDs and sum them
  wf.assign(len_fids, 0.0);

  for (auto val : gradient){

    pt_fid_->GetEntry(GetTreeIndex(val));

    for (int i = 0; i < wf.size(); i++){
      wf[i] += wf_[i] / gradient.size();
    }
  }
}

int GradientFidFactory::GetTreeIndex(double grad_strength)
{
  return (int)(std::nearbyint(grad_strength / d_grad_ + zero_idx_) + 0.5);
}

} // ::fid