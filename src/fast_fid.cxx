#include "fid/fast_fid.h"

namespace fid {

FastFid::FastFid(const std::string& fid_file)
{
  // Read and store the waveform and time from a .fid file.
  LoadTextData(fid_file);

  Init();
}


FastFid::FastFid(const char* fid_file)
{
  // Convert the char pointer toa string.
  std::string fid_string(fid_file);

  // Read and store the waveform and time from a .fid file.
  LoadTextData(fid_string);

  Init();
}


FastFid::FastFid(const std::vector<double>& wf, const std::vector<double>& tm)
{
  // Copy the waveform and time to member vectors.
  wf_ = wf;
  tm_ = tm;

  Init();
}

FastFid::FastFid(const std::vector<double>& wf)
{
  // Copy the waveform and construct a generic time range.
  wf_ = wf;
  tm_ = construct_range(0.0, (double)wf_.size(), 1.0);

  Init();
}


void FastFid::Init()
{
  // Initialize the health properly.
  health_ = 100.0;

  // Resize the temp array (maybe others?)
  temp_.reserve(wf_.size());

  // Initialize the FastFid for analysis
  LoadParams();
  CenterFid();
  CalcNoise();
  CalcMaxAmp();
  FindFidRange();
  CalcFreq();

  // Flag the FID as bad if it's negative.
  if (freq_ < 0.0) {

    health_ = 0.0;

  }
  
  // Else calculate a health based on signal to noise.
  if (max_amp_ < noise_ * snr_thresh_) {
    health_ *= max_amp_ / (noise_ * snr_thresh_);
  }

  // And factor fid duration into health.
  if (f_wf_ - i_wf_ < wf_.size() * len_thresh_) {
    health_ *= (f_wf_ - i_wf_) / (wf_.size() * len_thresh_);
  }
}


// Load all the current parameters in the fid::params namespace.
void FastFid::LoadParams()
{
  freq_method_ = params::freq_method;
  len_thresh_ = params::len_thresh;
  snr_thresh_ = params::snr_thresh;
  hyst_thresh_ = params::hyst_thresh; 
  centroid_thresh_ = params::centroid_thresh; 
  low_pass_freq_ = params::low_pass_freq; 
  max_phase_jump_ = params::max_phase_jump; 
  zc_alpha_ = params::zc_alpha; 
  start_thresh_ = params::start_thresh; 
  edge_ignore_ = params::edge_ignore; 
  zc_width_ = params::zc_width;
  fit_width_ = params::fit_width;
}


void FastFid::CenterFid()
{
  int w = zc_width_;
  double sum  = std::accumulate(wf_.begin(), wf_.begin() + w, 0.0);
  double avg = sum / w; // to pass into lambda
  mean_ = avg; // save to class

  std::for_each(wf_.begin(), wf_.end(), [avg](double& x){ x -= avg; });
}


void FastFid::CalcNoise()
{ 
  // Grab a new handle to the noise window width for aesthetics.
  int i = edge_ignore_;
  int f = zc_width_ + i;

  // Find the noise level in the head and tail.
  double head = stdev(wf_.begin() + i, wf_.begin() + f);
  double tail = stdev(wf_.rbegin() + i, wf_.rbegin() + f);

  // Take the smaller of the two.
  noise_ = (tail < head) ? (tail) : (head);
}


void FastFid::CalcMaxAmp() 
{
  auto mm = std::minmax_element(wf_.begin(), wf_.end());

  if (std::abs(*mm.first) > std::abs(*mm.second)) {

    max_amp_ = std::abs(*mm.first);

  } else {

    max_amp_ = std::abs(*mm.second);
  }
}


void FastFid::FindFidRange()
{
  // Find the starting and ending points
  double thresh = start_thresh_ * max_amp_;
  bool checks_out = false;

  // Find the first element with magnitude larger than thresh
  auto it_1 = wf_.begin() + edge_ignore_;

  while (!checks_out) {

    // Check if the point is above threshold.
    auto it_i = std::find_if(it_1, wf_.end(), 
      [thresh](double x) { 
        return std::abs(x) > thresh; 
    });

    // Make sure the point is not with one of the vector's end.
    if ((it_i != wf_.end()) && (it_i + 1 != wf_.end())) {

      // Check if the next point is also over threshold.
      checks_out = std::abs(*(it_i + 1)) > thresh;

      // Increase the comparison starting point.
      it_1 = it_i + 1;

      // Turn the iterator into an index
      if (checks_out) {
        i_wf_ = std::distance(wf_.begin(), it_i);
      }

    } else {

      // If we have reached the end, mark it as the last.
      i_wf_ = std::distance(wf_.begin(), wf_.end());
      break;
    }
  }

  // Find the next element with magnitude lower than thresh
  auto it_2 = wf_.begin() + i_wf_ + edge_ignore_;
  checks_out = false;

  while (!checks_out) {

    // Find the range around a peak.
    auto it_i = std::find_if(it_2, wf_.end(), 
      [thresh](double x) {
        return std::abs(x) > 0.8 * thresh;
    });

    auto it_f = std::find_if(it_i + 1, wf_.end(),
      [thresh](double x) {
        return std::abs(x) < 0.8 * thresh;
    });

    // Now check if the peak actually made it over threshold.
    if ((it_i != wf_.end()) && (it_f != wf_.end())) {

      auto mm = std::minmax_element(it_i, it_f);

      if ((*mm.first < -thresh) || (*mm.second > thresh)) {

        it_2 = it_f;

      } else {

        checks_out = true;
      }

      // Turn the iterator into an index
      if (checks_out) {
        f_wf_ = std::distance(wf_.begin(), it_f);
      }
    
    } else {

      f_wf_ = std::distance(wf_.begin(), wf_.end());
      break;
    }
  }

  // Gradients can cause a waist in the amplitude.
  // Mark the signal as bad if it didn't find signal above threshold.
  if (i_wf_ > wf_.size() * 0.95 || i_wf_ >= f_wf_) {

    health_ = 0.0;

    i_wf_ = 0;
    f_wf_ = wf_.size() * 0.01;
  } 
}


double FastFid::CalcFreq()
{
  // set up vectors to hold relevant stuff about the important part
  temp_.resize(f_wf_ - i_wf_);

  int nzeros = 0;
  bool pos = wf_[i_wf_] >= 0.0;
  bool hyst = false;
  
  int i_zero = -1;
  int f_zero = -1;
  double thresh = hyst_thresh_ * start_thresh_ * max_amp_; 

  // iterate over vector
  for (unsigned int i = i_wf_; i < f_wf_; i++){

    // hysteresis check
    if (hyst) {
      hyst = std::abs(wf_[i]) > thresh;
      continue;
    }

    // check for a sign change
    if ((wf_[i] >= 0.0) != pos){
      nzeros++;
      f_zero = i;
      if (i_zero == -1) i_zero = i;
      pos = !pos;
      hyst = true;
    }
  }

  // Use linear interpolation to find the zeros more accurately
  int i = i_zero;
  int f = f_zero;

  // do the interpolation
  double frac = std::abs(wf_[i] / (wf_[i-1] - wf_[i]));
  double ti = frac * tm_[i-1] + (1.0 - frac) * tm_[i];

  frac = std::abs(wf_[f] / (wf_[f-1] - wf_[f]));
  double tf = frac * tm_[f-1] + (1.0 - frac) * tm_[f];

  freq_ = 0.5 * (nzeros - 1) / (tf - ti);
  // todo: Fix this into a better error estimate. 
  freq_err_ = freq_ * sqrt(2) * (tm_[1] - tm_[0]) / (tf - ti);

  return freq_;
}


double FastFid::GetFreq()
{
  return freq_;
}


double FastFid::GetFreqError()
{
  return freq_err_;
}


void FastFid::PrintDiagnosticInfo()
{
  using std::cout;
  using std::endl;

  cout << endl << std::string(80, '<') << endl;
  cout << "Printing Diagostic Information for FastFid @ " << this << endl;
  cout << "noise level: " << noise_ << endl;
  cout << "waveform start, stop: " << i_wf_ << ", " << f_wf_ << endl;
  cout << std::string(80, '>') << endl;
}

void FastFid::PrintDiagnosticInfo(std::iostream out)
{
  using std::cout;
  using std::endl;

  out << std::string(80, '<') << endl;
  out << "Printing Diagostic Information for FastFid @ " << this << endl;
  out << "noise level: " << noise_ << endl;
  out << "waveform start, stop: " << i_wf_ << ", " << f_wf_ << endl;
  out << std::string(80, '>') << endl;
}


// Load FID data from a formatted text file.
void FastFid::LoadTextData(std::string filename)
{
  // open the file first
  std::ifstream in(filename);

  // shrink vectors
  wf_.resize(0);
  tm_.resize(0);

  double wf_temp;
  double tm_temp;

  while (in.good()) {
    in >> tm_temp >> wf_temp;
    tm_.push_back(tm_temp);
    wf_.push_back(wf_temp);
  } 
}

} // fid
