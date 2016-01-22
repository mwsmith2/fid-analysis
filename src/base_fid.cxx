#include "fid/base_fid.h"

namespace fid {

void BaseFid::Init()
{
  // Initialize the health properly.
  health_ = 100.0;

  // Resize the temp array (maybe others?)
  temp_.reserve(wf_.size());

  // Initialize the BaseFid for analysis
  LoadParams();
  CenterFid();
  CalcNoise();
  CalcMaxAmp();
  FindFidRange();
  InitHook();

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


// Load FID data from a formatted text file.
void BaseFid::LoadTextData(std::string filename)
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


// Load all the current parameters in the fid::params namespace.
void BaseFid::LoadParams()
{
  edge_width_ = params::edge_width;
  edge_ignore_ = params::edge_ignore; 
  start_amplitude_ = params::start_amplitude; 
  low_pass_freq_ = params::low_pass_freq; 
  fft_peak_width_ = params::fft_peak_width;
  centroid_thresh_ = params::centroid_thresh; 
  hyst_thresh_ = params::hyst_thresh; 
  snr_thresh_ = params::snr_thresh;
  len_thresh_ = params::len_thresh;
  freq_method_ = params::freq_method;
}


void BaseFid::CenterFid()
{
  int w = edge_width_;
  double sum  = std::accumulate(wf_.begin(), wf_.begin() + w, 0.0);
  double avg = sum / w; // to pass into lambda
  mean_ = avg; // save to class

  std::for_each(wf_.begin(), wf_.end(), [avg](double& x){ x -= avg; });
}


void BaseFid::CalcNoise()
{ 
  // Grab a new handle to the noise window width for aesthetics.
  int i = edge_ignore_;
  int f = edge_width_ + i;

  // Find the noise level in the head and tail.
  double head = stdev(wf_.begin() + i, wf_.begin() + f);
  double tail = stdev(wf_.rbegin() + i, wf_.rbegin() + f);

  // Take the smaller of the two.
  noise_ = (tail < head) ? (tail) : (head);
}


void BaseFid::CalcMaxAmp() 
{
  auto mm = std::minmax_element(wf_.begin(), wf_.end());

  if (std::abs(*mm.first) > std::abs(*mm.second)) {

    max_amp_ = std::abs(*mm.first);

  } else {

    max_amp_ = std::abs(*mm.second);
  }
}


void BaseFid::FindFidRange()
{
  // Find the starting and ending points
  double thresh = start_amplitude_ * max_amp_;
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
  auto it_2 = std::find_if(wf_.begin() + i_wf_, wf_.end(),
      [thresh](double x) {
        return std::abs(x) < 0.8 * thresh;
  });

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


// Save the interanl TGraph.
void BaseFid::SaveGraph(std::string filename, std::string title)
{ 
  gr_.SetTitle(title.c_str());
  gr_.Draw();
  c1_.Print(filename.c_str());
}


// Save a plot of FID waveform.
void BaseFid::SavePlot(std::string filename, std::string title)
{
  // If no title supplied give a reasonable default.
  if (title == "") {

    title = std::string("FID; time [ms]; amplitude [a.u.]");

  } else {

    // In case they didn't append x/y labels.
    title.append("; time [ms]; amplitude [a.u.]");
  }

  gr_ = TGraph(wf_.size(), &tm_[0], &wf_[0]);

  SaveGraph(filename, title);
}


// Print the time series fit from an FID.
void BaseFid::SaveTimeFit(std::string filename, std::string title)
{
  if (title == "") {

    title = std::string("Time Series Fit; time [ms]; amplitude [a.u.]");

  } else {

    // In case they didn't append x/y labels.
    title.append("; time [ms]; amplitude [a.u.]");
  }  

  // Copy the current time fit graph.
  gr_ = gr_time_series_;
  SaveGraph(filename, title);
}

// Print the time series fit from an FID.
void BaseFid::SaveFreqFit(std::string filename, std::string title)
{
  if (title == "") {

    title = std::string("Frequency Series Fit; time [ms]; amplitude [a.u.]");

  } else {

    // In case they didn't append x/y labels.
    title.append("; freq [kHz]; amplitude [a.u.]");
  }  

  // Copy the current time fit graph.
  gr_ = gr_freq_series_;
  SaveGraph(filename, title);
}

void BaseFid::SaveTimeRes(std::string filename, std::string title)
{
  if (title == "") {

    title = std::string("Time Series Fit Residuals; time [ms]; amplitude [a.u.]");

  } else {

    // In case they didn't append x/y labels.
    title.append("; time [ms]; amplitude [a.u.]");
  }  

  // Copy the current time fit.
  gr_ = gr_time_series_;

  // Set the points
  for (uint i = 0; i < res_.size(); ++i){
    static double x, y;

    gr_.GetPoint(i, x, y);
    gr_.SetPoint(i, x, res_[i]); 
  }

  SaveGraph(filename, title);
}


void BaseFid::SaveFreqRes(std::string filename, std::string title)
{
  if (title == "") {

    title = std::string("Freq Series Fit Residuals; time [ms]; amplitude [a.u.]");

  } else {

    // In case they didn't append x/y labels.
    title.append("; freq [kHz]; amplitude [a.u.]");
  }  

  // Copy the current time fit.
  gr_ = gr_freq_series_;

  // Set the points
  for (uint i = 0; i < res_.size(); ++i){
    static double x, y;

    gr_.GetPoint(i, x, y);
    gr_.SetPoint(i, x, res_[i]); 
  }

  SaveGraph(filename, title);
}


// Save the FID data to a text file as "<time> <amp>".
void BaseFid::SaveData(std::string filename)
{
  // open the file first
  std::ofstream out(filename);

  for (int i = 0; i < tm_.size(); ++i) {
    out << tm_[i] << " " << wf_[i] << std::endl;
  }
}


void BaseFid::DiagnosticInfo(std::ostream& out)
{
  using std::endl;

  // Save the flags, set them to defaults.
  auto flags = out.flags();
  std::ofstream testout;
  out.flags(testout.flags());

  out << std::string(80, '<') << endl << std::string(4, ' ');
  out << "Diagostic Information for Fid @ " << this << endl;
  out << std::string(80, '<') << endl;

  out << "    Fid Waveform Characteristics" << endl;
  out << "        mean:       " << mean_ << endl;
  out << "        amplitude:  " << max_amp_ << endl;
  out << "        noise:      " << noise_ << endl;
  out << "        start time: " << i_wf_;
  out << " (" << tm_[i_wf_] << " ms)" << endl;
  out << "        stop time:  " << f_wf_ - 1;
  out << " (" << tm_[f_wf_ - 1] << " ms)" << endl;
  out << "        health:     " << health_ << endl;
  out << std::string(80, '>') << endl << endl;
  
  // Restore set flags.
  out.flags(flags);
}


void BaseFid::DiagnosticPlot(std::string dirname, std::string filestub)
{
  boost::filesystem::path dir(dirname);
  boost::filesystem::create_directories(dir);

  SaveFreqFit(dir.string() + filestub + std::string("_freq_fit.png"));
  SaveTimeFit(dir.string() + filestub + std::string("_time_fit.png"));
  SaveFreqRes(dir.string() + filestub + std::string("_freq_res.png"));
  SaveTimeRes(dir.string() + filestub + std::string("_time_res.png"));
}


void BaseFid::DiagnosticDump(std::string dirname, std::string filestub)
{
  // Make the plots first, that will create the directory if needed.
  DiagnosticPlot(dirname, filestub);

  std::ofstream out;
  boost::filesystem::path dir(dirname);

  std::string str = dir.string() + std::string("libfid.log");
  out.open(str , std::ofstream::out | std::ofstream::app);

  DiagnosticInfo(out);
  out.close();
}

} // fid
