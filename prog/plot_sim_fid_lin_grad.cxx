/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   11/02/14

Detail: The program should plot the unwrapped phase, the linear phase fit, and
the residuals for each gradient strength.  It print the figures to the current
directory.

\*==========================================================================*/


//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <numeric>
#include <algorithm>
using std::vector;
using std::cout;
using std::endl;

//-- other includes ---------------------------------------------------------//
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"

//--- project includes ------------------------------------------------------//
#include "fid.h"
#include "fid_class.h"

// unnamed local global namespace
namespace {
  // declare variables
  int fid_length = 10000;
  int nfids = 1000;
  double ti = -1.0;
  double dt = 0.01;

  // gradient loop
  double gmin = 0.0;
  double gmax = 200.0;
  double dg = 10.0;

  int npoints = 21;

  int fig_w = 1200;
  int fig_h = 900;

  // Get our root data
  TFile *pf = new TFile("input/sim_fids.root");
  TTree *pt = (TTree *)pf->Get("t");
  vector<double> my_fid;
}

// Declare methods
void AddWhiteNoise(vector<double> &wf, double s2n=100.0);
void SetTimeVector(int ntimes, double t0, double dt, vector<double> &tm);
void ConstructGradient(int npoints, vector<double> &grad);
void ConstructFID(vector<double> &grad, vector<double> &wf);

inline int GetTreeIndex(double grad_strength){
  return (int)(std::nearbyint(grad_strength / (0.1) + 5000) + 0.5);
}

// Implement main function
int main(int argc, char** argv)
{
  // set precision
  cout.precision(10);
  cout.setf(std::ios::fixed, std:: ios::floatfield);

  vector<double> wf;
  vector<double> tm;
  vector<double> phase;
  vector<double> fit;
  vector<double> res;
  wf.reserve(fid_length);
  tm.reserve(fid_length);
  phase.reserve(fid_length);
  fit.reserve(fid_length);
  res.reserve(fid_length);
  
  SetTimeVector(fid_length, ti, dt, tm);

  vector<double> grad_strengths;
  vector<double> grad_norm;
  vector<double> gradient;

  ConstructGradient(npoints, grad_norm);
  gradient = grad_norm;

  for (double g = gmin; g <= gmax; g += dg){
    grad_strengths.push_back(g);
  }

  // root stuff
  my_fid.resize(fid_length);
  pt->SetBranchAddress("fid", &::my_fid[0]);
  pt->GetEntry(0);
  TCanvas c1("c1", "", fig_w, fig_h);

  for (auto val : grad_strengths){

    // Construct the gradient
    for (int i = 0; i < grad_norm.size(); i++){
      gradient[i] = val * grad_norm[i];
    }

    ConstructFID(gradient, wf);
    AddWhiteNoise(wf);

    fid::FID my_fid(wf, tm);
    my_fid.CalcPhaseFreq();

    // Get the phase
    phase = my_fid.phase();
    cout << phase.size() << ", " << tm.size() << endl;
    TGraph gr_phs(phase.size(), &tm[0], &phase[0]);

    // Get the phase fit
    TF1 fit_func(my_fid.f_fit());

    fit.resize(0);
    for (auto t : tm){
      fit.push_back(fit_func.Eval(t));
    }

    TGraph gr_fit(fit.size(), &tm[0], &fit[0]);

    res.resize(phase.size());
    std::transform(phase.begin(), phase.end(), fit.begin(), res.begin(),
      [](double x1, double x2) {return x1 - x2;});

    TGraph gr_res(res.size(), &tm[0], &res[0]);

    gr_phs.Draw();
    c1.Print("test1.pdf");
    gr_fit.Draw();
    c1.Print("test2.pdf");
    gr_res.Draw();
    c1.Print("test3.pdf");

  } // grad_strengths

  // clean up
  pf->Close();
  delete pf;

  return 0;
}

// Implement helper functions
void AddWhiteNoise(vector<double> &wf, double s2n){
  static std::default_random_engine gen;
  static std::normal_distribution<double> nrm(0.0, s2n);

  double max = *std::max_element(wf.begin(), wf.end());
  double min = *std::min_element(wf.begin(), wf.end());
  double scale = max > min ? max : min;

  for (auto x : wf){
    x += scale * nrm(gen);
  }

  return;
}

void SetTimeVector(int ntimes, double t0, double dt, vector<double> &tm)
{
  if (tm.size() != ntimes){
    tm.resize(ntimes);
  }

  for (int i = 0; i < ntimes; i++){
    tm[i] = dt * i - t0;
  }

  return;
}

void ConstructGradient(int npoints, vector<double> &grad)
{
  // construct a normalize linear gradient

  // first get the spacing right
  for (int i = 0; i < npoints; i++){
    grad.push_back((double)i);
  }

  // subtract off the average
  double avg = std::accumulate(grad.begin(), grad.end(), 0.0) / grad.size();
  for (int i = 0; i < grad.size(); i++){
    grad[i] -= avg;
  }

  // normalize by largest value
  double max = *std::max_element(grad.begin(), grad.end());
  for (int i = 0; i < grad.size(); i++){
    grad[i] /= max;
  }

  cout << endl;

}

void ConstructFID(vector<double> &grad, vector<double> &wf)
{
  // Find the appropriate FIDs and sum them
  wf.assign(fid_length, 0.0);

  for (auto val : grad){

    pt->GetEntry(GetTreeIndex(val));

    for (int i = 0; i < wf.size(); i++){
      wf[i] += my_fid[i] / grad.size();
    }
  }

  return;
}





