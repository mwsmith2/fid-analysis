/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   11/02/14

File:   prog/sim_fid_lin_grad.cxx
Detail: The program tests the effects of linear gradients on frequency extraction from free induction decays.  The waveforms are simulated.

\*===========================================================================*/


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

//--- project includes ------------------------------------------------------//
#include "fid.h"
#include "fid_class.h"

// unnamed local global namespace
namespace {
  // declare variables
  int fid_length = 10000;
  int nfids = 200;
  double ti = -1.0;
  double dt = 0.001;
  double dg = 50.0;

  int npoints = 5;

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
  return (int)(std::nearbyint(grad_strength / 0.1 + 5000) + 0.5);
}

// Implement main function
int main(int argc, char** argv)
{
  // set precision
  cout.precision(10);
  cout.setf(std::ios::fixed, std:: ios::floatfield);

  // Declare objects
  std::ofstream out;
  out.precision(10);
  out.setf(std::ios::fixed, std:: ios::floatfield);
  out.open("quad_grad_data.csv");

  vector<double> wf;
  vector<double> tm;
  wf.reserve(fid_length);
  tm.reserve(fid_length);
  
  SetTimeVector(fid_length, ti, dt, tm);

  vector<double> grad_strengths;
  vector<double> grad_norm;
  vector<double> gradient;

  ConstructGradient(npoints, grad_norm);
  gradient = grad_norm;

  for (double g = 0.0; g <= 500.0; g += dg){
    grad_strengths.push_back(g);
  }

  // root stuff

  my_fid.resize(fid_length);
  pt->SetBranchAddress("fid", &::my_fid[0]);
  pt->GetEntry(0);

  for (auto val : grad_strengths){

    for (int i = 0; i < grad_norm.size(); i++){
      gradient[i] = val * grad_norm[i];
    }

    cout << "Running for gradient strength " << val << " ppm.\n";

    for (int i = 0; i < nfids; i++){

        ConstructFID(gradient, wf);
        AddWhiteNoise(wf);

        fid::FID my_fid(wf, tm);

        out << val << ", ";
        out << my_fid.CalcZeroCountFreq() << ", ";
        out << my_fid.CalcCentroidFreq() << ", ";
        out << my_fid.CalcAnalyticalFreq() << ", ";
        out << my_fid.CalcLorentzianFreq() << ", ";
        out << my_fid.CalcSoftLorentzianFreq() << ", ";
        out << my_fid.CalcExponentialFreq() << ", ";
        out << my_fid.CalcPhaseFreq() << ", ";
        out << my_fid.CalcSinusoidFreq() << endl;
    }
  }

  out.close();
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
    grad.push_back((double)i * i);
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





