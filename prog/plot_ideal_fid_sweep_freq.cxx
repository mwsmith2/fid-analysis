/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   11/02/14

Detail: This is a new test program for my FID libraries

\*==========================================================================*/


//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <fstream>
#include <vector>
using std::vector;
using std::cout;
using std::endl;

//-- other inclues ----------------------------------------------------------//
#include "TCanvas.h"
#include "TString.h"
#include "TGraph.h"

//--- project includes ------------------------------------------------------//
#include "fid.h"
#include "fid_class.h"


int main(int argc, char** argv)
{
  // set precision
  cout.precision(10);
  cout.setf(std::ios::fixed, std:: ios::floatfield);

  // declare variables
  int fid_length = 10000;
  int nfids = 1;
  double ti = -1.0;
  double dt = 0.001;

  double fi = 23.0;
  double ff = 23.2;
  double df = 0.01;

  vector<double> wf;
  vector<double> tm;
  wf.reserve(fid_length);
  tm.reserve(fid_length);

  TString name;
  TGraph gr;
  TCanvas c1("c1", "", 1800, 600);

  for (int i = 0; i < fid_length; i++){
    tm.push_back(i * dt + ti);
  }

  for (double f = fi; f <= ff; f += df){

    cout << "Running for frequency " << f << ".\n";

    for (int i = 0; i < nfids; i++){

      fid::getIdealFID(wf, tm, f);

      fid::FID my_fid(wf, tm);

      name = "an";
      my_fid.CalcAnalyticalFreq();
      gr = my_fid.gr_freq_series();
      my_fid.gr_freq_series().Draw();
      my_fid.f_fit().Draw("Same");
      c1.Print(TString::Format("%s_%.2f_%d.C", name.Data(), f, i));
      c1.Print(TString::Format("%s_%.2f_%d.png", name.Data(), f, i));
      c1.Print(TString::Format("%s.gif+20", name.Data()));

      name = "lz";
      my_fid.CalcLorentzianFreq();
      my_fid.gr_freq_series().Draw();
      my_fid.f_fit().Draw("Same");
      c1.Print(TString::Format("%s_%.2f_%d.C", name.Data(), f, i));
      c1.Print(TString::Format("%s_%.2f_%d.png", name.Data(), f, i));
      c1.Print(TString::Format("%s.gif+20", name.Data()));

      name = "sl";
      my_fid.CalcSoftLorentzianFreq();
      gr = my_fid.gr_freq_series();
      my_fid.gr_freq_series().Draw();
      my_fid.f_fit().Draw("Same");
      c1.Print(TString::Format("%s_%.2f_%d.C", name.Data(), f, i));
      c1.Print(TString::Format("%s_%.2f_%d.png", name.Data(), f, i));
      c1.Print(TString::Format("%s.gif+20", name.Data()));

      name = "ex";
      my_fid.CalcExponentialFreq();
      gr = my_fid.gr_freq_series();
      my_fid.gr_freq_series().Draw();
      my_fid.f_fit().Draw("Same");
      c1.Print(TString::Format("%s_%.2f_%d.C", name.Data(), f, i));
      c1.Print(TString::Format("%s_%.2f_%d.png", name.Data(), f, i));
      c1.Print(TString::Format("%s.gif+20", name.Data()));

      name = "p1";
      my_fid.CalcPhaseFreq();
      gr = my_fid.gr_time_series();
      my_fid.gr_freq_series().Draw();
      my_fid.f_fit().Draw("Same");
      c1.Print(TString::Format("%s_%.2f_%d.C", name.Data(), f, i));
      c1.Print(TString::Format("%s_%.2f_%d.png", name.Data(), f, i));
      c1.Print(TString::Format("%s.gif+20", name.Data()));

      name = "sn";
      my_fid.CalcSinusoidFreq();
      gr = my_fid.gr_time_series();
      my_fid.gr_freq_series().Draw();
      my_fid.f_fit().Draw("Same");
      c1.Print(TString::Format("%s_%.2f_%d.C", name.Data(), f, i));
      c1.Print(TString::Format("%s_%.2f_%d.png", name.Data(), f, i));
      c1.Print(TString::Format("%s.gif+20", name.Data()));

    }
  }

  return 0;
}