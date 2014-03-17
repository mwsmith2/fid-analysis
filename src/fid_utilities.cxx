#include "fid_utilities.h"

void fid::DrawFID(const vec &wf, 
                  const vec &tm, 
                  const string filename,
                  const string title)
{
  // Get our own TCanvas
  TCanvas c1;

  // Set up the graph
  TGraph gr(wf.size(), &tm[0], &wf[0]);
  string new_title(title);
  new_title.append("; time [ms]; amplitude [a.u.]");
  gr.SetTitle(title.c_str());

  // Draw the waveform
  gr.Draw();
  c1.Print(filename.c_str());

  return;
}

void fid::DrawFID(fid::FID &my_fid, const string filename, const string title)
{
  // Get our own TCanvas
  TCanvas c1;

  // Set up the graph
  int N = my_fid.wf().size();
  TGraph gr(N, &my_fid.tm()[0], &my_fid.wf()[0]);
  string new_title(title);
  new_title.append("; time [ms]; amplitude [a.u.]");
  gr.SetTitle(title.c_str());

  // Draw the waveform
  gr.Draw();
  c1.Print(filename.c_str());

  return;
}