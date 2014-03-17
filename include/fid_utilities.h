#ifndef __FID_UTILITIES_H__
#define __FID_UTILITIES_H__

//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <vector>
using std::vector;
using std::string;

//--- root includes ---------------------------------------------------------//
#include "TCanvas.h"
#include "TGraph.h"

//--- project includes ------------------------------------------------------//
#include "fid_class.h"
#include "fid.h"

typedef vector<double> vec;

namespace fid{

  void DrawFID(const vec &wf, const vec &tm, 
    const string title, const string filename);

  void DrawFID(fid::FID &my_fid, const string title, const string filename);

} // fid

#endif