/*===========================================================================*\

author: Matthias W. Smith
email:  mwsmith2@uw.edu
date:   2015/01/09
file:   simple_fid_analyzer.cxx

notes: This program extracts the frequency of a single fid file and 
      prints the frequency and error.

usage:

./single_fid_analyzer <fid_data>

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <string>
#include <cassert>

//--- project includes ------------------------------------------------------//
#include "fid.h"

using namespace fid;

int main(int argc, char **argv)
{
  // Make sure a data file was specified and get it.
  assert(argc > 1);
  FID myfid(argv[1]);

  std::cout << "Frequency, Error" << std::endl;
  std::cout << myfid.GetFreq() << ", " << myfid.GetFreqError() << std::endl;

  return 0;
}
