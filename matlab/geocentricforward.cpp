/**
 * \file geocentricforward.cpp
 * \brief Matlab mex file for geographic to UTM/UPS conversions
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed under
 * the LGPL.  For more information, see http://geographiclib.sourceforge.net/
 **********************************************************************/

// Compile in Matlab with
// [Unix]
// mex -I/usr/local/include -L/usr/local/lib -Wl,-rpath=/usr/local/lib -lGeographic geocentricforward.cpp
// [Windows]
// mex -I../include -L../windows/Release -lGeographic geocentricforward.cpp

// $Id$

#include <algorithm>
#include <GeographicLib/Geocentric.hpp>
#include <mex.h>

using namespace std;
using namespace GeographicLib;

void mexFunction( int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[] ) {

  if (nrhs < 1)
    mexErrMsgTxt("One input argument required.");
  else if (nrhs > 3)
    mexErrMsgTxt("More than three input arguments specified.");
  else if (nrhs == 2)
    mexErrMsgTxt("Must specify flattening with the major radius.");
  else if (nlhs > 2)
    mexErrMsgTxt("More than two output arguments specified.");

  if (!( mxIsDouble(prhs[0]) && !mxIsComplex(prhs[0]) ))
    mexErrMsgTxt("geodetic coordinates are not of type double.");

  if (mxGetN(prhs[0]) != 3)
    mexErrMsgTxt("geodetic coordinates must be M x 3 matrix.");

  double a = Constants::WGS84_a<double>(), f = Constants::WGS84_f<double>();
  if (nrhs == 3) {
    if (!( mxIsDouble(prhs[1]) && !mxIsComplex(prhs[1]) &&
           mxGetNumberOfElements(prhs[1]) == 1 ))
      mexErrMsgTxt("Major radius is not a real scalar.");
    a = mxGetScalar(prhs[1]);
    if (!( mxIsDouble(prhs[2]) && !mxIsComplex(prhs[2]) &&
           mxGetNumberOfElements(prhs[2]) == 1 ))
      mexErrMsgTxt("Flattening is not a real scalar.");
    f = mxGetScalar(prhs[2]);
  }

  int m = mxGetM(prhs[0]);

  double* lat = mxGetPr(prhs[0]);
  double* lon = lat + m;
  double* h = lat + 2*m;

  plhs[0] = mxCreateDoubleMatrix(m, 3, mxREAL);
  double* x = mxGetPr(plhs[0]);
  std::fill(x, x + 3*m, Math::NaN<double>());
  double* y = x + m;
  double* z = x + 2*m;
  double* rot = NULL;
  bool rotp = nlhs == 2;

  if (rotp) {
    plhs[1] = mxCreateDoubleMatrix(m, 9, mxREAL);
    rot = mxGetPr(plhs[1]);
    std::fill(rot, rot + 9*m, Math::NaN<double>());
  }

  try {
    std::vector<double> rotv(rotp ? 9 : 0);
    const Geocentric c(a, f);
    for (int i = 0; i < m; ++i) {
      if (!(abs(lat[i]) > 90) && !(lon[i] < -180 || lon[i] > 360)) {
        c.Forward(lat[i], lon[i], h[i], x[i], y[i], z[i], rotv);
        if (rotp) {
          for (int k = 0; k < 9; ++k)
            rot[m * k + i] = rotv[k];
        }
      }
    }
  }
  catch (const std::exception& e) {
    mexErrMsgTxt(e.what());
  }
}