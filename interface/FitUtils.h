#ifndef __FIT_UTILS__
#define __FIT_UTILS__

#include <iostream>
#include <cmath>

#include "TH1F.h"
#include "TTree.h"



/*** double crystalBall ***/
double crystalBallLowHigh(double* x, double* par);

/*** find effective sigma ***/
void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction, const bool& verbosity);

#endif
