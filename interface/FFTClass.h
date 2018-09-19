#ifndef __FFTCLASS_H__
#define __FFTCLASS_H__

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <algorithm>

#include "TObject.h"
#include "TMath.h"

using namespace std;

class FFTClass : public TObject
{
public:
    //---ctors---
    FFTClass();
    FFTClass(int N, double* Re, double* Im);
    FFTClass(vector<double>& Re, vector<double>& Im);
    
    //---dtor---
    ~FFTClass() {};

    //---getters---
    const vector<double>* Get(string var);
    const vector<double>* GetRe() {return &vars_["re"];};
    const vector<double>* GetIm() {return &vars_["im"];};
    const vector<double>* GetAmplitudes();
    const vector<double>* GetPhases();
    
    //---setters---
    void SetPointsComplex(int N, double* Re, double* Im);
    void SetPointsComplex(vector<double>& Re, vector<double>& Im);

    //---utils---
    void Reset();
    
private:
    map<string, vector<double> > vars_;
};

#endif
