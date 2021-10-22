#ifndef __DIGI_TREE__
#define __DIGI_TREE__

#include <memory>
#include <iostream>
#include <string>
#include <vector>

#include "TTree.h"
#include "TString.h"

using namespace std;

typedef unsigned long int uint32;
typedef unsigned long long int uint64;
 
//****************************************************************************************

class DigiTree
{
public: 
    //---ctors---
    DigiTree(){};
    DigiTree(uint64* idx, TTree* tree=NULL, string prefix="");
    //---dtor---
    ~DigiTree();

    //---utils---
    void Init(vector<string>& names, vector<string>& timetypes);
    void Fill() {tree_->Fill();};
    void FillVoidChannel(int ch=0);
    
    TTree*        tree_; 
    string        prefix_;

    uint64*       index;
    unsigned int  n_channels;
    unsigned int  n_times;
    float*        gain;
    float*        pedestal;
    float*        b_charge;
    float*        b_slope;
    float*        b_rms;
    float*        time;
    float*        time_chi2;
    float*        time_error;    
    float*        time_slope;
    float*        period;    
    float*        maximum;
    float*        time_maximum;
    float*        amp_max;
    float*        time_max;
    float*        chi2_max;
    float*        charge_tot;
    float*        charge_sig;
    float*        fit_ampl;
    float*        fit_time;
    float*        fit_terr;    
    float*        fit_chi2;
    float*        fit_period;
    float*        fit_ampl_scint;
    float*        fit_time_scint;
    float*        fit_ampl_spike;
    float*        fit_time_spike;
    float*        fit_chi2_scint_plus_spike;
    bool*         fit_converged_scint_plus_spike;
    int*          channels;
    int*          time_types;
    float*        ampl_calib;
    float*        time_calib;
};

#endif
