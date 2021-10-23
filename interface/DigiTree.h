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

/**
 Reco format for pulse shape data.
 Structure:
 - Each branch is an array of size DigiTree#n_channels
 - Each channel name is stored in the DigiTree as a separate branch which serve as 
   a mapping between the mnemonic name and the index at which the channel information
   is stored in each reconstructed quantity branch.
   Example: `h4->Draw("amp_max[Ch0]") accesses `Ch0` value of `amp_max`.
 - The only exception is the time branch for which multiple values are stored for each
   channel following the DigiTree#time reconstruction methods specified in the WFAnalyzer plugin 
   options. For the DigiTree#time branch 
   Example: `h4->Draw("time[Ch0+CFD]") accesses the time computed with the `CFD` method of `Ch0`.
 */
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
    int*          channels;
    int*          time_types;
    /**
       The maximum gain among all samples recorded for a channel in each event.
     */
    float*        gain;
    /**
       Pedestal value (before baseline subtraction).
     */
    float*        pedestal;
    /**
       Integral over a window containing only baseline samples (the window is defined by the bIntWin option
       in the channel block of the config file)
     */
    float*        b_charge;
    /**
       Baseline slope from linear fit inside the baseline window (the window is defined by the bIntWin option
       in the channel block of the config file)
     */
    float*        b_slope;
    /**
       Noise RMS.
     */
    float*        b_rms;
    /**
       Time computed accordingly to all methods specified in the config file (supported: CFD, LED, TED)
     */
    float*        time;
    /**
       Chi2 of the linear interpolation used to estimate each time value.
     */
    float*        time_chi2;
    /**
       Uncertainty of the linear interpolation used to estimate each time value.
     */
    float*        time_error;    
    /**
       Slope of the linear interpolation used to estimate each time value.
     */
    float*        time_slope;
    /**
       For channel of type WFClassClock the measured period.
     */
    float*        period;    
    /**
       Maximum sample value (no fit)
     */
    float*        maximum;
    /**
       The time of the largest sample (no interpolation).
     */
    float*        time_maximum;
    /**
       The interpolated maximum amplitude (the interpolation funtion can either be polynomial or gaussian).
     */
    float*        amp_max;
    /**
       The interpolated time of the maximum (the interpolation funtion can either be polynomial or gaussian).
     */
    float*        time_max;
    /**
       The Chi2 of the interpolation to fit the maximum (the interpolation funtion can either be polynomial or gaussian).
     */
    float*        chi2_max;
    /**
       Sum over all samples
     */
    float*        charge_tot;
    /**
       Sum over all samples inside sIntWin
     */
    float*        charge_sig;
    /**
       The amplitude returned by the template fit.
     */
    float*        fit_ampl;
    /**
       The time returned by the template fit.
     */
    float*        fit_time;
    /**
       The time uncertainty returned by the template fit.
     */
    float*        fit_terr;    
    /**
       The Chi2 of the template fit.
     */
    float*        fit_chi2;
    /**
       The period returned by the template fit (for WFClassClock channels).
     */
    float*        fit_period;
    float*        fit_ampl_scint;
    float*        fit_time_scint;
    float*        fit_ampl_spike;
    float*        fit_time_spike;
    float*        fit_chi2_scint_plus_spike;
    bool*         fit_converged_scint_plus_spike;
    /**
       Amplitude calibration constant (multiplicative factor).
     */
    float*        ampl_calib;
    /**
       Time calibration constant (shift).
     */
    float*        time_calib;
};

#endif
