#ifndef __CALIB_DIGITIZER__
#define __CALIB_DIGITIZER__

#include "interface/PluginBase.h"
#include "interface/WFClass.h"
#include "interface/CalibDigiTree.h"
#include "TLinearFitter.h"
#include "TGraph.h"

using namespace std;


class CalibDigitizer: public PluginBase
{
public:

    //---ctors---
    CalibDigitizer() {};
  
    //---dtor---
    ~CalibDigitizer() 
        {
        };
   
    //---utils---
    bool Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(map<string, PluginBase*>& plugins, CfgManager& opts);

    bool BeginLoop(int iLoop, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool EndLoop(int iLoop, map<string, PluginBase*>& plugins, CfgManager& opts);

private:
    void   minimize();

    int                   nTotSamples_;
    string                srcInstance_;
    vector<string>        channelsNames_;

    bool                  fitDeltaV_;
    bool                  fitSlopeV_;
    bool                  fitQuadraticV_;
    bool                  fitDeltaT_;
    string                functionType_;

    map<string, WFClass*> WFs_;
    vector<TLinearFitter> fitters_;
    
    DigitizerCalibration  digiCalibration_;
    CalibDigiTree*        calibTree_;
    
    double                prevFit_;
    double                thisFit_;
};

DEFINE_PLUGIN(CalibDigitizer);

#endif
