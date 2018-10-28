#ifndef __CALIB_DIGITIZER__
#define __CALIB_DIGITIZER__

#include "interface/PluginBase.h"
#include "interface/WFClass.h"
#include "TLinearFitter.h"

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
    bool Begin(CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(CfgManager& opts);

    bool BeginLoop(int iLoop, CfgManager& opts);
    bool EndLoop(int iLoop, CfgManager& opts);

private:
    void   minimize();

    int                         nTotSamples_;
    string                      srcInstance_;
    vector<string>              channelsNames_;

    map<string, WFClass*>       WFs_;
    vector<TLinearFitter>       fitters_;
    
    DigitizerCalibration        digiCalibration_;
};

DEFINE_PLUGIN(CalibDigitizer);

#endif
