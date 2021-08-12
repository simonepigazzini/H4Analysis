#ifndef __DIGITIZER_RECO__
#define __DIGITIZER_RECO__

#include <iostream>

#include "interface/PluginBase.h"
#include "interface/DigiTree.h"
#include "interface/WFTree.h"
#include "interface/WFClass.h"
#include "interface/WFClassLiTEDTU.h"
#include "interface/WFClassNINO.h"
#include "interface/WFClassClock.h"
#include "interface/WFViewer.h"

class DigitizerReco: public PluginBase
{
public:
    //---ctors---
    DigitizerReco(){};

    //---dtor---
    ~DigitizerReco(){};

    //---utils---
    bool Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(map<string, PluginBase*>& plugins, CfgManager& opts) { return true; };
    
private:    
    //---internal data
    map<string, int>            nSamples_;
    vector<string>              channelsNames_;
    map<string, WFClass*>       WFs_;
    bool                        isLiTEDTU_;
    DigitizerCalibration        digitizerCalib_; 
};

DEFINE_PLUGIN(DigitizerReco);

#endif
