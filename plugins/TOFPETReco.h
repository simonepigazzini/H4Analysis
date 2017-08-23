#ifndef __TOFPET_RECO__
#define __TOFPET_RECO__

#include <iostream>

#include "interface/PluginBase.h"
#include "interface/TOFPETRawTree.h"
#include "interface/TOFPETRecoTree.h"

class TOFPETReco: public PluginBase
{
public:
    //---ctors---
    TOFPETReco():
        currentSpill_(-1), spillAdjust_(-1), h4daqRefTime_(-1), tofpetRefTime_(-1)
        {};

    //---dtor---
    ~TOFPETReco(){};

    //---utils---
    bool Begin(CfgManager& opts, uint64* index);
    bool ProcessEvent(const H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    
private:    
    //---internal data
    TFile*         inputFile_;
    TTree*         inputTree_;
    TOFPETRawTree* rawTree_;
    TOFPETRecoTree recoTree_;
    
    int            currentSpill_;
    double         spillAdjust_;
    double         h4daqRefTime_;
    double         tofpetRefTime_;
};

DEFINE_PLUGIN(TOFPETReco);

#endif
