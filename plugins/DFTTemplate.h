#ifndef __DFT_Template__
#define __DFT_Template__

#include <iostream>
#include <math.h>

#include "TVirtualFFT.h"

#include "interface/utils.h"
#include "interface/PluginBase.h"
#include "interface/WFClass.h"
#include "interface/FFTClass.h"

class DFTTemplate: public PluginBase
{
public:
    //---ctors---
    DFTTemplate(){};

    //---dtor---
    ~DFTTemplate(){};

    //---utils---
    bool Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(map<string, PluginBase*>& plugins, CfgManager& opts) {return true;};
    
private:    
    //---internal data
    uint64*                          index_;
    string                           srcInstance_;
    vector<string>                   channelsNames_;
    map<string, FFTClass*>           FFTs_;
    map<string, WFClass*>            WFs_;
    map<string, pair<float, float> > oversamplingMap_;
};

DEFINE_PLUGIN(DFTTemplate);

#endif
