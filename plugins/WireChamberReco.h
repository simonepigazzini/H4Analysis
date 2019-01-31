#ifndef __WIRE_CHMABER_RECO__
#define __WIRE_CHMABER_RECO__

#include "interface/PluginBase.h"
#include "interface/PositionTree.h"

using namespace std;

class WireChamberReco: public PluginBase
{
public:
    //---ctors---
    WireChamberReco() {};

    //---dtor---
    ~WireChamberReco() {};

    //---utils---
    bool Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(map<string, PluginBase*>& plugins, CfgManager& opts) { return true; };
    
private:
    PositionTree     wireTree_;
    int              chXl_;
    int              chXr_;
    int              chYu_;
    int              chYd_;
};

DEFINE_PLUGIN(WireChamberReco);

#endif
