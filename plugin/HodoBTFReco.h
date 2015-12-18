#ifndef __HODO_BTFRECO__
#define __HODO_BTFRECO__

#include "PluginBase.h"
#include "interface/PositionTree.h"

using namespace std;

class HodoBTFReco: public PluginBase
{
public:
    //---ctors---
    HodoBTFReco() {};

    //---dtor---
    ~HodoBTFReco() {};

    //---utils---
    bool Begin(CfgManager& opts, int* index);
    bool ProcessEvent(const H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    
private:
    map<int, int> ADC_to_PMT_map_;
    map<int, int> PMT_to_hodoX_map_;
    map<int, int> PMT_to_hodoY_map_;
    PositionTree  hodoTree_;

    //---datamember for registration
    PLUGIN_REGISTER(HodoBTFReco)
};

DEFINE_PLUGIN(HodoBTFReco);


#endif
