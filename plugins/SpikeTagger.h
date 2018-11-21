#ifndef __SPIKE_TAGGER__
#define __SPIKE_TAGGER__

#include <iostream>
#include <algorithm>
#include <iterator>

#include "interface/PluginBase.h"
#include "interface/WFTree.h"
#include "interface/SpikesTree.h"
#include "interface/WFClass.h"

class SpikeTagger: public PluginBase
{
public:
    //---ctors---
    SpikeTagger(){};

    //---dtor---
    ~SpikeTagger(){};

    //---utils---
    bool Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(map<string, PluginBase*>& plugins, CfgManager& opts) { return true; };
    
private:    
    //---internal data
    string                      srcInstance_;
    vector<string>              channelsNames_;
    map<string, vector<string>> channelsNamesSwissCross_;
    map<string, vector<string>> channelsNames3By3_;
    SpikesTree                  spikesTree_;
    WFTree                      outWFTree_;
    map<string, WFClass*>       WFs_;
};

DEFINE_PLUGIN(SpikeTagger);

#endif
