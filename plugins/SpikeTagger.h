#ifndef __SPIKE_TAGGER__
#define __SPIKE_TAGGER__

#include <iostream>

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
    bool Begin(CfgManager& opts, uint64* index);
    bool ProcessEvent(const H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(CfgManager& opts) { return true; };
    
private:    
    //---internal data
    string                      srcInstance_;
    vector<string>              channelsNames_;
    SpikesTree                  spikesTree_;
    WFTree                      outWFTree_;
    map<string, WFClass*>       WFs_;
};

DEFINE_PLUGIN(SpikeTagger);

#endif
