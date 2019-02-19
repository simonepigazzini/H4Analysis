#ifndef __INFO_TREE_MAKER__
#define __INFO_TREE_MAKER__

#include "interface/PluginBase.h"

class InfoTreeMaker: public PluginBase
{
public:
    //---ctors---
    InfoTreeMaker(){};

    //---dtor---
    ~InfoTreeMaker(){};

    //---utils---
    bool Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(map<string, PluginBase*>& plugins, CfgManager& opts) { return true; };
    
private:    
    //---internal data
    uint64* index_;
    TTree*  info_tree_;
    string                            trackedVariable_;
    map<string, float*>               mappedVariablesNum_;
    map<string, string*>              mappedVariablesStr_;    
    map<string, map<string, float> >  remapNum_;
    map<string, map<string, string> > remapStr_;    
};

DEFINE_PLUGIN(InfoTreeMaker);

#endif
