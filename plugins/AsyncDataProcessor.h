#ifndef __ASYNC_DATA_PROCESSOR__
#define __ASYNC_DATA_PROCESSOR__

#include "TFile.h"
#include "TTreeFormula.h"

#include "interface/PluginLoader.h"
#include "interface/PluginBase.h"
#include "interface/H4Tree.h"

class AsyncDataProcessor: public PluginBase
{
public:
    //---ctors---
    AsyncDataProcessor();

    //---dtor---
    ~AsyncDataProcessor(){};

    //---utils---
    bool Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(map<string, PluginBase*>& plugins, CfgManager& opts) { return true; };
    
private:    
    //---internal data
    PluginLoader<PluginBase>* loader_;
    vector<PluginLoader<PluginBase>* > pluginLoaders_;    
    map<string, PluginBase*> pluginMap_;
    vector<PluginBase*> pluginSequence_;
    vector<string> asyncPluginList_;

    TFile*        asyncDataFile_;
    H4Tree*       h4Tree_;
    TTreeFormula* dataSelector_;
    int           currentSpill_;
    double        deltaT_;
    double        syncTolerance_;
    int           maxForwardTries_;
};

DEFINE_PLUGIN(AsyncDataProcessor);

#endif
