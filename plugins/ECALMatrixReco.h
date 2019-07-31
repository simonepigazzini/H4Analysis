#ifndef __ECAL_MATRIX_RECO__
#define __ECAL_MATRIX_RECO__

#include <iostream>

#include "interface/PluginBase.h"
#include "interface/WFClass.h"
#include "interface/RecoEventAnalyzer.h"
#include "interface/ECALMatrixTree.h"

class ECALMatrixReco: public PluginBase
{
public:
    //---ctors---
    ECALMatrixReco(){};

    //---dtor---
    ~ECALMatrixReco(){};

    //---utils---
    bool Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(map<string, PluginBase*>& plugins, CfgManager& opts) { return true; };
    
private:    
    //---internal data
    vector<string>               channels_;
    map<string, pair<int, int> > crystalPosition_;
    ECALMatrixTree               outTree_;
    RecoEventAnalyzer            eventAnalyzer_;
};

DEFINE_PLUGIN(ECALMatrixReco);

#endif
