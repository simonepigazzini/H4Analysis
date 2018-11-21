#ifndef __MAKE_COVARIANCE_MATRIX__
#define __MAKE_COVARIANCE_MATRIX__

#include <iostream>

#include "TH2F.h"

#include "interface/WFClass.h"
#include "interface/PluginBase.h"

class MakeCovarianceMatrix: public PluginBase
{
public:
    //---ctors---
    MakeCovarianceMatrix(){};

    //---dtor---
    ~MakeCovarianceMatrix(){};

    //---utils---
    bool Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index);
    bool ProcessEvent(const H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool EndLoop(int iLoop, map<string, PluginBase*>& plugins, CfgManager& opts);
    
private:    
    //---internal data
    int                                    events_;
    string                                 digiInstance_;             
    vector<string>                         channelsNames_;
    map<string, TH2F>                      mapCovariances_;
    map<string, vector<float> >            sums_;
    map<string, vector<float> >            sum2s_;
    map<string, map<int, vector<float> > > values_;
};

DEFINE_PLUGIN(MakeCovarianceMatrix);

#endif
