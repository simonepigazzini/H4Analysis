#ifndef __TEMPLATE_MAKER__
#define __TEMPLATE_MAKER__

#include "TLinearFitter.h"
#include "TProfile.h"

#include "interface/PluginBase.h"
#include "interface/WFClass.h"
#include "interface/RecoEventAnalyzer.h"

using namespace std;


class TemplateMaker: public PluginBase
{
public:

    //---ctors---
    TemplateMaker() {};
  
    //---dtor---
    ~TemplateMaker() {};
   
    //---utils---
    bool Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    //bool End(map<string, PluginBase*>& plugins, CfgManager& opts);

    bool BeginLoop(int iLoop, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool EndLoop(int iLoop, map<string, PluginBase*>& plugins, CfgManager& opts);

private:
    string                srcInstance_;
    vector<string>        channelsNames_;
    bool                  storeRatios_;
    RecoEventAnalyzer     eventAnalyzer_;
    map<string, WFClass*> WFs_; 
    map<string, TProfile> hTemplates_;
    map<string, TH1D>     hRatios_;
    
};

DEFINE_PLUGIN(TemplateMaker);

#endif
