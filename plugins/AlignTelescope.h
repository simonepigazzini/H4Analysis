#ifndef __ALIGN_TELESCOPE__
#define __ALIGN_TELESCOPE__

#include "interface/PluginBase.h"
#include "interface/Track.h"

using namespace std;


class AlignTelescope: public PluginBase
{
public:

    //---ctors---
    AlignTelescope() {};
  
    //---dtor---
    ~AlignTelescope() 
        {
        };

   
    //---utils---
    bool Begin(CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(CfgManager& opts);

    bool BeginLoop(int iLoop, CfgManager& opts);
    bool EndLoop(int iLoop, CfgManager& opts);

private:

    double globalChi2(const double* par);
    void   minimize();

    string srcInstance_;
    bool   alignZ_;
    Tracking::TelescopeLayout* tLayout_; //pointer to current layout 
    Tracking::TrackContainer tracks_;
};

DEFINE_PLUGIN(AlignTelescope);

#endif
