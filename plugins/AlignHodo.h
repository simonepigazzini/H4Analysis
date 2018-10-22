#ifndef __TRACK_RECO__
#define __TRACK_RECO__

#include "interface/PluginBase.h"
#include "interface/Track.h"

using namespace std;


class AlignHodo: public PluginBase
{
public:

    //---ctors---
    AlignHodo() {};
  
    //---dtor---
    ~AlignHodo() 
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
    Tracking::TelescopeLayout* hodo_; //pointer to current layout 
    Tracking::TrackContainer tracks_;
};

DEFINE_PLUGIN(AlignHodo);

#endif
