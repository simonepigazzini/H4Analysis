#ifndef __TRACK_RECO__
#define __TRACK_RECO__

#include "interface/PluginBase.h"
#include "interface/TrackTree.h"
#include "interface/Track.h"

using namespace std;


class TrackReco: public PluginBase
{
public:

    //---ctors---
    TrackReco() {};
  
    //---dtor---
    ~TrackReco() {};

   
    //---utils---
    bool Begin(CfgManager& opts, uint64* index);
    bool BeginLoop(int iLoop, CfgManager& opts);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    
private:

    void buildTracks();
    void cleanTracks();

    Tracking::TrackContainer tracks_;
    std::vector<string> hitProducers_;
    std::map<string,Tracking::LayerMeasurements*> hits_;
    Tracking::TelescopeLayout hodo_;
    TrackTree*     trackTree_;
    float          maxChi2_;
    float          cleaningChi2Cut_;
};

DEFINE_PLUGIN(TrackReco);

#endif
