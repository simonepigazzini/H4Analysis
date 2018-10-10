#include "TrackReco.h"

//**********Utils*************************************************************************
//----------Begin-------------------------------------------------------------------------
bool TrackReco::Begin(CfgManager& opts, uint64* index)
{
    //---create a position tree
//    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
//        opts.GetOpt<bool>(instanceName_+".storeTree") : true;

    return true;
}

//----------ProcessEvent------------------------------------------------------------------
bool TrackReco::ProcessEvent(H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{

    tracks_.clear();
    //---fill output tree
    //fitpixTree_->Fill();

    return true;
}
