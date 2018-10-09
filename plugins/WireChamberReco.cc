#include "WireChamberReco.h"

//**********Utils*************************************************************************
//----------Begin-------------------------------------------------------------------------
bool WireChamberReco::Begin(CfgManager& opts, uint64* index)
{
    //---get wire chamber mapping
    chXl_ = opts.GetOpt<int>(instanceName_+".chXleft");
    chXr_ = opts.GetOpt<int>(instanceName_+".chXright");
    chYu_ = opts.GetOpt<int>(instanceName_+".chYup");
    chYd_ = opts.GetOpt<int>(instanceName_+".chYdown");
    
    //---create a position tree
    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;
    RegisterSharedData(new TTree("wire", "wire_tree"), "wire_tree", storeTree);
    wireTree_ = PositionTree(index, (TTree*)data_.back().obj);
    wireTree_.Init();
    
    return true;
}

//----------ProcessEvent------------------------------------------------------------------
bool WireChamberReco::ProcessEvent(H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    //---search the first (in time) hit for each channel
    vector<float> timeL, timeR, timeU, timeD;
    for(unsigned int iCh=0; iCh<h4Tree.nTdcChannels; ++iCh)
    {
        if(int(h4Tree.tdcChannel[iCh])==chXl_)
            timeL.push_back(h4Tree.tdcData[iCh]);
        if(int(h4Tree.tdcChannel[iCh])==chXr_)
            timeR.push_back(h4Tree.tdcData[iCh]);
        if(int(h4Tree.tdcChannel[iCh])==chYu_)
            timeU.push_back(h4Tree.tdcData[iCh]);
        if(int(h4Tree.tdcChannel[iCh])==chYd_)
            timeD.push_back(h4Tree.tdcData[iCh]);
    }

    //---compute X and Y from channels times
    // wireTree_.n_hitsX = timeR.size()+timeL.size();
    // wireTree_.n_hitsY = timeD.size()+timeU.size();
    if(timeR.size()!=0 && timeL.size()!=0)
    {
        wireTree_.n_clusters_X = timeR.size()+timeL.size();
        wireTree_.cluster_X_size.push_back(1);
        wireTree_.X.push_back((*min_element(timeR.begin(), timeR.begin()+timeR.size()) -
                               *min_element(timeL.begin(), timeL.begin()+timeL.size()))*0.005);
    }
    else
    {
        wireTree_.n_clusters_X = 0;
        wireTree_.cluster_X_size.push_back(0);
        wireTree_.X.push_back(-1000);
    }
    if(timeU.size()!=0 && timeD.size()!=0)
    {
        wireTree_.n_clusters_Y = timeD.size()+timeU.size();
        wireTree_.cluster_Y_size.push_back(1.);
        wireTree_.Y.push_back((*min_element(timeU.begin(), timeU.begin()+timeU.size()) -
                               *min_element(timeD.begin(), timeD.begin()+timeD.size()))*0.005);
    }
    else
    {
        wireTree_.n_clusters_Y = 0;
        wireTree_.cluster_Y_size.push_back(0);
        wireTree_.Y.push_back(-1000);
    }

    //---fill output tree
    wireTree_.Fill();

    //---clear tree
    wireTree_.X.clear();
    wireTree_.cluster_X_size.clear();
    wireTree_.n_clusters_X=0;
    wireTree_.Y.clear();
    wireTree_.cluster_Y_size.clear();
    wireTree_.n_clusters_Y=0;
    
    return true;
}
