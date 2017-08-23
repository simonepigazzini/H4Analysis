#include "TOFPETReco.h"

//**********Utils*************************************************************************
//----------Begin-------------------------------------------------------------------------
bool TOFPETReco::Begin(CfgManager& opts, uint64* index)
{
    //---inputs---
    auto fileName = opts.GetOpt<string>(instanceName_+".inputFile");
    inputFile_ = TFile::Open(fileName.c_str());
    inputTree_ = (TTree*)inputFile_->Get("data");
    rawTree_ = new TOFPETRawTree(inputTree_);
    rawTree_->NextEntry();
    while(rawTree_->channelID2 != 896)
        rawTree_->NextEntry();
    
    //---create a reco tree
    RegisterSharedData(new TTree("tofp", "tofpet_tree"), "tofpet_tree", true);
    recoTree_ = TOFPETRecoTree(index, (TTree*)data_.back().obj);
    recoTree_.Init();

    return true;
}

//----------ProcessEvent------------------------------------------------------------------
bool TOFPETReco::ProcessEvent(const H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    double h4daq_time = -1;
    for(int iT=0; iT<h4Tree.nEvtTimes; ++iT)
    { 
        if(h4daqRefTime_ == -1 && h4Tree.evtTimeBoard[iT] == 16842753)
            h4daqRefTime_ = h4Tree.evtTime[iT]; 
        if(h4Tree.evtTimeBoard[iT] == 16842753)
            h4daq_time = h4Tree.evtTime[iT]-h4daqRefTime_;
    }

    if(tofpetRefTime_ == -1)
        tofpetRefTime_ = rawTree_->time2/1e6;
    double tofpet_time = rawTree_->time2/1e6-tofpetRefTime_;   
    if(currentSpill_ != h4Tree.spillNumber)
    {
        if(fabs(tofpet_time - h4daq_time -spillAdjust_)<100)
        {
            currentSpill_ = h4Tree.spillNumber;
            spillAdjust_ = tofpet_time - h4daq_time;
        }
        else
        {
            recoTree_.t_h4daq = 9999;
            recoTree_.t_tofpet = 0;
            recoTree_.Fill();
            return false;                       
        }
    }

    bool matched=false;
    double time_diff = tofpet_time - h4daq_time - spillAdjust_;
    if(fabs(time_diff)<50)
    {
        recoTree_.t_sipm = rawTree_->time1-rawTree_->time2;
        recoTree_.tot = rawTree_->tot1;
        recoTree_.energy = rawTree_->energy1;
        recoTree_.t_h4daq = h4daq_time;
        recoTree_.t_tofpet = tofpet_time;
        matched = true;
        rawTree_->NextEntry();
        while(rawTree_->channelID2 != 896)            
            rawTree_->NextEntry();    
    }
    else
    {
        recoTree_.t_sipm = 9999;
        recoTree_.tot = 9999;
        recoTree_.energy = -1;
        recoTree_.t_h4daq = 9999;
        recoTree_.t_tofpet = 0;        
    }
    recoTree_.Fill();
    
    return matched;
}
