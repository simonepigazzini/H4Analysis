#include "SpikeTagger.h"

//----------Utils-------------------------------------------------------------------------
bool SpikeTagger::Begin(CfgManager& opts, uint64* index)
{
    //---inputs---
    if(!opts.OptExist(instanceName_+".srcInstanceName"))
    {
        cout << ">>> SpikeTagger ERROR: no source plugin specified" << endl;
        return false;
    }
    srcInstance_ = opts.GetOpt<string>(instanceName_+".srcInstanceName");
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");
    
    //---outputs---    
    string spikesTreeName = opts.OptExist(instanceName_+".spikeTreeName") ?
        opts.GetOpt<string>(instanceName_+".spikeTreeName") : "spikes";
    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;
    RegisterSharedData(new TTree(spikesTreeName.c_str(), "spikes_tree"), "spikes_tree", storeTree);
    spikesTree_ = SpikesTree(index, (TTree*)data_.back().obj);
    spikesTree_.Init(channelsNames_);
    if(opts.GetOpt<int>(instanceName_+".fillWFtree"))
    {
        int nSamples = 0;
        if(opts.OptExist(instanceName_+".storeNSampleAroundMax"))
            nSamples = channelsNames_.size() * opts.GetOpt<int>(instanceName_+".storeNSampleAroundMax");
        else
        {
            for(auto& channel : channelsNames_)
                nSamples += opts.GetOpt<int>(channel+".nSamples");
        }
        
        string wfTreeName = opts.OptExist(instanceName_+".wfTreeName") ?
            opts.GetOpt<string>(instanceName_+".wfTreeName") : "spikes_wf";
        RegisterSharedData(new TTree(wfTreeName.c_str(), "spikes_wf_tree"), "spikes_wf_tree", true);
        outWFTree_ = WFTree(nSamples, index, (TTree*)data_.back().obj);
        outWFTree_.Init();
    }

    return true;
}

bool SpikeTagger::ProcessEvent(const H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    //---setup output event 
    int outCh=0;
    bool fillWFtree=false;
    if(opts.GetOpt<int>(instanceName_+".fillWFtree"))
        fillWFtree = *spikesTree_.index % opts.GetOpt<int>(instanceName_+".WFtreePrescale") == 0;

    //---load WFs from source instance shared data
    for(auto& channel : channelsNames_)
    {
        auto shared_data = plugins[srcInstance_]->GetSharedData(srcInstance_+"_"+channel, "", false);
        if(shared_data.size() != 0)
            WFs_[channel] = (WFClass*)shared_data.at(0).obj;
        else
            cout << "[SpikeTagger::" << instanceName_ << "]: channels samples not found check DigiReco step" << endl; 
    }

    //---compute reco variables
    for(auto& channel : channelsNames_)
    {
        //---skip dead channels
        if(WFs_.find(channel) == WFs_.end())
        {
            ++outCh;
            continue;
        }

        //---Look for undershoot after maximum
        auto undershoot_window = opts.GetOpt<int>(instanceName_+".undershootFinderWindow");
        auto analizedWF = WFs_[channel]->GetSamples();
        int max_sample = WFs_[channel]->GetTimeCF(1).first/WFs_[channel]->GetTUnit();
        if(max_sample+undershoot_window < analizedWF->size())
        {
            auto undershoot_sample = std::min_element(analizedWF->begin() + max_sample,
                                                      analizedWF->begin() + max_sample + undershoot_window);        
            spikesTree_.undershoot[outCh] = *undershoot_sample;
        }
        else
            spikesTree_.undershoot[outCh] = 1e5;

        //---compute A/A_matrix locking the phase
        if(spikesTree_.max_hit == -1 || WFs_[channelsNames_[spikesTree_.max_hit]]->GetAmpMax() < WFs_[channel]->GetAmpMax())
            spikesTree_.max_hit = outCh;
        float matrix_A_sum = 0.;
        for(auto other_ch : channelsNames_)
            matrix_A_sum += WFs_[other_ch]->GetSamples()->at(max_sample);

        spikesTree_.amp_sum_matrix[outCh] = matrix_A_sum;

        //---increase output tree channel counter
        ++outCh;        
    }


    //---WFs---
    if(fillWFtree)
    {
        outCh = 0;
        //---fix WF window around maximum sample of largest hit in the event.
        //   also assuming all channels have the same sampling frequency (i.e. don't mix VFEs with V1742 channels)
        float tUnit = WFs_[channelsNames_[spikesTree_.max_hit]]->GetTUnit();
        int max_sample = WFs_[channelsNames_[spikesTree_.max_hit]]->GetTimeCF(1).first/tUnit;            
        for(auto& channel : channelsNames_)
        {
            //---skip dead channels
            if(WFs_.find(channel) == WFs_.end())
            {
                ++outCh;
                continue;
            }
            auto analizedWF = WFs_[channel]->GetSamples();
            //---symmetric window (FIXME: odd windows)
            unsigned int firstSample = max_sample - opts.GetOpt<int>(instanceName_+".storeNSampleAroundMax")/2;
            unsigned int lastSample = max_sample + opts.GetOpt<int>(instanceName_+".storeNSampleAroundMax")/2;
            for(unsigned int jSample=firstSample; jSample<lastSample; ++jSample)
            {
                outWFTree_.WF_ch.push_back(outCh);
                outWFTree_.WF_time.push_back(jSample*tUnit);
                if(jSample>=0 && jSample<analizedWF->size())
                    outWFTree_.WF_val.push_back(analizedWF->at(jSample));
                else
                    outWFTree_.WF_val.push_back(0.);
            }
            //---increase output tree channel counter
            ++outCh;        
        }
    }

    //---fill the output trees and clean
    //---reco var
    spikesTree_.Fill();
    spikesTree_.max_hit = -1;
    
    //---WFs
    if(fillWFtree)
        outWFTree_.Fill();

    return true;
}
