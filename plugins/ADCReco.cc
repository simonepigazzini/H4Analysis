#include "ADCReco.h"

//----------Begin-------------------------------------------------------------------------
bool ADCReco::Begin(CfgManager& opts, uint64* index)
{
    chNames_ = opts.GetOpt<vector<string> >(GetInstanceName()+".chNames");
    nChannels_ = chNames_.size();
    chNum_.resize(nChannels_);
    values_ = new vector<float>;
    values_->resize(nChannels_);

    //---book and instanciate tree
    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;
    RegisterSharedData(new TTree("adc", "adc_tree"), "adc_tree", storeTree);
    adcTree_ = (TTree*)data_.back().obj;

    adcTree_->Branch("index", index, "index/l");
    adcTree_->Branch("adc_data", "std::vector<float>", &values_);
    for(unsigned int iCh=0; iCh<chNames_.size(); ++iCh)
    {
        chNum_[iCh] = iCh;
        adcTree_->Branch(chNames_[iCh].c_str(), &chNum_[iCh], (chNames_[iCh]+"/I").c_str());
    }

    return true;
}

//----------Event loop--------------------------------------------------------------------
bool ADCReco::ProcessEvent(const H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    for(unsigned int iADC=0; iADC<event.nAdcChannels; ++iADC)
    {
        for(int iCh=0; iCh<nChannels_; ++iCh)
        {
          if(opts.GetOpt<int>(chNames_[iCh]+".boardNum") == int(event.adcBoard[iADC]) &&
             opts.GetOpt<int>(chNames_[iCh]+".channelNum") == int(event.adcChannel[iADC]))
            {
                values_->at(iCh) = event.adcData[iADC];
                break;
            }
        }
    }
    adcTree_->Fill();
    
    return true;
}
