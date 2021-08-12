#include "DigitizerReco.h"

//----------Utils-------------------------------------------------------------------------
bool DigitizerReco::Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index)
{
    //---inputs---
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");
    isLiTEDTU_ = opts.OptExist(instanceName_+".litedtu") ? opts.GetOpt<bool>(instanceName_+".litedtu") : false; 
    //---read dV and dT digitizer calibration from file and register it to the shared data
    if(opts.OptExist(instanceName_+".calibration"))
    {
        TFile* calibFile = TFile::Open(opts.GetOpt<string>(instanceName_+".calibration", 0).c_str(), "READ");
        if(calibFile)
        {
            auto calibPtr = calibFile->Get(opts.GetOpt<string>(instanceName_+".calibration", 1).c_str());
            if(calibPtr)
            {
                cout << ">>>DigiReco INFO: using calibration "
                     << opts.GetOpt<string>(instanceName_+".calibration", 1)
                     << " from " << opts.GetOpt<string>(instanceName_+".calibration", 0) << endl;
                
                digitizerCalib_ = *((DigitizerCalibration*)calibPtr->Clone("dVdTCalibration"));
                RegisterSharedData(&digitizerCalib_, "dVdTCalibration", false);
            }
            else
            {
                cout << ">>>DigiReco WARNING: calibration "
                     << opts.GetOpt<string>(instanceName_+".calibration", 1)
                     << " not found in " << opts.GetOpt<string>(instanceName_+".calibration", 0)
                     << ". No calibration will be applied" << endl;
            }
        }
        else
        {
            cout << ">>>DigiReco WARNING: impossible to open calibration file "
                 << opts.GetOpt<string>(instanceName_+".calibration", 0)
                 << ". No calibration will be applied" << endl;
        }
    }
    
    //---create channel container (WFClass)
    for(auto& channel : channelsNames_)
    {
        nSamples_[channel] = opts.GetOpt<int>(channel+".nSamples");
        auto tUnit = opts.GetOpt<float>(channel+".tUnit");
        if(opts.OptExist(channel+".type"))
        {
            if(opts.GetOpt<string>(channel+".type") == "NINO")
                WFs_[channel] = new WFClassNINO(opts.GetOpt<int>(channel+".polarity"), tUnit);
            else if(opts.GetOpt<string>(channel+".type") == "Clock")
                WFs_[channel] = new WFClassClock(tUnit);
            else if(isLiTEDTU_)
	      WFs_[channel] = new WFClassLiTEDTU(opts.GetOpt<int>(channel+".polarity"), tUnit);
        }
        else
	  WFs_[channel] = new WFClass(opts.GetOpt<int>(channel+".polarity"), tUnit);

        
        //---set channel calibration if available
        unsigned int digiBd = opts.GetOpt<unsigned int>(channel+".digiBoard");
        unsigned int digiGr = opts.GetOpt<unsigned int>(channel+".digiGroup");
        unsigned int digiCh = opts.GetOpt<unsigned int>(channel+".digiChannel");
        auto ch_key = make_tuple(digiBd, digiGr, digiCh);
        if(digitizerCalib_.find(ch_key) != digitizerCalib_.end())
            WFs_[channel]->SetCalibration(&digitizerCalib_[ch_key]);

        //---register WF in the shared data
        RegisterSharedData(WFs_[channel], channel, false);
    }
    
    return true;
}

bool DigitizerReco::ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{        
    //---read the digitizer
    //---set time reference from digitized trigger
    int trigRef=0;
    // FIXME
    // for(int iSample=nSamples_[channel]*8; iSample<nSamples_[channel]*9; ++iSample)
    // {
    //     if(event.digiSampleValue[iSample] < 1000)
    //     {
    //         trigRef = iSample-nSamples_[channel]*8;
    //         break;
    //     }
    // }
    
    //---user channels
    bool evtStatus = true;

    for(auto& channel : channelsNames_)
    {
        //---reset and read new WFs_
        WFs_[channel]->Reset();
        auto digiBd = opts.GetOpt<unsigned int>(channel+".digiBoard");
        auto digiGr = opts.GetOpt<unsigned int>(channel+".digiGroup");
        auto digiCh = opts.GetOpt<unsigned int>(channel+".digiChannel");
        auto offset = event.digiMap.at(make_tuple(digiBd, digiGr, digiCh));
        auto max_sample = offset+std::min(nSamples_[channel], event.digiNSamplesMap[make_tuple(digiBd, digiGr, digiCh)]); 
        auto iSample = offset;
        while(iSample < max_sample && event.digiBoard[iSample] != -1)
        {
            //Set the start index cell
            if (iSample==offset)
                WFs_[channel]->SetStartIndexCell(event.digiStartIndexCell[iSample]);

            //---H4DAQ bug: sometimes ADC value is out of bound.
            //---skip everything if one channel is bad
            if(event.digiSampleValue[iSample] > 1e6)
	      {
		evtStatus = false;
                WFs_[channel]->AddSample(4095);
	      }
            else if(isLiTEDTU_)
	      {	      
		WFs_[channel]->AddSample(event.digiSampleValue[iSample], event.digiSampleGain[iSample]);
	      }
	    else
	      WFs_[channel]->AddSample(event.digiSampleValue[iSample]); 

            iSample++;
        }
        if(opts.OptExist(channel+".useTrigRef") && opts.GetOpt<bool>(channel+".useTrigRef"))
            WFs_[channel]->SetTrigRef(trigRef);
        if(WFs_[channel]->GetCalibration())
            WFs_[channel]->ApplyCalibration();
    }

    if(!evtStatus)
        cout << ">>>DigiReco WARNING: bad amplitude detected" << endl;

    return evtStatus;
}
