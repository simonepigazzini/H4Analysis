#include "CalibDigitizer.h"
#include "TF1.h"
#include "TVectorD.h"

#include <time.h>
#include <regex>

//----------Utils-------------------------------------------------------------------------
bool CalibDigitizer::Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index)
{
    prevFit_ = 0.;
    thisFit_ = 0.;
    //---inputs---
    if(!opts.OptExist(instanceName_+".srcInstanceName"))
    {
        Log("no source plugin specified", ERR);
        return false;
    }

    srcInstance_ = opts.GetOpt<string>(instanceName_+".srcInstanceName");
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");

    nTotSamples_ = 0;
    for(auto& channel : channelsNames_)
        nTotSamples_ += opts.GetOpt<int>(channel+".nSamples");

    auto shared_calib = plugins[srcInstance_]->GetSharedData("dVdTCalibration", "", false);
    if(shared_calib.size() == 0 || !shared_calib[0].obj)
    {
        //---create empty calibration
        for(auto& channel : channelsNames_)
        {
            int nSamples = opts.GetOpt<int>(channel+".nSamples");
            bgc_key_t channelID(opts.GetOpt<int>(channel+".digiBoard"), opts.GetOpt<int>(channel+".digiGroup"), opts.GetOpt<int>(channel+".digiChannel"));
            for(int iSample=0; iSample<nSamples; ++iSample)
                digiCalibration_[channelID].push_back({0.,0.,0.});
        }
    }
    else
    {
        digiCalibration_ = *((DigitizerCalibration*)shared_calib[0].obj);
    }

    RegisterSharedData(&digiCalibration_, "digiCalib", true);

    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;

    string treeName = opts.OptExist(instanceName_+".treeName") ?
        opts.GetOpt<string>(instanceName_+".treeName") : "calib_tree";

    RegisterSharedData(new TTree(treeName.c_str(), treeName.c_str()), treeName.c_str(), storeTree);
  
    calibTree_ = new CalibDigiTree(index, (TTree*)data_.back().obj);
    calibTree_->Init();

    fitDeltaV_ = opts.OptExist(instanceName_+".fitDeltaV") ?
        opts.GetOpt<bool>(instanceName_+".fitDeltaV") : true;

    fitSlopeV_ = opts.OptExist(instanceName_+".fitSlopeV") ?
        opts.GetOpt<bool>(instanceName_+".fitSlopeV") : true;

    fitQuadraticV_ = opts.OptExist(instanceName_+".fitQuadraticV") ?
        opts.GetOpt<bool>(instanceName_+".fitQuadraticV") : true;

    fitDeltaT_ = opts.OptExist(instanceName_+".fitDeltaT") ?
        opts.GetOpt<bool>(instanceName_+".fitDeltaT") : true;

    functionType_ =  opts.GetOpt<string>(instanceName_+".functionType");

    return true;
}


//----------Begin-------------------------------------------------------------------------
bool CalibDigitizer::BeginLoop(int iLoop, map<string, PluginBase*>& plugins, CfgManager& opts)
{  
    prevFit_ = thisFit_;
    thisFit_ = 0;

    fitters_.clear();
    for (int i = 0;i<nTotSamples_;++i)
    {
        fitters_.push_back(TLinearFitter(4)); //using a linear fitter with 3 parameters for each cell i % DeltaV_i + DeltaT_i*Derivative + + SlopeV_i*Sample + QuadraticV_i*sample*sample %
        fitters_[i].SetFormula("hyp3"); 
    }
  
    return true;
}

//----------ProcessEvent------------------------------------------------------------------
bool CalibDigitizer::ProcessEvent(H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    calibTree_->Clear();
  
    //---load WFs from source instance shared data
    for(auto& channel : channelsNames_)
    {
        auto shared_data = plugins[srcInstance_]->GetSharedData(srcInstance_+"_"+channel, "", false);
        if(shared_data.size() != 0)
            WFs_[channel] = (WFClass*)shared_data.at(0).obj;
        else
            Log("channels samples not found check DigiReco step", WARN); 
    }

    //---compute reco variables
    int offset = 0;
    for(auto& channel : channelsNames_)
    {
        //---skip dead channels
        if(WFs_.find(channel) == WFs_.end())
            continue;

        bgc_key_t channelID(opts.GetOpt<int>(channel+".digiBoard"), opts.GetOpt<int>(channel+".digiGroup"), opts.GetOpt<int>(channel+".digiChannel"));
        WFs_[channel]->SetCalibration(&digiCalibration_[channelID]);
        WFs_[channel]->ApplyCalibration();

        auto analizedWF = WFs_[channel]->GetSamples();
        float tUnit = WFs_[channel]->GetTUnit();

        TF1* wFit = NULL;
        int ampParameter = -1;
        if (functionType_ == "sin")
        {
            //Set function and fit initial values and ranges
            wFit = new TF1("wFit","[0]+[1]*TMath::Sin([2]*(x+[3]))",0,analizedWF->size()*tUnit);
            wFit->SetParameter(0,2400);
            wFit->SetParLimits(0,1800.,2600.);
            wFit->SetParameter(1,1350); 
            wFit->SetParLimits(1,900,1900);
            wFit->SetParameter(2,0.62); //80MHz 
            wFit->SetParLimits(2,0.5,0.7);
            wFit->SetParameter(3,0);
            wFit->SetParLimits(3,-50.,50.);
            ampParameter = 1;
        }
        else if (functionType_ == "pol0")
        {
            //Set function and fit initial values and ranges
            wFit = new TF1("wFit","[0]",0,analizedWF->size()*tUnit);
            wFit->SetParameter(0,2000);
            wFit->SetParLimits(0,0.,4096.);
            ampParameter = 0;
        }

        if (!wFit)
        {
            Log("Unknown function", ERR);
            return false;
        }
	  
        double chi2 = WFs_[channel]->AnalyticFit(wFit,10,analizedWF->size()-50); //fit using most of the samples avoiding the last samples of the digitizer

        if(chi2>50000 || wFit->GetParameter(ampParameter)<100) //removing bad fits or empty channel
            continue;

        thisFit_ += chi2;

        //get residual and slope and add it to the linear minimizer
        auto samples = WFs_[channel]->GetSamples();
        auto times = WFs_[channel]->GetTimes();
        for (int iSample = 0;iSample<samples->size();++iSample)
        {
            double x[3];
            x[0] = wFit->Derivative((*times)[iSample]);
            x[1] = wFit->Eval((*times)[iSample]);
            x[2] = x[1]*x[1];
            double residual = wFit->Eval((*times)[iSample])-(*samples)[iSample];
            int cellIndex = (iSample + WFs_[channel]->GetStartIndexCell())%samples->size();
            fitters_[offset+cellIndex].AddPoint(x,residual);
        }	

        //---Fill calibTree
        // calibTree_->n_calibWf++;
        // calibTree_->wfChi2.push_back(chi2);

        delete wFit;
        offset += opts.GetOpt<int>(channel+".nSamples");
    }
  
    // calibTree_->Fill();
    return true;
}

void CalibDigitizer::minimize()
{
    Log("Starting Mininimization for #"+to_string(nTotSamples_)+" cells");
    //  clock_t begin = clock();
    for (int i = 0;i<nTotSamples_;++i)
    {
        if (!fitDeltaV_)
            fitters_[i].FixParameter(0,0.);
        if (!fitDeltaT_)
            fitters_[i].FixParameter(1,0.);
        if (!fitSlopeV_)
            fitters_[i].FixParameter(2,0.);
        if (!fitQuadraticV_)
            fitters_[i].FixParameter(3,0.);

        fitters_[i].Eval();
    }
    //clock_t end = clock();
}

bool CalibDigitizer::EndLoop(int iLoop, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    minimize();
    int offset = 0;
    bool ret=true;
    for(auto& channel : channelsNames_)
    {
        int nSamples = opts.GetOpt<int>(channel+".nSamples");
        bgc_key_t channelID(opts.GetOpt<int>(channel+".digiBoard"), opts.GetOpt<int>(channel+".digiGroup"), opts.GetOpt<int>(channel+".digiChannel"));
        double deltaV_mean = 0;
        double deltaT_mean = 0;
        for (int iSample = 0;iSample<nSamples;++iSample)
        {
            TVectorD params;

            fitters_[offset+iSample].GetParameters(params);
            // TVectorD errors;
            // fitters_[iSample].GetErrors(errors);
            digiCalibration_[channelID][iSample].deltaV += params(0);
            digiCalibration_[channelID][iSample].deltaT += params(1);
            digiCalibration_[channelID][iSample].slopeV += params(2);
            digiCalibration_[channelID][iSample].quadraticV += params(3);
            deltaV_mean += params(0);
            deltaT_mean += params(1);
        }

        //---force the calibrations not to shift the average V and T
        deltaV_mean = deltaV_mean/(double)nSamples;
        deltaT_mean = deltaT_mean/(double)nSamples;
        for (int iSample = 0;iSample<nSamples;++iSample)
        {
            if(!fitSlopeV_ && !fitQuadraticV_ ) //do this when only fitting delta
                digiCalibration_[channelID][iSample].deltaV -= deltaV_mean;
            digiCalibration_[channelID][iSample].deltaT -= deltaT_mean;
        }

        offset += nSamples;

        Log("Relative Improvement "+to_string(fabs(thisFit_-prevFit_)/thisFit_));
        
        calibTree_->board = (int)get<0>(channelID);
        calibTree_->group = (int)get<1>(channelID);
        calibTree_->channel = (int)get<2>(channelID);
        calibTree_->calib->assign(digiCalibration_[channelID].begin(), digiCalibration_[channelID].end());
        calibTree_->Fill();

        if(fabs(thisFit_-prevFit_)/thisFit_ < 1e-3)
            ret = false;
    }

    return ret;
}


bool CalibDigitizer::End(map<string, PluginBase*>& plugins, CfgManager& opts)
{
    return true;
}

