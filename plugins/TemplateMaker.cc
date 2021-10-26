#include "plugins/TemplateMaker.h"

//----------Utils-------------------------------------------------------------------------
bool TemplateMaker::Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index)
{
    gErrorIgnoreLevel = kError;
    
    //---inputs---
    if(!opts.OptExist(instanceName_+".srcInstanceName"))
    {
        Log("no source plugin specified", ERR);
        return false;
    }
    srcInstance_ = opts.GetOpt<string>(instanceName_+".srcInstanceName");
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");
    
    storeRatios_ = opts.OptExist(instanceName_+".storeRatios") ?
        opts.GetOpt<bool>(instanceName_+".storeRatios") : false;

    eventAnalyzer_ = RecoEventAnalyzer(plugins, index);
    
    //---outputs
    for(auto& channel : channelsNames_)
    {
        //---load WFs from source instance shared data
        auto shared_data = plugins[srcInstance_]->GetSharedData(srcInstance_+"_"+channel, "", false);
        if(shared_data.size() != 0)
            WFs_[channel] = (WFClass*)shared_data.at(0).obj;
        else
            Log("channels samples not found check DigiReco step", WARN);
        
        //---shape templates
        auto bins = opts.GetOpt<vector<double> >(channel+".templateMaker.timeBins");
        vector<double> range = {-0.1, 1.1};
        if(opts.OptExist(channel+".templateMaker.amplitudeRange"))
            range = opts.GetOpt<vector<double> >(channel+".templateMaker.amplitudeRange");
        // fixed width time bins
        if(bins.size() == 3)
            hTemplates_[channel] = TProfile(("tmpl_"+channel).c_str(), "", bins[0], bins[1], bins[2], range[0], range[1]);
        // variable width time bins
        else if(bins.size() > 3)
            hTemplates_[channel] = TProfile(("tmpl_"+channel).c_str(), "", bins.size()-1, bins.data(), range[0], range[1]);
        // wrong binning
        else
        {
            Log("incorrect number of binning parameters, 3 or more needed only "+to_string(bins.size())+" specified.", ERR);
            return false;
        }

        RegisterSharedData(&hTemplates_[channel], "tmpl_"+channel, true);
        
        //---shape template ratios
        if(storeRatios_)
        {
            hRatios_[channel] = *((TH1D*)hTemplates_[channel].Clone());
            RegisterSharedData(&hRatios_[channel], "tmpl_"+channel, true);
        }
    }

    return true;
}

bool TemplateMaker::BeginLoop(int iLoop, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    if(iLoop > 0)
    {
        for(auto& channel : channelsNames_)
        {
            if(storeRatios_)
                hRatios_[channel] = *((TH1D*)hTemplates_[channel].Clone());
            hRatios_[channel].Reset();
        }
    }
    
    return true;
}

bool TemplateMaker::ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    for(auto& channel : channelsNames_)
    {
        //---fill template after checking the event selections
        auto selection = opts.OptExist(channel+".templateMaker.eventSelection") ?
            opts.GetOpt<string>(channel+".templateMaker.eventSelection") : string("");
        if(eventAnalyzer_.EvaluateSelection(selection))
        {
            auto vars = opts.GetOpt<string>(channel+".templateMaker.expression");
            auto time_and_ampl = eventAnalyzer_.GetValues(vars, selection);
            for(unsigned int i=0; i<time_and_ampl.back().size(); ++i)
                hTemplates_[channel].Fill(time_and_ampl[1][i], time_and_ampl[0][i]);
        }            
    }

    return true;
}

bool TemplateMaker::EndLoop(int iLoop, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    for(auto& channel : channelsNames_)
    {
        //---Update template in WFClass
        WFs_[channel]->SetTemplate(&hTemplates_[channel]);
        //---compute ratio with template from previous step
        if(iLoop > 0)
            hRatios_[channel].Divide(&hTemplates_[channel]);
    }
    
    return true;
}
    
