#include "plugins/InfoTreeMaker.h"

//**********Utils*************************************************************************
//----------Begin-------------------------------------------------------------------------
bool InfoTreeMaker::Begin(CfgManager& opts, uint64* index)
{
    index_ = index;
    //---Register the output file
    RegisterSharedData(new TTree(opts.GetOpt<string>(instanceName_+".treeName").c_str(), "info_tree"),
                       "info", true);
    info_tree_ = (TTree*)data_.back().obj;
    info_tree_->Branch("index", index_, "index/l");

    //---get original variable and list of variables to which the first one is remapped to
    trackedVariable_ = opts.GetOpt<string>(instanceName_+".trackedVariable");
    vector<string> originValues = opts.GetOpt<vector<string> >(instanceName_+".originValues");
    vector<string> mappedVars = opts.GetOpt<vector<string> >(instanceName_+".mappedVars");
    for(auto& mappedVar : mappedVars)
    {
        //---crate a branch for each remapped variable
        mappedVariables_[mappedVar] = new float;
        info_tree_->Branch(mappedVar.c_str(), mappedVariables_[mappedVar], (mappedVar+"/F").c_str());
        
        //---get remapped values
        //---seg fault protection: fill the map for as many values as the min between
        //   the orginal variable values and the remapped values.
        vector<float> values = opts.GetOpt<vector<float> >(instanceName_+"."+mappedVar);
        for(int i=0; i<min(values.size(), originValues.size()); ++i)
            remap_[mappedVar][originValues[i]] = values[i];
    }
    
    return true;
}

//----------ProcessEvent------------------------------------------------------------------
bool InfoTreeMaker::ProcessEvent(const H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    for(auto& var : mappedVariables_)
        *var.second = remap_[var.first][opts.GetOpt<string>(trackedVariable_)];

    info_tree_->Fill();
    
    return true;
}
