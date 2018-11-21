#include "plugins/InfoTreeMaker.h"

//**********Utils*************************************************************************
//----------Begin-------------------------------------------------------------------------
bool InfoTreeMaker::Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index)
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

    //---mapped numerical variables
    vector<string> mappedVars = opts.GetOpt<vector<string> >(instanceName_+".mappedVarsNum");
    for(auto& mappedVar : mappedVars)
    {
        //---crate a branch for each remapped variable
        mappedVariablesNum_[mappedVar] = new float;
        info_tree_->Branch(mappedVar.c_str(), mappedVariablesNum_[mappedVar], (mappedVar+"/F").c_str());
        
        //---get remapped values
        //---seg fault protection: fill the map for as many values as the min between
        //   the orginal variable values and the remapped values.
        vector<float> values = opts.GetOpt<vector<float> >(instanceName_+"."+mappedVar);
        for(unsigned int i=0; i<min(values.size(), originValues.size()); ++i)
            remapNum_[mappedVar][originValues[i]] = values[i];
    }

    //---mapped string variables
    mappedVars = opts.GetOpt<vector<string> >(instanceName_+".mappedVarsStr");
    for(auto& mappedVar : mappedVars)
    {
        //---crate a branch for each remapped variable
        mappedVariablesStr_[mappedVar] = new string;
        info_tree_->Branch(mappedVar.c_str(), &mappedVariablesStr_[mappedVar]);
        
        //---get remapped values
        //---seg fault protection: fill the map for as many values as the min between
        //   the orginal variable values and the remapped values.
        vector<string> values = opts.GetOpt<vector<string> >(instanceName_+"."+mappedVar);
        for(unsigned int i=0; i<min(values.size(), originValues.size()); ++i)
            remapStr_[mappedVar][originValues[i]] = values[i];
    }
    
    return true;
}

//----------ProcessEvent------------------------------------------------------------------
bool InfoTreeMaker::ProcessEvent(H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    //---numeric variables
    for(auto& var : mappedVariablesNum_)
        *var.second = remapNum_[var.first][opts.GetOpt<string>(trackedVariable_)];

    //---string variables
    for(auto& var : mappedVariablesStr_)
        *var.second = remapStr_[var.first][opts.GetOpt<string>(trackedVariable_)];

    info_tree_->Fill();
    
    return true;
}
