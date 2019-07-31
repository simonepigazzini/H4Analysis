#include "ECALMatrixReco.h"

//----------Utils-------------------------------------------------------------------------
bool ECALMatrixReco::Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index)
{
    channels_ = opts.GetOpt<vector<string> >(instanceName_+".channels");
    for(auto& ch : channels_)
    {
        auto pos = opts.GetOpt<vector<float> >(ch+".matrixPos");
        crystalPosition_[ch] = make_pair(pos[0], pos[1]);
    }

    //---output tree
    outTree_ = ECALMatrixTree("ecal_tree", "ECAL matrix tree", index);
    RegisterSharedData(outTree_.GetTTreePtr(), "ecal_tree", true);
    // RegisterSharedData(new TTree("ecal_tree", "ECAL matrix tree"), "ecal_tree", true);
    // outTree_ = ECALMatrixTree((TTree*)data_.back().obj, index);

    //---EventAnalyzer to retrieve channel reco info
    eventAnalyzer_ = RecoEventAnalyzer(plugins, index);
    
    return true;

}

bool ECALMatrixReco::ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    outTree_.Reset();
    
    //---Get the energies
    map<string, float> energies;
    string seed_ch;
    for(auto& ch : channels_)
    {
        auto expr = "fit_ampl["+ch+"]";
        energies[ch] = eventAnalyzer_.GetValues(expr).back()[0];
        if(seed_ch == "" || energies[ch]>energies[seed_ch])
            seed_ch = ch;
    }
            
    //---Energy sums
    float e3x3=0, e5x5=0;
    auto& seed_pos = crystalPosition_[seed_ch]; 
    for(auto& ch : channels_)
    {
        auto& pos = crystalPosition_[ch];
        if(std::abs(pos.first-seed_pos.first)<2 && std::abs(pos.second-seed_pos.second)<2)
            e3x3 += energies[ch];

        e5x5 += energies[ch];
    }

    //---Hit position
    float sum_w=0, ecal_x=0, ecal_y=0;
    for(auto& ch : channels_)
    {        
        auto& pos = crystalPosition_[ch];        
        if(std::abs(pos.first-seed_pos.first)<2 && std::abs(pos.second-seed_pos.second)<2)
        {
            float w = std::max(0., 4.2 + log(energies[ch]/e3x3));
            ecal_x += (pos.first*22-11)*w;
            ecal_y += (pos.second*22-11)*w;
            sum_w += w;
        }
    }

    //---Fill output tree
    outTree_.seed = eventAnalyzer_.GetValues(seed_ch).back()[0];
    outTree_.e3x3 = e3x3;
    outTree_.e5x5 = e5x5;
    outTree_.ecal_x = ecal_x/sum_w;
    outTree_.ecal_y = ecal_y/sum_w; 
    outTree_.GetTTreePtr()->Fill();

    return true;
}

