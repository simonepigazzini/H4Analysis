#include "HodoBTFReco.h"

//**********Utils*************************************************************************
//----------Begin*************************************************************************
bool HodoBTFReco::Begin(CfgManager& opts, uint64* index)
{
    // PMT In - PNT Out  [both 0->63]
    ADC_to_PMT_map.insert(pair<int,int>(0,35));
    ADC_to_PMT_map.insert(pair<int,int>(1,41));
    ADC_to_PMT_map.insert(pair<int,int>(2,48));
    ADC_to_PMT_map.insert(pair<int,int>(3,56));
    ADC_to_PMT_map.insert(pair<int,int>(4,57));
    ADC_to_PMT_map.insert(pair<int,int>(5,49));
    ADC_to_PMT_map.insert(pair<int,int>(6,58));
    ADC_to_PMT_map.insert(pair<int,int>(7,50));
    ADC_to_PMT_map.insert(pair<int,int>(8,32));
    ADC_to_PMT_map.insert(pair<int,int>(9,59));
    ADC_to_PMT_map.insert(pair<int,int>(10,33));
    ADC_to_PMT_map.insert(pair<int,int>(11,51));
    ADC_to_PMT_map.insert(pair<int,int>(12,40));
    ADC_to_PMT_map.insert(pair<int,int>(13,60));
    ADC_to_PMT_map.insert(pair<int,int>(14,34));
    ADC_to_PMT_map.insert(pair<int,int>(15,52));
    ADC_to_PMT_map.insert(pair<int,int>(16,42));
    ADC_to_PMT_map.insert(pair<int,int>(17,61));
    ADC_to_PMT_map.insert(pair<int,int>(18,43));
    ADC_to_PMT_map.insert(pair<int,int>(19,53));
    ADC_to_PMT_map.insert(pair<int,int>(20,62));
    ADC_to_PMT_map.insert(pair<int,int>(21,44));
    ADC_to_PMT_map.insert(pair<int,int>(22,63));
    ADC_to_PMT_map.insert(pair<int,int>(23,54));
    ADC_to_PMT_map.insert(pair<int,int>(24,36));
    ADC_to_PMT_map.insert(pair<int,int>(25,55));
    ADC_to_PMT_map.insert(pair<int,int>(26,46));
    ADC_to_PMT_map.insert(pair<int,int>(27,45));
    ADC_to_PMT_map.insert(pair<int,int>(28,47));
    ADC_to_PMT_map.insert(pair<int,int>(29,38));
    ADC_to_PMT_map.insert(pair<int,int>(30,37));
    ADC_to_PMT_map.insert(pair<int,int>(31,39));
    ADC_to_PMT_map.insert(pair<int,int>(32,29));
    ADC_to_PMT_map.insert(pair<int,int>(33,31));
    ADC_to_PMT_map.insert(pair<int,int>(34,21));
    ADC_to_PMT_map.insert(pair<int,int>(35,30));
    ADC_to_PMT_map.insert(pair<int,int>(36,23));
    ADC_to_PMT_map.insert(pair<int,int>(37,28));
    ADC_to_PMT_map.insert(pair<int,int>(38,22));
    ADC_to_PMT_map.insert(pair<int,int>(39,15));
    ADC_to_PMT_map.insert(pair<int,int>(40,20));
    ADC_to_PMT_map.insert(pair<int,int>(41,14));
    ADC_to_PMT_map.insert(pair<int,int>(42,7));
    ADC_to_PMT_map.insert(pair<int,int>(43,27));
    ADC_to_PMT_map.insert(pair<int,int>(44,6));
    ADC_to_PMT_map.insert(pair<int,int>(45,13));
    ADC_to_PMT_map.insert(pair<int,int>(46,3));
    ADC_to_PMT_map.insert(pair<int,int>(47,5));
    ADC_to_PMT_map.insert(pair<int,int>(48,12));
    ADC_to_PMT_map.insert(pair<int,int>(49,4));
    ADC_to_PMT_map.insert(pair<int,int>(50,11));
    ADC_to_PMT_map.insert(pair<int,int>(51,19));
    ADC_to_PMT_map.insert(pair<int,int>(52,10));
    ADC_to_PMT_map.insert(pair<int,int>(53,18));
    ADC_to_PMT_map.insert(pair<int,int>(54,2));
    ADC_to_PMT_map.insert(pair<int,int>(55,17));
    ADC_to_PMT_map.insert(pair<int,int>(56,9));
    ADC_to_PMT_map.insert(pair<int,int>(57,16));
    ADC_to_PMT_map.insert(pair<int,int>(58,1));
    ADC_to_PMT_map.insert(pair<int,int>(59,0));
    ADC_to_PMT_map.insert(pair<int,int>(60,25));
    ADC_to_PMT_map.insert(pair<int,int>(61,8));
    ADC_to_PMT_map.insert(pair<int,int>(62,26));
    ADC_to_PMT_map.insert(pair<int,int>(63,24));

    // PMT Out - Hodo X  [both 0->63]
    PMT_to_hodoX_map.insert(pair<int,int>(0,28));
    PMT_to_hodoX_map.insert(pair<int,int>(1,17));
    PMT_to_hodoX_map.insert(pair<int,int>(2,16));
    PMT_to_hodoX_map.insert(pair<int,int>(3,15));
    PMT_to_hodoX_map.insert(pair<int,int>(4,14));
    PMT_to_hodoX_map.insert(pair<int,int>(5,13));
    PMT_to_hodoX_map.insert(pair<int,int>(6,12));
    PMT_to_hodoX_map.insert(pair<int,int>(7,0));
    PMT_to_hodoX_map.insert(pair<int,int>(8,29));
    PMT_to_hodoX_map.insert(pair<int,int>(9,18));
    PMT_to_hodoX_map.insert(pair<int,int>(10,19));
    PMT_to_hodoX_map.insert(pair<int,int>(11,20));
    PMT_to_hodoX_map.insert(pair<int,int>(12,9));
    PMT_to_hodoX_map.insert(pair<int,int>(13,10));
    PMT_to_hodoX_map.insert(pair<int,int>(14,11));
    PMT_to_hodoX_map.insert(pair<int,int>(15,1));
    PMT_to_hodoX_map.insert(pair<int,int>(16,30));
    PMT_to_hodoX_map.insert(pair<int,int>(17,23));
    PMT_to_hodoX_map.insert(pair<int,int>(18,22));
    PMT_to_hodoX_map.insert(pair<int,int>(19,21));
    PMT_to_hodoX_map.insert(pair<int,int>(20,8));
    PMT_to_hodoX_map.insert(pair<int,int>(21,7));
    PMT_to_hodoX_map.insert(pair<int,int>(22,6));
    PMT_to_hodoX_map.insert(pair<int,int>(23,2));
    PMT_to_hodoX_map.insert(pair<int,int>(24,31));
    PMT_to_hodoX_map.insert(pair<int,int>(25,24));
    PMT_to_hodoX_map.insert(pair<int,int>(26,25));
    PMT_to_hodoX_map.insert(pair<int,int>(27,26));
    PMT_to_hodoX_map.insert(pair<int,int>(28,27));
    PMT_to_hodoX_map.insert(pair<int,int>(29,5));
    PMT_to_hodoX_map.insert(pair<int,int>(30,4));
    PMT_to_hodoX_map.insert(pair<int,int>(31,3));

    // PMT Out - Hodo Y  [both 0->63]
    PMT_to_hodoY_map.insert(pair<int,int>(32,3));
    PMT_to_hodoY_map.insert(pair<int,int>(33,4));
    PMT_to_hodoY_map.insert(pair<int,int>(34,5));
    PMT_to_hodoY_map.insert(pair<int,int>(35,27));
    PMT_to_hodoY_map.insert(pair<int,int>(36,26));
    PMT_to_hodoY_map.insert(pair<int,int>(37,25));
    PMT_to_hodoY_map.insert(pair<int,int>(38,24));
    PMT_to_hodoY_map.insert(pair<int,int>(39,31));
    PMT_to_hodoY_map.insert(pair<int,int>(40,2));
    PMT_to_hodoY_map.insert(pair<int,int>(41,6));
    PMT_to_hodoY_map.insert(pair<int,int>(42,7));
    PMT_to_hodoY_map.insert(pair<int,int>(43,8));
    PMT_to_hodoY_map.insert(pair<int,int>(44,21));
    PMT_to_hodoY_map.insert(pair<int,int>(45,22));
    PMT_to_hodoY_map.insert(pair<int,int>(46,23));
    PMT_to_hodoY_map.insert(pair<int,int>(47,30));
    PMT_to_hodoY_map.insert(pair<int,int>(48,1));
    PMT_to_hodoY_map.insert(pair<int,int>(49,11));
    PMT_to_hodoY_map.insert(pair<int,int>(50,10));
    PMT_to_hodoY_map.insert(pair<int,int>(51,9));
    PMT_to_hodoY_map.insert(pair<int,int>(52,20));
    PMT_to_hodoY_map.insert(pair<int,int>(53,19));
    PMT_to_hodoY_map.insert(pair<int,int>(54,18));
    PMT_to_hodoY_map.insert(pair<int,int>(55,29));
    PMT_to_hodoY_map.insert(pair<int,int>(56,0));
    PMT_to_hodoY_map.insert(pair<int,int>(57,12));
    PMT_to_hodoY_map.insert(pair<int,int>(58,13));
    PMT_to_hodoY_map.insert(pair<int,int>(59,14));
    PMT_to_hodoY_map.insert(pair<int,int>(60,15));
    PMT_to_hodoY_map.insert(pair<int,int>(61,16));
    PMT_to_hodoY_map.insert(pair<int,int>(62,17));
    PMT_to_hodoY_map.insert(pair<int,int>(63,28));

    //---create a position tree
    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;
    RegisterSharedData(new TTree("hX", "hodoX_tree"), "hodoX_tree", storeTree);
    hodoTreeX_ = PositionTree(index, (TTree*)data_.back().obj);
    hodoTreeX_.Init();
    RegisterSharedData(new TTree("hY", "hodoY_tree"), "hodoY_tree", storeTree);
    hodoTreeY_ = PositionTree(index, (TTree*)data_.back().obj);
    hodoTreeY_.Init();

    return true;
}

bool HodoBTFReco::ProcessEvent(H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    hodoTreeX_.Clear();
    hodoTreeY_.Clear();    
    hodoXpos.clear();
    hodoYpos.clear();
    for(unsigned int iCh=0; iCh<h4Tree.nAdcChannels; iCh++)
    {
        if(h4Tree.adcBoard[iCh] == 201392129)
        {
            int tmpPAD = ADC_to_PMT_map[h4Tree.adcChannel[iCh]];
            if(tmpPAD >= 0 && tmpPAD < 32)
            {
                if(int(h4Tree.adcData[iCh]) > opts.GetOpt<int>(instanceName_+".threshold"))
                    hodoXpos.push_back(PMT_to_hodoX_map[tmpPAD]);
            }
            else if(tmpPAD >= 32 && tmpPAD < 64)
            {
                if(int(h4Tree.adcData[iCh]) > opts.GetOpt<int>(instanceName_+".threshold"))
                    hodoYpos.push_back(PMT_to_hodoY_map[tmpPAD]);
            }
        }	
    }
    //--fill output trees with default values if needed
    hodoTreeX_.n_clusters = hodoXpos.size();
    hodoTreeX_.clusters.emplace_back(hodoXpos.size(), hodoXpos.size() > 0 ? hodoXpos[0] : -999, -999);
    hodoTreeY_.n_clusters = hodoYpos.size();
    hodoTreeY_.clusters.emplace_back(hodoYpos.size(), -999, hodoYpos.size() > 0 ? hodoYpos[0] : -999);
    hodoTreeX_.Fill();
    hodoTreeY_.Fill();

    //---clear trees
    hodoTreeX_.Clear();
    hodoTreeY_.Clear();
    
    return true;
}
