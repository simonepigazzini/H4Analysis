#include "interface/DigiTree.h"

DigiTree::DigiTree(uint64* idx, TTree* tree, string prefix)
{
    prefix_=prefix;
    tree_ = tree ? tree : new TTree();

    index=idx;
}

void DigiTree::Init(vector<string>& names, vector<string>& timetypes)
{
    //---allocate enough space for all channels
    n_channels = names.size();
    n_times = timetypes.size()*n_channels;
    channels = new int[n_channels];
    time_types = new int[n_times];    
    pedestal = new float[n_channels];
    b_charge = new float[n_channels];
    b_slope = new float[n_channels];
    b_rms = new float[n_channels];
    time = new float[n_channels*n_times];
    time_chi2 = new float[n_channels*n_times];
    maximum = new float[n_channels];
    time_maximum = new float[n_channels];
    amp_max = new float[n_channels];
    time_max = new float[n_channels];
    chi2_max = new float[n_channels];
    charge_tot = new float[n_channels];
    charge_sig = new float[n_channels];
    fit_ampl = new float[n_channels];
    fit_time = new float[n_channels];
    fit_chi2 = new float[n_channels];
    calibration = new float[n_channels];

    //---channels branches
    for(unsigned int iCh=0; iCh<n_channels; iCh++)
    {
        channels[iCh]=iCh;
        tree_->Branch(names[iCh].c_str(), &(channels[iCh]), (names[iCh]+"/I").c_str());
    }
    //---time types branches
    for(unsigned int iT=0; iT<timetypes.size(); iT++)
    {
        time_types[iT]=iT*n_channels;
        tree_->Branch(timetypes[iT].c_str(), &(time_types[iT]), (timetypes[iT]+"/I").c_str());
    }
    //---global branches    
    string size_var = "n_"+prefix_+"channels";
    string size_time_var = "n_"+prefix_+"timetypes";
    tree_->Branch("index", index, "index/l");
    tree_->Branch(size_var.c_str(), &n_channels, (size_var+"/i").c_str());
    tree_->Branch(size_time_var.c_str(), &n_times, (size_time_var+"/i").c_str());
    tree_->Branch((prefix_+"pedestal").c_str(), pedestal, (prefix_+"pedestal["+size_var+"]/F").c_str());    
    tree_->Branch((prefix_+"b_charge").c_str(), b_charge, (prefix_+"b_charge["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"b_slope").c_str(), b_slope, (prefix_+"b_slope["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"b_rms").c_str(), b_rms, (prefix_+"b_rms["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"time").c_str(), time, (prefix_+"time["+size_time_var+"]/F").c_str());
    tree_->Branch((prefix_+"time_chi2").c_str(), time_chi2, (prefix_+"time_chi2["+size_time_var+"]/F").c_str());
    tree_->Branch((prefix_+"maximum").c_str(), maximum, (prefix_+"maximum["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"time_maximum").c_str(), time_maximum, (prefix_+"time_maximum["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"amp_max").c_str(), amp_max, (prefix_+"amp_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"time_max").c_str(), time_max, (prefix_+"time_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"chi2_max").c_str(), chi2_max, (prefix_+"chi2_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"charge_tot").c_str(), charge_tot, (prefix_+"charge_tot["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"charge_sig").c_str(), charge_sig, (prefix_+"charge_sig["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"fit_ampl").c_str(), fit_ampl, (prefix_+"fit_ampl["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"fit_time").c_str(), fit_time, (prefix_+"fit_time["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"fit_chi2").c_str(), fit_chi2, (prefix_+"fit_chi2["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"calibration").c_str(), calibration, (prefix_+"calibration["+size_var+"]/F").c_str());
}

DigiTree::~DigiTree()
{}
