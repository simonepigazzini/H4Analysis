#include "RecoTree.h"

RecoTree::RecoTree(int nCh, int* idx, string prefix)
{
    prefix_=prefix;
    tree_ = new TTree("reco_tree", "reco_tree");
    tree_->SetMaxVirtualSize(10000000);

    index=idx;
    time_stamp=0;
    //---allocate enough space for all channels
    n_channels = nCh;
    b_charge = new float[nCh];
    b_slope = new float[nCh];
    b_rms = new float[nCh];
    time = new float[nCh];
    time_chi2 = new float[nCh];
    maximum = new float[nCh];
    amp_max = new float[nCh];
    charge_tot = new float[nCh];
    charge_sig = new float[nCh];
    fit_ampl = new float[nCh];
    fit_time = new float[nCh];
    fit_chi2 = new float[nCh];

    //---global branches
    string size_var = "n_"+prefix_+"channels";
    tree_->Branch("index", index, "index/I");
    tree_->Branch("time_stamp", &time_stamp, "time_stamp/l");
    tree_->Branch(size_var.c_str(), &n_channels, (size_var+"/i").c_str());
    tree_->Branch((prefix_+"b_charge").c_str(), b_charge, ("b_charge["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"b_slope").c_str(), b_slope, ("b_slope["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"b_rms").c_str(), b_rms, ("b_rms["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"time").c_str(), time, ("time["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"time_chi2").c_str(), time_chi2, ("time_chi2["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"maximum").c_str(), maximum, ("maximum["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"amp_max").c_str(), amp_max, ("amp_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"charge_tot").c_str(), charge_tot, ("charge_tot["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"charge_sig").c_str(), charge_sig, ("charge_sig["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"fit_ampl").c_str(), fit_ampl, ("fit_ampl["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"fit_time").c_str(), fit_time, ("fit_time["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"fit_chi2").c_str(), fit_chi2, ("fit_chi2["+size_var+"]/F").c_str());
}

RecoTree::~RecoTree()
{
    delete[] b_charge;
    delete[] b_slope;
    delete[] b_rms;
    delete[] time;
    delete[] time_chi2;
    delete[] maximum;
    delete[] amp_max;
    delete[] charge_tot;
    delete[] charge_sig;
    delete[] fit_ampl;
    delete[] fit_time;
    delete[] fit_chi2;

    delete tree_;
}
