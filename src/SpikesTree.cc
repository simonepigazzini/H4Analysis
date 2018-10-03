#include "interface/SpikesTree.h"

SpikesTree::SpikesTree(uint64* idx, TTree* tree, string prefix)
{
    prefix_=prefix;
    tree_ = tree ? tree : new TTree();

    index=idx;
}

void SpikesTree::Init(vector<string>& names)
{
    //---allocate enough space for all channels
    n_channels = names.size();
    max_hit = -1;
    channels = new int[n_channels];
    undershoot = new float[n_channels];
    amp_sum_matrix = new float[n_channels];
    n_swiss_cross_neighbours = new unsigned int[n_channels];
    swiss_cross = new float[n_channels];
    n_channels_3by3 = new unsigned int[n_channels];
    amp_sum_3by3 = new float[n_channels];
    n_samples_above_25perc_max = new unsigned int[n_channels];
    n_samples_above_50perc_max = new unsigned int[n_channels];
    n_samples_above_75perc_max = new unsigned int[n_channels];
    tot_25perc_max = new float[n_channels];
    tot_50perc_max = new float[n_channels];
    tot_75perc_max = new float[n_channels];
    sample_max_minus1_over_sample_max = new float[n_channels];
    sample_max_minus2_over_sample_max = new float[n_channels];
    sample_max_minus3_over_sample_max = new float[n_channels];
    sample_max_plus1_over_sample_max = new float[n_channels];
    sample_max_plus2_over_sample_max = new float[n_channels];
    sample_max_plus3_over_sample_max = new float[n_channels];
    t_undershoot_minus_t_sample_max = new float[n_channels];
    t_3sigma_noise_minus_t_sample_max = new float[n_channels];
    
    //---channels branches
    for(unsigned int iCh=0; iCh<n_channels; iCh++)
    {
        channels[iCh]=iCh;
        tree_->Branch(names[iCh].c_str(), &(channels[iCh]), (names[iCh]+"/I").c_str());
    }

    //---global branches    
    string size_var = "n_"+prefix_+"channels";
    string size_time_var = "n_"+prefix_+"timetypes";
    tree_->Branch("index", index, "index/l");
    tree_->Branch("max_hit", &max_hit, "max_hit/I");
    tree_->Branch(size_var.c_str(), &n_channels, (size_var+"/i").c_str());
    tree_->Branch((prefix_+"undershoot").c_str(), undershoot, (prefix_+"undershoot["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"amp_sum_matrix").c_str(), amp_sum_matrix, (prefix_+"amp_sum_matrix["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"n_swiss_cross_neighbours").c_str(), n_swiss_cross_neighbours, (prefix_+"n_swiss_cross_neighbours["+size_var+"]/i").c_str());
    tree_->Branch((prefix_+"swiss_cross").c_str(), swiss_cross, (prefix_+"swiss_cross["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"n_channels_3by3").c_str(), n_channels_3by3, (prefix_+"n_channels_3by3["+size_var+"]/i").c_str());
    tree_->Branch((prefix_+"amp_sum_3by3").c_str(), amp_sum_3by3, (prefix_+"amp_sum_3by3["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"n_samples_above_25perc_max").c_str(), n_samples_above_25perc_max, (prefix_+"n_samples_above_25perc_max["+size_var+"]/i").c_str());
    tree_->Branch((prefix_+"n_samples_above_50perc_max").c_str(), n_samples_above_50perc_max, (prefix_+"n_samples_above_50perc_max["+size_var+"]/i").c_str());
    tree_->Branch((prefix_+"n_samples_above_75perc_max").c_str(), n_samples_above_75perc_max, (prefix_+"n_samples_above_75perc_max["+size_var+"]/i").c_str());
    tree_->Branch((prefix_+"tot_25perc_max").c_str(), tot_25perc_max, (prefix_+"tot_25perc_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"tot_50perc_max").c_str(), tot_50perc_max, (prefix_+"tot_50perc_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"tot_75perc_max").c_str(), tot_75perc_max, (prefix_+"tot_75perc_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"sample_max_minus1_over_sample_max").c_str(), sample_max_minus1_over_sample_max, (prefix_+"sample_max_minus1_over_sample_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"sample_max_minus2_over_sample_max").c_str(), sample_max_minus2_over_sample_max, (prefix_+"sample_max_minus2_over_sample_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"sample_max_minus3_over_sample_max").c_str(), sample_max_minus3_over_sample_max, (prefix_+"sample_max_minus3_over_sample_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"sample_max_plus1_over_sample_max").c_str(), sample_max_plus1_over_sample_max, (prefix_+"sample_max_plus1_over_sample_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"sample_max_plus2_over_sample_max").c_str(), sample_max_plus2_over_sample_max, (prefix_+"sample_max_plus2_over_sample_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"sample_max_plus3_over_sample_max").c_str(), sample_max_plus3_over_sample_max, (prefix_+"sample_max_plus3_over_sample_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"t_undershoot_minus_t_sample_max").c_str(), t_undershoot_minus_t_sample_max, (prefix_+"t_undershoot_minus_t_sample_max["+size_var+"]/F").c_str());
    tree_->Branch((prefix_+"t_3sigma_noise_minus_t_sample_max").c_str(), t_3sigma_noise_minus_t_sample_max, (prefix_+"t_3sigma_noise_minus_t_sample_max["+size_var+"]/F").c_str());
}

SpikesTree::~SpikesTree()
{}
