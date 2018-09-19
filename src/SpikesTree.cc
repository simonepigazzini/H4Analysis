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
}

SpikesTree::~SpikesTree()
{}
