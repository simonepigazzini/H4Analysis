#ifndef __SPIKES_TREE__
#define __SPIKES_TREE__

#include <memory>
#include <iostream>
#include <string>
#include <vector>

#include "TTree.h"
#include "TString.h"

using namespace std;

typedef unsigned long int uint32;
typedef unsigned long long int uint64;
 
//****************************************************************************************

class SpikesTree
{
public: 
    //---ctors---
    SpikesTree(){};
    SpikesTree(uint64* idx, TTree* tree=NULL, string prefix="");
    //---dtor---
    ~SpikesTree();

    //---utils---
    void Init(vector<string>& names);
    void Fill() {tree_->Fill();};
    
    TTree*        tree_; 
    string        prefix_;

    uint64*       index;
    unsigned int  n_channels;
    int           max_hit;
    int*          channels;
    float*        undershoot;
    float*        amp_sum_matrix;
    unsigned int* n_swiss_cross_neighbours;
    float*        swiss_cross;
    unsigned int* n_channels_3by3;
    float*        amp_sum_3by3;
    unsigned int* n_samples_above_25perc_max;
    unsigned int* n_samples_above_50perc_max;
    unsigned int* n_samples_above_75perc_max;
    float*        tot_25perc_max;
    float*        tot_50perc_max;
    float*        tot_75perc_max;
    float*        sample_max_minus1_over_sample_max;
    float*        sample_max_minus2_over_sample_max;
    float*        sample_max_minus3_over_sample_max;
    float*        sample_max_plus1_over_sample_max;
    float*        sample_max_plus2_over_sample_max;
    float*        sample_max_plus3_over_sample_max;
    float*        t_undershoot_minus_t_sample_max;
    float*        t_3sigma_noise_minus_t_sample_max;
    float*        ld;
};

#endif
