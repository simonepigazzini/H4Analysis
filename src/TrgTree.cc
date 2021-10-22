#include "interface/TrgTree.h"

TrgTree::TrgTree(uint64* idx, TTree* tree)
{
    tree_ = tree ? tree : new TTree();

    index=idx;
    trg_type_=0;
}

void TrgTree::Init(map<int, string>& triggers)
{
    //---trigger branches
    trgs_.reserve(triggers.size());
    for(auto const& [mask, name] : triggers)
    {
        trgs_.push_back(mask);
        tree_->Branch(name.c_str(), &(trgs_.back()), (name+"/I").c_str());
    }
    tree_->Branch("index", index, "index/l");
    tree_->Branch("trg", &trg_type_, "trg/i");

    return;
}
