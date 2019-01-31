#include "interface/RecoEventAnalyzer.h"

//**********Constructors******************************************************************
RecoEventAnalyzer::RecoEventAnalyzer(map<string, PluginBase*>& plugins, uint64* index)    
{
    index_ = index;
    currentEntry_ = 0;
    
    sTree_ = new TTree();
    sTree_->Branch("index", index, "index/l");
    
    //---look for all the trees registered by the plugins
    for(auto& plugin : plugins)
    {
       //---get permanent data from each plugin and store them in the out file
        for(auto& shared : plugin.second->GetSharedData())
        {
            if(shared.obj->IsA()->GetName() == string("TTree"))
            {
                auto tree = (TTree*)shared.obj;
                sTree_->AddFriend(tree, tree->GetName());
                indexCheckStr_ += std::string(tree->GetName())+".index == index && ";
            }
        }
    }
}

//**********Utils*************************************************************************
//----------Load Trees--------------------------------------------------------------------
//---Load all trees that has last entry equal to the current entry
void RecoEventAnalyzer::UpdateEntry()
{
    //---Load trees only once pre event
    if(*index_ > currentEntry_)
    {        
        sTree_->Fill();
        currentEntry_ = *index_;
    }
}

//----------EvaluateSelection--------------------------------------------------------------
//---Evaluate expression on the last event, return passed/failed
bool RecoEventAnalyzer::EvaluateSelection(std::string& expr)
{
    UpdateEntry();    
    
    return sTree_->Draw("1", (indexCheckStr_+expr).c_str(), "goff", 1, sTree_->GetEntriesFast()-1);
}

//----------GetValues----------------------------------------------------------------------
//---Return the values of the expression specified in expr (column separated as in TTree::Draw)
const std::vector<std::vector<double> >& RecoEventAnalyzer::GetValues(std::string& expr, std::string cuts)
{
    UpdateEntry();

    //---prepare return value
    auto n_vals = std::count(expr.begin(), expr.end(), ':')+1;
    lastValues_.clear();
    lastValues_.resize(n_vals);

    //---get values
    auto size = sTree_->Draw(expr.c_str(), (indexCheckStr_+cuts).c_str(), "goff", 1, sTree_->GetEntriesFast()-1);
    for(unsigned int i=0; i<n_vals; ++i)
        lastValues_[i].assign(sTree_->GetVal(i), sTree_->GetVal(i)+size);

    return lastValues_;
}
    
