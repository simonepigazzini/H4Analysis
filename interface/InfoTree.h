#ifndef __INFO_TREE__
#define __INFO_TREE__

#include <string>
#include <vector>

#include "TTree.h"
#include "TString.h"

using namespace std;

typedef unsigned long int uint32;
typedef unsigned long long int uint64;
 
//****************************************************************************************

class InfoTree
{
public: 
    //---ctors---
    InfoTree(){};
    InfoTree(vector<string> names, int nSamples, bool H4hodo, int nWireCh, int* idx);
    //---dtor---
    ~InfoTree();

    //---utils---
    void Fill() {tree_->Fill();};
    void Write(const char* name="info_tree", const char* title="info_tree")
        {tree_->BuildIndex("index"); tree_->SetTitle(title); tree_->Write(name);};

    TTree* tree_; 

    int*         index;
    uint64       start_time;
    uint64       time_stamp;
    unsigned int run;
    unsigned int spill;
    unsigned int event;   
    int*         channels;
    int          n_wchamber;
    float*       wireX;
    float*       wireY;
    float        hodoX1;
    float        hodoY1;
    float        hodoX2;
    float        hodoY2;
};

InfoTree::InfoTree(vector<string> names, int nSamples, bool H4hodo, int nWireCh, int* idx)
{
    tree_ = new TTree("info_tree", "info_tree");
    tree_->SetMaxVirtualSize(10000000);
    unsigned int nCh=names.size();

    index=idx;
    start_time=0;
    time_stamp=0;
    run=0;
    spill=0;
    event=0;   
    //---allocate enough space for all channels
    channels = new int[nCh];
    //---global branches
    tree_->Branch("index", index, "index/I");
    tree_->Branch("start_time", &start_time, "start_time/l");
    tree_->Branch("time_stamp", &time_stamp, "time_stamp/l");
    tree_->Branch("run", &run, "run/i");
    tree_->Branch("spill", &spill, "spill/i");
    tree_->Branch("event", &event, "event/i");
    //---channels branches
    for(int iCh=0; iCh<nCh; iCh++)
    {
        channels[iCh]=iCh;
        tree_->Branch(names[iCh].c_str(), &(channels[iCh]), (names[iCh]+"/I").c_str());
    }
    //---wire chamber branches
    if(nWireCh>0)
    {
        n_wchamber = nWireCh;
        wireX = new float[nWireCh];
        wireY = new float[nWireCh];
        tree_->Branch("n_wchamber", &n_wchamber, "n_wchamber/I");
        tree_->Branch("wireX", wireX, "wireX[n_wchamber]/F");
        tree_->Branch("wireY", wireY, "wireY[n_wchamber]/F");
    }
    //---hodoscope branches
    hodoX1=0;
    hodoY1=0; 
    hodoX2=0;
    hodoY2=0;
    if(H4hodo)
    {
        tree_->Branch("hodoX1", &hodoX1, "hodoX1/F");
        tree_->Branch("hodoY1", &hodoY1, "hodoY1/F");
        tree_->Branch("hodoX2", &hodoX2, "hodoX2/F");
        tree_->Branch("hodoY2", &hodoY2, "hodoY2/F");
    }
}

InfoTree::~InfoTree()
{
    delete tree_;
}

#endif
