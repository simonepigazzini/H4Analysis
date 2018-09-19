#ifndef __WFVIEWER_H__
#define __WFVIEWER_H__

#include <math.h>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"
#include "TObject.h"
#include "Math/Interpolator.h"
#include "Math/Factory.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"

#include "WFClass.h"

using namespace std;

class InterpolatorFunc: public ROOT::Math::Interpolator
{
public:
    InterpolatorFunc(unsigned int ndata, ROOT::Math::Interpolation::Type type, float ampl, float time):
        ROOT::Math::Interpolator(ndata, type), ampl_(ampl), time_(time) {};
    
    double operator() (double* x, double*) const
        {
            double value=ampl_*ROOT::Math::Interpolator::Eval(x[0]-time_);
            return std::isnan(value) ? 0 : value;
        };

    float ampl_;
    float time_;
    
    //---ROOT dictionary
    ClassDef(InterpolatorFunc, 0)  
};

class WFViewer: public TObject
{
public:
    //---ctors---
    WFViewer();
    WFViewer(const char* tree_name);
    WFViewer(TTree* tree);
    WFViewer(string name, TH1F* h_template);
    //---dtor---
    ~WFViewer();

    //---setters---
    void         SetTemplate(TH1F* h_template);
    void         SetTree(const char* digi_tree="digi_tree", const char* wf_tree="wf");
    void         SetTree(TTree* tree) {tree_ = tree;};
    void         SetChannelName(string name) {name_ = name;};

    //---getters---
    string       GetChannelName() {return name_;};
    
    //---utils---
    virtual void Draw(string ref, const char* cut="", Long64_t max_entries=std::numeric_limits<Long64_t>::max());
    virtual void Draw(unsigned int iEntry, const char* wf_tree="wf");

private:
    string            name_;
    int               channel_;
    TH1F              h_pull_;
    TH1F              h_template_;
    TF1*              f_template_;
    TTree*            tree_;
    
    //---ROOT dictionary
    ClassDef(WFViewer, 1)
};

#endif
