#ifndef __MCPANALYZER_H__
#define __MCPANALYZER_H__

#include <math.h>
#include <iostream>

#include "TROOT.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"

using namespace std;

class MCPAnalyzer: public TObject
{
public:
    //---ctors---
    MCPAnalyzer();
    MCPAnalyzer(string channel, string tree="h4");
    //---dtor---
    ~MCPAnalyzer(){};

    //---getters---
    TTree* GetTree(){return tree_;};
    string GetChannelName(){return channel_;};
    
    //---utils---
    TGraphAsymmErrors* Efficiency(string scan_var, string eff_var="charge_sig", string cut="");
    
private:
    string channel_;
    TTree* tree_;
    
    //---ROOT dictionary
    ClassDef(MCPAnalyzer, 1)
        };

#endif
