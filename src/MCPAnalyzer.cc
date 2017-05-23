#include "interface/MCPAnalyzer.h"

//**********constructors******************************************************************
MCPAnalyzer::MCPAnalyzer():
    channel_(""), tree_((TTree*)gDirectory->GetFile()->Get("h4"))
{}

MCPAnalyzer::MCPAnalyzer(string channel, string tree):
    channel_(channel), tree_((TTree*)gDirectory->GetFile()->Get(tree.c_str()))
{}
    
TGraphAsymmErrors* MCPAnalyzer::Efficiency(string scan_var, string eff_var, string cut)
{
    //---consistency check
    if(!tree_ || channel_=="")
    {
        cout << "ERROR: MCPAnalyzer has not been set up properly!" << endl;
        return NULL;
    }

    //---delete old histo
    if(gDirectory->GetFile()->Get("h2_eff"))
        gDirectory->GetFile()->Get("h2_eff")->Delete();

    //---creat 2D histo: eff[0,1] vs scan_var
    //---if X binning is not specified creat a TH2F with automatic binning
    TH2F* h2_eff = NULL;
    if(scan_var.find_last_of(">") != string::npos)
    {
        string binning = scan_var.substr(scan_var.find_last_of("(")+1, scan_var.find_last_of(")"));
        string bins = binning.substr(0, binning.find_first_of(","));
        string min = binning.substr(binning.find_first_of(",")+1, binning.find_last_of(","));
        string max = binning.substr(binning.find_last_of(",")+1);
        h2_eff = new TH2F("h2_eff", "", stoi(binning), stof(min), stof(max), 2, 0, 1.1);
    }
    else
        h2_eff = new TH2F("h2_eff", "", 1, 1, 1, 2, 0, 1.1);

    string draw_var = eff_var+":"+scan_var.substr(0, scan_var.find_first_of(">"))+">>h2_eff";
    cout << "Drawing: " << draw_var << endl;
    tree_->Draw(draw_var.c_str(), cut.c_str(), "goff");

    //---fed the TEfficiency for efficiency computation
    TH1D* h_total = h2_eff->ProjectionX("_total", 1, 2);
    TH1D* h_pass = h2_eff->ProjectionX("_pass", 2, 2);
    TEfficiency* eff = NULL;
    if(TEfficiency::CheckConsistency(*h_pass, *h_total))
        eff = new TEfficiency(*h_pass, *h_total);

    eff->Draw();

    return eff->CreateGraph();
}
    
        
