#include "interface/WFViewer.h"
#include "interface/DigiTree.h"

#include "TStyle.h"
#include "TLine.h"
#include "TCanvas.h"

//**********constructors******************************************************************
WFViewer::WFViewer():
    name_(), channel_(0), h_pull_("h_pull_", "", 1024, 0, 204.8), h_template_(), f_template_()
{
    TFile* currentFile = gROOT->GetFile();    
    if(currentFile)
        tree_ = (TTree*)currentFile->Get("h4");
}

WFViewer::WFViewer(const char* tree_name):
    name_(), channel_(0), h_pull_("h_pull_", "", 1024, 0, 204.8), h_template_(), f_template_()
{
    TFile* currentFile = gROOT->GetFile();
    if(!tree_name)
        tree_name = "h4";
    if(currentFile)
        tree_ = (TTree*)currentFile->Get(tree_name);
    else
        cout << "ERROR WFViewer(): TTree '" << tree_name << "' not found" << endl;
}

WFViewer::WFViewer(TTree* tree):
    name_(), channel_(0), h_pull_("h_pull_", "", 1024, 0, 204.8), h_template_(), f_template_()
{
    tree_ = tree;
}

WFViewer::WFViewer(string name, TH1F* h_template):
    name_(name), channel_(0), h_pull_(("h_pull_"+name).c_str(), "", 1024, 0, 204.8)
{
    h_template_ = *h_template;
    f_template_ = new TF1();
    tree_ = new TTree();
}

//**********destructors*******************************************************************
WFViewer::~WFViewer()
{}

//**********Setters***********************************************************************
//----------Set the DigiTree--------------------------------------------------------------
void WFViewer::SetTemplate(TH1F* h_template)
{
    h_template_ = *h_template;
    
    return;
}    

//----------Set the DigiTree--------------------------------------------------------------
void WFViewer::SetTree(const char* digi_tree, const char* wf_tree)
{
    if(!tree_ || tree_->GetName() != digi_tree)
    {
        TFile* currentFile = gROOT->GetFile();
        if(currentFile)
        {            
            tree_ = (TTree*)currentFile->Get(digi_tree);
            h_pull_.SetAxisRange(-1, 1, "Y");
            h_pull_.SetFillColor(kBlack);
            h_pull_.SetFillStyle(3244);
            tree_->AddFriend(wf_tree);
            gStyle->SetOptStat("");
        }
    }
    if(!tree_)
        cout << "ERROR SetTree(): TTree '" << digi_tree << "' not found" << endl;
    
    return;
}

//**********Utils*************************************************************************
//----------Draw template with parameters taken from fit----------------------------------
//----------all the events up to *max_entries*,
//----------wf are shifted accordingly to the time of *ref* channel
void WFViewer::Draw(string ref, const char* cut, Long64_t max_entries)
{
    //---if a tree hasn't been loaded yet try to get the default one
    if(tree_->GetNbranches() == 0)
        SetTree();
    
    //---fill interpolator data (only once)
    if(f_template_)
        f_template_->Delete();
    InterpolatorFunc* interpolator_ = new InterpolatorFunc(MAX_INTERPOLATOR_POINTS,
                                                           ROOT::Math::Interpolation::kCSPLINE,
                                                           1., 0.);
    float offset = h_template_.GetBinCenter(h_template_.GetMaximumBin());
    vector<double> x, y;
    for(int iBin=1; iBin<=h_template_.GetNbinsX(); ++iBin)
    {
        x.push_back(h_template_.GetBinCenter(iBin));
        y.push_back(h_template_.GetBinContent(iBin));
    }
    interpolator_->SetData(x, y);
    f_template_ = new TF1(("f_"+name_).c_str(), interpolator_, 0, 180, 0);
    f_template_->SetNpx(10000);
    f_template_->SetLineColor(kRed+1);

    //---Draw event WF + fit result + residuals
    TCanvas* cnv = new TCanvas();
    TPad* p1 = new TPad("wf", "", 0.0, 0.4, 1.0, 1.0, 21);
    TPad* p2 = new TPad("pull", "", 0.0, 0.0, 1.0, 0.4, 21);
    p1->SetFillColor(kWhite);
    p2->SetFillColor(kWhite);
    p1->Draw();
    p2->Draw();
    p1->cd();
    tree_->Draw(("WF_val/fit_ampl["+name_+"]:WF_time-time["+ref+"]>>h(1500,-100,200,1000,-0.1,1.1)").c_str(),
                ("WF_ch=="+name_+"&&"+cut).c_str(), "", max_entries);
    tree_->Draw(("WF_val/fit_ampl["+name_+"]:WF_time-time["+ref+"]>>prof(1500,-100,200,-0.1,1.1)").c_str(),
                ("WF_ch=="+name_+"&&"+cut).c_str(), "PROFgoff", max_entries);
    f_template_->Draw("same");

    p2->cd();
    TH1F* h_tmp = (TH1F*)gDirectory->Get("prof");
    TH1F* h_pull = new TH1F("h_pull", "", 1500, -100, 200);
    h_pull->SetAxisRange(-1, 1, "Y");
    for(int i=1; i<h_tmp->GetNbinsX(); ++i)
    {
        if(f_template_->Eval(h_tmp->GetBinCenter(i))!=0)
            h_pull->SetBinContent(i, (h_tmp->GetBinContent(i)-f_template_->Eval(h_tmp->GetBinCenter(i)))/
                                  f_template_->Eval(h_tmp->GetBinCenter(i)));
        else
        {
            h_pull->SetBinContent(i, 0);
            h_pull->SetBinError(i+1, 1);
        }
    }    
    h_pull->Draw("E3");        

    return;
}    

//----------Draw template with parameters taken from fit----------------------------------
//----------only one selected event
void WFViewer::Draw(unsigned int iEntry, const char* wf_tree)
{
    //---if a tree hasn't been loaded yet try to get the default one
    if(tree_->GetNbranches() == 0)
        SetTree();

    //---set relevant branches digi_tree
    unsigned int n_channels = 0;
    tree_->SetBranchAddress("n_channels", &n_channels);    
    tree_->GetEntry(iEntry);
    float fit_ampl[n_channels];
    float fit_time[n_channels];
    float fit_chi2[n_channels];
    tree_->SetBranchAddress("fit_ampl", fit_ampl);
    tree_->SetBranchAddress("fit_time", fit_time);
    tree_->SetBranchAddress("fit_chi2", fit_chi2);
    tree_->SetBranchAddress(name_.c_str(), &channel_);
    tree_->GetEntry(iEntry);
    //---set relevant branches wf_tree
    int wf_samples=0;
    tree_->SetBranchAddress("WF_samples", &wf_samples);    
    tree_->GetEntry(iEntry);
    float* wf_val = new float[wf_samples];
    float* wf_time = new float[wf_samples];
    tree_->SetBranchAddress("WF_val", wf_val);
    tree_->SetBranchAddress("WF_time", wf_time);
    tree_->GetEntry(iEntry);

    //---fill interpolator data (only once)
    if(f_template_)
        f_template_->Delete();
    InterpolatorFunc* interpolator_ = new InterpolatorFunc(MAX_INTERPOLATOR_POINTS,
                                                           ROOT::Math::Interpolation::kCSPLINE,
                                                           fit_ampl[channel_], fit_time[channel_]);
    float offset = h_template_.GetBinCenter(h_template_.GetMaximumBin());
    vector<double> x, y;
    for(int iBin=1; iBin<=h_template_.GetNbinsX(); ++iBin)
    {
        x.push_back(h_template_.GetBinCenter(iBin)-offset);
        y.push_back(h_template_.GetBinContent(iBin));
    }
    interpolator_->SetData(x, y);
    f_template_ = new TF1(("f_"+name_).c_str(), interpolator_, 0, 180, 0);
    f_template_->SetNpx(10000);
    f_template_->SetLineColor(kRed+1);

    //---Draw event WF + fit result + residuals
    TCanvas* cnv = new TCanvas();
    TPad* p1 = new TPad("wf", "", 0.0, 0.4, 1.0, 1.0, 21);
    TPad* p2 = new TPad("pull", "", 0.0, 0.0, 1.0, 0.4, 21);
    p1->SetFillColor(kWhite);
    p2->SetFillColor(kWhite);
    p1->Draw();
    p2->Draw();
    p1->cd();
    tree_->SetMarkerStyle(7);
    tree_->Draw("WF_val:WF_time>>h(1024,0,204.8)", ("WF_ch=="+name_+"&&index=="+to_string(iEntry)).c_str());
    f_template_->Draw("same");
    p2->cd();
    int samples=wf_samples/n_channels;
    for(int i=0; i<samples; ++i)
    {
        if(f_template_->Eval(wf_time[i+samples*channel_])!=0)            
            h_pull_.SetBinContent(i+1,
                                  (wf_val[i+samples*channel_]-f_template_->Eval(wf_time[i+samples*channel_]))/
                                  f_template_->Eval(wf_time[i+samples*channel_]));
        else
        {
            h_pull_.SetBinContent(i+1, 0);
            h_pull_.SetBinError(i+1, 1);
        }
    }    
    h_pull_.Draw("E3");
        
    return;
}
