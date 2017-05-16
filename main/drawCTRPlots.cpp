#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <map>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TApplication.h"


int ch1;
int ch2;

/*** tree variables ***/
struct TreeVars
{
  float t_beamX;
  float t_beamY;
  float* t_Vbias;
  float t_NINOthr;
  float* t_ped;
  float* t_amp;
  float* t_dur;
  float* t_time;
  int* t_isOk;
};

void InitTreeVars(TTree* chain1, TTree* chain2, TTree* chain3,
                  TreeVars& treeVars);

bool AcceptEvent(TreeVars treeVars);
bool AcceptEventTh(TreeVars treeVars, const float& thMin, const float& thMax);
bool AcceptEventAmp(TreeVars treeVars, const float& ampMin1, const float& ampMax1, const float& ampMin2, const float& ampMax2);
bool AcceptEventDur(TreeVars treeVars);



int main(int argc, char** argv)
{
  gSystem -> Load("CfgManager/lib/libCFGMan.so");
  
  if( argc < 2 )
  {
    std::cerr << ">>>>> drawCTRPlots.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return 1;
  }
  
  
  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  std::vector<std::string> inputFiles = opts.GetOpt<std::vector<std::string> >("Input.inputFiles");
  
  ch1 = opts.GetOpt<int>("Input.ch1");
  ch2 = opts.GetOpt<int>("Input.ch2");
  int VbiasIndex1 = opts.GetOpt<int>("Input.VbiasIndex1");
  int VbiasIndex2 = opts.GetOpt<int>("Input.VbiasIndex2");
  
  int cut_batch = opts.GetOpt<int>("Cuts.batch");
  
  float cut_NINOthrMin = opts.GetOpt<float>("Cuts.minThreshold");
  float cut_NINOthrMax = opts.GetOpt<float>("Cuts.maxThreshold");
  
  std::vector<float> cut_VbiasValues = opts.GetOpt<std::vector<float> >("Cuts.VbiasValues");
  
  std::vector<float> cut_minAmplitudes1 = opts.GetOpt<std::vector<float> >("Cuts.minAmplitudes1");
  std::vector<float> cut_maxAmplitudes1 = opts.GetOpt<std::vector<float> >("Cuts.maxAmplitudes1");
  std::vector<float> cut_minAmplitudes2 = opts.GetOpt<std::vector<float> >("Cuts.minAmplitudes2");
  std::vector<float> cut_maxAmplitudes2 = opts.GetOpt<std::vector<float> >("Cuts.maxAmplitudes2");
  std::map<float,float> cut_ampMin1;
  std::map<float,float> cut_ampMax1;
  std::map<float,float> cut_ampMin2;
  std::map<float,float> cut_ampMax2;
  for(unsigned int it = 0; it < cut_VbiasValues.size(); ++it)
  {
    cut_ampMin1[cut_VbiasValues.at(it)] = cut_minAmplitudes1.at(it);
    cut_ampMax1[cut_VbiasValues.at(it)] = cut_maxAmplitudes1.at(it);
    cut_ampMin2[cut_VbiasValues.at(it)] = cut_minAmplitudes2.at(it);
    cut_ampMax2[cut_VbiasValues.at(it)] = cut_maxAmplitudes2.at(it);
  }
  
  std::vector<float> cut_minDurations1 = opts.GetOpt<std::vector<float> >("Cuts.minDurations1");
  std::vector<float> cut_maxDurations1 = opts.GetOpt<std::vector<float> >("Cuts.maxDurations1");
  std::vector<float> cut_minDurations2 = opts.GetOpt<std::vector<float> >("Cuts.minDurations2");
  std::vector<float> cut_maxDurations2 = opts.GetOpt<std::vector<float> >("Cuts.maxDurations2");
  std::map<float,float> cut_durMin1;
  std::map<float,float> cut_durMax1;
  std::map<float,float> cut_durMin2;
  std::map<float,float> cut_durMax2;
  for(unsigned int it = 0; it < cut_VbiasValues.size(); ++it)
  {
    cut_durMin1[cut_VbiasValues.at(it)] = cut_minDurations1.at(it);
    cut_durMax1[cut_VbiasValues.at(it)] = cut_maxDurations1.at(it);
    cut_durMin2[cut_VbiasValues.at(it)] = cut_minDurations2.at(it);
    cut_durMax2[cut_VbiasValues.at(it)] = cut_maxDurations2.at(it);
  }
  
  int rebin = opts.GetOpt<int>("Plots.rebin");
  std::string label1 = opts.GetOpt<std::string>("Plots.label1");
  std::string label2 = opts.GetOpt<std::string>("Plots.label2");
  
  
  //------------------------
  // labels and canvas style
  setTDRStyle();
  TApplication* theApp = new TApplication("App", &argc, argv);

  TLatex* latexLabel1 = new TLatex(0.13,0.97,Form("%s",label1.c_str()));
  TLatex* latexLabel2 = new TLatex(0.13,0.97,Form("%s",label2.c_str()));
  TLatex* latexLabel12 = new TLatex(0.13,0.97,Form("%s -- %s",label1.c_str(),label2.c_str()));
  latexLabel1 -> SetNDC();
  latexLabel1 -> SetTextFont(42);
  latexLabel1 -> SetTextSize(0.03);
  latexLabel2 -> SetNDC();
  latexLabel2 -> SetTextFont(42);
  latexLabel2 -> SetTextSize(0.03);
  latexLabel12 -> SetNDC();
  latexLabel12 -> SetTextFont(42);
  latexLabel12 -> SetTextSize(0.03);  
  
  
  //---------------------------
  // open input files and trees
  TChain* chain1 = new TChain("info","info");
  TChain* chain2 = new TChain("wire","wire");
  TChain* chain3 = new TChain("digi","digi");
  for(unsigned int fileIt = 0; fileIt < inputFiles.size(); ++fileIt)
  {
    chain1 -> Add(inputFiles.at(fileIt).c_str());
    chain2 -> Add(inputFiles.at(fileIt).c_str());
    chain3 -> Add(inputFiles.at(fileIt).c_str());
  }
  chain2 -> BuildIndex("index");
  chain1 -> AddFriend("wire");
  chain3 -> BuildIndex("index");
  chain1 -> AddFriend("digi");
  chain1 -> BuildIndex("index");
  std::cout << " Read " << chain1->GetEntries() << " total events in tree " << chain1->GetName() << std::endl;
  std::cout << " Read " << chain2->GetEntries() << " total events in tree " << chain2->GetName() << std::endl;
  std::cout << " Read " << chain3->GetEntries() << " total events in tree " << chain3->GetName() << std::endl;
  
  // set branches
  TreeVars treeVars;
  InitTreeVars(chain1,chain2,chain1,treeVars);
  
  
  //------------------
  // Define histograms
  std::map<std::string,TH1F*> h_amp1;
  std::map<std::string,TH1F*> h_amp2;
  std::map<std::string,TH1F*> h_ampRatio;
  
  std::map<std::string,TH1F*> h_dur1;
  std::map<std::string,TH1F*> h_dur2;
  
  TH1F* h_CTR = new TH1F("h_CTR","",20000,-10.,10.);
  TH1F* h_CTR_corr12 = new TH1F("h_CTR_corr12","",20000,-10.,10.);
  TH1F* h_CTR_corrRatio = new TH1F("h_CTR_corrRatio","",20000,-10.,10.);
  std::map<std::string,TH1F*> map_CTR_vs_Vbias_th;
  std::map<std::string,TH1F*> map_CTR_corr12_vs_Vbias_th;
  std::map<std::string,TH1F*> map_CTR_corrRatio_vs_Vbias_th;
  
  TProfile* p_CTR_vs_amp1 = new TProfile("p_CTR_vs_amp1","",100,0.,1.00);
  TProfile* p_CTR_vs_amp2 = new TProfile("p_CTR_vs_amp2","",100,0.,1.00);
  TProfile* p_CTR_vs_ampRatio = new TProfile("p_CTR_vs_ampRatio","",100,0.,5.);  
  std::map<std::string,TProfile*> map_CTR_vs_amp1_vs_Vbias_th;
  std::map<std::string,TProfile*> map_CTR_vs_amp2_vs_Vbias_th;
  std::map<std::string,TProfile*> map_CTR_vs_ampRatio_vs_Vbias_th;
  
  TH2F* h2_beam_Y_vs_X = new TH2F("h2_beam_Y_vs_X","",200,-50.,50.,200,-50.,50.);
  TH2F* h2_beam_Y_vs_X_cut = new TH2F("h2_beam_Y_vs_X_cut","",200,-50.,50.,200,-50.,50.);
  
  
  //-----------------------
  // first loop over events
  int nEntries = chain1 -> GetEntries();
  std::vector<std::pair<float,float> > pairs_Vbias;
  std::vector<std::pair<std::pair<float,float>,float> > pairs_Vbias_th;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> loop 1/4: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    chain1 -> GetEntry(entry);
    
    if( !AcceptEvent(treeVars) ) continue;
    
    float Vbias1 = treeVars.t_Vbias[VbiasIndex1];
    float Vbias2 = treeVars.t_Vbias[VbiasIndex2];
    float amp1 = treeVars.t_amp[8+ch1] * 0.25;
    float amp2 = treeVars.t_amp[8+ch2] * 0.25;
    float dur1 = treeVars.t_dur[ch1] * 0.2;
    float dur2 = treeVars.t_dur[ch2] * 0.2;
    float time1 = treeVars.t_time[ch1];
    float time2 = treeVars.t_time[ch2];
    float CTR = time2 - time1;
    
    std::pair<float,float> pair_Vbias(Vbias1,Vbias2);
    std::string label_Vbias = std::string(Form("%.0fV-%.0fV",Vbias1,Vbias2));
    if( h_amp1[label_Vbias] == NULL )
    {
      pairs_Vbias.push_back(pair_Vbias);
      
      h_amp1[label_Vbias] = new TH1F(Form("h_amp1_%s",label_Vbias.c_str()),"",1000,0.,1000.);
      h_amp2[label_Vbias] = new TH1F(Form("h_amp2_%s",label_Vbias.c_str()),"",1000,0.,1000.);
      h_ampRatio[label_Vbias] = new TH1F(Form("h_ampRatio_%s",label_Vbias.c_str()),"",1000,0.,3.);
      
      h_dur1[label_Vbias] = new TH1F(Form("h_dur1_%s",label_Vbias.c_str()),"",1000,0.,500.);
      h_dur2[label_Vbias] = new TH1F(Form("h_dur2_%s",label_Vbias.c_str()),"",1000,0.,500.);
    }
    
    h_amp1[label_Vbias] -> Fill( amp1 );
    h_amp2[label_Vbias] -> Fill( amp2 );
    h_ampRatio[label_Vbias] -> Fill( amp2/amp1 );
    
    if( !AcceptEventAmp(treeVars,cut_ampMin1[Vbias1],cut_ampMax1[Vbias1],cut_ampMin2[Vbias2],cut_ampMax2[Vbias2]) ) continue;
    
    h_dur1[label_Vbias] -> Fill( dur1 );
    h_dur2[label_Vbias] -> Fill( dur2 );
    
    if( !AcceptEventDur(treeVars) ) continue;
    
    std::string label_th = std::string(Form("%.0fmV",treeVars.t_NINOthr));
    std::string label_Vbias_th = std::string(Form("%s_%s",label_Vbias.c_str(),label_th.c_str()));
    if( map_CTR_vs_Vbias_th[label_Vbias_th] == NULL )
    {
      std::pair<std::pair<float,float>,float> pair_Vbias_th(pair_Vbias,treeVars.t_NINOthr);
      pairs_Vbias_th.push_back(pair_Vbias_th);
      map_CTR_vs_Vbias_th[label_Vbias_th] = new TH1F(Form("h_CTR_%s",label_Vbias_th.c_str()),"",20000,-10.,10.);
    }
    map_CTR_vs_Vbias_th[label_Vbias_th] -> Fill( CTR );
  }
  std::cout << std::endl;
  
  
  
  // draw plots
  TCanvas* c1;
  
  for(unsigned int it = 0; it < pairs_Vbias.size(); ++it)
  {
    float Vbias1 = pairs_Vbias.at(it).first;
    float Vbias2 = pairs_Vbias.at(it).second;
    std::string label_Vbias = std::string(Form("%.0fV-%.0fV",Vbias1,Vbias2));
    
    c1 = new TCanvas(Form("c1_amplitudes_%s",label_Vbias.c_str()),Form("Vbias %s",label_Vbias.c_str()),900,700);
    c1 -> Divide(2,2);
    c1 -> cd(1);
    gPad -> SetLogy();
    h_amp1[label_Vbias] -> GetXaxis() -> SetRangeUser(0.,1000.);
    h_amp1[label_Vbias] -> SetTitle(";max. amplitude (mV); events");
    h_amp1[label_Vbias] -> Draw();
    TLine* line_cutAmpMin1 = new TLine(cut_ampMin1[Vbias1],0.,cut_ampMin1[Vbias1],1.05*h_amp1[label_Vbias]->GetMaximum());
    TLine* line_cutAmpMax1 = new TLine(cut_ampMax1[Vbias1],0.,cut_ampMax1[Vbias1],1.05*h_amp1[label_Vbias]->GetMaximum());
    line_cutAmpMin1 -> SetLineColor(kRed);
    line_cutAmpMax1 -> SetLineColor(kRed);
    line_cutAmpMin1 -> Draw("same");
    line_cutAmpMax1 -> Draw("same");
    latexLabel1 -> Draw("same");
    c1 -> cd(2);
    gPad -> SetLogy();
    h_amp2[label_Vbias] -> GetXaxis() -> SetRangeUser(0.,1000.);
    h_amp2[label_Vbias] -> SetTitle(";max. amplitude (mV); events");
    h_amp2[label_Vbias] -> Draw();
    TLine* line_cutAmpMin2 = new TLine(cut_ampMin2[Vbias2],0.,cut_ampMin2[Vbias2],1.05*h_amp2[label_Vbias]->GetMaximum());
    TLine* line_cutAmpMax2 = new TLine(cut_ampMax2[Vbias2],0.,cut_ampMax2[Vbias2],1.05*h_amp2[label_Vbias]->GetMaximum());
    line_cutAmpMin2 -> SetLineColor(kRed);
    line_cutAmpMax2 -> SetLineColor(kRed);
    line_cutAmpMin2 -> Draw("same");
    line_cutAmpMax2 -> Draw("same");
    latexLabel2 -> Draw("same");
    gPad -> Update();
    
    c1 -> cd(3);
    gPad -> SetLogy();
    h_dur1[label_Vbias] -> GetXaxis() -> SetRangeUser(0.,200.);
    h_dur1[label_Vbias] -> SetTitle(";NINO pulse length (ns); events");
    h_dur1[label_Vbias] -> Draw();
    TLine* line_cutDurMin1 = new TLine(cut_durMin1[Vbias1],0.,cut_durMin1[Vbias1],1.05*h_dur1[label_Vbias]->GetMaximum());
    TLine* line_cutDurMax1 = new TLine(cut_durMax1[Vbias1],0.,cut_durMax1[Vbias1],1.05*h_dur1[label_Vbias]->GetMaximum());
    line_cutDurMin1 -> SetLineColor(kRed);
    line_cutDurMax1 -> SetLineColor(kRed);
    line_cutDurMin1 -> Draw("same");
    line_cutDurMax1 -> Draw("same");
    latexLabel1 -> Draw("same");
    c1 -> cd(4);
    gPad -> SetLogy();
    h_dur2[label_Vbias] -> GetXaxis() -> SetRangeUser(0.,200.);
    h_dur2[label_Vbias] -> SetTitle(";NINO pulse length (ns); events");
    h_dur2[label_Vbias] -> Draw();
    TLine* line_cutDurMin2 = new TLine(cut_durMin2[Vbias2],0.,cut_durMin2[Vbias2],1.05*h_dur2[label_Vbias]->GetMaximum());
    TLine* line_cutDurMax2 = new TLine(cut_durMax2[Vbias2],0.,cut_durMax2[Vbias2],1.05*h_dur2[label_Vbias]->GetMaximum());
    line_cutDurMin2 -> SetLineColor(kRed);
    line_cutDurMax2 -> SetLineColor(kRed);
    line_cutDurMin2 -> Draw("same");
    line_cutDurMax2 -> Draw("same");
    latexLabel2 -> Draw("same");
    gPad -> Update();
  }
  
  
  
  //-------------------
  // CTR threshold scan
  
  std::map<std::string,float> map_mean_vs_Vbias_th;
  std::map<std::string,TGraph*> g_mean_vs_Vbias;
  std::map<std::string,TGraph*> g_CTR_vs_Vbias;
  std::map<std::string,TGraph*> g_mean_vs_th;
  std::map<std::string,TGraph*> g_CTR_vs_th;
  
  std::cout << "size: " << pairs_Vbias_th.size() << std::endl;
  for(unsigned int it = 0; it < pairs_Vbias_th.size(); ++it)
  {
    float Vbias1 = pairs_Vbias_th.at(it).first.first;
    float Vbias2 = pairs_Vbias_th.at(it).first.second;
    float th = pairs_Vbias_th.at(it).second;
    std::string label_Vbias = std::string(Form("%.0fV-%.0fV"),Vbias1,Vbias2);
    std::string label_th    = std::string(Form("%.0fmV"),th);
    std::string label_Vbias_th = label_Vbias + "_" + label_th;
    std::cout << label_Vbias_th << std::endl;
    
    // if( mapIt->first > maxThVal) maxThVal = mapIt -> first;
    
    TH1F* histo = map_CTR_vs_Vbias_th[label_Vbias_th];
    float* vals = new float[4];
    FindSmallestInterval(vals,histo,0.68,true); 
    
    float mean = vals[0];      
    float min = vals[2];
    float max = vals[3];
    float delta = max-min;
    float sigma = 0.5*delta;
    float effSigma = sigma;
    
    map_mean_vs_Vbias_th[label_Vbias_th] = mean;
    g_mean_vs_Vbias[label_th] -> SetPoint(it,Vbias1,mean);
    g_CTR_vs_Vbias[label_th] -> SetPoint(it,Vbias1,effSigma);
    g_mean_vs_th[label_Vbias] -> SetPoint(it,th,mean);
    g_CTR_vs_th[label_Vbias] -> SetPoint(it,th,effSigma);
  }
  
  
  /*
  //------------------------
  // second loop over events
  
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> loop 2/4: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    chain1 -> GetEntry(entry);
    
    if( !AcceptEvent(treeVars) ) continue;
    
    h2_beam_Y_vs_X -> Fill(treeVars.t_beamX,treeVars.t_beamY);
    
    if( !AcceptEventAmp(treeVars,cut_ampMin1,cut_ampMax1,cut_ampMin2,cut_ampMax2) ) continue;
    if( !AcceptEventDur(treeVars) ) continue;
    
    h2_beam_Y_vs_X_cut -> Fill(treeVars.t_beamX,treeVars.t_beamY);
    
    float amp1 = treeVars.t_amp[8+ch1] / 4000.;
    float amp2 = treeVars.t_amp[8+ch2] / 4000.;
    float time1 = treeVars.t_time[ch1];
    float time2 = treeVars.t_time[ch2];
    float CTR = (time2 - time1) - map_mean_vs_th[treeVars.t_NINOthr] + map_mean_vs_th[maxThVal];
    
    if( CTR > -1. && CTR < 1. )
    {
      if( map_CTR_vs_amp1_vs_th[treeVars.t_NINOthr] == NULL )
      {
        map_CTR_vs_amp1_vs_th[treeVars.t_NINOthr] = new TProfile(Form("p_CTR_vs_amp1_th%f",treeVars.t_NINOthr),"",100,0.,0.50);
        map_CTR_vs_amp2_vs_th[treeVars.t_NINOthr] = new TProfile(Form("p_CTR_vs_amp2_th%f",treeVars.t_NINOthr),"",100,0.,0.50);
        map_CTR_vs_ampRatio_vs_th[treeVars.t_NINOthr] = new TProfile(Form("p_CTR_vs_ampRatio_th%f",treeVars.t_NINOthr),"",100,0.,5.);
      }
      map_CTR_vs_amp1_vs_th[treeVars.t_NINOthr] -> Fill( amp1,CTR );
      map_CTR_vs_amp2_vs_th[treeVars.t_NINOthr] -> Fill( amp2,CTR );
      map_CTR_vs_ampRatio_vs_th[treeVars.t_NINOthr] -> Fill( amp2/amp1,CTR );
    }
    
    if( !AcceptEventTh(treeVars,cut_NINOthrMin,cut_NINOthrMax) ) continue;
    
    h_CTR -> Fill( CTR );
    if( CTR > -1. && CTR < 1. )
    {
      p_CTR_vs_amp1 -> Fill( amp1,CTR );
      p_CTR_vs_amp2 -> Fill( amp2,CTR );
      p_CTR_vs_ampRatio -> Fill( amp2/amp1,CTR );
    }
  }
  std::cout << std::endl;
  
  
  
  
  TCanvas* c1_beam_Y_vs_X = new TCanvas("c1_beam_Y_vs_X","beam profile",1050,450);
  c1_beam_Y_vs_X -> Divide(2,1);
  c1_beam_Y_vs_X -> cd(1);
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(-50.,-50.,50.,50.) );
  hPad -> SetTitle(";beam x (mm); beam y (mm)");
  hPad -> Draw();
  h2_beam_Y_vs_X -> Draw("colz");
  c1_beam_Y_vs_X -> cd(2);
  hPad = (TH1F*)( gPad->DrawFrame(-50.,-50.,50.,50.) );
  hPad -> SetTitle(";beam x (mm); beam y (mm)");
  hPad -> Draw();
  h2_beam_Y_vs_X_cut -> Draw("colz");
  gPad -> Update();
  c1_beam_Y_vs_X -> Print(Form("plots/c_beamProfile__%s__%s.png",label1.c_str(),label2.c_str()));
  
  
  
  //------------------
  // fit CTR histogram
  
  float* vals = new float[4];
  FindSmallestInterval(vals,h_CTR,0.68,true); 
  h_CTR -> Rebin(rebin);
  
  float norm = h_CTR -> GetMaximum();
  float mean = vals[0];
  float min = vals[2];
  float max = vals[3];
  float delta = max-min;
  float sigma = 0.5*delta;
  float effSigma = sigma;
  min = min - 2.*delta;
  max = max + 2.*delta;

  TCanvas* c1_CTR = new TCanvas("c1_CTR","CTR",700,600);
  h_CTR -> SetTitle(";#Deltat = time_{xtal2} #minus time_{xtal1} (ns);events");
  h_CTR -> SetMarkerStyle(24);
  h_CTR -> SetMarkerColor(kRed);
  h_CTR -> SetLineColor(kRed);
  h_CTR -> GetXaxis() -> SetRangeUser(min,max);
  h_CTR -> Draw("PE");
  latexLabel12 -> Draw("same");
  
  // gaus fit
  std::string gausName = Form("fitFunc_gaus");
  TF1* fitFunc_gaus = new TF1(gausName.c_str(),"[0]*exp(-1.*(x-[1])*(x-[1])/(2.*[2]*[2]))",mean-sigma,mean+sigma);
  fitFunc_gaus -> SetNpx(10000);
  fitFunc_gaus -> SetParameters(norm,mean,sigma);
  h_CTR -> Fit(gausName.c_str(),"QNRSL");
  norm = fitFunc_gaus -> GetParameter(0);
  mean = fitFunc_gaus -> GetParameter(1);
  sigma = fitFunc_gaus -> GetParameter(2);
  float sigmaErr = fitFunc_gaus -> GetParError(2);
  
  TLatex* latex22 = new TLatex(0.16,0.55,Form("#sigma_{CTR}^{gaus} = (%.1f #pm %.1f) ps",fabs(sigma*1000),fabs(sigmaErr*1000)));
  latex22 -> SetNDC();
  latex22 -> SetTextFont(42);
  latex22 -> SetTextSize(0.03);
  latex22 -> SetTextColor(kRed);
  latex22 -> Draw("same");
  
  // crystal ball fit
  std::string cbName = Form("fitFunc_cb");
  TF1* fitFunc_cb = new TF1(cbName.c_str(),crystalBallLowHigh,mean-4.*sigma,mean+4.*sigma,8);
  fitFunc_cb -> SetNpx(10000);
  fitFunc_cb -> SetParameters(norm,mean,sigma,1.5,10.,1.5,10.);
  h_CTR -> Fit(cbName.c_str(),"QNRSL");
  norm = fitFunc_cb -> GetParameter(0);
  mean = fitFunc_cb -> GetParameter(1);
  sigma = fitFunc_cb -> GetParameter(2);
  sigmaErr = fitFunc_cb -> GetParError(2);
  fitFunc_cb -> Draw("same");
  
  TLatex* latex2 = new TLatex(0.16,0.89,Form("eff. #sigma_{CTR} = %.1f ps",fabs(effSigma*1000)));
  latex2 -> SetNDC();
  latex2 -> SetTextFont(42);
  latex2 -> SetTextSize(0.05);
  latex2 -> SetTextColor(kRed);
  latex2 -> Draw("same");

  TLatex* latex21 = new TLatex(0.16,0.40,Form("#sigma_{CTR}^{c.b.} = (%.1f #pm %.1f) ps",fabs(sigma*1000),fabs(sigmaErr*1000)));
  latex21 -> SetNDC();
  latex21 -> SetTextFont(42);
  latex21 -> SetTextSize(0.03);
  latex21 -> SetTextColor(kRed);
  latex21 -> Draw("same");
  
  gPad -> Update();
  
  
  TCanvas* c1_thScan = new TCanvas("c1_CTR_vs_th","threshold scan",1050,450);
  c1_thScan -> Divide(2,1);
  c1_thScan -> cd(1);
  hPad = (TH1F*)( gPad->DrawFrame(0.,min,2100.,max) );
  hPad -> SetTitle(";NINO th. (mV); #LT #Deltat #GT (ns)");
  hPad -> Draw();
  g_mean_vs_th -> SetMarkerStyle(24);
  g_mean_vs_th -> SetMarkerColor(kRed);
  g_mean_vs_th -> SetLineColor(kRed);
  g_mean_vs_th -> Draw("PL,same");
  c1_thScan -> cd(2);
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,2100.,0.2) );
  hPad -> SetTitle(";NINO th. (mV); CTR (ns)");
  hPad -> Draw();
  g_CTR_vs_th -> SetMarkerStyle(24);
  g_CTR_vs_th -> SetMarkerColor(kRed);
  g_CTR_vs_th -> SetLineColor(kRed);
  g_CTR_vs_th -> Draw("PL,same");
  
  gPad -> Update();
  
  
  
  
  //---------------------
  // time walk correction
  
  c1 = new TCanvas("c1_ampCorrection","time walk",1050,450);
  c1 -> Divide(2,1);
  c1 -> cd(1);
  p_CTR_vs_amp1 -> SetMarkerStyle(24);
  p_CTR_vs_amp1 -> SetMarkerColor(kRed);
  p_CTR_vs_amp1 -> SetLineColor(kRed);
  p_CTR_vs_amp1 -> SetTitle(";amplitude_{xtal1 or xtal2};#Deltat = time_{xtal2} #minus time_{xtal1} (ns)");
  p_CTR_vs_amp1 -> SetMinimum(min-2.*delta);
  p_CTR_vs_amp1 -> SetMaximum(max+2.*delta);
  p_CTR_vs_amp1 -> Draw();
  p_CTR_vs_amp2 -> Draw("same");
  p_CTR_vs_amp2 -> SetMarkerStyle(21);
  p_CTR_vs_amp2 -> SetMarkerColor(kBlue);
  p_CTR_vs_amp2 -> SetLineColor(kBlue);
  p_CTR_vs_amp2 -> SetTitle(";amplitude_{xtal1 or xtal2};#Deltat = time_{xtal2} #minus time_{xtal1} (ns)");
  p_CTR_vs_amp2 -> SetMinimum(min-2.*delta);
  p_CTR_vs_amp2 -> SetMaximum(max+2.*delta);
  p_CTR_vs_amp2 -> Draw("same");
  
  TF1* fitFunc_corr1 = new TF1("fitFunc_corr1","[0]*log([1]*x)+[2]",0.,1.0);
  fitFunc_corr1 -> SetParameters(0.05,0.0001,0.);
  fitFunc_corr1 -> SetLineColor(kRed);
  p_CTR_vs_amp1 -> Fit("fitFunc_corr1","QNS+","",
                       h_amp1->GetMean()-4.*h_amp1->GetRMS(),
                       h_amp1->GetMean()+2.*h_amp1->GetRMS());
  fitFunc_corr1 -> Draw("same");
  
  TF1* fitFunc_corr2 = new TF1("fitFunc_corr2","[0]*log([1]*x)+[2]",0.,1.0);
  fitFunc_corr2 -> SetParameters(-1.,20.,0.);
  fitFunc_corr2 -> SetLineColor(kBlue);
  p_CTR_vs_amp2 -> Fit("fitFunc_corr2","QNS+","",
                       h_amp2->GetMean()-2.*h_amp2->GetRMS(),
                       h_amp2->GetMean()+4.*h_amp2->GetRMS());
  fitFunc_corr2 -> Draw("same");
  
  c1 -> cd(2);
  p_CTR_vs_ampRatio -> SetTitle(";amplitude_{xtal2} / amplitude_{xtal1};#Deltat = time_{xtal2} #minus time_{xtal1} (ns)");
  p_CTR_vs_ampRatio -> SetMinimum(min-2.*delta);
  p_CTR_vs_ampRatio -> SetMaximum(max+2.*delta);
  p_CTR_vs_ampRatio -> Draw();
    
  TF1* fitFunc_corrRatio = new TF1("fitFunc_corrRatio","[0]*log([1]*x)+[2]",0.,5.);
  fitFunc_corrRatio -> SetParameters(-0.2,0.1,0.);
  p_CTR_vs_ampRatio -> Fit("fitFunc_corrRatio","QNS+","",
                           h_ampRatio->GetMean()-2.*h_ampRatio->GetRMS(),
                           h_ampRatio->GetMean()+2.*h_ampRatio->GetRMS());
  fitFunc_corrRatio -> Draw("same");

  std::map<float,TF1*> map_fitFunc_corrRatio_vs_th;
  for(std::map<float,TProfile*>::const_iterator mapIt = map_CTR_vs_amp1_vs_th.begin(); mapIt != map_CTR_vs_amp1_vs_th.end(); ++mapIt)
  {
    map_fitFunc_corrRatio_vs_th[mapIt->first] = new TF1(Form("fitFunc_corrRatio_th%f",mapIt->first),"[0]*log([1]*x)+[2]",0.,5.);
    map_fitFunc_corrRatio_vs_th[mapIt->first] -> SetParameters(-0.3,1.,0.);
    map_CTR_vs_ampRatio_vs_th[mapIt->first] -> Fit(Form("fitFunc_corrRatio_th%f",mapIt->first),"QNS+","",
                                                   h_ampRatio->GetMean()-2.*h_ampRatio->GetRMS(),
                                                   h_ampRatio->GetMean()+2.*h_ampRatio->GetRMS());
  }
  gPad -> Update();
  
  
  
  
  //------------------------
  // third loop over events
  
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> loop 3/4: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    chain1 -> GetEntry(entry);
    
    if( !AcceptEvent(treeVars) ) continue;
    if( !AcceptEventAmp(treeVars,cut_ampMin1,cut_ampMax1,cut_ampMin2,cut_ampMax2) ) continue;
    if( !AcceptEventDur(treeVars) ) continue;
    
    float amp1 = treeVars.t_amp[8+ch1] / 4000.;
    float amp2 = treeVars.t_amp[8+ch2] / 4000.;
    float time1 = treeVars.t_time[ch1];
    float time2 = treeVars.t_time[ch2];
    float CTR = (time2 - time1);// - map_mean_vs_th[treeVars.t_NINOthr] + map_mean_vs_th[maxThVal];
    float CTR_corr12 = CTR -
      fitFunc_corr2->Eval(amp2) + fitFunc_corr2->Eval(h_amp2->GetMean()) -
      fitFunc_corr1->Eval(amp1) + fitFunc_corr1->Eval(h_amp1->GetMean());
    float CTR_corrRatio = CTR -
      map_fitFunc_corrRatio_vs_th[treeVars.t_NINOthr]->Eval(amp2/amp1) + map_fitFunc_corrRatio_vs_th[treeVars.t_NINOthr]->Eval(h_ampRatio->GetMean());
    
    if( map_CTR_corr12_vs_th[treeVars.t_NINOthr] == NULL )
    {
      map_CTR_corr12_vs_th[treeVars.t_NINOthr] = new TH1F(Form("h_CTR_corr12_th%f",treeVars.t_NINOthr),"",10000,-5.,5.);
    }
    map_CTR_corr12_vs_th[treeVars.t_NINOthr] -> Fill( CTR_corr12 );

    if( map_CTR_corrRatio_vs_th[treeVars.t_NINOthr] == NULL )
    {
      map_CTR_corrRatio_vs_th[treeVars.t_NINOthr] = new TH1F(Form("h_CTR_corrRatio_th%f",treeVars.t_NINOthr),"",10000,-5.,5.);
    }
    map_CTR_corrRatio_vs_th[treeVars.t_NINOthr] -> Fill( CTR_corrRatio );
  }
  std::cout << std::endl;
  

  
  //------------------------
  // CTR corr threshold scan

  std::map<float,float> map_mean_corr12_vs_th;
  std::map<float,float> map_mean_corrRatio_vs_th;
  TGraph* g_mean_corr12_vs_th = new TGraph();
  TGraph* g_mean_corrRatio_vs_th = new TGraph();
  TGraph* g_CTR_corr12_vs_th = new TGraph();
  TGraph* g_CTR_corrRatio_vs_th = new TGraph();
  
  it = 0;
  maxThVal = -999.;
  for(std::map<float,TH1F*>::const_iterator mapIt = map_CTR_corrRatio_vs_th.begin(); mapIt != map_CTR_corrRatio_vs_th.end(); ++mapIt)
  {
    std::cout << "th: " << mapIt->first << std::endl;
    if( mapIt->first > maxThVal) maxThVal = mapIt -> first;
    
    vals = new float[4];
    FindSmallestInterval(vals,mapIt->second,0.68,true); 
    
    mean = vals[0];      
    min = vals[2];
    max = vals[3];
    delta = max-min;
    sigma = 0.5*delta;
    effSigma = sigma;
    
    map_mean_corrRatio_vs_th[mapIt->first] = mean;
    g_mean_corrRatio_vs_th -> SetPoint(it,mapIt->first,mean);
    g_CTR_corrRatio_vs_th -> SetPoint(it,mapIt->first,effSigma);
    ++it;
  }

  c1_thScan -> cd(1);
  g_mean_corrRatio_vs_th -> SetMarkerStyle(21);
  g_mean_corrRatio_vs_th -> SetMarkerColor(kBlue);
  g_mean_corrRatio_vs_th -> SetLineColor(kBlue);
  g_mean_corrRatio_vs_th -> Draw("PL,same");
  c1_thScan -> cd(2);
  g_CTR_corrRatio_vs_th -> SetMarkerStyle(21);
  g_CTR_corrRatio_vs_th -> SetMarkerColor(kBlue);
  g_CTR_corrRatio_vs_th -> SetLineColor(kBlue);
  g_CTR_corrRatio_vs_th -> Draw("PL,same");
  
  gPad -> Update();
  
  
  

  //------------------------
  // fourth loop over events
  
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> loop 4/4: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    chain1 -> GetEntry(entry);
    
    if( !AcceptEvent(treeVars) ) continue;
    if( !AcceptEventAmp(treeVars,cut_ampMin1,cut_ampMax1,cut_ampMin2,cut_ampMax2) ) continue;
    if( !AcceptEventDur(treeVars) ) continue;
    if( !AcceptEventTh(treeVars,cut_NINOthrMin,cut_NINOthrMax) ) continue;
    
    float amp1 = treeVars.t_amp[8+ch1] / 4000.;
    float amp2 = treeVars.t_amp[8+ch2] / 4000.;
    float time1 = treeVars.t_time[ch1];
    float time2 = treeVars.t_time[ch2];
    float CTR = (time2 - time1) - map_mean_corrRatio_vs_th[treeVars.t_NINOthr] + map_mean_vs_th[maxThVal];
    
    float CTR_corr12 = CTR -
      fitFunc_corr2->Eval(amp2) + fitFunc_corr2->Eval(h_amp2->GetMean()) -
      fitFunc_corr1->Eval(amp1) + fitFunc_corr1->Eval(h_amp1->GetMean());
    float CTR_corrRatio = CTR -
      map_fitFunc_corrRatio_vs_th[treeVars.t_NINOthr]->Eval(amp2/amp1) + map_fitFunc_corrRatio_vs_th[treeVars.t_NINOthr]->Eval(h_ampRatio->GetMean());
    
    h_CTR_corr12 -> Fill( CTR_corr12 );
    h_CTR_corrRatio -> Fill( CTR_corrRatio );
    
  }
  std::cout << std::endl;

  
  c1_CTR -> cd();
  
  h_CTR_corr12 -> Rebin(rebin);
  h_CTR_corr12 -> SetMarkerStyle(21);
  h_CTR_corr12 -> SetMarkerColor(kBlue);
  h_CTR_corr12 -> SetLineColor(kBlue);
  
  FindSmallestInterval(vals,h_CTR_corrRatio,0.68,true); 
  h_CTR_corrRatio -> Rebin(rebin);
  h_CTR_corrRatio -> SetMarkerStyle(21);
  h_CTR_corrRatio -> SetMarkerColor(kBlue);
  h_CTR_corrRatio -> SetLineColor(kBlue);
  
  h_CTR -> SetMaximum(1.5*h_CTR_corrRatio->GetMaximum());
  h_CTR_corrRatio -> Draw("PE,same");
  
  gPad -> Update();  
  
  
  norm = h_CTR_corrRatio -> GetMaximum();  
  mean = vals[0];
  min = vals[2];
  max = vals[3];
  delta = max-min;
  sigma = 0.5*delta;
  effSigma = sigma;
  min = min - 2.*delta;
  max = max + 2.*delta;
  
  h_CTR_corrRatio -> GetXaxis() -> SetRangeUser(min,max);
  
  // gaus fit
  gausName = Form("fitFunc_gaus_corrRatio");
  TF1* fitFunc_gaus_corrRatio = new TF1(gausName.c_str(),"[0]*exp(-1.*(x-[1])*(x-[1])/(2.*[2]*[2]))",mean-sigma,mean+sigma);
  fitFunc_gaus_corrRatio -> SetNpx(10000);
  fitFunc_gaus_corrRatio -> SetParameters(norm,mean,sigma);
  fitFunc_gaus_corrRatio -> SetLineColor(kBlue);
  h_CTR_corrRatio -> Fit(gausName.c_str(),"QNRSL");
  norm = fitFunc_gaus_corrRatio -> GetParameter(0);
  mean = fitFunc_gaus_corrRatio -> GetParameter(1);
  sigma = fitFunc_gaus_corrRatio -> GetParameter(2);
  sigmaErr = fitFunc_gaus_corrRatio -> GetParError(2);
 
  TLatex* latex22_corr = new TLatex(0.16,0.50,Form("#sigma_{CTR}^{gaus} = (%.1f #pm %.1f) ps",fabs(sigma*1000),fabs(sigmaErr*1000)));
  latex22_corr -> SetNDC();
  latex22_corr -> SetTextFont(42);
  latex22_corr -> SetTextSize(0.03);
  latex22_corr -> SetTextColor(kBlue);
  latex22_corr -> Draw("same");
  
  // crystal ball fit
  cbName = Form("fitFunc_cb_corrRatio");
  TF1* fitFunc_cb_corrRatio = new TF1(cbName.c_str(),crystalBallLowHigh,mean-4.*sigma,mean+4.*sigma,8);
  fitFunc_cb_corrRatio -> SetNpx(10000);
  fitFunc_cb_corrRatio -> SetParameters(norm,mean,sigma,1.5,10.,1.5,10.);
  fitFunc_cb_corrRatio -> SetLineColor(kBlue);
  h_CTR_corrRatio -> Fit(cbName.c_str(),"QNRSL");
  norm = fitFunc_cb_corrRatio -> GetParameter(0);
  mean = fitFunc_cb_corrRatio -> GetParameter(1);
  sigma = fitFunc_cb_corrRatio -> GetParameter(2);
  sigmaErr = fitFunc_cb_corrRatio -> GetParError(2);
  
  fitFunc_cb_corrRatio -> Draw("same");
  
  TLatex* latex2_corr = new TLatex(0.16,0.80,Form("eff. #sigma_{CTR} = %.1f ps",fabs(effSigma*1000)));
  latex2_corr -> SetNDC();
  latex2_corr -> SetTextFont(42);
  latex2_corr -> SetTextSize(0.05);
  latex2_corr -> SetTextColor(kBlue);
  latex2_corr -> Draw("same");

  TLatex* latex21_corr = new TLatex(0.16,0.35,Form("#sigma_{CTR}^{c.b.} = (%.1f #pm %.1f) ps",fabs(sigma*1000),fabs(sigmaErr*1000)));
  latex21_corr -> SetNDC();
  latex21_corr -> SetTextFont(42);
  latex21_corr -> SetTextSize(0.03);
  latex21_corr -> SetTextColor(kBlue);
  latex21_corr -> Draw("same");
  
  gPad -> Update();
  c1_CTR -> Print(Form("plots/c_CTR__%s__%s.png",label1.c_str(),label2.c_str()));
  */
  
  
  theApp -> Run();
  return 0;
}



bool AcceptEvent(TreeVars treeVars)
{
  //if( treeVars.t_batch != cut_batch ) return false;
  
  //if( !(treeVars.t_isOk[ch1] && treeVars.t_isOk[ch2]) ) return false;
  
  // if( fabs(treeVars.t_beamX+22.) > 5. ) return false;
  // if( fabs(treeVars.t_beamY-6.) > 5. ) return false;
  
  return true;
}

bool AcceptEventTh(TreeVars treeVars, const float& thMin, const float& thMax)
{
  if( (treeVars.t_NINOthr < thMin) || (treeVars.t_NINOthr > thMax) ) return false;
  
  return true;
}

bool AcceptEventAmp(TreeVars treeVars, const float& ampMin1, const float& ampMax1, const float& ampMin2, const float& ampMax2)
{
  if( (treeVars.t_amp[8+ch1] < ampMin1) || (treeVars.t_amp[8+ch1] > ampMax1) ) return false;
  if( (treeVars.t_amp[8+ch2] < ampMin2) || (treeVars.t_amp[8+ch2] > ampMax2) ) return false;
  
  return true;
}

bool AcceptEventDur(TreeVars treeVars)
{
  // if( (treeVars.t_dur[ch1] < durMin1) || (treeVars.t_dur[ch1] > durMax1) ) return false;
  // if( (treeVars.t_dur[ch2] < durMin2) || (treeVars.t_dur[ch2] > durMax2) ) return false;
  
  return true;
}



void InitTreeVars(TTree* chain1, TTree* chain2, TTree* chain3,
                  TreeVars& treeVars)
{
  treeVars.t_Vbias = new float[3];
  treeVars.t_ped = new float[16];
  treeVars.t_amp = new float[16];
  treeVars.t_dur = new float[16];
  treeVars.t_time = new float[16];
  treeVars.t_isOk = new int[16];
  
  //chain1 -> SetBranchStatus("*",0);
  chain1 -> SetBranchStatus("NINOthr",1); chain1 -> SetBranchAddress("NINOthr",&treeVars.t_NINOthr);
  chain1 -> SetBranchStatus("Vbias1" ,1); chain1 -> SetBranchAddress("Vbias1", &treeVars.t_Vbias[0]);
  chain1 -> SetBranchStatus("Vbias2" ,1); chain1 -> SetBranchAddress("Vbias2", &treeVars.t_Vbias[1]);
  chain1 -> SetBranchStatus("Vbias3" ,1); chain1 -> SetBranchAddress("Vbias3", &treeVars.t_Vbias[2]);
  
  //chain2 -> SetBranchStatus("*",0);
  chain2 -> SetBranchStatus("X",1); chain2 -> SetBranchAddress("X",&treeVars.t_beamX);
  chain2 -> SetBranchStatus("Y",1); chain2 -> SetBranchAddress("Y",&treeVars.t_beamY);
  
  //chain3 -> SetBranchStatus("*",0);
  chain3 -> SetBranchStatus("amp_max",   1); chain3 -> SetBranchAddress("amp_max",   treeVars.t_amp);
  chain3 -> SetBranchStatus("charge_sig",1); chain3 -> SetBranchAddress("charge_sig",treeVars.t_dur);
  chain3 -> SetBranchStatus("time",      1); chain3 -> SetBranchAddress("time",      treeVars.t_time);
}
