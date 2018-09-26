#ifndef __DRAW_EVENT__
#define __DRAW_EVENT__

#include <unistd.h>
#include <csignal>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TBox.h"
#include "TApplication.h"
#include "TLatex.h"

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include "interface/PluginLoader.h"
#include "interface/RecoTree.h"
#include "interface/PluginBase.h"
#include "plugins/WFAnalyzer.h"

#include <netdb.h>
#include <unistd.h>

bool popupPlots = true;

string getMachineDomain()
{
  char hn[254];
  char *dn;
  struct hostent *hp;

  gethostname(hn, 254);
  hp = gethostbyname(hn);
  dn = strchr(hp->h_name, '.');
  if ( dn != NULL )
  {
    return std::string(++dn);
  }
  else
  {
    return "";
  }
}

//----------Get input files---------------------------------------------------------------
void ReadInputFiles(CfgManager& opts, TChain* inTree)
{
    int nFiles=0;
    string ls_command;
    string file;
    string path=opts.GetOpt<string>("h4reco.path2data");
    string run=opts.GetOpt<string>("h4reco.run");

    //---Get file list searching in specified path (eos or locally)
    if(path.find("/eos/cms") != string::npos)
    {
	if ( getMachineDomain() != "cern.ch" )
            ls_command = string("gfal-ls root://eoscms/"+path+run+" | grep 'root' > tmp/"+run+".list");
	else
            ls_command = string("ls "+path+run+" | grep 'root' > tmp/"+run+".list");
    }
    else if(path.find("srm://") != string::npos)
        ls_command = string("echo "+path+run+"/`gfal-ls "+path+run+
                            "` | sed -e 's:^.*\\/cms\\/:root\\:\\/\\/xrootd-cms.infn.it\\/\\/:g' | grep 'root' > tmp/"+run+".list");
    else
        ls_command = string("ls "+path+run+" | grep 'root' > tmp/"+run+".list");
    system(ls_command.c_str());
    ifstream waveList(string("tmp/"+run+".list").c_str(), ios::in);
    while(waveList >> file && (opts.GetOpt<int>("h4reco.maxFiles")<0 || nFiles<opts.GetOpt<int>("h4reco.maxFiles")) )
    {
        if(path.find("/eos/cms") != string::npos)
        {
            if ( getMachineDomain() != "cern.ch" )
            {
                std::cout << "+++ Adding file " << ("root://eoscms/"+path+run+"/"+file).c_str() << std::endl;
                inTree->AddFile(("root://eoscms/"+path+run+"/"+file).c_str());
            }
            else
            {
                std::cout << "+++ Adding file " << (path+run+"/"+file).c_str() << std::endl;
                inTree->AddFile((path+run+"/"+file).c_str());
            }
        }
        else if(path.find("srm://") != string::npos)
        {
            std::cout << "+++ Adding file " << file << std::endl;
            inTree->AddFile((file).c_str());
        }
        else
        {
            std::cout << "+++ Adding file " << (path+run+"/"+file).c_str() << std::endl;
            inTree->AddFile((path+run+"/"+file).c_str());
        }
        ++nFiles;
    }
    std::cout << "+++ Added " << nFiles << " files with " << inTree->GetEntries() << " events" << std::endl;
    return;
}

//**********MAIN**************************************************************************
int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        cout << argv[0] << " cfg file " << "[run] [event]" << endl; 
        return -1;
    }
    
    //---load options---    
    CfgManager opts;
    opts.ParseConfigFile(argv[1]);

    //-----input setup-----    
    if(argc > 2)
    {
        vector<string> run(1, argv[2]);
        opts.SetOpt("h4reco.run", run);
    }
    string run = opts.GetOpt<string>("h4reco.run");
    if(argc > 3)
    {
        vector<string> event(1, argv[3]);
        opts.SetOpt("h4reco.event", event);
    }
    string event = opts.GetOpt<string>("h4reco.event");
    
    int maxEvents = opts.OptExist("h4reco.maxEvents") ? opts.GetOpt<int>("h4reco.maxEvents") : -1;
    int spill     = opts.OptExist("h4reco.spill")     ? opts.GetOpt<int>("h4reco.spill")     : -1;
    
    std::cout << ">>> Getting tree" << std::endl;
    TChain* inTree = new TChain("H4tree");
    ReadInputFiles(opts, inTree);
    H4Tree h4Tree(inTree);
    
    //-----output setup-----
    uint64 index=stoul(run)*1e9;
    // TFile* outROOT = new TFile("H4draw"+TString(run)+"_"+TString(event)+".root", "RECREATE");
    // outROOT->cd();
    RecoTree mainTree(&index);
    
    //---Get plugin sequence---
    PluginLoader<PluginBase>* loader;
    vector<PluginLoader<PluginBase>* > pluginLoaders;
    map<string, PluginBase*> pluginMap;
    vector<PluginBase*> pluginSequence;
    vector<string> pluginList = opts.GetOpt<vector<string> >("h4reco.pluginList");
    //---plugin creation
    pluginLoaders.reserve(pluginList.size());
    for(auto& plugin : pluginList)
    {
      cout << ">>> Loading plugin <" << plugin << ">" << endl;
      //---create loader
      loader = new PluginLoader<PluginBase>(opts.GetOpt<string>(plugin+".pluginType"));
      pluginLoaders.push_back(loader);
      pluginLoaders.back()->Create();
      //---get instance and put it in the plugin sequence
      PluginBase* newPlugin = pluginLoaders.back()->CreateInstance();
      if(newPlugin)
      {
        pluginSequence.push_back(newPlugin);
        pluginSequence.back()->SetInstanceName(plugin);
        pluginMap[plugin] = pluginSequence.back();
      }
      else
      {
        cout << ">>> ERROR: plugin type " << opts.GetOpt<string>(plugin+".pluginType") << " is not defined." << endl;
        return 0;
      }
    }
    
    
    //---begin
    for(auto& plugin : pluginSequence)
    {
      //---call Begin() methods and check the return status
      bool r_status = plugin->Begin(opts, &index);
      if(!r_status)
      {
        cout << ">>> ERROR: plugin returned bad flag from Begin() call: " << plugin->GetInstanceName() << endl;
        exit(-1);
      }
      //---Get plugin shared data
      for(auto& shared : plugin->GetSharedData("", "TTree", true))
      {
        TTree* tree = (TTree*)shared.obj;
        tree->SetMaxVirtualSize(10000);
        // tree->SetDirectory(outROOT);
      }
    }
    
    
    //---interactive plots
    TApplication* theApp;
    if( popupPlots )
      theApp = new TApplication("App", &argc, argv);
    
    
    //---events loop
    cout << ">>> Processing H4DAQ run #" << run << " event # " << event << " / " << h4Tree.GetEntries() << " <<<" << endl;
    
    if( (atoi(event.c_str()) != -1) && (atoi(event.c_str()) < int(h4Tree.GetEntries())) )
    {
      h4Tree.NextEntry(atoi(event.c_str())-1);
      
      //---call ProcessEvent for each plugin and check the return status
      bool status=true;
      for(auto& plugin : pluginSequence)
      {
        if(status)
        {
          status = plugin->ProcessEvent(h4Tree, pluginMap, opts);
        }
      }
      
      //---fill the main tree with info variables and increase event counter
      if(status)
      {
        mainTree.time_stamp = h4Tree.evtTimeStart;
        mainTree.run = h4Tree.runNumber;
        mainTree.spill = h4Tree.spillNumber;
        mainTree.event = h4Tree.evtNumber;
        mainTree.Fill();
        ++index;
      }
      
      //---draw waveforms
      vector<string> channels = opts.GetOpt<vector<string> >("DigiReco.channelsNames");
      for(auto& channel : channels)
      {
        TCanvas* c1 = new TCanvas(Form("%s",channel.c_str()),Form("%s",channel.c_str()));
        c1 -> cd();
        
        auto shared_data = pluginMap["DigiReco"]->GetSharedData(std::string("DigiReco")+"_"+channel, "", false);
        WFClass* WF = (WFClass*)shared_data.at(0).obj;
        
        auto analizedWF = WF->GetSamples();
        float tUnit = WF->GetTUnit();
        int bWinMin = WF->GetBWinMin();
        int bWinMax = WF->GetBWinMax();
        int bIntWinMin = WF->GetBIntWinMin();
        int bIntWinMax = WF->GetBIntWinMax();
        int sWinMin = WF->GetSWinMin();
        int sWinMax = WF->GetSWinMax();
        int sIntWinMin = WF->GetSIntWinMin();
        int sIntWinMax = WF->GetSIntWinMax();
        float baselineRMS = WF->GetBaselineRMS();
        float maxSample = WF->GetMaxSample();
        float fitAmpMax = WF->GetFitAmpMax();
        float fitTimeMax = WF->GetFitTimeMax();
        float cfFrac = WF->GetCFFrac();
        float cfTime = WF->GetCFTime();
        float cfSlope = WF->GetCFSlope();
        float leThr = WF->GetLEThr();
        float leTime = WF->GetLETime();
        float leSlope = WF->GetLESlope();
        float teThr = WF->GetTEThr();
        float teTime = WF->GetTETime();
        float teSlope = WF->GetTESlope();
        TF1* funcAmp = WF->GetAmpFunc();
        
        std::cout << "\n\n\n***CHANNEL: " << channel << std::endl;
        std::cout << ">>> baseline RMS: " << baselineRMS << std::endl;
        
        TGraphErrors* g = new TGraphErrors();
        for(unsigned int jSample=0; jSample<analizedWF->size(); ++jSample)
        {
          g -> SetPoint(jSample,jSample,analizedWF->at(jSample));
          g -> SetPointError(jSample,0,baselineRMS);
        }
        
        g -> SetTitle(";sample;ADC");
        g -> SetMarkerStyle(20);
        g -> SetMarkerSize(0.7);
        g -> Draw("AL");
        g -> Draw("P,same");
        
        TLine* line_baseline = new TLine(bWinMin,0.,bWinMax,0.);
        line_baseline -> SetLineColor(kRed+2);
        line_baseline -> SetLineWidth(3);
        line_baseline -> Draw("same");
        
        TBox* box_baselineInt = new TBox(std::max(float(0.),float(bIntWinMin)),0.-fitAmpMax/10.,std::min(float(1023.),float(bIntWinMax)),0+fitAmpMax/10.);
        box_baselineInt -> SetFillColor(kRed);
        box_baselineInt -> SetFillStyle(3001);
        box_baselineInt -> Draw("same");
        
        TLine* line_signal = new TLine(sWinMin,0.,sWinMax,0.);
        line_signal -> SetLineColor(kGreen+2);
        line_signal -> SetLineWidth(2);
        line_signal -> Draw("same");
        
        TBox* box_signalInt = new TBox(std::max(float(0.),float(maxSample-sIntWinMin)),0.-fitAmpMax/10.,std::min(float(1023.),float(maxSample+sIntWinMax)),0+fitAmpMax/10.);
        box_signalInt -> SetFillColor(kGreen);
        box_signalInt -> SetFillStyle(3001);
        box_signalInt -> Draw("same");
        
        if( funcAmp )
        {
          funcAmp->SetLineColor(kRed);
          funcAmp->Draw("same");
          for(int param = 0; param < funcAmp->GetNpar(); ++param)
            std::cout << "[" << param << "] = " << funcAmp->GetParameter(param) << " +/- " << funcAmp->GetParError(param) << std::endl;
          std::cout << "chi2/NDF: " << funcAmp->GetChisquare()/funcAmp->GetNDF() << std::endl;
          TLine* line_maxTime = new TLine(fitTimeMax/tUnit,0.,fitTimeMax/tUnit,0.+fitAmpMax);
          line_maxTime -> SetLineColor(kRed);
          line_maxTime -> SetLineStyle(2);
          line_maxTime -> Draw("same");
        }
        if( leSlope > 0. )
        {
          // funcTimeLE -> SetLineColor(kBlue);
          // funcTimeLE -> Draw("same");
          
          TLine* line_thr = new TLine(leTime/tUnit-10.,leThr,leTime/tUnit+10.,leThr);
          line_thr -> SetLineColor(kBlue);
          line_thr -> SetLineStyle(2);
          line_thr -> Draw("same");
          
          TLine* line_leTime = new TLine(leTime/tUnit,0.,leTime/tUnit,0.+fitAmpMax);
          line_leTime -> SetLineColor(kBlue);
          line_leTime -> SetLineStyle(2);
          line_leTime -> Draw("same");
          
          TLatex* latexLabel = new TLatex(0.50,0.80,Form("#splitline{LED}{#splitline{#sigma_{V}: %.1f ADC   dV/dt = %.1f ADC/ns}{#sigma_{V} / (dV/dt) = %.0f ps}}",
                                                         baselineRMS,leSlope,baselineRMS/leSlope*1000.));
          latexLabel -> SetNDC();
          latexLabel -> SetTextFont(42);
          latexLabel -> SetTextSize(0.03);
          latexLabel -> SetTextColor(kBlue);
          latexLabel -> Draw("same");
        }
        if( teSlope > 0. )
        {
          // funcTimeTE -> SetLineColor(kBlue);
          // funcTimeTE -> Draw("same");
          
          TLine* line_thr = new TLine(teTime/tUnit-10.,teThr,teTime/tUnit+10.,teThr);
          line_thr -> SetLineColor(kBlue);
          line_thr -> SetLineStyle(2);
          line_thr -> Draw("same");
          
          TLine* line_teTime = new TLine(teTime/tUnit,0.,teTime/tUnit,0.+fitAmpMax);
          line_teTime -> SetLineColor(kBlue);
          line_teTime -> SetLineStyle(2);
          line_teTime -> Draw("same");
          
          TLatex* latexLabel = new TLatex(0.50,0.60,Form("#splitline{TED}{#splitline{#sigma_{V}: %.1f ADC   dV/dt = %.1f ADC/ns}{#sigma_{V} / (dV/dt) = %.0f ps}}",
                                                         baselineRMS,teSlope,baselineRMS/teSlope*1000.));
          latexLabel -> SetNDC();
          latexLabel -> SetTextFont(42);
          latexLabel -> SetTextSize(0.03);
          latexLabel -> SetTextColor(kBlue);
          latexLabel -> Draw("same");
        }
        if( cfSlope > 0. )
        {
          // funcTimeCF -> SetLineColor(kTeal);
          // funcTimeCF -> Draw("same");
          
          TLine* line_cfFrac = new TLine(cfTime/tUnit-10.,0.+cfFrac*fitAmpMax,cfTime/tUnit+10.,0.+cfFrac*fitAmpMax);
          line_cfFrac -> SetLineColor(kTeal);
          line_cfFrac -> SetLineStyle(2);
          line_cfFrac -> Draw("same");
          
          TLine* line_cfTime = new TLine(cfTime/tUnit,0.,cfTime/tUnit,0.+fitAmpMax);
          line_cfTime -> SetLineColor(kTeal);
          line_cfTime -> SetLineStyle(2);
          line_cfTime -> Draw("same");
          
          TLatex* latexLabel = new TLatex(0.50,0.40,Form("#splitline{CFD}{#splitline{#sigma_{V}: %.1f ADC   dV/dt = %.1f ADC/ns}{#sigma_{V} / (dV/dt) = %.0f ps}}",
                                                         baselineRMS,cfSlope,baselineRMS/cfSlope*1000.));
          latexLabel -> SetNDC();
          latexLabel -> SetTextFont(42);
          latexLabel -> SetTextSize(0.03);
          latexLabel -> SetTextColor(kTeal);
          latexLabel -> Draw("same");
        }
        
        gPad -> Update();
      } 
    }
    
    
    if( atoi(event.c_str()) == -1 )
    {
      vector<string> channels = opts.GetOpt<vector<string> >("DigiReco.channelsNames");
      std::map<std::string,TCanvas*> c1;
      for(auto& channel : channels)
      {
        c1[channel] = new TCanvas(Form("%s",channel.c_str()),Form("%s",channel.c_str()));
        TH1F* hPad = (TH1F*)( gPad->DrawFrame(-10.,0.,210.,1000.) );
        hPad -> SetTitle(";time [ns]; amplitude [mV]");
        hPad -> Draw();
      }
      
      
      while( h4Tree.NextEntry() && (index-stoul(run)*1e9<maxEvents || maxEvents==-1) )
      {
        if( index%1 == 0 ) std::cout << ">>> processing event " << index-stoul(run)*1e9 << " / " << h4Tree.GetEntries() << "\r" << std::flush;
        
        if( spill != -1 && int(h4Tree.spillNumber) != spill ) continue;
        
        //---call ProcessEvent for each plugin and check the return status
        bool status=true;
        for(auto& plugin : pluginSequence)
        {
          if(status)
          {
            status = plugin->ProcessEvent(h4Tree, pluginMap, opts);
          }
        }
        
        //---fill the main tree with info variables and increase event counter
        if(status)
        {
          mainTree.time_stamp = h4Tree.evtTimeStart;
          mainTree.run = h4Tree.runNumber;
          mainTree.spill = h4Tree.spillNumber;
          mainTree.event = h4Tree.evtNumber;
          mainTree.Fill();
          ++index;
        }
        
        //---draw waveforms
        vector<string> channels = opts.GetOpt<vector<string> >("DigiReco.channelsNames");
        for(auto& channel : channels)
        {
          auto shared_data = pluginMap["DigiReco"]->GetSharedData(std::string("DigiReco")+"_"+channel, "", false);
          WFClass* WF = (WFClass*)shared_data.at(0).obj;
          
          auto analizedWF = WF->GetSamples();
          float tUnit = WF->GetTUnit();
          
          TGraphErrors* g = new TGraphErrors();
          for(unsigned int jSample=0; jSample<analizedWF->size(); ++jSample)
          {
            g -> SetPoint(jSample,jSample*tUnit,analizedWF->at(jSample)*0.25);
          }
          
          c1[channel] -> cd();
          g -> Draw("L,same");
          
          // gPad -> Update();
        } 
      }
    }
    
    //---end
    // for(auto& plugin : pluginSequence)
    // {
    //     //---call endjob for each plugin        
    //     // bool r_status = plugin->End(opts);
    //     // if(!r_status)
    //     // {
    //     //     cout << ">>> ERROR: plugin returned bad flag from End() call: " << plugin->GetInstanceName() << endl;
    //     //     exit(-1);
    //     // }

    //     //---get permanent data from each plugin and store them in the out file
    //     for(auto& shared : plugin->GetSharedData())
    //     {
    //         if(shared.obj->IsA()->GetName() == string("TTree"))
    //         {
    //             TTree* currentTree = (TTree*)shared.obj;
    //             outROOT->cd();
    //             currentTree->BuildIndex("index");
    //             currentTree->Write(currentTree->GetName(), TObject::kOverwrite);
    //             mainTree.AddFriend(currentTree->GetName());
    //         }
    //         else
    //         {
    //             outROOT->cd();
    //             shared.obj->Write(shared.tag.c_str(), TObject::kOverwrite);
    //         }
    //     }
    // }
    
    // for(auto& loader : pluginLoaders)
    //   loader->Destroy();
    
    if( popupPlots ) { std::cout << ">>> popping up the plots" << std::endl; theApp -> Run(); }
    
    exit(0);
}    

#endif
