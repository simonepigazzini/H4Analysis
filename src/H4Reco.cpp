#ifndef __MAIN__
#define __MAIN__

#include <string>
#include <sstream>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include "interface/CfgManager.h"
#include "interface/CfgManagerT.h"
#include "interface/H4Tree.h"
#include "interface/RecoTree.h"
#include "interface/WFTree.h"
#include "interface/InfoTree.h"
#include "interface/WFClass.h"
#include "interface/WFClassNINO.h"
#include "interface/WFViewer.h"
#include "interface/WireChamber.h"
#include "interface/HodoUtils.h"

using namespace std;

//**********MAIN**************************************************************************

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        cout << argv[0] << " cfg file " << "[run]" << endl; 
        return -1;
    }

    //---load options---
    CfgManager opts;
    opts.ParseConfigFile(argv[1]);
    //cout << opts << endl;
    //---data opts
    string path=opts.GetOpt<string>("global.path2data");
    string run=opts.GetOpt<string>("global.run");

    if(argc > 2)
        run = argv[2];


    int maxEvents = opts.GetOpt<int>("global.maxEvents");

    //---setup run tag
    ifstream runList(opts.GetOpt<string>("global.runList").c_str(), ios::in);
    std::string line;
    std::map<int, std::string> runTags;
    while(std::getline(runList,line))
    {
      std::istringstream iss(line);
      int run;
      string tag;
      if (!(iss >> run >> tag)) { break; }
      runTags[run]=tag;
    }
    std::string tag = runTags[atoi(run.c_str())];
    std::cout << ">>> RUN# " << run << " TAG " << tag << " <<<" <<std::endl;

    //---channels setup
    int nSamples = opts.GetOpt<int>("global.nSamples");
    float tUnit = opts.GetOpt<float>("global.tUnit");
    bool useTrigRef = opts.GetOpt<int>("global.useTrigRef");
    bool fillToys = opts.GetOpt<int>("global.fillToys");
    bool FFTCleanWF = opts.GetOpt<int>("global.FFTCleanWF");

    vector<string> channelsNames = opts.GetOpt<vector<string>& >("global.channelsNames");
    map<string, WFClass*> WFs;
    map<string, WFViewer> WFViewers;
    map<string, WFClass*> emulatedWFs;
    map<string, WFClass*> cleanedWFs;
    map<string, TH1*> templates;
    map<string, vector<float> > timeOpts;
    for(auto& channel : channelsNames)
    {
        if(opts.GetOpt<bool>(channel+".type") && opts.GetOpt<string>(channel+".type") == "NINO")
        {
            WFs[channel] = new WFClassNINO(opts.GetOpt<int>(channel+".polarity"), tUnit);
            emulatedWFs[channel] = new WFClassNINO(1, tUnit);
            cleanedWFs[channel] = new WFClassNINO(1, tUnit);
        }
        else
        {
            WFs[channel] = new WFClass(opts.GetOpt<int>(channel+".polarity"), tUnit);
            emulatedWFs[channel] = new WFClass(1, tUnit);
            cleanedWFs[channel] = new WFClass(1, tUnit);
        }
        timeOpts[channel] = opts.GetOpt<vector<float> >(channel+".timeOpts");
        if(opts.GetOpt<bool>(channel+".templateFit"))
        {
            TFile* templateFile = TFile::Open(opts.GetOpt<string>(channel+".templateFit.file", 0).c_str(), ".READ");
            TH1* wfTemplate=(TH1*)templateFile->Get(Form("%s_%s_%s",channel.c_str(),tag.c_str(),opts.GetOpt<string>(channel+".templateFit.file", 1).c_str()));
            WFs[channel]->SetTemplate(wfTemplate);
	    emulatedWFs[channel]->SetTemplate(wfTemplate);
	    cleanedWFs[channel]->SetTemplate(wfTemplate);
            WFViewers[channel]=WFViewer(channel, (TH1F*)wfTemplate);
            templateFile->Close();
        }
    }
  
    //-----output setup-----
    int iEvent=-1;
    TFile* outROOT = TFile::Open("ntuples/"+TString(run)+".root", "RECREATE");
    InfoTree outInfo(channelsNames, nSamples, opts.GetOpt<bool>("global.H4hodo"),
                     opts.GetOpt<int>("global.nWireChambers"), &iEvent);
    RecoTree outTree(channelsNames.size(), &iEvent);
    RecoTree outToysTree(channelsNames.size(), &iEvent);
    RecoTree outCleanedRecoTree(channelsNames.size(), &iEvent);
    WFTree outWFTree(channelsNames.size(), nSamples, &iEvent);
    WFTree outCleanedWFTree(channelsNames.size(), nSamples, &iEvent);
  
    //-----input setup-----
    TChain* inTree = new TChain("H4tree");
    string ls_command;
    if(path.find("/eos/cms") != string::npos)
        ls_command = string("gfal-ls root://eoscms/"+path+run+" | grep 'root' > tmp/"+run+".list");
    else
        ls_command = string("ls "+path+run+" | grep 'root' > tmp/"+run+".list");
    system(ls_command.c_str());
    ifstream waveList(string("tmp/"+run+".list").c_str(), ios::in);
    string file;

    int nFiles=0;
    while(waveList >> file && (opts.GetOpt<int>("global.maxFiles")<0 || nFiles<opts.GetOpt<int>("global.maxFiles")) )
    {
        if(path.find("/eos/cms") != string::npos)
        {
	    std::cout << "+++ Adding file " << ("root://eoscms/"+path+run+"/"+file).c_str() << std::endl;
            inTree->AddFile(("root://eoscms/"+path+run+"/"+file).c_str());
        }
        else
        {
	    std::cout << "+++ Adding file " << (path+run+"/"+file).c_str() << std::endl;
            inTree->AddFile((path+run+"/"+file).c_str());
        }
	++nFiles;
    }

    H4Tree h4Tree(inTree);
    WireChamber wires;
    if(opts.GetOpt<int>("global.nWireChambers") > 0)
    {
        wires = WireChamber(&h4Tree,
                            opts.GetOpt<int>("WireChamber.chXleft"),
                            opts.GetOpt<int>("WireChamber.chXright"),
                            opts.GetOpt<int>("WireChamber.chYup"),
                            opts.GetOpt<int>("WireChamber.chYdown"));
    }


    //---setup H4hodo
    vector<int> hodoFiberOrderA;
    vector<int> hodoFiberOrderB;
    if(opts.GetOpt<bool>("global.H4hodo"))
        FillHodoFiberOrder(hodoFiberOrderA,hodoFiberOrderB);

    //---process WFs---
    cout << ">>> Processing H4DAQ run #" << run << " <<<" << endl;
    while(h4Tree.NextEvt() && (iEvent<maxEvents || maxEvents==-1))
    {
        //---setup output event 
        int outCh=0;
        bool badEvent=false;
        bool fillWFtree=false;
        if(opts.GetOpt<int>("global.fillWFtree"))
            fillWFtree = iEvent % opts.GetOpt<int>("global.WFtreePrescale") == 0;

        //---read the wire chamber
        if(opts.GetOpt<int>("global.nWireChambers")>0)
        {
            wires.Unpack();
            outInfo.wireX[0] = wires.GetX();
            outInfo.wireY[0] = wires.GetY();
        }
        
        //---read the hodoscopes
        if(opts.GetOpt<bool>("global.H4hodo"))
        {
            map<int,std::map<int,bool> > hodoFiberOn;
            map<int,float> beamXY_hodo12;
            FillHodo(h4Tree, hodoFiberOn, hodoFiberOrderA, hodoFiberOrderB, beamXY_hodo12);
            outInfo.hodoX1 = beamXY_hodo12[0];
            outInfo.hodoY1 = beamXY_hodo12[1];
            outInfo.hodoX2 = beamXY_hodo12[2];
            outInfo.hodoY2 = beamXY_hodo12[3]; 
        }
        
        //---read the digitizer
        //---set time reference from digitized trigger
        int timeRef=0;
        if(useTrigRef)
        {
            for(int iSample=nSamples*8; iSample<nSamples*9; ++iSample)
            {
                if(h4Tree.digiSampleValue[iSample] < 1000)
                {
                    timeRef = iSample-nSamples*8;
                    break;
                }
            }
        }
        //---user channels
        for(auto& channel : channelsNames)
        {
            //---reset and read new WFs
            WFs[channel]->Reset();
            int digiGr = opts.GetOpt<int>(channel+".digiGroup");
            int digiCh = opts.GetOpt<int>(channel+".digiChannel");
            int offset = h4Tree.digiMap[make_pair(digiGr, digiCh)];
            for(int iSample=offset; iSample<offset+nSamples; ++iSample)
            {
                //---H4DAQ bug: sometimes ADC value is out of bound.
                if(h4Tree.digiSampleValue[iSample] > 4096)
                {
                    badEvent = true;
                    break;
                }
                WFs[channel]->AddSample(h4Tree.digiSampleValue[iSample]);
            }
            //---skip everything if one channel has trigger the badEvent flag 
            if(badEvent)
                break;
            //---compute reco variables
            WFs[channel]->SetBaselineWindow(opts.GetOpt<int>(channel+".baselineWin", 0), 
                                            opts.GetOpt<int>(channel+".baselineWin", 1));
            WFs[channel]->SetSignalWindow(opts.GetOpt<int>(channel+".signalWin", 0)+timeRef, 
                                          opts.GetOpt<int>(channel+".signalWin", 1)+timeRef);
            WFBaseline baselineInfo = WFs[channel]->SubtractBaseline();
            pair<float, float> timeInfo = WFs[channel]->GetTime(opts.GetOpt<string>(channel+".timeType"), timeOpts[channel]);
            outTree.b_charge[outCh] = WFs[channel]->GetIntegral(opts.GetOpt<int>(channel+".baselineInt", 0), 
                                                                opts.GetOpt<int>(channel+".baselineInt", 1));
            outTree.b_slope[outCh] = baselineInfo.slope;
            outTree.b_rms[outCh] = baselineInfo.rms;
            outTree.time[outCh] = timeInfo.first;
            outTree.time_chi2[outCh] = timeInfo.second;
            outTree.maximum[outCh] = WFs[channel]->GetAmpMax();
            outTree.amp_max[outCh] = WFs[channel]->GetInterpolatedAmpMax(-1,-1,opts.GetOpt<int>(channel+".signalWin", 2));
            outTree.charge_tot[outCh] = WFs[channel]->GetModIntegral(opts.GetOpt<int>(channel+".baselineInt", 1), 
                                                                     nSamples);
            outTree.charge_sig[outCh] = WFs[channel]->GetSignalIntegral(opts.GetOpt<int>(channel+".signalInt", 0), 
                                                                        opts.GetOpt<int>(channel+".signalInt", 1));
            //---template fit (only specified channels)
            WFFitResults fitResults{-1, -1000, -1};
            if(opts.GetOpt<bool>(channel+".templateFit"))
            {
                fitResults = WFs[channel]->TemplateFit(opts.GetOpt<float>(channel+".templateFit.fitWin", 0),
                                                       opts.GetOpt<int>(channel+".templateFit.fitWin", 1),
                                                       opts.GetOpt<int>(channel+".templateFit.fitWin", 2));
                outTree.fit_ampl[outCh] = fitResults.ampl;
                outTree.fit_time[outCh] = fitResults.time;
                outTree.fit_chi2[outCh] = fitResults.chi2;

                // --- Throw toys
                if (fillToys)
		{
                    emulatedWFs[channel]->Reset();
                    WFs[channel]->EmulatedWF(*(emulatedWFs[channel]), baselineInfo.rms, fitResults.ampl, fitResults.time);
                    emulatedWFs[channel]->SetBaselineWindow(opts.GetOpt<int>(channel+".baselineWin", 0), 
                                                            opts.GetOpt<int>(channel+".baselineWin", 1));
                    emulatedWFs[channel]->SetSignalWindow(opts.GetOpt<int>(channel+".signalWin", 0)+timeRef, 
                                                          opts.GetOpt<int>(channel+".signalWin", 1)+timeRef);
                    WFBaseline baselineInfo = emulatedWFs[channel]->SubtractBaseline();
                    double toy_maximum=emulatedWFs[channel]->GetAmpMax();
                    double toy_amp_max=emulatedWFs[channel]->GetInterpolatedAmpMax(-1,-1,opts.GetOpt<int>(channel+".signalWin", 2));	       
                    pair<float, float> toy_timeInfo = emulatedWFs[channel]->GetTime(opts.GetOpt<string>(channel+".timeType"), timeOpts[channel]);
                    outToysTree.time[outCh] = toy_timeInfo.first;
                    outToysTree.time_chi2[outCh] = toy_timeInfo.second;
                    outToysTree.maximum[outCh] = toy_maximum;
                    outToysTree.amp_max[outCh] = toy_amp_max;
		  
                    outToysTree.charge_tot[outCh] = emulatedWFs[channel]->GetModIntegral(opts.GetOpt<int>(channel+".baselineInt", 1), 
                                                                                         nSamples);
                    outToysTree.charge_sig[outCh] = emulatedWFs[channel]->GetSignalIntegral(opts.GetOpt<int>(channel+".signalInt", 0), 
                                                                                            opts.GetOpt<int>(channel+".signalInt", 1));
		  
                    WFFitResults toy_fitResults{-1, -1000, -1};
                    toy_fitResults = emulatedWFs[channel]->TemplateFit(opts.GetOpt<float>(channel+".templateFit.fitWin", 0), opts.GetOpt<int>(channel+".templateFit.fitWin", 1), opts.GetOpt<int>(channel+".templateFit.fitWin", 2));
		  
                    outToysTree.fit_ampl[outCh] = toy_fitResults.ampl;
                    outToysTree.fit_time[outCh] = toy_fitResults.time;
                    outToysTree.fit_chi2[outCh] = toy_fitResults.chi2;
		}
            }         
	    if (fillToys)
                outToysTree.time_stamp = h4Tree.evtTime;        

	    // ---- Cleaned (FFT) WF -----
	    if (FFTCleanWF)
            {
		cleanedWFs[channel]->Reset();
		WFs[channel]->FFT(*(cleanedWFs[channel]),opts.GetOpt<int>(channel+".FFTCuts", 0),
                                  opts.GetOpt<int>(channel+".FFTCuts", 1)); 
		cleanedWFs[channel]->SetBaselineWindow(opts.GetOpt<int>(channel+".baselineWin", 0), 
						       opts.GetOpt<int>(channel+".baselineWin", 1));
		cleanedWFs[channel]->SetSignalWindow(opts.GetOpt<int>(channel+".signalWin", 0)+timeRef, 
						     opts.GetOpt<int>(channel+".signalWin", 1)+timeRef);
		double cleaned_maximum=cleanedWFs[channel]->GetAmpMax();
		double cleaned_amp_max=cleanedWFs[channel]->GetInterpolatedAmpMax(-1,-1,opts.GetOpt<int>(channel+".signalWin", 2));
		pair<float, float> cleaned_timeInfo = cleanedWFs[channel]->GetTime(opts.GetOpt<string>(channel+".timeType"), timeOpts[channel]);
		outCleanedRecoTree.time[outCh] = cleaned_timeInfo.first;
		outCleanedRecoTree.time_chi2[outCh] = cleaned_timeInfo.second;
		outCleanedRecoTree.maximum[outCh] = cleaned_maximum;
		outCleanedRecoTree.amp_max[outCh] = cleaned_amp_max;
		
		outCleanedRecoTree.charge_tot[outCh] = cleanedWFs[channel]->GetModIntegral(opts.GetOpt<int>(channel+".baselineInt", 1), 
                                                                                           nSamples);
		outCleanedRecoTree.charge_sig[outCh] = cleanedWFs[channel]->GetSignalIntegral(opts.GetOpt<int>(channel+".signalInt", 0), 
                                                                                              opts.GetOpt<int>(channel+".signalInt", 1));
		
		WFFitResults cleaned_fitResults{-1, -1000, -1};
		if(opts.GetOpt<bool>(channel+".templateFit"))
                {
		    cleaned_fitResults = cleanedWFs[channel]->TemplateFit(opts.GetOpt<float>(channel+".templateFit.fitWin", 0), opts.GetOpt<int>(channel+".templateFit.fitWin", 1), opts.GetOpt<int>(channel+".templateFit.fitWin", 2));
		    outCleanedRecoTree.fit_ampl[outCh] = cleaned_fitResults.ampl;
		    outCleanedRecoTree.fit_time[outCh] = cleaned_fitResults.time;
		    outCleanedRecoTree.fit_chi2[outCh] = cleaned_fitResults.chi2;
                }
            }
  
	    //---WFs---
            if(fillWFtree)
            {
                vector<float>* analizedWF = WFs[channel]->GetSamples();
                for(int jSample=0; jSample<analizedWF->size(); ++jSample)
                {
                    outWFTree.time_stamp = h4Tree.evtTime;
                    outWFTree.WF_ch[jSample+outCh*nSamples] = outCh;
                    outWFTree.WF_time[jSample+outCh*nSamples] = jSample*tUnit;
                    outWFTree.WF_val[jSample+outCh*nSamples] = analizedWF->at(jSample);
                }
            }


            if(fillWFtree && FFTCleanWF)
            {
                vector<float>* analizedWF = cleanedWFs[channel]->GetSamples();
                for(int jSample=0; jSample<analizedWF->size(); ++jSample)
                {
                    outCleanedWFTree.time_stamp = h4Tree.evtTime;
                    outCleanedWFTree.WF_ch[jSample+outCh*nSamples] = outCh;
                    outCleanedWFTree.WF_time[jSample+outCh*nSamples] = jSample*tUnit;
                    outCleanedWFTree.WF_val[jSample+outCh*nSamples] = analizedWF->at(jSample);
                }
            }

            //---increase output tree channel counter
            ++outCh;
        }

        //---fill the output trees if the event is good
        if(!badEvent)
        {
            ++iEvent;
            if(iEvent % 1000 == 0)
                cout << ">>>Processed good events: " << iEvent << "/" << h4Tree.GetEntries() << endl;
            //---reco var
            outTree.time_stamp = h4Tree.evtTime;
            outTree.Fill();
            //---info var
            outInfo.time_stamp = h4Tree.evtTime;
            outInfo.run = h4Tree.runNumber;
            outInfo.spill = h4Tree.spillNumber;
            outInfo.event = h4Tree.evtNumber;
            outInfo.Fill();

	    //---Alternative WFs
	    if(fillToys)
                outToysTree.Fill();
	    if(FFTCleanWF)
                outCleanedRecoTree.Fill();

            //---WFs
            if(fillWFtree)
            {
	      outWFTree.Fill();
	      outWFTree.tree_->AutoSave("FlushBaskets");
		if (FFTCleanWF)
		  {
                    outCleanedWFTree.Fill();
                    outCleanedWFTree.tree_->AutoSave("FlushBaskets");
		  }
            }
        }
    }

    //---write output file 
    outROOT->cd();
    //---trees
    if(opts.GetOpt<int>("global.fillWFtree"))
    {
        outWFTree.Write("wf", "wf_tree");
        outTree.AddFriend("wf");
	if(FFTCleanWF)
        {
	    outCleanedWFTree.Write("wf_c", "cleaned_wf_tree");
	    outTree.AddFriend("wf_c");
        }
    }

    if(fillToys)
    {
        outToysTree.Write("toy", "toys_tree");
        outTree.AddFriend("toy");
    }

    if(FFTCleanWF)
    {
	outCleanedRecoTree.Write("reco_c", "cleaned_reco_tree");
	outTree.AddFriend("reco_c");
    }
    outInfo.Write();
    outTree.AddFriend("info_tree");
    outTree.Write();
    //---templates for re-fit
    for(auto& viewer : WFViewers)
        viewer.second.Write(viewer.first.c_str());
    //---cfg
    opts.Write("cfg");
    outROOT->Close();

    return 0;
}

#endif

