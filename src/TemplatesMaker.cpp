#ifndef __MAIN__
#define __MAIN__

#include <string>

#include "TString.h"
#include "TFile.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TROOT.h"

#include "interface/CfgManager.h"
#include "interface/CfgManagerT.h"
#include "interface/H4Tree.h"
#include "interface/WFClass.h"

using namespace std;

//----------DFT---------------------------------------------------------------------------
pair<TH1F, TH1F> DFT_cut(TProfile* inWave, string name, float fCut=800)
{
    float *Re, *Im, *mod;
    int nbins = inWave->GetNbinsX();
    Re = (float*)malloc(sizeof(float)*nbins);
    Im = (float*)malloc(sizeof(float)*nbins);
    mod = (float*)malloc(sizeof(float)*nbins);
    string nameDFT = "DFT_"+name;
    name = "cleaned_"+name;
    TH1F outDFT(nameDFT.c_str(), nameDFT.c_str(), nbins/2, 0, nbins/2);
    TH1F outWave(name.c_str(), name.c_str(), nbins, -40, 160); 
    for(int k=0; k<nbins/2; k++) 
    { 
        Re[k]=0;
        Im[k]=0; 
	mod[k]=0;
        for(int t=1; t<=nbins; t++) 
        { 
            float angle = 2 * TMath::Pi() * t * k / (nbins); 
            Re[k] += 2*inWave->GetBinContent(t) * cos(angle)/(nbins); 
            Im[k] += -2*inWave->GetBinContent(t) * sin(angle)/(nbins);
        }	    
	if(k == 0 || k == nbins-1)
	{
	    Re[k] = Re[k]/2;
	    Im[k] = Im[k]/2;
	}
	mod[k] = TMath::Sqrt(TMath::Abs(pow(Re[k],2) - pow(Im[k],2)));
        outDFT.SetBinContent(k, mod[k]);
    }
    for(int t=1; t<=nbins; t++) 
    { 
        float recoWave=0;
        for(int k=0; k<nbins/2; k++) 
        { 
	    if(k > fCut)
	    	continue;
            float angle = 2 * TMath::Pi() * t * k / nbins; 
	    recoWave += Re[k]*cos(angle) - Im[k]*sin(angle);
	}  
	outWave.SetBinContent(t, recoWave);
    }
    outWave.Scale(inWave->GetMaximum()/outWave.GetMaximum());

    return make_pair(outWave, outDFT);
}

TH1F* getMeanProfile(TH2F* waveForm)
{
  TH1F* prof=new TH1F(TString(waveForm->GetName())+"_prof",TString(waveForm->GetName())+"_prof",waveForm->GetNbinsX(),waveForm->GetXaxis()->GetXmin(),waveForm->GetXaxis()->GetXmax());
  for (int ibin=1;ibin<waveForm->GetNbinsX()+1;++ibin)
    {
      TH1F* h=new TH1F(Form("%s_xbin_%d",waveForm->GetName(),ibin),Form("%s_xbin_%d",waveForm->GetName(),ibin),waveForm->GetNbinsY(),waveForm->GetYaxis()->GetXmin(),waveForm->GetYaxis()->GetXmax());
      for (int ybin=1;ybin<waveForm->GetNbinsY()+1;++ybin)
	{
	  h->SetBinContent(ybin,waveForm->GetBinContent(ibin,ybin));
	  h->SetBinError(ybin,waveForm->GetBinError(ibin,ybin));
	}
      float deltaMean=9999;
      float deltaRMS=9999;
      float oldMean=h->GetMean();
      float oldRMS=h->GetRMS();
      while (deltaMean>1E-4 || deltaRMS>1E-5)
	{
	  h->GetXaxis()->SetRangeUser(oldMean-3*oldRMS,oldMean+3*oldRMS);
	  float newMean=h->GetMean();
	  float newRMS=h->GetRMS();
	  deltaMean=abs(newMean-oldMean);
	  deltaRMS=abs(newRMS-oldRMS);
	  oldMean=newMean;
	  oldRMS=newRMS;
	}
      prof->SetBinContent(ibin,h->GetMean());
      prof->SetBinError(ibin,h->GetMeanError());

      delete h;
    }
  return prof;
}
//**********MAIN**************************************************************************
int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        cout << argv[0] << " cfg file " << "[run]" << endl; 
        return -1;
    }

    gROOT->SetBatch(kTRUE);
    //---load options---
    CfgManager opts;
    opts.ParseConfigFile(argv[1]);

    //---data opts
    string path=opts.GetOpt<string>("global.path2data");
    string run=opts.GetOpt<string>("global.run");
    string outSuffix=opts.GetOpt<string>("global.outFileSuffix");
    string tag=opts.GetOpt<string>("global.tag");

    if(argc > 2)
        run = argv[2];
    int maxEvents = opts.GetOpt<int>("global.maxEvents");
    //---channels opts
    int nCh = opts.GetOpt<int>("global.nCh");
    int nSamples = opts.GetOpt<int>("global.nSamples");
    float tUnit = opts.GetOpt<float>("global.tUnit");
    string refChannel = opts.GetOpt<string>("global.refChannel");
    vector<string> channelsNames = opts.GetOpt<vector<string>& >("global.channelsNames");
    map<string, vector<float> > timeOpts;
    timeOpts[refChannel] = opts.GetOpt<vector<float> >(refChannel+".timeOpts");
    for(auto& channel : channelsNames)
         timeOpts[channel] = opts.GetOpt<vector<float> >(channel+".timeOpts");

    //---definitions---
    int iEvent=0;
  
    //-----output setup-----
    TString outF="ntuples/Templates_"+TString(outSuffix)+"_"+TString(run)+".root";
    TFile* outROOT = TFile::Open(outF, "RECREATE");
    outROOT->cd();
    map<string, TH2F*> templates;
    for(auto channel : channelsNames)
      templates[channel] = new TH2F(Form("%s_%s",channel.c_str(),tag.c_str()), Form("%s_%s",channel.c_str(),tag.c_str()),
				      16000, -40, 160, 1200, -0.1, 1.1);
  
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

    //---process WFs---
    long nentries=h4Tree.GetEntries();
    cout << ">>> Processing H4DAQ run #" << run << " events #" << nentries << " <<<" << endl;
    while(h4Tree.NextEvt() && (iEvent<maxEvents || maxEvents==-1))
    {        
        ++iEvent;
	if (iEvent%1000==0) std::cout << "Processing event " << iEvent << "/" << nentries << std::endl;
        //---setup output event 
        int outCh=0;
        bool badEvent=false;
        
        //---read the digitizer
        //---read time raference channel
        float refTime=0, refAmpl=0;
        WFClass WF(opts.GetOpt<int>(refChannel+".polarity"), tUnit);
        int digiGr = opts.GetOpt<int>(refChannel+".digiGroup");
        int digiCh = opts.GetOpt<int>(refChannel+".digiChannel");
        int offset = digiGr*9*nSamples + digiCh*nSamples;
        for(int iSample=offset; iSample<offset+nSamples; ++iSample)
        {
            //---H4DAQ bug: sometimes ADC value is out of bound.
            if(h4Tree.digiSampleValue[iSample] > 4096)
            {
                badEvent = true;
                break;
            }
            WF.AddSample(h4Tree.digiSampleValue[iSample]);
        }
        //---skip bad events
        if(badEvent)
            continue;
        WF.SetBaselineWindow(opts.GetOpt<int>(refChannel+".baselineWin", 0), 
                             opts.GetOpt<int>(refChannel+".baselineWin", 1));
        WF.SetSignalWindow(opts.GetOpt<int>(refChannel+".signalWin", 0), 
                           opts.GetOpt<int>(refChannel+".signalWin", 1));
        WFBaseline refBaseline=WF.SubtractBaseline();
        refAmpl = WF.GetInterpolatedAmpMax();            
        refTime = WF.GetTime(opts.GetOpt<string>(refChannel+".timeType"), timeOpts[refChannel]).first;
        //---require reference channel to be good
        if(refTime/tUnit < opts.GetOpt<int>(refChannel+".signalWin", 0) ||
           refTime/tUnit > opts.GetOpt<int>(refChannel+".signalWin", 1) ||
	   refBaseline.rms > opts.GetOpt<float>(refChannel+".noiseThreshold") ||
	   refAmpl < opts.GetOpt<int>(refChannel+".amplitudeThreshold"))
            continue;

        //---template channels
        for(auto& channel : channelsNames)
        {
            //---read WFs
            WFClass WF(opts.GetOpt<int>(channel+".polarity"), tUnit);
            int digiGr = opts.GetOpt<int>(channel+".digiGroup");
            int digiCh = opts.GetOpt<int>(channel+".digiChannel");
            int offset = digiGr*9*nSamples + digiCh*nSamples;
            for(int iSample=offset; iSample<offset+nSamples; ++iSample)
            {
                //---H4DAQ bug: sometimes ADC value is out of bound.
                if(h4Tree.digiSampleValue[iSample] > 4096)
                {
                    badEvent = true;
                    break;
                }
                WF.AddSample(h4Tree.digiSampleValue[iSample]);
            }
            //---skip bad events
            if(badEvent)
                continue;
            //---compute reco variables
            float channelTime=0, channelAmpl=0;            
            WF.SetBaselineWindow(opts.GetOpt<int>(channel+".baselineWin", 0), 
                                 opts.GetOpt<int>(channel+".baselineWin", 1));
            WF.SetSignalWindow(opts.GetOpt<int>(channel+".signalWin", 0), 
                               opts.GetOpt<int>(channel+".signalWin", 1));
            WFBaseline channelBaseline=WF.SubtractBaseline();
	    channelAmpl = WF.GetInterpolatedAmpMax(-1,-1,opts.GetOpt<int>(channel+".signalWin", 2));
            channelTime = WF.GetTime(opts.GetOpt<string>(channel+".timeType"), timeOpts[channel]).first;
            //---skip bad events or events with no signal
            if(channelTime/tUnit > opts.GetOpt<int>(channel+".signalWin", 0) &&
               channelTime/tUnit < opts.GetOpt<int>(channel+".signalWin", 1) &&
	       channelBaseline.rms < opts.GetOpt<float>(channel+".noiseThreshold") &&
               channelAmpl > opts.GetOpt<int>(channel+".amplitudeThreshold") &&
               channelAmpl < 4000)
            {                
                vector<float>* analizedWF = WF.GetSamples();
                for(int iSample=0; iSample<analizedWF->size(); ++iSample)
		  templates[channel]->Fill(iSample*tUnit-refTime, analizedWF->at(iSample)/channelAmpl);
            }
        }
    }   


    cout << ">>> Writing output " << outF << " <<<" << endl;    

    outROOT->cd();

    for(auto& channel : channelsNames)
    {
        // pair<TH1F, TH1F> dft = DFT_cut(templates[channel], channel, 800);
        // dft.first.Write();
        // dft.second.Write();
      if (templates[channel]->GetEntries()>1024*100)
	{
	  TH1F* prof=getMeanProfile(templates[channel]);
	  templates[channel]->Write();
	  prof->Write();
	}
    }
    opts.Write("cfg");
    outROOT->Close();
}

#endif
