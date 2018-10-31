#ifndef __H4_RECO__
#define __H4_RECO__

#include <unistd.h>
#include <csignal>
#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TChain.h"

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include "interface/PluginLoader.h"
#include "interface/PluginBase.h"
#include "interface/DataLoader.h"
#include "interface/RecoTree.h"

#include <netdb.h>
#include <unistd.h>

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

//----------Simple function to track memory and CPU usage---------------------------------
void TrackProcess(float* cpu, float* mem, float* vsz, float* rss)
{
    string dummy1, dummy2, dummy3, dummy4, dummy5;
    string time;
    //---get cpu/mem info
    int pid = getpid();
    string ps_command = "ps up "+to_string(pid)+" >.H4Reco_proc.tmp";
    system(ps_command.c_str());
    ifstream proc_tmp(".H4Reco_proc.tmp", ios::in);
    getline(proc_tmp, dummy1);
    proc_tmp >> dummy1 >> dummy2 >> cpu[0] >> mem[0] >> vsz[0] >> rss[0]
             >> dummy3 >> dummy4 >> dummy5 >> time;
    vsz[0] = vsz[0]/1000;
    rss[0] = rss[0]/1000;
    proc_tmp.close();
    if(cpu[0]>cpu[1])
        cpu[1] = cpu[0];
    if(mem[0]>mem[1])
        mem[1] = mem[0];
    if(vsz[0]>vsz[1])
        vsz[1] = vsz[0];
    if(rss[0]>rss[1])
        rss[1] = rss[0];
    //---print statistics
    cout << "-----Machine stats---current/max-----" << endl
         << "CPU(%): " << cpu[0] << "/" << cpu[1] << endl
         << "MEM(%): " << mem[0] << "/" << mem[1] << endl
         << "VSZ(M): " << vsz[0] << "/" << vsz[1] << endl
         << "RSS(M): " << rss[0] << "/" << rss[1] << endl
         << "time lasted: " << time << endl;
}
                  
//**********MAIN**************************************************************************
int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        cout << argv[0] << " cfg file " << "[run] " << "[first spill] " << "[number of spills] " << endl;
        return -1;
    }

    //---memory consumption tracking---
    float cpu[2]{0}, mem[2]={0}, vsz[2]={0}, rss[2]={0};

    //---load options---    
    CfgManager opts;
    opts.ParseConfigFile(argv[1]);

    //-----input setup-----    
    int spill=-1;
    int nspills=-1;
    if(argc > 2)
    {
        vector<string> run(1, argv[2]);
        opts.SetOpt("h4reco.run", run);
    }
    if(argc > 3) {
        spill = atoi(argv[3]);
        vector<string> files(1, "1");
        opts.SetOpt("h4reco.maxFiles", files);
    }
    if(argc > 4)
    {
        nspills = atoi(argv[4]);
        vector<string> files(1, argv[4]);
        opts.SetOpt("h4reco.maxFiles", files);
    }
    string outSuffix = opts.GetOpt<string>("h4reco.outNameSuffix");
    bool storeTree = opts.OptExist("h4reco.storeTree") ?
        opts.GetOpt<bool>("h4reco.storeTree") : true;

    string run = opts.GetOpt<string>("h4reco.run");
    int totLoops= opts.OptExist("h4reco.totLoops") ? opts.GetOpt<int>("h4reco.totLoops") : 1;

    //-----Load raw data-----
    vector<string> spillOpt(1, to_string(spill));
    opts.SetOpt("h4reco.firstSpill", spillOpt);
    DataLoader dataLoader(opts);

    //-----output setup-----
    uint64 index=0;    
    TFile* outROOT;
    if(spill == -1)
        outROOT = new TFile(outSuffix+"run"+TString(run)+".root", "RECREATE");
    else if(nspills == -1)
        outROOT = new TFile(outSuffix+"run"+TString(run)+"_spill"+spillOpt.back()+".root", "RECREATE");
    else {
        // TODO: generate correct file name if less than the requested number of spills are available
        outROOT = new TFile(outSuffix+"run"+TString(run)+"_spills"+spillOpt.back()+"-"+to_string(spill+nspills-1)+".root", "RECREATE");
    }
    outROOT->cd();
    TDirectory* outDIR=outROOT;

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
    }
            

    int iLoop=0;
    int maxEvents = opts.OptExist("h4reco.maxEvents") ? opts.GetOpt<int>("h4reco.maxEvents") : -1;
    bool isSim = opts.OptExist("h4reco.generateEvents") ? opts.GetOpt<bool>("h4reco.generateEvents") : false;

    if(isSim)
        cout << ">>> Processing H4DAQ simulation <<<" << endl;
    else
        cout << ">>> Processing H4DAQ run #" << run << " <<<" << endl;

    while(iLoop<totLoops)
    {
	int nEvents = 0;

	if (totLoops>1)
        {
	    cout << ">>> Starting iteration #" << iLoop << " <<<" << endl;
	    outDIR=outROOT->mkdir(Form("Iter_%d",iLoop));
        }

	//---begin loop
	mainTree.tree_->Reset();
	for(auto& plugin : pluginSequence)
        {
	    //---call Begin() methods and check the return status
	    bool r_status = plugin->BeginLoop(iLoop, opts);
	    if(!r_status)
            {
		cout << ">>> ERROR: plugin returned bad flag from BeginLoop() call: " << plugin->GetInstanceName() << endl;
		exit(-1);
            }
	    
	    //---Get plugin shared data
	    for(auto& shared : plugin->GetSharedData("", "TTree", true))
            {
	    	TTree* tree = (TTree*)shared.obj;
		tree->Reset();
	    	tree->SetDirectory(outDIR);
            }
        }


	//---events loop
	while((dataLoader.NextEvent() && (nEvents < maxEvents || maxEvents == -1)) || (isSim && (nEvents < maxEvents)))
        {
	    if(dataLoader.FirstEventInSpill())
            {
                cout << ">>> Processed spills: " << dataLoader.GetNFilesProcessed() << "/" << dataLoader.GetNFiles() << endl;
                cout << ">>> Processed events: " << nEvents << endl;
		TrackProcess(cpu, mem, vsz, rss);
            }
        
	    //---set index value run*1e10+spill*1e4+event
	    index = dataLoader.GetTree().runNumber*1e9 + dataLoader.GetTree().spillNumber*1e5 + dataLoader.GetTree().evtNumber;
	    
	    //---call ProcessEvent for each plugin and check the return status
	    bool status=true;
	    for(auto& plugin : pluginSequence)
            {
	    	status &= plugin->Clear(); 
	    	status &= plugin->ProcessEvent(dataLoader.GetTree(), pluginMap, opts);
            }
        
	    //---Fill the main tree with info variables and increase event counter
            mainTree.time_stamps.clear();
            for(int iT=0; iT<dataLoader.GetTree().nEvtTimes; ++iT)
                mainTree.time_stamps.push_back(dataLoader.GetTree().evtTime[iT]);
	    mainTree.evt_flag = status;
	    mainTree.run = dataLoader.GetTree().runNumber;
	    mainTree.spill = dataLoader.GetTree().spillNumber;
	    mainTree.event = dataLoader.GetTree().evtNumber;
	    mainTree.Fill();
	    ++nEvents;
        }

	//---end
	for(auto& plugin : pluginSequence)
        {
	    //---call endjob for each plugin        
	    bool r_status = plugin->EndLoop(iLoop,opts);
	    if(!r_status)
            {
		cout << ">>> ERROR: plugin returned bad flag from EndLoop() call: " << plugin->GetInstanceName() << endl;
		exit(-1);
            }

	    //---get permanent data from each plugin and store them in the out file
	    for(auto& shared : plugin->GetSharedData())
            {
	    	if(shared.obj->IsA()->GetName() == string("TTree") && storeTree)
                {
	    	    TTree* currentTree = (TTree*)shared.obj;
	    	    outDIR->cd();
	    	    currentTree->BuildIndex("index");
	    	    currentTree->Write(currentTree->GetName(), TObject::kOverwrite);
	    	    mainTree.AddFriend(currentTree->GetName());
                }
	    	else
                {
	    	    outDIR->cd();
	    	    shared.obj->Write(shared.tag.c_str(), TObject::kOverwrite);
                }
            }
        }

	outDIR->cd();
	if(storeTree)
            mainTree.Write();
	++iLoop;
    }

    //---end
    for(auto& plugin : pluginSequence)
    {
        //---call endjob for each plugin        
        bool r_status = plugin->End(opts);
        if(!r_status)
        {
            cout << ">>> ERROR: plugin returned bad flag from End() call: " << plugin->GetInstanceName() << endl;
            exit(-1);
        }
    }
    
    //---close
    outROOT->cd();
    opts.Write("cfg");
    outROOT->Close();
    for(auto& loader : pluginLoaders)
        loader->Destroy();

    //---info
    TrackProcess(cpu, mem, vsz, rss);

    exit(0);
}    

#endif
