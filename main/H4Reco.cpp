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
#include "interface/RecoTree.h"
#include "interface/PluginBase.h"

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
                  
//----------Get input files---------------------------------------------------------------
void ReadInputFiles(CfgManager& opts, int& firstSpill, TChain* inTree)
{
    int nFiles=0;
    string ls_command;
    string file;
    string path=opts.GetOpt<string>("h4reco.path2data");
    string run=opts.GetOpt<string>("h4reco.run");    

    //---Get file list searching in specified path (eos or locally)
    if(path.find("/eos/cms") != string::npos)
    {
	// if ( getMachineDomain() != "cern.ch" )
        //     ls_command = string("gfal-ls root://eoscms/"+path+run+" | grep 'root' > /tmp/"+run+".list");
	// else
        ls_command = string("ls "+path+run+" | grep 'root' > /tmp/"+run+".list");
    }
    else if(path.find("srm://") != string::npos)
        ls_command = string("echo "+path+run+"/`gfal-ls "+path+run+
                            "` | sed -e 's:^.*\\/cms\\/:root\\:\\/\\/xrootd-cms.infn.it\\/\\/:g' | grep 'root' > /tmp/"+run+".list");
    else
        ls_command = string("ls "+path+run+" | grep 'root' > /tmp/"+run+".list");
    system(ls_command.c_str());

    ifstream waveList(string("/tmp/"+run+".list").c_str(), ios::in);
    while(waveList >> file && (opts.GetOpt<int>("h4reco.maxFiles")<0 || nFiles<opts.GetOpt<int>("h4reco.maxFiles")) )
    {
        //---skip files before specified spill
        auto currentSpill = std::stoi(file.substr(0, file.size()-4));
        if(firstSpill == -1 || (currentSpill >= firstSpill && currentSpill < firstSpill + opts.GetOpt<int>("h4reco.maxFiles")))
        {
            if(path.find("/eos/cms") != string::npos)
            {
                // if ( getMachineDomain() != "cern.ch" )
                // {
                //     std::cout << "+++ Adding file " << ("root://eoscms/"+path+run+"/"+file).c_str() << std::endl;
                //     inTree->AddFile(("root://eoscms/"+path+run+"/"+file).c_str());
                // }
                // else
                // {
                std::cout << "+++ Adding file " << (path+run+"/"+file).c_str() << std::endl;
                inTree->AddFile((path+run+"/"+file).c_str());
                // }
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
    }
    std::cout << "+++ Added " << nFiles << " files with " << inTree->GetEntries() << " events" << std::endl;
    return;
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
    if(argc > 3)
    {
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

    TChain* inTree = new TChain("H4tree");
    ReadInputFiles(opts, spill, inTree);
    H4Tree h4Tree(inTree);

    //-----output setup-----
    uint64 index=0;    
    TFile* outROOT;
    if(spill == -1)
        outROOT = new TFile(outSuffix+"run"+TString(run)+".root", "RECREATE");
    else if (nspills == -1)
        outROOT = new TFile(outSuffix+"run"+TString(run)+"_spill"+to_string(spill).c_str()+".root", "RECREATE");
    else
        // TODO: generate correct file name if less than the requested number of spills are available
        outROOT = new TFile(outSuffix+"run"+TString(run)+"_spills"+to_string(spill)+"-"+to_string(spill+nspills-1).c_str()+".root", "RECREATE");
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
	while((h4Tree.NextEntry() && (nEvents < maxEvents || maxEvents == -1)) || (isSim && (nEvents < maxEvents)))
        {
	    if(nEvents % 1000 == 0)
            {
		if(isSim)
                    cout << ">>> Generated events: " << nEvents << "/"
                         << maxEvents 
                         << endl;
		else                
                    cout << ">>> Processed events: " << nEvents << "/"
                         << (maxEvents<0 ? h4Tree.GetEntries() : min(h4Tree.GetEntries(), (uint64)maxEvents))
                         << endl;
		TrackProcess(cpu, mem, vsz, rss);
            }
        
	    //---set index value run*1e10+spill*1e4+event
	    index = h4Tree.runNumber*1e9 + h4Tree.spillNumber*1e5 + h4Tree.evtNumber;
	    
	    //---call ProcessEvent for each plugin and check the return status
	    bool status=true;
	    for(auto& plugin : pluginSequence)
            {
	    	status &= plugin->Clear(); 
	    	status &= plugin->ProcessEvent(h4Tree, pluginMap, opts);
            }
        
	    //---Fill the main tree with info variables and increase event counter
	    mainTree.time_stamp = h4Tree.evtTimeStart;
	    mainTree.evt_flag = status;
	    mainTree.run = h4Tree.runNumber;
	    mainTree.spill = h4Tree.spillNumber;
	    mainTree.event = h4Tree.evtNumber;
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
