#include "plugins/AsyncDataProcessor.h"

//---------Default constructor------------------------------------------------------------
AsyncDataProcessor::AsyncDataProcessor():
    asyncDataFile_(NULL),
    h4Tree_(NULL),
    dataSelector_(NULL),
    currentSpill_(0),
    syncTolerance_(0),
    maxForwardTries_(0)
{}
    
//---------Begin--------------------------------------------------------------------------
bool AsyncDataProcessor::Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index)
{
    if(!opts.OptExist(instanceName_+".srcPath"))
    {
        Log("no data source path specified", ERR);
        return false;
    }
   
    
    syncTolerance_ = opts.OptExist(instanceName_+".syncTolerance_us") ? opts.GetOpt<double>(instanceName_+".syncTolerance_us") : 1e7;
    maxForwardTries_ = opts.OptExist(instanceName_+".maxForwardTries") ? opts.GetOpt<int>(instanceName_+".maxForwardTries") : 1;

    asyncPluginList_ = opts.GetOpt<vector<string> >(instanceName_+".asyncPluginList"); 

    //---create async data plugins
    //   - each plugin ProcessEvent method will be called by AsyncDataProcessor instead of
    //     H4Reco main loop.
    //   - Data acquired asynchronously are fetch from the H4tree in srcPath/#run/#spill.root
    //   - Any data registered by an asynchronously plugin is shared through the AsyncDataProcessor
    //     that handles it.
    pluginLoaders_.reserve(asyncPluginList_.size());
    for(auto& plugin : asyncPluginList_)
    {
        Log("Loading asynchronous plugin <"+plugin+">");
        //---create loader 
        loader_ = new PluginLoader<PluginBase>(opts.GetOpt<string>(plugin+".pluginType"));
        pluginLoaders_.push_back(loader_);
        pluginLoaders_.back()->Create();
        //---get instance and put it in the plugin sequence   
        PluginBase* newPlugin = pluginLoaders_.back()->CreateInstance(plugin);
        if(newPlugin)
        {
            pluginSequence_.push_back(newPlugin);
            pluginMap_[plugin] = pluginSequence_.back();
        }
        else
        {
            Log("plugin type "+opts.GetOpt<string>(plugin+".pluginType")+" is not defined.", ERR);
            return 0;
        }
    }

    //---begin
    for(auto& plugin : pluginSequence_)
    {
        //---call Begin() methods and check the return status
        bool r_status = plugin->Begin(plugins, opts, index);
        if(!r_status)
        {
            Log("plugin returned bad flag from Begin() call: "+plugin->GetInstanceName(), ERR);
            return r_status;
        }
        //---Get plugin permanent shared data
        for(auto& shared : plugin->GetSharedData("", "TTree", true))
            RegisterSharedData(shared.obj, shared.tag, shared.permanent); 

        //---Get plugin transient shared data
        for(auto& shared : plugin->GetSharedData("", "", false))
            RegisterSharedData(shared.obj, shared.tag, shared.permanent); 
    }
    
    return true;
}

//---------ProcessEvent-------------------------------------------------------------------
bool AsyncDataProcessor::ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{

    //---check if spill number has changed, if so open new async data file
    if(currentSpill_ != event.spillNumber)
    {
        Log("opening new file");
        
        //---clean data from previous spill
        if(h4Tree_)
        {
            h4Tree_->GetTTreePtr()->Delete();
            delete h4Tree_;
            h4Tree_ = NULL;
        }
        if(asyncDataFile_ && asyncDataFile_->IsOpen())
            asyncDataFile_->Close();
        if(dataSelector_)
        {
            dataSelector_->Delete();
            dataSelector_ = NULL;
        }

        deltaT_ = 1e7;
        
        //---get new file and data tree
        //   if the file is not found no event the spill is skipped by the async plugins
        currentSpill_ = event.spillNumber;
        char formatted_spill_number[5];
        std::sprintf(formatted_spill_number, "%04d", currentSpill_);
        asyncDataFile_ = TFile::Open((opts.GetOpt<string>(instanceName_+".srcPath")+"/"+
                                      to_string(event.runNumber)+"/"+
                                      formatted_spill_number+".root").c_str(), "READ");
        if(asyncDataFile_)
        {
            h4Tree_ = new H4Tree((TTree*)asyncDataFile_->Get("H4tree"));

            //---Initialize TTreeFormula
            //   The TTreeFormula specified by 'asyncEventSelection' is used
            //   to select events in the RC tree that are matched to the ones
            //   collected by the asynchrohous DR.
            //   The whole logic somehow assume that the number of events recorded
            //   by RC is always >= than those recorded by any DR
            auto selection = opts.OptExist(instanceName_+".asyncEventSelection") ?
                opts.GetOpt<string>(instanceName_+".asyncEventSelection") : "1";
            dataSelector_ = new TTreeFormula((instanceName_+"_selector").c_str(), selection.c_str(), event.GetTTreePtr());
        }
    }

    //---Event loop. Match RC event with asynchronous DR events
    if(dataSelector_ && dataSelector_->EvalInstance())
    {
        bool status=true;
        int retries=0;
        while(h4Tree_->NextEntry())
        {
            auto diff = double(h4Tree_->evtTime[0])-event.evtTime[0];
            //---Get time offset from first event in spill
            if(deltaT_ == 1e7)
                deltaT_ = diff;

            diff -= deltaT_;

            //---If 
            if(fabs(diff) > syncTolerance_)
                if(retries < maxForwardTries_)
                {
                    ++retries;
                    continue;
                }
                else
                {
                    h4Tree_->NextEntry(h4Tree_->GetCurrentEntry()-retries-1);
                    return true;
                }                    
            else
                break;
        }
       
        for(auto& plugin : pluginSequence_)
            status &= plugin->ProcessEvent(*h4Tree_, plugins, opts);

        return status;
    }

    return true;
}

