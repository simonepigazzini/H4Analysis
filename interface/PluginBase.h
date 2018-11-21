#ifndef __PLUGIN_BASE__
#define __PLUGIN_BASE__

#include <string>
#include <map>

#include "TObject.h"

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/H4Tree.h"

//**********Helper macros*****************************************************************
//---Helper macro to keep track of the plugin running: must be inserted at the
//   beginning of each method
#define CHECKPOINT()                            \
    SetCurrentMethod(__FUNCTION__);

//**********SHARED DATAFORMAT*************************************************************
struct SharedData
{    
    bool        permanent;
    std::string tag;
    TObject*    obj;
};

//**********PLUGIN BASE CLASS*************************************************************
class PluginBase
{
public:
    //---ctors---
    PluginBase(){};

    //---dtor---
    virtual ~PluginBase(){};    

    //---setters---
    void SetInstanceName(const std::string& instance) { instanceName_=instance; };
    void SetCurrentMethod(std::string method) { currentMethod_ = method; };

    //---getters---
    std::string        GetInstanceName() { return instanceName_; };
    std::string        GetCurrentMethod() { return currentMethod_; };
    vector<SharedData> GetSharedData(std::string tag="", std::string type="", bool permanent=true);
    
    //---utils---
    virtual bool Begin(map<std::string, PluginBase*>& plugins, CfgManager& opts, uint64* index)        { CHECKPOINT(); return true; };
    virtual bool ProcessEvent(H4Tree& event, map<std::string, PluginBase*>& plugins, CfgManager& opts) { CHECKPOINT(); return true; };
    virtual bool End(map<std::string, PluginBase*>& plugins, CfgManager& opts)                         { CHECKPOINT(); return true; };
    virtual bool BeginLoop(int iLoop, map<std::string, PluginBase*>& plugins, CfgManager& opts)        { CHECKPOINT(); return true; };
    virtual bool EndLoop(int iLoop, map<std::string, PluginBase*>& plugins, CfgManager& opts)          { CHECKPOINT(); return true; };

protected:
    //---utils---
    void RegisterSharedData(TObject* obj, std::string tag, bool isPermanent);

protected:
    //---keep track of the plugin name as defined in the cfg
    std::string instanceName_;
    //---keep track of the current running plugin method. Must be called manually in each function;
    std::string currentMethod_;
    //---shared data
    vector<SharedData> data_;
};

//---To be inserted at the end of each plugin definition
#define DEFINE_PLUGIN(NAME)                                         \
    extern "C" PluginBase* create() { return new NAME; }            \
    extern "C" void destroy(PluginBase* plugin) { delete plugin; }

#endif

