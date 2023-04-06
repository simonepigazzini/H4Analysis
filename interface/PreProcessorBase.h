#ifndef __PREPROCESSOR_BASE__
#define __PREPROCESSOR_BASE__

#include <string>
#include <map>

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/H4Tree.h"

//**********Helper macros*****************************************************************
//---Helper macro to keep track of the preProcessor running: must be inserted at the
//   beginning of each method
#define CHECKPOINT()                            \
    SetCurrentMethod(__FUNCTION__);


//**********PREPROCESSOR BASE CLASS*************************************************************
class PreProcessorBase
{
public:

  //**********LOGGER LEVELS*****************************************************************
  enum LoggerLevel {
    INFO = 0,
    WARN = 1,
    ERR  = 2
  };
  
    //---ctors---
    PreProcessorBase(){};

    //---dtor---
    virtual ~PreProcessorBase(){};    

    //---setters---
    void SetPreProcessorType(const std::string& preProcessor) { preProcessorType_ = preProcessor; };
    void SetInstanceName(const std::string& instance) { instanceName_ = instance; };
    void SetCurrentMethod(std::string method) { currentMethod_ = method; };

    //---getters---
    std::string        GetPreProcessorType() { return preProcessorType_; };
    std::string        GetInstanceName() { return instanceName_; };
    std::string        GetCurrentMethod() { return currentMethod_; };
    
    //---utils---
    virtual H4Tree* ProcessEvent(DynamicTTreeBase* event)
    { CHECKPOINT(); return 0; };

    void Log(std::string message, LoggerLevel lv=INFO);

protected:
    //---keep track of the preProcessor type as defined in the cfg
    std::string preProcessorType_;
    //---keep track of the preProcessor name as defined in the cfg
    std::string instanceName_;
    //---keep track of the current running preProcessor method. Must be called manually in each function;
    std::string currentMethod_;
};

#undef CHECKPOINT

//---To be inserted at the end of each preProcessor definition
#define DEFINE_PREPROCESSOR(NAME)                                         \
    extern "C" PreProcessorBase* create() { return new NAME; }            \
    extern "C" void destroy(PreProcessorBase* preProcessor) { delete preProcessor; }

#endif

