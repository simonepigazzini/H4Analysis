#ifndef __H4DAQ_PREPROCESSOR__
#define __H4DAQ_PREPROCESSOR__

#include <iostream>
#include <algorithm>

#include "interface/PreProcessorBase.h"
#include "interface/H4Tree.h"

class H4DAQPreProcessor: public PreProcessorBase
{
public: 
    //---ctors---
  H4DAQPreProcessor(){};

  //---dtor---
  ~H4DAQPreProcessor(){};

  //---utils---
  bool Begin(CfgManager& opts);
  H4Tree* ProcessEvent(DynamicTTreeBase* event, CfgManager& opts);
  
private:    
};

DEFINE_PREPROCESSOR(H4DAQPreProcessor);

#endif
