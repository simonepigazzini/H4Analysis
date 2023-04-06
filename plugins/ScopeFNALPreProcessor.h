#ifndef __SCOPEFNAL_PREPROCESSOR__
#define __SCOPEFNAL_PREPROCESSOR__

#include <iostream>
#include <algorithm>

#include "interface/PreProcessorBase.h"
#include "interface/H4Tree.h"
#include "interface/scopeFNALTree.h"

class ScopeFNALPreProcessor: public PreProcessorBase
{
public: 
    //---ctors---
  ScopeFNALPreProcessor(): h4Tree_("H4Tree","H4TreeBase"){};

  //---dtor---
  ~ScopeFNALPreProcessor(){};

  //---utils---
  bool Begin(CfgManager& opts);
  H4Tree* ProcessEvent(DynamicTTreeBase* event, CfgManager& opts);
  
private:    
    //---internal data
  int runNumber_;
  H4Tree h4Tree_;
};

DEFINE_PREPROCESSOR(ScopeFNALPreProcessor);

#endif
