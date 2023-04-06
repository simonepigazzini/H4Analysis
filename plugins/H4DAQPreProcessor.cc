#include "H4DAQPreProcessor.h"

bool H4DAQPreProcessor::Begin(CfgManager& opts)
{
  return true;
}

H4Tree* H4DAQPreProcessor::ProcessEvent(DynamicTTreeBase* event, CfgManager& opts)
{
  H4Tree* thisEvent=dynamic_cast<H4Tree*>(event);
  return thisEvent;
}
