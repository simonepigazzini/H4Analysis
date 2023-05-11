#include "ScopeFNALPreProcessor.h"

bool ScopeFNALPreProcessor::Begin(CfgManager& opts)
{
  runNumber_=std::stoi(opts.GetOpt<std::string>("h4reco.run"));

  for (int ich=0;ich<NCHANNELS;++ich)
    {
      h4Tree_.digiMap[make_tuple(0, 0, ich)] = ich*NDIGIS;
      h4Tree_.digiNSamplesMap[make_tuple(0, 0, ich)] = NDIGIS;
    }
  return true;
}

H4Tree* ScopeFNALPreProcessor::ProcessEvent(DynamicTTreeBase* event, CfgManager& opts)
{
  scopeFNALTree* thisEvent=dynamic_cast<scopeFNALTree*>(event);
  h4Tree_.runNumber=runNumber_;
  h4Tree_.spillNumber=1;
  h4Tree_.evtNumber=thisEvent->i_evt;

  //filling digis
  h4Tree_.nDigiSamples=NCHANNELS*NDIGIS;
  for (int ich=0;ich<NCHANNELS;++ich)
    for (int is=0;is<NDIGIS;++is)
      {
	h4Tree_.digiBoard[ich*NDIGIS+is]=0;
	h4Tree_.digiGroup[ich*NDIGIS+is]=0;
	h4Tree_.digiChannel[ich*NDIGIS+is]=ich;
	h4Tree_.digiStartIndexCell[ich*NDIGIS+is]=0;
	h4Tree_.digiSampleGain[ich*NDIGIS+is]=0;
	h4Tree_.digiSampleValue[ich*NDIGIS+is]=(*thisEvent).channel[ich*NDIGIS+is];
	h4Tree_.digiSampleTime[ich*NDIGIS+is]=(*thisEvent).time[is]*1e9;
      }

  //now fill the digis map
  return &h4Tree_;
}
