#include "CalibDigitizer.h"
#include "TF1.h"
#include "TVectorD.h"

//----------Utils-------------------------------------------------------------------------
bool CalibDigitizer::Begin(CfgManager& opts, uint64* index)
{
    //---inputs---
    if(!opts.OptExist(instanceName_+".srcInstanceName"))
    {
        cout << ">>> CalibDigitizer ERROR: no source plugin specified" << endl;
        return false;
    }
    srcInstance_ = opts.GetOpt<string>(instanceName_+".srcInstanceName");
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");

    nTotSamples_ = 0;
    for(auto& channel : channelsNames_)
        nTotSamples_ += opts.GetOpt<int>(channel+".nSamples");

    for (int i=0;i<nTotSamples_;++i)
      {
	fitters_.push_back(TLinearFitter(2)); //using a linear fitter with 2 parameters for each cell i % DeltaV_i + DeltaT_i*Slope %
	fitters_[i].SetFormula("hyp1"); 
      }

    for(auto& channel : channelsNames_)
      {
	int nSamples=opts.GetOpt<int>(channel+".nSamples");
	pair<int,int> channelID(opts.GetOpt<int>(channel+".digiGroup"),opts.GetOpt<int>(channel+".digiChannel"));
	for (int iSample=0;iSample<nSamples;++iSample)
	  digiCalibration_.channelCalibrations_[channelID].calibrations_.push_back({0.,0.});
      }

    RegisterSharedData(&digiCalibration_, "digiCalib", true);

    return true;
}


//----------Begin-------------------------------------------------------------------------
bool CalibDigitizer::BeginLoop(int iLoop, CfgManager& opts)
{  
  if (iLoop>0)
    {
      for (int i=0;i<nTotSamples_;++i)
	fitters_[i].ClearPoints();
    }
     
    return true;
}

//----------ProcessEvent------------------------------------------------------------------
bool CalibDigitizer::ProcessEvent(H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    //---load WFs from source instance shared data
    for(auto& channel : channelsNames_)
    {
        auto shared_data = plugins[srcInstance_]->GetSharedData(srcInstance_+"_"+channel, "", false);
        if(shared_data.size() != 0)
            WFs_[channel] = (WFClass*)shared_data.at(0).obj;
        else
            cout << "[CalibDigitizer::" << instanceName_ << "]: channels samples not found check DigiReco step" << endl; 
    }

    //---compute reco variables
    int offset=0;
    for(auto& channel : channelsNames_)
    {
        //---skip dead channels
        if(WFs_.find(channel) == WFs_.end())
            continue;

	auto analizedWF = WFs_[channel]->GetSamples();
	float tUnit = WFs_[channel]->GetTUnit();

	//Set function and fit initial values and ranges
	TF1* wFit=new TF1("wFit","[0]+[1]*TMath::Sin([2]*(x+[3]))",0,analizedWF->size()*tUnit);
	wFit->SetParameter(0,2400);
	wFit->SetParLimits(0,1800.,2600.);
	wFit->SetParameter(1,1350); 
	wFit->SetParLimits(1,900,1900);
	wFit->SetParameter(2,0.62); //80MHz 
	wFit->SetParLimits(2,0.5,0.7);
	wFit->SetParameter(3,0);
	wFit->SetParLimits(3,-50.,50.);

	WFs_[channel]->AnalyticFit(wFit,0,analizedWF->size()-50); //fit using most of the samples avoiding the last samples of the digitizer

	//should add a check on the fit status
	//get residual and slope and add it to the linear minimizer
	const std::vector<double>* samples = WFs_[channel]->GetSamples();
	for (int iSample=0;iSample<samples->size();++iSample)
	  {
	    double slope=wFit->Derivative(iSample*tUnit);
	    double residual=wFit->Eval(iSample*tUnit)-(*samples)[iSample];
	    int cellIndex=(iSample + WFs_[channel]->GetStartIndexCell())%samples->size();
	    //	    std::cout << Form("cellIndex %d %d:%f,%f",WFs_[channel]->GetStartIndexCell(),cellIndex,slope,residual) << std::endl; 
	    fitters_[offset+cellIndex].AddPoint(&slope,residual);
	  }	

	delete wFit;
	offset+=opts.GetOpt<int>(channel+".nSamples");
    }
  
    return true;
}

void CalibDigitizer::minimize()
{
    for (int i=0;i<nTotSamples_;++i)
	fitters_[i].Eval();
}

bool CalibDigitizer::EndLoop(int iLoop, CfgManager& opts)
{
  minimize();
  for(auto& channel : channelsNames_)
    {
      int nSamples=opts.GetOpt<int>(channel+".nSamples");
      pair<int,int> channelID(opts.GetOpt<int>(channel+".digiGroup"),opts.GetOpt<int>(channel+".digiChannel"));
      for (int iSample=0;iSample<nSamples;++iSample)
	{
	  TVectorD params;
	  TVectorD errors;
	  fitters_[iSample].GetParameters(params);
	  //	  fitters_[iSample].GetErrors(errors);

	  digiCalibration_.channelCalibrations_[channelID].calibrations_[iSample]=DigiChannelCalibration::calibrationParameters{params(0),params(1)};
	  // for (int i=0; i<2; i++)
	  //   printf("par[%d]=%f+-%f\n", i, params(i), errors(i));
	  //printf("chisquare=%f\n", fitters_[i].GetChisquare());
	}
    }

  return true;
}


bool CalibDigitizer::End(CfgManager& opts)
{
    return true;
}

