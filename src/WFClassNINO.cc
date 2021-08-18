#include "interface/WFClassNINO.h"

//**********Constructors******************************************************************

WFClassNINO::WFClassNINO(int polarity, float tUnit):
    WFClass(polarity, tUnit)
{}

//**********Getters***********************************************************************

//----------Get the signal integral around the the max-------------------------------------
float WFClassNINO::GetSignalIntegral(int thr, int min=-1)
{
    //---compute position of the max
    if(min==-1 && sWinMin_!=-1)
        min = sWinMin_;

    //---find pulse borders
    SubtractBaseline();
    int begin=-1, end=-1;
    for(unsigned int iSample=min; iSample<samples_.size(); ++iSample)
    {
        if(begin==-1 && samples_[iSample]>thr)
            begin = iSample;
        else if(begin!=-1 && samples_[iSample]<thr)
        {
            end = iSample;
            break;
        }
    }            
    
    // //---compute integral
    // float integral=0;
    // for(int iSample=begin; iSample<end; ++iSample)
    // {
    //     //---if signal window goes out of bound return a bad value
    //     if(iSample > samples_.size() || iSample < 0)
    //         return -1000;        
    //     integral += samples_.at(iSample);
    // }
    
    return 1.*(end-begin);
}

//---------Add waveform sample to the list of uncalibrated samples------------------------
//---sample is inserted at the end of uncalibSamples_
//---the times vector is filled with the uncalibrated sample time computed from the time unit
void WFClassNINO::AddSample(float sample)
{
  uncalibSamples_.push_back(polarity_*sample); 
  times_.push_back( (samples_.size()-1.)*tUnit_ );
  samples_ = uncalibSamples_;
};
