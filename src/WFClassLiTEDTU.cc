#include "interface/WFClassLiTEDTU.h"

//**********Constructors******************************************************************

WFClassLiTEDTU::WFClassLiTEDTU(int polarity, float tUnit, DigiChannelCalibration* calibration):
    WFClass(polarity, tUnit)
{}

//**********Getters***********************************************************************

//---------Add waveform sample to the list of uncalibrated samples------------------------
//---sample is inserted at the end of uncalibSamples_
//---the times vector is filled with the uncalibrated sample time computed from the time unit
//---a gain for each sample can be added. This is stored in a separate vector as well multiplied to the sample value.
void WFClassLiTEDTU::AddSample(float sample, float gain) 
{
    uncalibSamples_.push_back(polarity_*sample*gain);
    gain_.push_back(gain);
    times_.push_back( (samples_.size()-1.)*tUnit_ );
    samples_ = uncalibSamples_;
};

//----------template fit to the WF--------------------------------------------------------
WFFitResults WFClassLiTEDTU::TemplateFit(float amp_threshold, float offset, int lW, int hW)
{
    double tmplFitChi2=0;
    if(tmplFitAmp_ == -1)
    {
        if(samples_[maxSample_]>amp_threshold)
        {
            //---set template fit window around maximum, [min, max)
            BaselineRMS();
            GetAmpMax();    
            fWinMin_ = maxSample_ + int(offset/tUnit_) - lW;
            fWinMax_ = maxSample_ + int(offset/tUnit_) + hW;
            //---setup minimization
            auto t0 = GetInterpolatedAmpMax().time;
            ROOT::Math::Functor chi2(this, &WFClassLiTEDTU::TemplateChi2, 2);
            ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
            minimizer->SetMaxFunctionCalls(100000);
            minimizer->SetMaxIterations(1000);
            minimizer->SetTolerance(1e-4);
            minimizer->SetPrintLevel(0);
            minimizer->SetFunction(chi2);
            minimizer->SetLimitedVariable(0, "amplitude", GetAmpMax(), 1e-2, 0., GetAmpMax()*2.);
            minimizer->SetLimitedVariable(1, "deltaT", t0, 1e-3, times_[fWinMin_], times_[fWinMax_]);        
            //---fit
            minimizer->Minimize();
            tmplFitAmp_ = minimizer->X()[0];
            tmplFitTime_ = minimizer->X()[1];
            tmplFitTimeErr_ = minimizer->Errors()[1];

            delete minimizer;
        }
        else
            return GetInterpolatedAmpMax();
    }    

    return WFFitResults{tmplFitAmp_, tmplFitTime_, tmplFitTimeErr_, TemplateChi2()/(fWinMax_-fWinMin_+1-2), 0};
}
