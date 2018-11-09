#include "interface/WFClass.h"

#include "TRandom3.h"
#include "TVirtualFFT.h"
#include "TMath.h"

//**********Constructors******************************************************************
WFClass::WFClass(int polarity, float tUnit):
    tUnit_(tUnit), polarity_(polarity), trigRef_(0), sWinMin_(-1), sWinMax_(-1), 
    bWinMin_(-1), bWinMax_(-1),  maxSample_(-1), fitAmpMax_(-1), fitTimeMax_(-1),
    fitChi2Max_(-1), baseline_(-1), bRMS_(-1), cfSample_(-1), cfFrac_(-1), cfTime_(-1),
    leSample_(-1), leTime_(-1), chi2cf_(-1), chi2le_(-1),
    fWinMin_(-1), fWinMax_(-1), tempFitTime_(-1), tempFitAmp_(-1),
    tempFitTimeScint_(-1), tempFitAmpScint_(-1), tempFitTimeSpike_(-1), tempFitAmpSpike_(-1),
    tempFitConverged_(false), f_max_(NULL),
    interpolator_(NULL), interpolatorScint_(NULL), interpolatorSpike_(NULL)
{}
//**********Getters***********************************************************************

//----------Get the max/min amplitude wrt polarity----------------------------------------
float WFClass::GetAmpMax(int min, int max)
{
    //---check if signal window is valid
    if(min==max && max==-1 && sWinMin_==sWinMax_ && sWinMax_==-1)
        return -1;
    //---setup signal window
    if(min!=-1 && max!=-1)
        SetSignalWindow(min, max);
    //---return the max if already computed
    else if(maxSample_ != -1)
        return samples_.at(maxSample_);

    //---find the max
    maxSample_=sWinMin_;
    for(int iSample=sWinMin_; iSample<sWinMax_; ++iSample)
    {
        if(iSample < 0)
            continue;
        if(iSample >= samples_.size())
            break;
        if(samples_.at(iSample) > samples_.at(maxSample_)) 
            maxSample_ = iSample;
    }    
    return samples_.at(maxSample_);
}

//----------Get the interpolated max/min amplitude wrt polarity---------------------------
WFFitResults WFClass::GetInterpolatedAmpMax(int min, int max, int nmFitSamples, int npFitSamples, string function, vector<float> params)
{
    //---check if already computed
    if(min==-1 && max==-1 && fitAmpMax_!=-1)
        return WFFitResults{fitAmpMax_, fitTimeMax_*tUnit_, fitChi2Max_, 0};
    //---check if signal window is valid
    if(min==max && max==-1 && sWinMin_==sWinMax_ && sWinMax_==-1)
        return {-1, -1000, -1, 0};
    //---setup signal window
    if(min!=-1 && max!=-1)
        SetSignalWindow(min, max);
    //---return the max if already computed
    else if(maxSample_ == -1) 
        GetAmpMax(min, max); 

    //---fit the max
    TH1F h_max("h_max", "", nmFitSamples+npFitSamples+1, maxSample_-nmFitSamples, maxSample_+npFitSamples);
    if(f_max_ == NULL)
        f_max_ = new TF1("f_max", function.c_str(), maxSample_-nmFitSamples-0.5, maxSample_+npFitSamples+0.5, TF1::EAddToList::kNo);
    else
        f_max_->SetRange(maxSample_-nmFitSamples-0.5, maxSample_+npFitSamples+0.5);
    for(unsigned int ipar = 0; ipar < params.size(); ++ipar)
        f_max_->SetParameter(ipar, params.at(ipar));
    
    int bin=1;
    for(int iSample=maxSample_-nmFitSamples; iSample<=maxSample_+npFitSamples; ++iSample)
    {
        h_max.SetBinContent(bin, samples_[iSample]);
        h_max.SetBinError(bin, BaselineRMS());
        ++bin;
    }
    
    if(h_max.GetMaximum() != 0)
    {
        auto fit_result = h_max.Fit(f_max_, "QRSO");
        fitTimeMax_ = f_max_->GetMaximumX();
        fitAmpMax_ = f_max_->GetMaximum();
        fitChi2Max_ = f_max_->GetChisquare() / f_max_->GetNDF();
    }
    else
    {
        fitTimeMax_ = -1;
        fitAmpMax_ = 1000;
        fitChi2Max_ = -1;
    }
    
    return WFFitResults{fitAmpMax_, fitTimeMax_*tUnit_, fitChi2Max_, 0};
}

//----------Get time with the specified method--------------------------------------------
WFFitResults WFClass::GetTime(string method, vector<float>& params)
{
    //---CFD
    if(method.find("CFD") != string::npos)
    {
        if(params.size()<1)
            cout << ">>>ERROR, WFClass: to few arguments passed for CFD time computation" << endl;
        else if(params.size()<2)
            return GetTimeCF(params[0]);
        else if(params.size()<3)
            return GetTimeCF(params[0], params[1]);
        else
            return GetTimeCF(params[0], params[1], params[2], params[3]);
    }
    
    //---LED
    else if(method.find("LED") != string::npos)
    {
        if(params.size()<1)
            cout << ">>>ERROR, WFClass: to few arguments passed for LED time computation" << endl;
        else if(params.size()<2)
            return GetTimeLE(params[0]);
        else if(params.size()<4)
            return GetTimeLE(params[0], params[1], params[2]);
        else
            return GetTimeLE(params[0], params[1], params[2], params[3], params[4]);
    }
    
    //---TED
    else if(method.find("TED") != string::npos)
    {
        if(params.size()<1)
            cout << ">>>ERROR, WFClass: to few arguments passed for TED time computation" << endl;
        else if(params.size()<2)
            return GetTimeTE(params[0]);
        else if(params.size()<4)
            return GetTimeTE(params[0], params[1], params[2]);
        else
            return GetTimeTE(params[0], params[1], params[2], params[3], params[4]);
    }
    
    cout << ">>>ERROR, WFClass: time reconstruction method <" << method << "> not supported" << endl;
    return WFFitResults{-1, -1000, -1, 0};
}

//----------Get CF time for a given fraction and in a given range-------------------------
WFFitResults WFClass::GetTimeCF(float frac, int nFitSamples, int min, int max)
{
    float A=0, B=0;
    if(frac != cfFrac_ || cfSample_ != -1)
    {
        //---setups---
        int tStart=min;
        if(tStart == -1)
            tStart=sWinMin_ == -1 ? 0 : sWinMin_;
        cfSample_ = tStart;
        cfFrac_ = frac;
        if(fitAmpMax_ == -1)
            GetInterpolatedAmpMax(min, max);
        if(frac == 1) 
            return WFFitResults{fitAmpMax_, maxSample_*tUnit_, 1, 0};
        
        //---find first sample above Amax*frac
        for(int iSample=maxSample_; iSample>tStart; --iSample)
        {
            if(samples_.at(iSample) < fitAmpMax_*frac) 
            {
                cfSample_ = iSample;
                break;
            }
        }
        //---interpolate -- A+Bx = frac * amp
        chi2cf_ = LinearInterpolation(A, B, cfSample_-(nFitSamples-1)/2, cfSample_+(nFitSamples-1)/2);
        cfTime_ = (fitAmpMax_ * frac - A) / B;
        cfSlope_ = B;
    }

    return WFFitResults{fitAmpMax_*frac, cfTime_, chi2cf_, cfSlope_};
}

//----------Get leading edge time at a given threshold and in a given range---------------
WFFitResults WFClass::GetTimeLE(float thr, int nmFitSamples, int npFitSamples, int min, int max)
{
    //---check if signal window is valid
    if(min==max && max==-1 && sWinMin_==sWinMax_ && sWinMax_==-1)
        return WFFitResults{leThr_, -1000, -1, 0};
    //---setup signal window
    if(min!=-1 && max!=-1)
        SetSignalWindow(min, max);
    //---compute LED time value 
    float A=0, B=0;
    if(thr != leThr_ || leSample_ == -1)
    {
        //---find first sample above thr
        leThr_ = thr;
        leSample_ = -1;
        for(int iSample=sWinMin_; iSample<sWinMax_; ++iSample)
        {
            if(samples_.at(iSample) > leThr_) 
            {
                leSample_ = iSample;
                break;
            }
        }
        if(leSample_>0)
        {
            //---interpolate -- A+Bx = amp
            chi2le_ = LinearInterpolation(A, B, leSample_-nmFitSamples, leSample_+npFitSamples);
            leTime_ = (leThr_ - A) / B;
            leSlope_ = B;
        }
        else
        {
            leTime_ = -1000;
            chi2le_ = -1;
            leSlope_ = 0;
        }
    }

    return WFFitResults{leThr_, leTime_, chi2le_, leSlope_};
}

//----------Get trailing edge time at a given threshold and in a given range---------------
WFFitResults WFClass::GetTimeTE(float thr, int nmFitSamples, int npFitSamples, int min, int max)
{
    //---check if signal window is valid
    if(min==max && max==-1 && sWinMin_==sWinMax_ && sWinMax_==-1)
        return WFFitResults{teThr_, -1000, -1, 0};
    //---setup signal window
    if(min!=-1 && max!=-1)
        SetSignalWindow(min, max);
    //---compute TED time value 
    float A=0, B=0;
    if(thr != teThr_ || teSample_ == -1)
    {
        //---find first sample above thr
        teThr_ = thr;
        teSample_ = -1;
        for(int iSample=maxSample_; iSample<samples_.size(); ++iSample)
        {
            if(samples_.at(iSample) < teThr_) 
            {
                teSample_ = iSample;
                break;
            }
        }
        if(teSample_>0)
        {
            //---interpolate -- A+Bx = amp
            chi2te_ = LinearInterpolation(A, B, teSample_-nmFitSamples, teSample_+npFitSamples);
            teTime_ = (teThr_ - A) / B;
            teSlope_ = B;
        }
        else
        {
            leTime_ = -1000;
            chi2le_ = -1;
            teSlope_ = 0;
        }
    }

    return WFFitResults{teThr_, teTime_, chi2te_, teSlope_};
}

//----------Get the waveform integral in the given range----------------------------------
float WFClass::GetIntegral(int min, int max)
{
    //---compute integral
    float integral=0;
    for(int iSample=min; iSample<max; ++iSample)
    {
        if(iSample < 0)
            continue;
        if(iSample >= samples_.size())
            break;
        integral += samples_.at(iSample);
    }

    return integral;
}

//----------Get the signal integral around the the max-------------------------------------
float WFClass::GetSignalIntegral(int riseWin, int fallWin)
{
    //---compute position of the max
    if(maxSample_ == -1)
        GetAmpMax();

    //---compute integral
    float integral=0;
    for(int iSample=maxSample_-riseWin; iSample<maxSample_+fallWin; ++iSample)
    {
        //---if signal window goes out of bound return a bad value
        if(iSample < 0)
            continue;
        if(iSample >= samples_.size())
            break;
        integral += samples_.at(iSample);
    }

    return integral;
}


//----------Get the integral of Abs(WF) over the given range------------------------------
float WFClass::GetModIntegral(int min, int max)
{   
    float integral=0;
    for(int iSample=min; iSample<max; ++iSample)
    {
        if(iSample < 0)
            continue;
        if(iSample >= samples_.size())
            break;
        if(samples_.at(iSample) < 0)
            integral -= samples_.at(iSample);
        else
            integral += samples_.at(iSample);
    }
    return integral;
}

//**********Setters***********************************************************************

//----------Set the signal windows---------------------------------------------------------
void WFClass::SetSignalWindow(int min, int max)
{
    sWinMin_ = std::max(int(min + trigRef_), 0);
    sWinMax_ = std::min(int(max + trigRef_), int(samples_.size()));
}

void WFClass::SetSignalIntegralWindow(int min, int max)
{
    sIntWinMin_ = min + trigRef_;
    sIntWinMax_ = max + trigRef_;
}

//----------Set the baseline windows-------------------------------------------------------
void WFClass::SetBaselineWindow(int min, int max)
{
    bWinMin_ = std::max(min, 0);
    bWinMax_ = std::min(max, int(samples_.size()));
}

void WFClass::SetBaselineIntegralWindow(int min, int max)
{
    bIntWinMin_ = min;
    bIntWinMax_ = max;
}

//----------Set the fit template----------------------------------------------------------
void WFClass::SetTemplate(TH1* templateWF)
{
    //---check input
    if(!templateWF)
    {
        cout << ">>>ERROR: template passed as input does not exist" << endl;
        return;
    }

    //---reset template fit variables
    if(interpolator_)
        return;

    interpolator_ = new ROOT::Math::Interpolator(0, ROOT::Math::Interpolation::kCSPLINE);
    tempFitTime_ = templateWF->GetBinCenter(templateWF->GetMaximumBin());
    tempFitAmp_ = -1;

    //---fill interpolator data
    vector<double> x, y;
    for(int iBin=1; iBin<=templateWF->GetNbinsX(); ++iBin)
    {
        x.push_back(templateWF->GetBinCenter(iBin)-tempFitTime_);
        y.push_back(templateWF->GetBinContent(iBin));
    }
    interpolator_->SetData(x, y);

    return;
}

//----------Set the scintillation plus spike fit templates--------------------------------
void WFClass::SetTemplateScint(TH1* templateWF)
{
    //---check input
    if(!templateWF)
    {
        cout << ">>>ERROR: scintillation template passed as input does not exist" << endl;
        return;
    }

    //---reset template fit variables
    if(interpolatorScint_)
        return;

    interpolatorScint_ = new ROOT::Math::Interpolator(0, ROOT::Math::Interpolation::kCSPLINE);
    tempFitTimeScint_ = templateWF->GetBinCenter(templateWF->GetMaximumBin());
    tempFitAmpScint_ = -1;

    //---fill interpolator data
    vector<double> x, y;
    for(int iBin=1; iBin<=templateWF->GetNbinsX(); ++iBin)
    {
        x.push_back(templateWF->GetBinCenter(iBin)-tempFitTimeScint_);
        y.push_back(templateWF->GetBinContent(iBin));
    }
    interpolatorScint_->SetData(x, y);

    return;
}

void WFClass::SetTemplateSpike(TH1* templateWF)
{
    //---check input
    if(!templateWF)
    {
        cout << ">>>ERROR: spike template passed as input does not exist" << endl;
        return;
    }

    //---reset template fit variables
    if(interpolatorSpike_)
        return;

    interpolatorSpike_ = new ROOT::Math::Interpolator(0, ROOT::Math::Interpolation::kCSPLINE);
    tempFitTimeSpike_ = templateWF->GetBinCenter(templateWF->GetMaximumBin());
    tempFitAmpSpike_ = -1;

    //---fill interpolator data
    vector<double> x, y;
    for(int iBin=1; iBin<=templateWF->GetNbinsX(); ++iBin)
    {
        x.push_back(templateWF->GetBinCenter(iBin)-tempFitTimeSpike_);
        y.push_back(templateWF->GetBinContent(iBin));
    }
    interpolatorSpike_->SetData(x, y);

    return;
}

//**********Utils*************************************************************************

//----------Reset: get new set of sample, keep interpolator-------------------------------
void WFClass::Reset()
{
    sWinMin_ = -1;
    sWinMax_ = -1;
    bWinMin_ = -1;
    bWinMax_ = -1;
    maxSample_ = -1;
    fitAmpMax_ = -1;
    fitTimeMax_ = -1;
    fitChi2Max_ = -1;
    baseline_ = -1;
    bRMS_ = -1;
    cfSample_ = -1;
    cfFrac_ = -1;
    cfTime_ = -1;
    leSample_ = -1;
    leTime_ = -1;
    chi2cf_ = -1;
    chi2le_ = -1;
    fWinMin_ = -1;
    fWinMax_ = -1;
    tempFitTime_ = -1;
    tempFitAmp_ = -1;
    tempFitTimeScint_ = -1;
    tempFitAmpScint_ = -1;
    tempFitTimeSpike_ = -1;
    tempFitAmpSpike_ = -1;
    tempFitConverged_ = false;
    samples_.clear();
} 

//---------estimate the baseline in a given range and then subtract it from the signal----
WFBaseline WFClass::SubtractBaseline(int min, int max)
{
    if(min!=-1 && max==-1)
    {
        bWinMin_=min;
        bWinMax_=max;
    }
    //---compute baseline
    float baseline_=0;
    for(int iSample=bWinMin_; iSample<bWinMax_; ++iSample)
    {
        if(iSample < 0)
            continue;
        if(iSample >= samples_.size())
            break;
        baseline_ += samples_.at(iSample);
    }
    baseline_ = baseline_/((float)(bWinMax_-bWinMin_));
    //---subtract baseline
    for(unsigned int iSample=0; iSample<samples_.size(); ++iSample)
        samples_.at(iSample) = (samples_.at(iSample) - baseline_);    
    //---interpolate baseline
    BaselineRMS();
    float A=0, B=0;
    float chi2 = LinearInterpolation(A, B, bWinMin_, bWinMax_);
    
    return WFBaseline{baseline_, bRMS_, A, B, chi2};
}

//----------template fit to the WF--------------------------------------------------------
WFFitResults WFClass::TemplateFit(float offset, int lW, int hW)
{
    if(tempFitAmp_ == -1)
    {
        //---set template fit window around maximum, [min, max)
        BaselineRMS();
        GetAmpMax();
        fWinMin_ = maxSample_ + int(offset/tUnit_) - lW;
        fWinMax_ = maxSample_ + int(offset/tUnit_) + hW;
        //---setup minimization
        ROOT::Math::Functor chi2(this, &WFClass::TemplateChi2, 2);
        ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimizer->SetMaxFunctionCalls(100000);
        minimizer->SetMaxIterations(1000);
        minimizer->SetTolerance(1e-3);
        minimizer->SetPrintLevel(0);
        minimizer->SetFunction(chi2);
        minimizer->SetLimitedVariable(0, "amplitude", GetAmpMax(), 1e-2, 0., GetAmpMax()*2.);
        minimizer->SetLimitedVariable(1, "deltaT", maxSample_*tUnit_, 1e-2, fWinMin_*tUnit_, fWinMax_*tUnit_);
        //---fit
        minimizer->Minimize();
        tempFitAmp_ = minimizer->X()[0];
        tempFitTime_ = minimizer->X()[1];

        delete minimizer;
    }

    return WFFitResults{tempFitAmp_, tempFitTime_, TemplateChi2()/(fWinMax_-fWinMin_-2), 0};
}

WFFitResultsScintPlusSpike WFClass::TemplateFitScintPlusSpike(float offset, int lW, int hW)
{
    if(tempFitAmpScint_ == -1 && tempFitAmpSpike_ == -1)
    {
        //---set template fit window around maximum, [min, max)
        BaselineRMS();
        GetAmpMax();
        fWinMin_ = maxSample_ + int(offset/tUnit_) - lW;
        fWinMax_ = maxSample_ + int(offset/tUnit_) + hW;
        //---setup minimization
        ROOT::Math::Functor chi2(this, &WFClass::TemplatesChi2, 4);
        ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimizer->SetMaxFunctionCalls(100000);
        minimizer->SetMaxIterations(1000);
        minimizer->SetTolerance(1e-3);
        minimizer->SetPrintLevel(0);
        minimizer->SetFunction(chi2);
        minimizer->SetLimitedVariable(0, "amplitude_scint", GetAmpMax(), 1e-2, 0., GetAmpMax()*2.);
        minimizer->SetLimitedVariable(1, "deltaT_scint", maxSample_*tUnit_, 1e-2, fWinMin_*tUnit_, fWinMax_*tUnit_);
        minimizer->SetLimitedVariable(2, "amplitude_spike", GetAmpMax(), 1e-2, 0., GetAmpMax()*2.);
        minimizer->SetLimitedVariable(3, "deltaT_spike", maxSample_*tUnit_, 1e-2, fWinMin_*tUnit_, fWinMax_*tUnit_);
        //minimizer->SetFixedVariable(2, "amplitude_spike", 0.);
        //minimizer->SetFixedVariable(3, "deltaT_spike", 0.);
        //---fit
        tempFitConverged_ = minimizer->Minimize();
        tempFitAmpScint_ = minimizer->X()[0];
        tempFitTimeScint_ = minimizer->X()[1];
        tempFitAmpSpike_ = minimizer->X()[2];
        tempFitTimeSpike_ = minimizer->X()[3];

        delete minimizer;
    }

    return WFFitResultsScintPlusSpike{tempFitAmpScint_, tempFitTimeScint_, tempFitAmpSpike_, tempFitTimeSpike_, TemplatesChi2()/(fWinMax_-fWinMin_-4), tempFitConverged_};
}

void WFClass::EmulatedWF(WFClass& wf,float rms, float amplitude, float time)
{
    TRandom3 rnd(0);

    wf.Reset();

    if (tempFitTime_ == -1)
    {
        std::cout << "ERROR: no TEMPLATE for this WF" << std::endl;
        return;
    }

    for (unsigned int i=0; i<samples_.size();++i)
    {
        float emulatedSample=amplitude*interpolator_->Eval(i*tUnit_-tempFitTime_-(time-tempFitTime_));
        emulatedSample+=rnd.Gaus(0,rms);
        wf.AddSample(emulatedSample);
    }
}


void WFClass::FFT(WFClass& wf, float tau, int cut)
{
    if(samples_.size() == 0)
    {
        std::cout << "ERROR: EMPTY WF" << std::endl;
        return;
    }

    wf.Reset();

    int n=samples_.size();
    TVirtualFFT *vfft = TVirtualFFT::FFT(1,&n,"C2CFORWARD");

    Double_t orig_re[n],orig_im[n];
    for(int i=0;i<n;i++) 
    {
        orig_re[i]=samples_[i];
        if(i>1000) orig_re[i]=orig_re[999];// DIGI CAENV1742 NOT USABLE
        orig_im[i]=0;
    }
    vfft->SetPointsComplex(orig_re,orig_im);
    vfft->Transform();
    Double_t re[n],im[n];
    vfft->GetPointsComplex(re,im);

    TVirtualFFT *vinvfft = TVirtualFFT::FFT(1,&n,"C2CBACKWARD M K");
    Double_t cut_re[n],cut_im[n];

    for(int i=0;i<n;i++) 
    {
        if( i> cut-1 && i<n-cut) 
        {
            int delta = TMath::Min(i-cut-1,n-cut-i); 
            double dump=TMath::Exp(-delta/tau);
            cut_im[i]=im[i]*dump;
            cut_re[i]=re[i]*dump;
            continue;
        }
        cut_re[i]= re[i];
        cut_im[i]= im[i];
    }

    vinvfft->SetPointsComplex(cut_re,cut_im);
    vinvfft->Transform();
    Double_t inv_re[n],inv_im[n];
    vinvfft->GetPointsComplex(inv_re,inv_im);

    for(int i=0;i<n ;i++)
        wf.AddSample(inv_re[i]/n);

    delete vinvfft;
    delete vfft;

    return;
}

//----------compute baseline RMS (noise)--------------------------------------------------
float WFClass::BaselineRMS()
{
    if(bRMS_ != -1)
        return bRMS_;

    int nSample=0;
    float sum=0, sum2=0;
    for(int iSample=bWinMin_; iSample<bWinMax_; ++iSample)
    {
        ++nSample;
        sum += samples_[iSample];
        sum2 += samples_[iSample]*samples_[iSample];
    }
    
    bRMS_=sqrt(sum2/nSample - pow(sum/nSample, 2));
    return bRMS_;
}

//----------Linear interpolation util-----------------------------------------------------
float WFClass::LinearInterpolation(float& A, float& B, const int& min, const int& max)
{
    //---definitions---
    float xx= 0.;
    float xy= 0.;
    float Sx = 0.;
    float Sy = 0.;
    float Sxx = 0.;
    float Sxy = 0.;

    //---compute sums
    int usedSamples=0;
    for(int iSample=min; iSample<=max; ++iSample)
    {
        if(iSample<0 || iSample>=int(samples_.size())) 
            continue;
        xx = iSample*iSample*tUnit_*tUnit_;
        xy = iSample*tUnit_*samples_[iSample];
        Sx = Sx + (iSample)*tUnit_;
        Sy = Sy + samples_[iSample];
        Sxx = Sxx + xx;
        Sxy = Sxy + xy;
        ++usedSamples;
    }
    
    float Delta = usedSamples*Sxx - Sx*Sx;
    A = (Sxx*Sy - Sx*Sxy) / Delta;
    B = (usedSamples*Sxy - Sx*Sy) / Delta;

    //---compute chi2---
    float chi2=0;
    float sigma2 = pow(bRMS_, 2);
    for(int iSample=min; iSample<=max; ++iSample)
    {
        if(iSample<0 || iSample>=int(samples_.size())) 
            continue;
        chi2 = chi2 + pow(samples_[iSample] - A - B*iSample*tUnit_, 2)/sigma2;
    } 

    return chi2/(usedSamples-2);
}

//----------chi2 for template fit---------------------------------------------------------
double WFClass::TemplateChi2(const double* par)
{
    double chi2 = 0;
    double delta = 0;
    for(int iSample=fWinMin_; iSample<=fWinMax_; ++iSample)
    {
        if(iSample < 0 || iSample >= int(samples_.size()))
        {
            //cout << ">>>WARNING: template fit out of samples rage (chi2 set to -1)" << endl;
            chi2 += 9999;
        }
        else
        {
            //---fit: par[0]*ref_shape(t-par[1]) par[0]=amplitude, par[1]=DeltaT
            //---if not fitting return chi2 value of best fit
            if(par)
                delta = (samples_[iSample] - par[0]*interpolator_->Eval(iSample*tUnit_-par[1]))/bRMS_;
            else
                delta = (samples_[iSample] - tempFitAmp_*interpolator_->Eval(iSample*tUnit_-tempFitTime_))/bRMS_;
            chi2 += delta*delta;
        }
    }

    return chi2;
}

//----------chi2 for scintillation plus spike template fit---------------------------------------------------------
double WFClass::TemplatesChi2(const double* par)
{
    double chi2 = 0;
    double delta = 0;
    for(int iSample=fWinMin_; iSample<=fWinMax_; ++iSample)
    {
        if(iSample < 0 || iSample >= int(samples_.size()))
        {
            //cout << ">>>WARNING: template fit out of samples rage (chi2 set to -1)" << endl;
            chi2 += 9999;
        }
        else
        {
            //---fit: par[0]*ref_shape_scint(t-par[1]) + par[2]*ref_shape_spike(t-par[3])
            //---with par[0]=amplitude scintillation, par[1]=DeltaT scintillation
            //---and  par[2]=amplitude spike, par[3]=DeltaT spike
            //---if not fitting return chi2 value of best fit
            const auto iSampleTime = iSample*tUnit_;
            if(par)
                delta = (samples_[iSample] - par[0]*interpolatorScint_->Eval(iSampleTime-par[1]) - par[2]*interpolatorSpike_->Eval(iSampleTime-par[3]))/bRMS_;
            else
                delta = (samples_[iSample] - tempFitAmpScint_*interpolatorScint_->Eval(iSampleTime-tempFitTimeScint_) - tempFitAmpSpike_*interpolatorSpike_->Eval(iSampleTime-tempFitTimeSpike_))/bRMS_;
            chi2 += delta*delta;
        }
    }

    return chi2;
}

void WFClass::Print()
{
    std::cout << "+++ DUMP WF +++" << std::endl;
    for (unsigned int i=0; i<samples_.size(); ++i)
        std::cout << "SAMPLE " << i << ": " << samples_[i] << std::endl;
}

//**********operators*********************************************************************
//----------assignment--------------------------------------------------------------------
WFClass& WFClass::operator=(const WFClass& origin)
{
    samples_ = origin.samples_;
    tUnit_ = origin.tUnit_;
    polarity_ = origin.polarity_;
    sWinMin_ = origin.sWinMin_;
    sWinMax_ = origin.sWinMax_;
    bWinMin_ = origin.bWinMin_;
    bWinMax_ = origin.bWinMax_;
    maxSample_ = origin.maxSample_;
    fitAmpMax_ = origin.fitAmpMax_;
    baseline_ = origin.baseline_;
    bRMS_ = origin.bRMS_;
    cfSample_ = origin.cfSample_;
    cfFrac_ = origin.cfFrac_;
    cfTime_ = origin.cfTime_;
    leSample_ = origin.leSample_;
    leThr_ = origin.leThr_;
    leTime_ = origin.leTime_;
    chi2cf_ = origin.chi2cf_;
    chi2le_ = origin.chi2le_;
    fWinMin_ = origin.fWinMin_;
    fWinMax_ = origin.fWinMax_;
    tempFitTime_ = origin.tempFitTime_;
    tempFitAmp_ = origin.tempFitAmp_;
    tempFitTimeScint_ = origin.tempFitTimeScint_;
    tempFitAmpScint_ = origin.tempFitAmpScint_;
    tempFitTimeSpike_ = origin.tempFitTimeSpike_;
    tempFitAmpSpike_ = origin.tempFitAmpSpike_;
    tempFitConverged_ = origin.tempFitConverged_;
    interpolator_ = NULL;
    interpolatorScint_ = NULL;
    interpolatorSpike_ = NULL;

    return *this;
}

//----------subtraction-------------------------------------------------------------------
WFClass WFClass::operator-(const WFClass& sub)
{
    if(tUnit_ != sub.tUnit_)
        return *this;

    WFClass diff(1, tUnit_);
    for(unsigned int iSample=0; iSample<min(samples_.size(), sub.samples_.size()); ++iSample)
        diff.AddSample(samples_[iSample] - sub.samples_[iSample]);

    return diff;
}

//----------addition----------------------------------------------------------------------
WFClass WFClass::operator+(const WFClass& sub)
{
    if(tUnit_ != sub.tUnit_)
        return *this;

    WFClass sum(1, tUnit_);
    for(unsigned int iSample=0; iSample<min(samples_.size(), sub.samples_.size()); ++iSample)
        sum.AddSample(samples_[iSample] + sub.samples_[iSample]);

    return sum;
}

//----------subtraction and assignment----------------------------------------------------
WFClass& WFClass::operator-=(const WFClass& sub)
{
    if(tUnit_ == sub.tUnit_)
        for(unsigned int iSample=0; iSample<min(samples_.size(), sub.samples_.size()); ++iSample)
            samples_[iSample] -= sub.samples_[iSample];
    
    return *this;
}

//----------addition and assignmet--------------------------------------------------------
WFClass& WFClass::operator+=(const WFClass& sub)
{
    if(tUnit_ == sub.tUnit_)
        for(unsigned int iSample=0; iSample<min(samples_.size(), sub.samples_.size()); ++iSample)
            samples_[iSample] += sub.samples_[iSample];
    
    return *this;
}
