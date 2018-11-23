#include "interface/WFClassClock.h"

//**********Constructors******************************************************************

WFClassClock::WFClassClock(float tUnit):
    WFClass(1., tUnit), clkPeriod_(-1), tmplFitPeriod_(-1)
{}

//**********Getters***********************************************************************

//----------Get time of a clock pulse-----------------------------------------------------
WFFitResults WFClassClock::GetTime(string method, vector<float>& params)
{
    //---CLK
    if(method.find("CLK") != string::npos)
    {
        if(params.size()==0)
            return GetTimeCLK();
        else if(params.size()==2)
            return GetTimeCLK(params[0], params[1]);
        else if(params.size()==3)
            return GetTimeCLK(params[0], params[1], params[2]);
        else
            cout << ">>>ERROR: wrong number of arguments passed for CLK time computation: " << params.size()
                 << ". usage: CLK [min max] [+/- 1]" << endl;

    }
    //---LED
    else if(method.find("LED") != string::npos)
    {
        if(params.size()<1)
            cout << ">>>ERROR: to few arguments passed for LED time computation" << endl;
        else if(params.size()<2)
            return GetTimeLE(params[0]);
        else if(params.size()<4)
            return GetTimeLE(params[0], params[1], params[2]);
        else
            return GetTimeLE(params[0], params[1], params[2], params[3], params[4]);

    }
   
    cout << ">>>ERROR, WFClassClock: time reconstruction method <" << method << "> not supported" << endl;
    return WFFitResults{-1, -1000, -1, 0};
}

//----------Get time of a clock pulse-----------------------------------------------------
WFFitResults WFClassClock::GetTimeLE(float thr, int nmFitSamples, int npFitSamples, int min, int max)
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
        //---find first sample above thr on RISING edge
        leThr_ = thr;
        for(int iSample=sWinMin_; iSample<sWinMax_; ++iSample)
        {
            if(samples_.at(iSample) > leThr_ && samples_.at(iSample) > samples_.at(iSample-1)) 
            {
                leSample_ = iSample;
                break;
            }
        }
        //---interpolate -- A+Bx = amp
        chi2le_ = LinearInterpolation(A, B, leSample_-nmFitSamples, leSample_+npFitSamples);
        leTime_ = (leThr_ - A) / B;
    }

    return WFFitResults{leThr_, leTime_, chi2le_, B};
}

//----------Get time of a clock pulse-----------------------------------------------------
//---Each rising edge zero crossing is measured, the clock period is
//---then extracted from a heuristic fit (cit Marc)
WFFitResults WFClassClock::GetTimeCLK(float wleft, float wright, int min, int max)
{
    //---check if signal window is valid
    if(min==max && max==-1 && sWinMin_==sWinMax_ && sWinMax_==-1)
        return WFFitResults{-1, -1000, -1, 0};
    //---setup signal window
    if(min!=-1 && max!=-1)
        SetSignalWindow(min, max);

    double S=0, S2=0, T=0, ST=0; 
    vector<double> clkZeros;
    TGraph gClockEdge(samples_.size(), times_.data(), samples_.data());
    TF1 fClockEdge("fClockEdge", "[0]*tanh((x-[1])/[2])+[3]", times_[min], times_[max]);
    
    for(int iSample=min; iSample<max; ++iSample)
    {
        if(samples_.at(iSample) < 0 && samples_.at(iSample+1) > 0)
        {
            //---edge fit
            auto t0 = times_[iSample];
            if(t0+wleft > times_[min] && t0+wright < times_[max])
            {
                fClockEdge.SetParameters(300., t0, 0.5, 0.);
                gClockEdge.Fit(&fClockEdge, "QRSO", "goff", t0+wleft, t0+wright);
                clkZeros.push_back(fClockEdge.GetParameter(1));            
                //---fill variables
                auto N = clkZeros.size();
                S += N;
                S2 += N*N;
                T += clkZeros.back();
                ST += N*clkZeros.back();
            
                ++iSample;
            }
        }
    }

    //---compute clock period, phase and phase error
    auto N = clkZeros.size();
    clkPeriod_ = (ST-S*T/N)/(S2-S*S/N);
    double phase = (T-S*clkPeriod_)/N;
    TH1F h_phase_err("h_phase_err", "", 75, -1, 1);
    for(unsigned int i=1; i<N; ++i)
    {
        auto t_clk = clkZeros[i]-clkZeros[0];
        h_phase_err.Fill(i*clkPeriod_-t_clk);
        //phase_err += (phase+i*period-t_clk)*(phase+i*period-t_clk);
    }
    //auto phase_err = h_phase_err.Fit("gaus", "QRSO", "goff")->Parameter(2);    
    auto phase_err = h_phase_err.GetRMS();
    
//    return WFFitResults{0, phase, std::sqrt(phase_err/(N-1)), 0};
    return WFFitResults{0, phase, phase_err, -1, 0};
}

//----------template fit to the WF--------------------------------------------------------
//---Template fit of each single clock cycle
WFFitResults WFClassClock::TemplateFit(float offset, int lW, int hW)
{
    if(tmplFitAmp_ == -1)
    {
        BaselineRMS();
        GetAmpMax();
        bRMS_ = 10.;
        int N=0;
        double S=0, S2=0, T=0, ST=0; 
        vector<double> clkZeros;

        //---cycles loop
        for(int iSample=sWinMin_; iSample<sWinMax_; ++iSample)
        {
            //---single clock cycle fit
            if(samples_.at(iSample) < 0 && samples_.at(iSample+1) > 0)
            {
                auto t0 = (times_.at(iSample+1)+times_.at(iSample))/2.;
                //---set template fit window around zero crossing, [min, max)
                fWinMin_ = iSample - lW;
                fWinMax_ = iSample + hW;
                //---setup minimization
                std::default_random_engine generator;
                std::uniform_real_distribution<float> rnd_t0(times_.at(iSample), times_.at(iSample+1));
                ROOT::Math::Functor chi2(this, &WFClassClock::TemplateChi2, 2);
                ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
                minimizer->SetMaxFunctionCalls(100000);
                minimizer->SetMaxIterations(1000);
                minimizer->SetTolerance(1e-2);
                minimizer->SetPrintLevel(0);
                minimizer->SetFunction(chi2);
                minimizer->SetLimitedVariable(0, "amplitude", GetAmpMax(), 1e-2, GetAmpMax()*0.8, GetAmpMax()*1.2);
                minimizer->SetLimitedVariable(1, "deltaT", t0, 1e-3, times_[fWinMin_], times_[fWinMax_]);
                //---fit
                minimizer->Minimize();
                int tries=0;
                while(minimizer->Status() != 0 && tries < 1000)
                {
                    minimizer->SetVariableValue(1, rnd_t0(generator));
                    minimizer->Minimize();
                    ++tries;
                }    
                if(tmplFitAmp_ == -1)
                    tmplFitAmp_ = minimizer->X()[0];
                tmplFitTime_ = minimizer->X()[1];
                clkZeros.push_back(minimizer->X()[1]);
                //---fill variables
                ++N;
                S += N;
                S2 += N*N;
                T += clkZeros.back();
                ST += N*clkZeros.back();
                
                ++iSample;
                
                delete minimizer;        
            }
        }

        //---compute clock period, phase and phase error
        tmplFitPeriod_ = (ST-S*T/N)/(S2-S*S/N);
        tmplFitTime_ = (T-S*tmplFitPeriod_)/N;
        TH1F h_phase_err("h_phase_err", "", 75, -1, 1);
        for(unsigned int i=1; i<N; ++i)
        {
            auto t_clk = clkZeros[i]-clkZeros[0];
            h_phase_err.Fill(i*tmplFitPeriod_-t_clk);
        }
        tmplFitTimeErr_ = h_phase_err.GetRMS();
        //tmplFitTimeErr_ = h_phase_err.Fit("gaus", "QRSO", "goff")->Parameter(2);
    }

    return WFFitResults{tmplFitAmp_, tmplFitTime_, tmplFitTimeErr_, TemplateChi2()/(hW-lW), 0};
}
