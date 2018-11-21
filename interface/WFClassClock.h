#ifndef __WFCLASS_CLOCK_H__
#define __WFCLASS_CLOCK_H__

#include <random>

#include "WFClass.h"
#include "TGraph.h"

using namespace std;

class WFClassClock : public WFClass
{
public:
    //---ctors---
    WFClassClock() {};
    WFClassClock(float tUnit);

    //---getters---
    WFFitResults GetTime(string method, vector<float>& params) override;
    WFFitResults GetTimeLE(float thr = 0, int nmFitSamples=2, int npFitSamples=2, int min=-1, int max=-1) override;
    WFFitResults GetTimeCLK(float wleft=-1.3, float wright=1.3, int min=130, int max=900);
    WFFitResults TemplateFit(float offset=0., int lW=0, int hW=0) override;
};

#endif
