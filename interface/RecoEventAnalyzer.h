#ifndef _RECO_EVENT_ANALYZER_
#define _RECO_EVENT_ANALYZER_

#include <map>
#include <string>
#include <vector>
#include <algorithm>

#include "TTree.h"

#include "interface/PluginBase.h"

class RecoEventAnalyzer
{
public:
    //---ctors---    
    RecoEventAnalyzer() {};
    RecoEventAnalyzer(std::map<std::string, PluginBase*>& plugins, uint64* index);

    //---utils---
    // Evaluate simple selection
    bool EvaluateSelection(const char* expr)
        {
            std::string s(expr);
            return EvaluateSelection(s);
        };
    bool EvaluateSelection(std::string& expr);
    // Evaluate expression and selections and return expression values for the current entry
    const std::vector<std::vector<double> >& GetValues(const char* expr, const char* cuts)
        {
            std::string s(expr), c(cuts);
            return GetValues(s, c);
        };
    const std::vector<std::vector<double> >& GetValues(std::string& expr, std::string cuts="1");

protected:
    void                UpdateEntry();
    
private:
    uint64*                 index_;
    uint64                  currentEntry_;
    string                  indexCheckStr_;
    TTree*                  sTree_;
    vector<vector<double> > lastValues_;
};

#endif
