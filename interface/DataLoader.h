#ifndef __DATA_LOADER__
#define __DATA_LOADER__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "CfgManager/interface/CfgManager.h"

#include "interface/H4Tree.h"

class DataLoader
{
public:
    //---ctors---
    DataLoader();
    DataLoader(CfgManager& opts);
    
    //---dtor---
    ~DataLoader() {};

    //----getters---
    inline H4Tree&  GetTree() {return *inTree_;};
    inline int      GetNFiles() {return fileList_.size();};
    inline int      GetNFilesProcessed() {return iFile_;};
    inline long int GetCurrentEntry() {return inTree_->GetCurrentEntry();};
    
    //---utils---    
    bool NextEvent();
    inline bool FirstEventInSpill() {return firstEventInSpill_;};

protected:
    bool ReadInputFiles();
    bool LoadNextFile();

private:
    CfgManager     opts_;
    vector<string> fileList_;
    int            iFile_;
    TFile*         currentFile_;
    H4Tree*        inTree_;
    bool           firstEventInSpill_;    
    
};

#endif
