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
    inline H4Tree& GetTree() {return *inTree_;};
    inline int&    GetNFiles() {return nFiles_;};
    inline int     GetNFilesProcessed() {return nFiles_-fileList_.size();};
    
    //---utils---    
    bool NextEvent();
    inline bool FirstEventInSpill() {return firstEventInSpill_;};

protected:
    bool ReadInputFiles();
    bool LoadNextFile();

private:
    CfgManager     opts_;
    vector<string> fileList_;
    int            nFiles_;
    TFile*         currentFile_;
    H4Tree*        inTree_;
    bool           firstEventInSpill_;    
    
};

#endif
