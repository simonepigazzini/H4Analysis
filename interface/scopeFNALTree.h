#ifndef __H4_TREE__
#define __H4_TREE__

#include <iostream>
#include <map>
#include <unordered_map>
#include <tuple>
#include <string>

#include "DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

#define NCHANNELS 4
#define NDIGIS 800

typedef unsigned long int uint32;
typedef unsigned long long int uint64;
 
//****************************************************************************************
//----------Tree reader class-------------------------------------------------------------
#define DYNAMIC_TREE_NAME scopeFNALTreeBase

#define DATA_TABLE                              \
    DATA(unsigned int,  i_evt)           \

#define DATA_VECT_TABLE                                     \
  DATA(float,channel, NCHANNELS*NDIGIS)		      \
  DATA(float,time, 1*800)       

#include "DynamicTTree/interface/DynamicTTreeInterface.h"

#undef DYNAMIC_TREE_NAME
#undef DATA_TABLE
#undef DATA_VECT_TABLE

class scopeFNALTree : public scopeFNALTreeBase
{
public:
    //---ctors
    scopeFNALTree(TChain* t):
        scopeFNALTreeBase(t)
        {
            Init();
        }
    scopeFNALTree(TTree* t):
        scopeFNALTreeBase(t)
        {
            Init();
        }    
    //---dtor
    ~scopeFNALTree();
    
    //---members
    void Init();
    uint64 GetEntries(){ return tree_->GetEntriesFast(); };
    
};
   
#endif 
