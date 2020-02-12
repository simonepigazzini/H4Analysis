#ifndef __ECAL_MATRIX_TREE__
#define __ECAL_MATRIX_TREE__

typedef unsigned long long int uint64;

#include "DynamicTTree/interface/DynamicTTreeBase.h"

//----------Tree reader class-------------------------------------------------------------
#define DYNAMIC_TREE_NAME ECALMatrixTreeBase

#define DATA_TABLE                              \
    DATA(int,   seed)                           \
    DATA(float, e3x3)                           \
    DATA(float, e5x5)                           \
    DATA(float, ecal_x)                         \
    DATA(float, ecal_y)

#include "DynamicTTree/interface/DynamicTTreeInterface.h"

#undef DYNAMIC_TREE_NAME
#undef DATA_TABLE

class ECALMatrixTree : public ECALMatrixTreeBase
{
public:
    //---ctors
    ECALMatrixTree():
        ECALMatrixTreeBase()
        {};
    ECALMatrixTree(TTree* t, uint64* index):
        ECALMatrixTreeBase(t)
        {
            index_ = index;
            Init();
        }    
    ECALMatrixTree(const char* name, const char* title, uint64* index):
        ECALMatrixTreeBase(name, title)
        {
            index_ = index;
            Init();
        }    
    
    //---members
    void Init();

    uint64* index_;
};

#endif
