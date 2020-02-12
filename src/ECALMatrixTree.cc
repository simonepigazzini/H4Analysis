#include "interface/ECALMatrixTree.h"

void ECALMatrixTree::Init()
{
    GetTTreePtr()->Branch("index", index_, "index/l");
}
