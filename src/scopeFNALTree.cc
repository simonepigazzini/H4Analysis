#include "interface/scopeFNALTree.h"

void scopeFNALTree::Init()
{
}

scopeFNALTree::~scopeFNALTree()
{
  delete[] channel;
  delete[] time;
}
