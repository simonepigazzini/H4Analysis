#!/bin/python
# import ROOT in batch mode
import sys
import re
import argparse
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

class customAction(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        setattr(args, self.dest, values.split(','))

def get_comma_separated_args(self, arg_line):
    return arg_line.split(',')

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

# Load the CalibrationFile format
ROOT.gSystem.Load("lib/H4Dict.so")

trees={}
mergedFile = ROOT.TFile(sys.argv[1], "RECREATE")
headFile = ROOT.TFile(sys.argv[2], "READ")
headFile.cd()
trees["h4"]=ROOT.TChain("h4", "main_h4reco_tree")        
nextKey = ROOT.TIter(ROOT.gDirectory.GetListOfKeys())
    
while True:
    key = nextKey()
    if key == nextKey.End()():
        break
    ## keep only the highest cycle number for each key and only TTree keys    
    obj = key.ReadObj()
    if obj.IsA().GetName() == "TTree":
        trees[key.GetName()]= ROOT.TChain(key.GetName(), key.GetTitle())
        trees[key.GetName()].SetDirectory(0)
    else:
        mergedFile.cd()
        obj.Write(key.GetName())        

for filename in sys.argv[2:]:
    for tree_name in trees:
        trees[tree_name].AddFile(filename)

for tree_name in trees:
    if tree_name == "h4" or tree_name == "wf":
        continue
    mergedFile.cd()
    trees[tree_name].Merge(mergedFile, 0, "keep")
    trees[tree_name].BuildIndex("index")
    trees["h4"].AddFriend(tree_name)

mergedFile.cd()
trees["h4"].BuildIndex("index")
trees["h4"].Merge(mergedFile, 0, "keep")

mergedFile.Close()
