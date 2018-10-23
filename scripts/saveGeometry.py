import ROOT
from parser_utils import *

if __name__ == '__main__':

    parser = argparse.ArgumentParser (description = 'Save Geometry into a separate root file')
    parser.add_argument('-i', '--input' , dest='input', default='', help='input file')
    parser.add_argument('-o', '--output' , dest='output', default='', help='output file')
    parser.add_argument('-n', '--name' , dest='name', default='', help='input object name')
    parser.add_argument('-w', '--write' , dest='write', default='', help='output object name')
    args = parser.parse_args ()

    ROOT.gSystem.Load("lib/libH4Analysis.so")

    _file0=ROOT.TFile.Open(args.input)
    tt=_file0.Get(args.name)
    tt.SetName(args.write)
    tt.Print()

    f=ROOT.TFile.Open(args.output,"RECREATE")
    f.cd()
    tt.Write(args.write)
    f.ls()
    f.Close()


