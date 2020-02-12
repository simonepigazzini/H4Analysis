import array
import numpy as np

import ROOT
import operations

def ecal_time_res_vs_effampl(trees=None, dt_var=None, cut=None, amplbins=None, mcpres=None):
    """
    Compute time resolution of single ECAL channel comparing ECAL and MCP times. Inputs:
    [0] - h4 tree
    [1] - ECAL-MCP time difference vs A/n expression
    [2] - event selections
    [3] - binning, comma separated
    [3] - nx_quantiles,y-binning. comma separated
    [5] - MCP time resolution expression (optional)
    """

    #---inputs
    h4trees = trees

    #---helper histograms
    bins = [float(v) for v in amplbins.split(",")]
    dt = {}
    ampl = {}
    noise = {}    
    #h_An = ROOT.TH1F("h_An", "", int(bins[0]), bins[1], bins[2])
    res_mcp = {}

    #---event loop to fill the get dt, ECAL_ampl and MCP resolution
    for e, h4 in h4trees.items():
        dt_tmp = []
        ampl_tmp = []
        noise_tmp = []        
        res_mcp_tmp = []
        for ientry in range(h4.GetEntriesFast()):
            vars =  dt_var+":"+mcpres if mcpres != None else dt_var
            ret = h4.Draw(vars, cut, "goff", 1, ientry)        
            if ret > 0:
                dt_tmp.append(h4.GetVal(0)[0])
                ampl_tmp.append(h4.GetVal(1)[0])
                noise_tmp.append(h4.GetVal(2)[0])
                #h_An.Fill(h4.GetVal(1)[0])
                if mcpres != None:
                    res_mcp_tmp.append(h4.GetVal(3)[0])

        dt[e] = np.array(dt_tmp)
        ampl[e] = np.array(ampl_tmp)
        noise[e] = np.array(noise_tmp)        
        res_mcp[e] = np.array(res_mcp_tmp)        

    return ampl, noise,
#dt, res_mcp

# def quantile_binning(binning=[], xd=[], yd=[]):

#     # nqx = int(binning[0]+1)
#     # probs_x = array.array('d', [x/(nqx-1) for x in range(0, nqx)])
#     # quantiles_x = array.array('d', [0 for x in range(0, nqx)])
#     # histo.GetQuantiles(nqx, quantiles_x, probs_x)
#     # print(quantiles_x)
#     nqx = int(binning[0]+1)
#     quantiles_x = np.array([])
#     if energies == None:
#         xvals = np.sort(xd)
#         np.append(quantiles_x, np.percentile(xvals, [x/(nqx-1)*100 for x in range(0, nqx)]))
#     else:
#         nqx = int(int(binning[0])/len(energies)+1)
#         for i in range(len(energies)):
#             xvals = np.sort(xd[energies[i]])
#             quantiles_x = np.append(quantiles_x, np.percentile(xvals, [x/(nqx-1)*100 for x in range(0, nqx)]))
#             if i < len(energies)-1:
#                 quantiles_x = quantiles_x[:-1]
#     print(quantiles_x, len(energies))

#     ret = ROOT.TH2F("tmp", "", int(len(quantiles_x)-1), quantiles_x, int(binning[1]), binning[2], binning[3])
#     if energies == None:
#         for i in range(len(xd)):
#             ret.Fill(xd[i], yd[i])
#     else:
#         histo.Reset()
#         for energy in energies:
#             xvals, yvals = xd[energy], yd[energy]
#             print(xvals.min(), xvals.max())
#             for i in range(len(xvals)):
#                 ret.Fill(xvals[i], yvals[i])
#                 histo.Fill(xvals[i])
            
#     return ret

def quantile_binning(binning=[], xd=[], yd=[]):

    nqx = binning[0]+1
    quantiles_x = np.array([])
    quantiles_x = np.append(quantiles_x, np.percentile(xd, [x/(nqx-1)*100 for x in range(0, nqx)]))
    print(quantiles_x)
    
    ret = ROOT.TH2F("tmp", "", int(len(quantiles_x)-1), quantiles_x, int(binning[1]), binning[2], binning[3])
    for i in range(len(xd)):
        ret.Fill(xd[i], yd[i])
            
    return ret


def make_resolution_vs_ampl_graph(h_An=None, h_dt_vs_An=None):


    #---compute resolution and fill graph
    xval = array.array('d')
    yval = array.array('d')
    xerrl = array.array('d')
    xerrh = array.array('d')    
    yerr = array.array('d')
    h_res_tmp = []
    for ibin in range(1, h_dt_vs_An.GetNbinsX()+1):
        # x-axis
        axis = h_dt_vs_An.GetXaxis()
        h_An.GetXaxis().SetRange(h_An.FindBin(axis.GetBinCenter(ibin)-axis.GetBinWidth(ibin)/2),
                                 h_An.FindBin(axis.GetBinCenter(ibin)+axis.GetBinWidth(ibin)/2))
        xval.append(h_An.GetMean())
        xerrl.append(xval[-1]-(axis.GetBinCenter(ibin)-axis.GetBinWidth(ibin)/2))
        xerrh.append((axis.GetBinCenter(ibin)+axis.GetBinWidth(ibin)/2)-xval[-1])
        # y-axis
        h_res_tmp.append(h_dt_vs_An.ProjectionY("_py_"+str(ibin), ibin, ibin+1, "eo"))
        fcb = ROOT.TF1("fcb"+str(ibin), "crystalball", h_res_tmp[-1].GetMean()-5*h_res_tmp[-1].GetRMS(), h_res_tmp[-1].GetMean()+5*h_res_tmp[-1].GetRMS())
        fcb.SetParameters(1000., h_res_tmp[-1].GetMean(), h_res_tmp[-1].GetRMS(), 2, 0.1)
        fcb.SetParLimits(1, -2*h_res_tmp[-1].GetRMS(), 2*h_res_tmp[-1].GetRMS())
        fcb.SetParLimits(2, 0, h_res_tmp[-1].GetRMS()*2)
        h_res_tmp[-1].Fit(fcb, "BR")
        yval.append(fcb.GetParameter(2))
        yerr.append(fcb.GetParError(2))

    tmp = ROOT.TGraphAsymmErrors(len(xval), xval, yval, xerrl, xerrh, yerr, yerr)

    return tmp, h_res_tmp

def ECAL_energy_res_from_txt(srcs, tree=None, var=None, err=None):

    tree = srcs[tree]
    
    x = array.array('d')
    y = array.array('d')
    xerr = array.array('d')
    yerr = array.array('d')
    
    for i in range(tree.GetEntriesFast()):
        tree.Draw(var+":"+err, "", "goff", 1, i)
        x.append(tree.GetVal(0)[0])
        y.append(tree.GetVal(1)[0])
        xerr.append(tree.GetVal(2)[0])
        yerr.append(tree.GetVal(3)[0])

    tmp = ROOT.TGraphErrors(len(x), x, y, xerr, yerr)

    return tmp


###--------------------------------------------------------
dictionary = dict(ECALvsMCPtimeres=ecal_time_res_vs_effampl,
                  ECALEnergyResFromTxt=ECAL_energy_res_from_txt
)

