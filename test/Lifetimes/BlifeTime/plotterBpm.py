#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys
import numpy as np


histo = r.TH1D("histo","Bpm lifetime",200,0.0,0.06)

data = np.loadtxt('times.txt')
print(data.shape)
print(data)
hdata = (1/3)*data[:17764] #from the cut on dR_min

print(hdata)

N = len(hdata)
print("N: ",N)
weights = np.full(N,1.0)

histo.FillN(N,hdata,weights)

outfile = r.TFile("BpmLfTm.root","recreate")
histo.Write()
outfile.Close()

canvas = r.TCanvas("canvas")
canvas.cd()
histo.Draw("h")
canvas.Print("lifetime.pdf")

input('press enter to exit')
