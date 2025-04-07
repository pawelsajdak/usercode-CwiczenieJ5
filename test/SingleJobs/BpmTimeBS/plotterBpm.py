#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys
import numpy as np


histo = r.TH1D("histo","dist Bpm",1000,0.0,0.2)

data = np.loadtxt('values2.txt',delimiter=' ')
data = data[:5282][:] #dR<0.02

N = data.shape[0]

weights = np.full(N,1.0)

histo.FillN(N,data[:][1],weights)

canvas = r.TCanvas("canvas")
canvas.cd()
histo.Draw("h")
canvas.Print("dist02.pdf")

input('press enter to exit')
