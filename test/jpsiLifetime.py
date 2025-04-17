#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys
import numpy as np

tupleFile = r.TFile("jpsiLifetime_wdM.root","READ")
tLifetime = tupleFile.Get("tLifetime")
# crashes when one closes the file with TNtuple(D)

histo = r.TH1D("histo","temp",1000,0.0,0.1)
tLifetime.Project("histo","properTime","deltaM > 0.04 && dR_min < 0.2")

canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetLogy(1)
histo.Draw()

canvas.Print("temp.pdf")
input("press enter to exit")