#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys
import numpy as np

tupleFile = r.TFile("jpsiLifetime.root","READ")
tLifetime = tupleFile.Get("tLifetime")
# crashes when one closes the file with TNtuple(D)

histo = r.TH1D("histo","dR_min < 0.15;properTime;Counts",1000,0.0,0.3)
tLifetime.Project("histo","properTime","dR_min < 0.15")

canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetLogy(1)
histo.Draw()

canvas.Print("temp.pdf")
input("press enter to exit")