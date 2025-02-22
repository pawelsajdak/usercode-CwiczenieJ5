#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

jpsiMass = 3.096900

histfilename = "fullhistogram.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histo")
histo.SetDirectory(0)
histfile.Close()

canvas = r.TCanvas("canvas")
canvas.cd()

histo.SetAxisRange(2.8, 3.4, "X")

histo.Draw()
canvas.Print("dev.pdf")
input('press enter to exit')

# 3.06 - 3.13 FWHM