#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys
import numpy as np

tupleFile = r.TFile("XX_B0s.root","READ")
tup = tupleFile.Get("tupKK")

cut_min = 5.39
cut_max = 5.47


histo = r.TH1D("histo","mJKK in ["+str(cut_min)+","+str(cut_max)+"]",500,0.9,2.5)
tup.Project("histo","mKK","mJKK > "+str(cut_min)+" && mJKK < "+str(cut_max))

canvas = r.TCanvas("canvas")
canvas.cd()
#canvas.SetLogy(1)
histo.Draw()

canvas.Print("temp.pdf")
input("press enter to exit")