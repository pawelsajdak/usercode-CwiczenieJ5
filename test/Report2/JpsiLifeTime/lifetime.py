#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys
import numpy as np

tupleFile = r.TFile("jpsiLifetime_wdM.root","READ")
tLifetime = tupleFile.Get("tLifetime")
# crashes when one closes the file with TNtuple(D)

'''
tLifetime.Print()
histo = r.TH1D("histo","temp",100,0.0,1.0)
tLifetime.Project("histo","dR_min")
'''

# Cuts
dR_cut =            r.TCut("dR_min < 0.15")
deltaMJ_psi_cut =   r.TCut("deltaM < 0.02")
cuts = r.TCut(dR_cut + deltaMJ_psi_cut)
cuts.Print()

histo = r.TH1D("histo","temp",100,0.0,8.e-12)
#tLifetime.Project("histo","deltaMJ_psi")
tLifetime.Project("histo","properTime/(3*TMath::Power(10.,10))",str(cuts))

# File for fitting
outfile = r.TFile("hJpsiLT.root","RECREATE")
histo.Write()
dR_cut.Write()
deltaMJ_psi_cut.Write()


canvas = r.TCanvas("canvas")
canvas.cd()
#canvas.SetLogy(1)
histo.Draw()

canvas.Print("temp.pdf")
input("press enter to exit")