#!/cvmfs/cms.cern.ch/slc7_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/slc7_amd64_gcc12/bin/python3

import sys
import math
import ROOT as r


print ("Hello ROOT")

fileName = "histos_000.root"
histo = r.TH1D("histo")
r.gROOT.Reset()
f = r.TFile(fileName)
histo = r.TH1D(r.gROOT.FindObject('histo'))


for j in range(1, 216):
    fileName = 'histos_{:03d}.root'.format(j)
    print ('Read data from: ', fileName)
    r.gROOT.Reset()
    f = r.TFile(fileName);
    histo2 = r.TH1D(r.gROOT.FindObject('histo'))
    histo.Add(histo2)


c1 = r.TCanvas('cHisto','cHisto',600,600)
histo.Draw()
c1.Print("histo.pdf")
input('press enter to exit')
