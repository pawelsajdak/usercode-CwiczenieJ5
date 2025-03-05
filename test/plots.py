#!/cvmfs/cms.cern.ch/slc7_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/slc7_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
fileName = "psi2S.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName)
f.ls()

c1 = TCanvas('cHisto','cHisto',600,600)
c1.SetLogy(1)
histo = gROOT.FindObject('histoPr')
histo.SetTitle("#psi(2S) + #it{p}^{#pm}; #it{M}_{inv} (GeV);Counts")
histo.GetXaxis().CenterTitle(True)
histo.SetAxisRange(4.5,11.)
histo.SetAxisRange(300.,1500.,"Y")
histo.Draw()
c1.Print("psi2SPr.pdf")
input('press enter to exit')
