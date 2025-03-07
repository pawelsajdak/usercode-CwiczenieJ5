#!/cvmfs/cms.cern.ch/slc7_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/slc7_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
fileName = "JXnewdR.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName)
f.ls()

c1 = TCanvas('cHisto','cHisto',600,600)
#c1.SetLogy(1)
histo = gROOT.FindObject('histoPr')
histo.SetTitle("J/#psi + #it{p}^{#pm} (#DeltaR > 0.0003); M_{inv};Counts")
histo.GetXaxis().CenterTitle(True)
histo.SetAxisRange(4.2,7.9)
#histo.SetAxisRange(0.,50.,"Y")
#histo.SetFillColor(19)
histo.Draw()
c1.Print("JPrnewdR.pdf")
input('press enter to exit')
