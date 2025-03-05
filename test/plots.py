#!/cvmfs/cms.cern.ch/slc7_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/slc7_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
fileName = "probvBXhill.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName)
f.ls()

c1 = TCanvas('cHisto','cHisto',600,600)
#c1.SetLogy(1)
histo = gROOT.FindObject('histo_probvBX_comp')
histo.SetTitle("Common vertex of J/#psi and K^{#pm} (#it{M}_{#mu#muK} #in [4.8, 5.0]);Probability;Counts")
histo.GetXaxis().CenterTitle(True)
#histo.SetAxisRange(4.2,6.)
#histo.SetAxisRange(200.,1500.,"Y")
histo.Draw()
c1.Print("probvBX_comp.pdf")
input('press enter to exit')
