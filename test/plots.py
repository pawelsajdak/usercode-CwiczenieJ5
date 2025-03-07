#!/cvmfs/cms.cern.ch/slc7_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/slc7_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
fileName = "dRdistribution.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName)
f.ls()

c1 = TCanvas('cHisto','cHisto',600,600)
#c1.SetLogy(1)
histo = gROOT.FindObject('histodR')
histo.SetTitle("J/#psi + #it{K}^{#pm}\t M_{mmK} #in [3.8, 6.0]; min(#DeltaR);Counts")
histo.GetXaxis().CenterTitle(True)
histo.SetAxisRange(0.,0.0008)
#histo.SetAxisRange(0.,50.,"Y")
histo.SetFillColor(19)
histo.Draw()
c1.Print("dRdistribution_0008F.pdf")
input('press enter to exit')
