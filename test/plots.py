#!/cvmfs/cms.cern.ch/slc7_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/slc7_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
fileName = "BpmTime.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName)
f.ls()

c1 = TCanvas('cHisto','cHisto',600,600)
c1.SetLogy(1)
histo = gROOT.FindObject('hproperTime')
histo.SetTitle("Lifetime of B^{#pm};t;Counts")
histo.GetXaxis().CenterTitle(True)
histo.SetAxisRange(0.0,0.03)
#histo.SetAxisRange(25.e2,85.e2,"Y")
#histo.SetFillColor(19)
histo.Draw()
c1.Print("BpmTime.pdf")
input('press enter to exit')
