#!/cvmfs/cms.cern.ch/slc7_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/slc7_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
fileName = "Xs_wdRcheck.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName)
f.ls()

c1 = TCanvas('cHisto','cHisto',600,600)
#c1.SetLogy(1)
histo = gROOT.FindObject('histoPr')
histo.SetTitle("J/#psi + pr (dR>0.01) ;M_{inv};Counts")
histo.GetXaxis().CenterTitle(True)
histo.SetAxisRange(4.2,6.)
#histo.SetAxisRange(200.,1500.,"Y")
histo.Draw()
c1.Print("JPrdRcheck.pdf")
input('press enter to exit')
