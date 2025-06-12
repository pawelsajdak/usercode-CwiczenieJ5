#!/cvmfs/cms.cern.ch/slc7_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/slc7_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
fileName = "Jpsi_P_K_merged.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName)
f.ls()

c1 = TCanvas('cHisto','cHisto',600,600)
#c1.SetLogy(1)
ihisto = gROOT.FindObject('histoPi')

histo = ihisto.Rebin(20,"histo")

#histo.SetTitle("PP")
histo.GetXaxis().CenterTitle(True)
histo.SetAxisRange(4.5,6.0)   # 601-1001    10 000 bins overall


#histo.SetAxisRange(1.e3,10.e3,"Y")
#histo.SetFillColor(19)
histo.Draw()
c1.Print("temp.pdf")
input('press enter to exit')
