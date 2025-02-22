#!/cvmfs/cms.cern.ch/slc7_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/slc7_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
fileName = "fullHistoJxCorrected.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName)
f.ls()

c1 = TCanvas('cHisto','cHisto',600,600)
#c1.SetLogy(1)
histo = gROOT.FindObject('histoK')
histo.SetTitle("kaon")
#histo.SetAxisRange(3.,10.)
#histo.SetAxisRange(2.e3,8.e3,"Y")
histo.Draw()
c1.Print("oldJxCorrNK.pdf")
input('press enter to exit')
