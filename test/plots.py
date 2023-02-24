#!/cvmfs/cms.cern.ch/slc7_amd64_gcc10/cms/cmssw/CMSSW_12_6_3/external/slc7_amd64_gcc10/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
fileName = "histos.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName);
f.ls();

c1 = TCanvas('cHisto','cHisto',600,600)
histo = gROOT.FindObject('histo')
histo.Draw()
c1.Print("histo.pdf")
input('press enter to exit')
