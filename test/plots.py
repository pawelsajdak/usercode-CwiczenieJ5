#!/cvmfs/cms.cern.ch/slc7_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/slc7_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
fileName = "dhisto.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName);
f.ls();

c1 = TCanvas('cHisto','cHisto',600,600)
histo = gROOT.FindObject('histo')
histo.Draw()
c1.Print("dhisto.pdf")
input('press enter to exit')
