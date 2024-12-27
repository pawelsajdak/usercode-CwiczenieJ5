#!/cvmfs/cms.cern.ch/slc7_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/slc7_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")

#first file
fileName = "histos_000.root"
print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName);
histo = gROOT.FindObject('histo')

#next files
for j in range(1, 216):
    fileName = 'histos_{:03d}.root'.format(j)
    print ('Read data from: ', fileName)
    gROOT.Reset()
    f = TFile(fileName);
    histo.Add(gROOT.FindObject('histo'))
    f.Close()

c1 = TCanvas('cHisto','cHisto',600,600)
histo.Draw()
c1.Print("histo.pdf")
input('press enter to exit')
