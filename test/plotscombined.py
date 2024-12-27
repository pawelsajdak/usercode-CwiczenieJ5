#!/cvmfs/cms.cern.ch/slc7_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/slc7_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
c1 = TCanvas('cHisto','cHisto',600,600)


for j in range(0, 216):
    fileName = 'histos_{:03d}.root'.format(j)
    print ('Read data from: ', fileName)
    gROOT.Reset()
    f = TFile(fileName);
    histo = gROOT.FindObject('histo')
    histo.Draw("same")
    
c1.Print("histo.pdf")
input('press enter to exit')
