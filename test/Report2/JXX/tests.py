#!/cvmfs/cms.cern.ch/slc7_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/slc7_amd64_gcc12/bin/python3

import sys
import math
from ROOT import *


print ("Hello ROOT")
fileName = "3partMinv.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName)
f.ls()

c1 = TCanvas('cHisto','cHisto',600,600)
#c1.SetLogy(1)
histo = gROOT.FindObject('histoK')


rbhisto = histo.Rebin(5,"rbhisto")

print(histo.GetNbinsX())
print(histo.GetBinContent(5000))
print(histo.GetBinError(5000))
print(TMath.Sqrt(histo.GetBinContent(5000)))
print("Rebinned histo")
print(rbhisto.GetNbinsX())
print(rbhisto.GetBinContent(5000))
print(rbhisto.GetBinError(5000))
print(TMath.Sqrt(rbhisto.GetBinContent(5000)))


histo.Draw()
c1.Print("temp.pdf")
input('press enter to exit')
