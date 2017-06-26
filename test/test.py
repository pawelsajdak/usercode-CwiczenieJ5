#!/usr/bin/python

import sys
import math
from ROOT import *

print "Hello ROOT"
fileName = "test.root"

print 'Read data from: ', fileName
gROOT.Reset()
f = TFile(fileName);
f.ls();

c1 = TCanvas('cHisto','cHisto',600,600)
histo.Draw()
c1.Print("histo.pdf")
raw_input('press enter to exit')

