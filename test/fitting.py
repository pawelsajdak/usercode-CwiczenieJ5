import ROOT as r
import sys
import os

import numpy as np

class FitFunc:
    def __call__(self, arr, par):
        return par[0]+arr[0]*par[1]+arr[0]*par[2]**2 + par[3]*np.exp(-(arr[0]-par[4])**2/(2*par[5]**2))
    # polynomial background + Gaussian peak
    # par[4] - peak position

t = FitFunc()
func = r.TF1('func',t,-2.,3.,6)
func.SetParameters(5.,0.,0.,0.,5.,6.)

c = r.TCanvas()
func.Draw()
c.Print('canvas.pdf')

