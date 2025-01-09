import ROOT as r

import numpy as np

class FitFunc:
    def __call__(self, arr, par):
        return par[0]+arr[0]*par[1]+arr[0]*par[2]**2 + par[3]*np.exp(-(arr[0]-par[4])**2/(2*par[5]**2))
    # polynomial background + Gaussian peak
    # par[4] - peak position

file = r.TFile("fullhistogram.root")
histo = file.Get('histo')

t = FitFunc()
func = r.TF1('func',t,0.,20.,6)
func.SetParameters(30000.,1.,1.,700000.,3.,0.01)

histo.Fit(func,"","",2.5,3.5)
print(histo.GetArray())

c = r.TCanvas()
histo.Draw()
c.Print('canvas2.pdf')
input('press enter to exit')