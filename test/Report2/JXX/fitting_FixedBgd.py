#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

histname = "histoPi20"
xmin = 4.6
xmax = 5.8
##########################################
histfilename = "rebin3partMinv.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get(histname)
histo.SetDirectory(0)
histfile.Close()

# Background function with fitted parameters
bgdfilename = histname+"_bgdFunc.root"
bgdfile = r.TFile.Open(bgdfilename)
bgd = bgdfile.Get("bgd")

# Fitting function
expression = "[0]+[1]*x+[2]*x*x + [3]*exp((-(x-[4])**2)/(2*[5]**2))"# + [6]*exp((-(x-[7])**2)/(2*[8]**2))"
fitFunc = r.TF1("fitFunc",expression,xmin,xmax,6)
fitFunc.SetParameters(1.,1.,1.,1000.,5.15,0.1)
fitFunc.FixParameter(0,bgd.GetParameter(0))
fitFunc.FixParameter(1,bgd.GetParameter(1))
fitFunc.FixParameter(2,bgd.GetParameter(2))

results = histo.Fit(fitFunc,"ERSLB")

'''
funcFile = r.TFile.Open(histname+"_fitFunc.root","UPDATE")
fitFunc.Write()
funcFile.Close()

with open(histname+'_fitResults.txt','a') as of:
    print(results, file=of)
'''
    
canvas = r.TCanvas("canvas")
canvas.cd()
#canvas.SetLogy(True)

histo.SetAxisRange(4.0,7.0)
#histo.SetAxisRange(3.5, 6., "X")
#histo.SetAxisRange(1500, 3.e3, "Y")
#histo.SetTitle("Lifetime of B^{#pm};t;Counts")
histo.SetStats(0)
histo.Draw("h")
fitFunc.Draw("same")


canvas.Print(histname+"_Fit.pdf")
input('press enter to exit')
