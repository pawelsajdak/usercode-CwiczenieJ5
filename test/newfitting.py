#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

peakname = "Bpm"
xmin = 4.9
xmax = 5.5
par0 = 100.
#axmin = 3.5
#axmax = 6.

##########################################
histfilename = "JxCorrN.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histoPi")
histo.SetDirectory(0)
histfile.Close()

expression = "[0]*exp((-(x-[1])**2)/(2*[2]**2)) + [3]*exp((-(x-[4])**2)/(2*[5]**2))+[6]"
fitFunc = r.TF1("fitFunc",expression,xmin,xmax)
fitFunc.SetParameters(par0,(xmin+xmax)/2,0.1,1.,1.,1.,1.)
#Background parameters (from "NPibgd.txt")
fitFunc.FixParameter(3,729052.)
fitFunc.FixParameter(4,-2.79952)
fitFunc.FixParameter(5,1.74473)
fitFunc.FixParameter(6,314.618)


results = histo.Fit(fitFunc,"ERSL")
funcFile = r.TFile.Open("NPifuncs.root","UPDATE")
fitFunc.Write("Bfunc")
#funcFile.Close()

with open('NPiresults.txt','a') as of:
    print(peakname,"\t",fitFunc.GetParameter(1),"\n", results, file=of)

canvas = r.TCanvas("canvas")
canvas.cd()
#canvas.SetLogy(True)

#histo.SetAxisRange(axmin, axmax)
#histo.SetAxisRange(3.5, 6., "X")
#histo.SetAxisRange(2.e3, 7.e3, "Y")
histo.SetTitle(peakname+"\t {:.3f}".format(fitFunc.GetParameter(1))+"; Minv; #events")
histo.SetStats(0)
histo.Draw("h")


canvas.Print("NPi_"+peakname+".pdf")
input('press enter to exit')