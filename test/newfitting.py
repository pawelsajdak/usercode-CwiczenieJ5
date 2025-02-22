#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

peakname = "Xk"
xmin = 4.0
xmax = 5.0
par0 = 600.
#axmin = 3.5
#axmax = 6.

##########################################
histfilename = "JxCorrN.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histoK")
histo.SetDirectory(0)
histfile.Close()
'''
expression = "[0]+x*[1]+x*x*[2]"
fitFunc = r.TF1("fitFunc",expression,xmin,xmax)
'''
expression = "[0]*exp((-(x-[1])**2)/(2*[2]**2)) + [3]+x*[4]+x*x*[5]"
fitFunc = r.TF1("fitFunc",expression,xmin,xmax)
fitFunc.SetParameters(par0,(xmin+xmax)/2,0.1,45808.7,-20854.5,2406.45)


results = histo.Fit(fitFunc,"ERS")
funcFile = r.TFile.Open("NKfuncs.root","RECREATE")
fitFunc.Write("Xfunc")
#funcFile.Close()

with open('NKresults.txt','a') as of:
    print(peakname,"\t",fitFunc.GetParameter(1),"\n", results, file=of)

canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetLogy(True)

#histo.SetAxisRange(axmin, axmax)
#histo.SetAxisRange(3.5, 6., "X")
#histo.SetAxisRange(2.e3, 7.e3, "Y")
histo.SetTitle(peakname+"\t {:.3f}".format(fitFunc.GetParameter(1))+"; Minv; #events")
histo.SetStats(0)
histo.Draw("h")


canvas.Print("NK_"+peakname+".pdf")
input('press enter to exit')