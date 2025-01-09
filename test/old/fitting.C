#include "TMath.h"
//#include "TFile.h"
//#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"

#include <iostream>
using namespace std;

// Background
Double_t background(Double_t *x, Double_t *par){
    return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}

// Gaussian peak
Double_t peak(Double_t *x, Double_t *par){
    return par[0]*TMath::Gaus(x[0],par[1],par[2]);
}

// Fitting function
Double_t fitFunction(Double_t *x, Double_t *par){
    return background(x, par) + peak(x, &par[3]);
}

int fitting(TH1D* histo){
    Double_t par[6];
    TF1 *fitFcn = new TF1("fitFcn",fitFunction,6.,10.,6);
    fitFcn->SetParameter(1,9.5);
    histo->Fit("fitFcn");

    return 0;
}