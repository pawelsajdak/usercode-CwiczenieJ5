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

int macro_fit(Double_t xmin, Double_t xmax){
    cout << "I am working" << endl;
    TFile f("fullhistogram.root", "UPDATE");
    cout << "file opened" << endl;
    TH1D* myhisto = (TH1D*)f.Get("histo");
    Double_t par[6];
    TF1* fitFun = new TF1("fitFun",fitFunction,xmin,xmax,6);
    fitFun->SetParameter(4,3.0);
    myhisto->Fit("fitFun");
    return 0;
}