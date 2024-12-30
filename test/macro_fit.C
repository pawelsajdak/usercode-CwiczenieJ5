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

void macro_fit(){
    cout << "I am working" << endl;
}
