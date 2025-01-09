int macro_function () 
{
  TH1F h1 ("hist1", "", 50, -5., 5.);
  TH1F* h2 = new TH1F ("hist2", "", 50, -5., 5.);

  TRandom3 r;  r.SetSeed ();

  for (int i = 0; i < 1e5; i++) {
    h1.Fill  ( r.Gaus() );
    h2->Fill ( r.Gaus() );
  }

TCanvas can1 ("c1", "", 640, 480);

  h1.Draw();

can1.Update();
cin.ignore();

  return 0;
}
