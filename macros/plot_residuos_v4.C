void plot_residuos_v4() {

  TFile* fmap = TFile::Open("/tmp/HitMap_O2.root");
  TGraph* gx = (TGraph*)fmap->Get("gx");
  TGraph* gy = (TGraph*)fmap->Get("gy");
  TGraph* gz = (TGraph*)fmap->Get("gz");

  const int   nBarsL0        = 3680;
  const int   nBarsL1        = 3200;
  const int   barsPerStaveL1 = 200;
  const float barLenL0       = 49.9f;

  TFile* fhits = TFile::Open("/tmp/SimALICE3/o2sim_HitsMI3.root");
  TFile* fmc   = TFile::Open("/tmp/SimALICE3/o2sim_Kine.root");
  TTree* thits = (TTree*)fhits->Get("o2sim");
  TTree* tmc   = (TTree*)fmc->Get("o2sim");

  std::vector<o2::MCTrack>* tracks = nullptr;
  std::vector<o2::itsmft::Hit>* hits = nullptr;
  tmc->SetBranchAddress("MCTrack", &tracks);
  thits->SetBranchAddress("MI3Hit", &hits);

  TH1F* hDX_raw = new TH1F("hDX_raw","",100,-10,10);
  TH1F* hDY_raw = new TH1F("hDY_raw","",100,-10,10);
  TH1F* hDZ_raw = new TH1F("hDZ_raw","",100,-10,10);
  TH1F* hDX_dig = new TH1F("hDX_dig","",100,-10,10);
  TH1F* hDY_dig = new TH1F("hDY_dig","",100,-10,10);
  TH1F* hDZ_dig = new TH1F("hDZ_dig","",100,-10,10);

  const double massparticle   = 0.10566;
  const double massparticleKG = 1.883531627e-28;
  const double clight         = 2.99792458e8;
  const double chargeparticle = -1.602176634e-19;
  const double Bz             = -1.619;

  int nMatched=0, nWindow_raw=0, nWindow_dig=0;
  double sumX_raw=0,sumY_raw=0,sumZ_raw=0;
  double sumX_dig=0,sumY_dig=0,sumZ_dig=0;

  for(int iev=0; iev<thits->GetEntries(); iev++){
    thits->GetEntry(iev); tmc->GetEntry(iev);

    std::vector<int> hitsL0, hitsL1;
    for(int i=0; i<(int)hits->size(); i++){
      int layer=(*hits)[i].GetDetectorID();
      if(layer==0) hitsL0.push_back(i);
      else if(layer==1) hitsL1.push_back(i);
    }
    if(hitsL0.empty() || hitsL1.empty()) continue;

    auto& h0=(*hits)[hitsL0[0]];
    int tid=h0.GetTrackID();
    if(tid<0 || tid>=(int)tracks->size()) continue;
    auto& tr=(*tracks)[tid];
    if(abs(tr.GetPdgCode())!=13) continue;

    double x0_mc=h0.GetStartX(), y0_mc=h0.GetStartY(), z0_mc=h0.GetStartZ();
    double phi_hit0=TMath::ATan2(y0_mc,x0_mc);

    double minXY0=1e9; int bestBar0=-1;
    for(int ib=0; ib<nBarsL0; ib++){
      double bx=gx->GetY()[ib], by=gy->GetY()[ib], bz=gz->GetY()[ib];
      double phi_bar=TMath::ATan2(by,bx);
      double dphi=TMath::Abs(phi_hit0-phi_bar);
      if(dphi>TMath::Pi()) dphi=TMath::TwoPi()-dphi;
      if(dphi>TMath::Pi()/6.) continue;
      if(z0_mc < bz-barLenL0 || z0_mc > bz+barLenL0) continue;
      double dXY=sqrt((bx-x0_mc)*(bx-x0_mc)+(by-y0_mc)*(by-y0_mc));
      if(dXY<minXY0){ minXY0=dXY; bestBar0=ib; }
    }
    if(bestBar0<0) continue;

    double x0_dig=gx->GetY()[bestBar0];
    double y0_dig=gy->GetY()[bestBar0];

    double px_ini=h0.GetPx(), py_ini=h0.GetPy(), pz_ini=h0.GetPz();
    double p_tot=sqrt(px_ini*px_ini+py_ini*py_ini+pz_ini*pz_ini);
    double E=sqrt(p_tot*p_tot+massparticle*massparticle);
    double gamma=E/massparticle;
    double vx_0=px_ini*clight/(massparticle*gamma);
    double vy_0=py_ini*clight/(massparticle*gamma);
    double vz_0=pz_ini*clight/(massparticle*gamma);
    double omega_c=chargeparticle*Bz/(massparticleKG*gamma);

    double x=x0_mc/100.,y=y0_mc/100.,z=z0_mc/100.,dt=0.; int step=0;
    while(sqrt(x*x+y*y)<=3.11 && step<100000){
      x=x0_mc/100.+(1./omega_c)*(vx_0*sin(omega_c*dt)+vy_0*(1.-cos(omega_c*dt)));
      y=y0_mc/100.+(1./omega_c)*(vx_0*(cos(omega_c*dt)-1.)+vy_0*sin(omega_c*dt));
      z=z0_mc/100.+vz_0*dt;
      dt+=(2.*TMath::Pi()/abs(omega_c))/7000.;
      step++; if(abs(z)>6.) break;
    }
    double xex_raw=x*100., yex_raw=y*100., zex_raw=z*100.;

    x=x0_dig/100.; y=y0_dig/100.; z=z0_mc/100.; dt=0.; step=0;
    while(sqrt(x*x+y*y)<=3.11 && step<100000){
      x=x0_dig/100.+(1./omega_c)*(vx_0*sin(omega_c*dt)+vy_0*(1.-cos(omega_c*dt)));
      y=y0_dig/100.+(1./omega_c)*(vx_0*(cos(omega_c*dt)-1.)+vy_0*sin(omega_c*dt));
      z=z0_mc/100.+vz_0*dt;
      dt+=(2.*TMath::Pi()/abs(omega_c))/7000.;
      step++; if(abs(z)>6.) break;
    }
    double xex_dig=x*100., yex_dig=y*100.;

    double minDist1=1e9, bx1=0, by1=0, bz1=0;
    for(int idx : hitsL1){
      auto& h1=(*hits)[idx];
      double x2=h1.GetStartX(), y2=h1.GetStartY(), z2=h1.GetStartZ();
      double d=sqrt((x2-xex_raw)*(x2-xex_raw)+(y2-yex_raw)*(y2-yex_raw)+(z2-zex_raw)*(z2-zex_raw));
      if(d<minDist1){ minDist1=d; bx1=x2; by1=y2; bz1=z2; }
    }
    if(minDist1>50.) continue;

    double phi1=TMath::ATan2(by1,bx1);
    double minDphi=1e9; int bestStave=-1;
    for(int istave=0; istave<16; istave++){
      int idx=nBarsL0+istave*barsPerStaveL1;
      double bx=gx->GetY()[idx], by=gy->GetY()[idx];
      double phi_bar=TMath::ATan2(by,bx);
      double dphi=TMath::Abs(phi1-phi_bar);
      if(dphi>TMath::Pi()) dphi=TMath::TwoPi()-dphi;
      if(dphi<minDphi){ minDphi=dphi; bestStave=istave; }
    }
    if(bestStave<0) continue;

    double minDZ=1e9; int bestBar1=-1;
    int staveStart=nBarsL0+bestStave*barsPerStaveL1;
    for(int ib=staveStart; ib<staveStart+barsPerStaveL1; ib++){
      double dZ=TMath::Abs(gz->GetY()[ib]-bz1);
      if(dZ<minDZ){ minDZ=dZ; bestBar1=ib; }
    }
    if(bestBar1<0) continue;

    double bz1_dig=gz->GetY()[bestBar1];
    double dX_raw=bx1-xex_raw, dY_raw=by1-yex_raw, dZ_raw=bz1-zex_raw;
    double dX_dig=bx1-xex_dig, dY_dig=by1-yex_dig, dZ_dig=bz1_dig-zex_raw;

    hDX_raw->Fill(dX_raw); hDY_raw->Fill(dY_raw); hDZ_raw->Fill(dZ_raw);
    hDX_dig->Fill(dX_dig); hDY_dig->Fill(dY_dig); hDZ_dig->Fill(dZ_dig);
    sumX_raw+=dX_raw; sumY_raw+=dY_raw; sumZ_raw+=dZ_raw;
    sumX_dig+=dX_dig; sumY_dig+=dY_dig; sumZ_dig+=dZ_dig;

    if(abs(dX_raw)<5. && abs(dY_raw)<5. && abs(dZ_raw)<5.) nWindow_raw++;
    if(abs(dX_dig)<5. && abs(dY_dig)<5. && abs(dZ_dig)<5.) nWindow_dig++;
    nMatched++;
  }

  double n=nMatched;
  printf("Coincidencias: %d\n",nMatched);
  printf("Dentro +-5cm (3D) sin digit.: %d (%.1f%%)\n",nWindow_raw,100.*nWindow_raw/n);
  printf("Dentro +-5cm (3D) con digit.: %d (%.1f%%)\n",nWindow_dig,100.*nWindow_dig/n);

  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(2);

  hDX_raw->SetLineColor(kAzure+2); hDX_raw->SetLineWidth(2); hDX_raw->SetFillColorAlpha(kAzure+2,0.15);
  hDY_raw->SetLineColor(kAzure+2); hDY_raw->SetLineWidth(2); hDY_raw->SetFillColorAlpha(kAzure+2,0.15);
  hDZ_raw->SetLineColor(kAzure+2); hDZ_raw->SetLineWidth(2); hDZ_raw->SetFillColorAlpha(kAzure+2,0.15);
  hDX_dig->SetLineColor(kRed+1);   hDX_dig->SetLineWidth(2); hDX_dig->SetFillColorAlpha(kRed+1,0.15);
  hDY_dig->SetLineColor(kRed+1);   hDY_dig->SetLineWidth(2); hDY_dig->SetFillColorAlpha(kRed+1,0.15);
  hDZ_dig->SetLineColor(kRed+1);   hDZ_dig->SetLineWidth(2); hDZ_dig->SetFillColorAlpha(kRed+1,0.15);

  // Canvas con titulo general
  TCanvas* c = new TCanvas("c","",1600,620);
  c->Divide(3,1);

  // Titulo general
  TPaveText* title = new TPaveText(0.01,0.94,0.99,0.99,"NDC");
  title->SetBorderSize(0); title->SetFillStyle(0);
  title->SetTextSize(0.045); title->SetTextFont(62);
  title->AddText(Form("MID Track Extrapolation: digitizing costs 1.8 pp (%.1f%% to %.1f%%) -- #pm5 cm window survives",
    100.*nWindow_raw/n, 100.*nWindow_dig/n));
  title->Draw();

  auto makeLine = [](double xval, double ymax) -> TLine* {
    TLine* l = new TLine(xval,0,xval,ymax);
    l->SetLineColor(kGreen+2); l->SetLineWidth(2); l->SetLineStyle(2);
    return l;
  };

  // Leyenda corta — solo etiquetas
  TLegend* leg = new TLegend(0.10,0.75,0.92,0.90);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.042);
  leg->AddEntry(hDX_raw,"MC truth (exact hit position)","lf");
  leg->AddEntry(hDX_dig,"Digitized: L0 bar center (X,Y) + L1 bar center (Z)","lf");
  leg->SetNColumns(2);

  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.040);

  // Panel DX
  c->cd(1); gPad->SetLeftMargin(0.14); gPad->SetTopMargin(0.08);
  double ymx=TMath::Max(hDX_raw->GetMaximum(),hDX_dig->GetMaximum())*1.50;
  hDX_raw->SetTitle(";#DeltaX = X^{L1}_{hit} - X_{extrap} (cm);Entries / 0.2 cm");
  hDX_raw->GetYaxis()->SetRangeUser(0,ymx);
  hDX_raw->Draw("hist"); hDX_dig->Draw("hist same");
  makeLine(-5,ymx)->Draw(); makeLine(5,ymx)->Draw();
  tx.DrawLatex(0.16,0.88,Form("Mean: %.3f / %.3f cm",sumX_raw/n,sumX_dig/n));
  tx.DrawLatex(0.16,0.82,Form("RMS:  %.2f / %.2f cm",hDX_raw->GetRMS(),hDX_dig->GetRMS()));
  tx.DrawLatex(0.16,0.74,"995 #mu^{-}, 5-50 GeV, |#eta|<1.2");
  tx.DrawLatex(0.16,0.68,"Bar pitch = 5.2 cm");
  tx.DrawLatex(0.16,0.62,Form("Eff: %.1f%% / %.1f%%",100.*nWindow_raw/n,100.*nWindow_dig/n));
  leg->Draw();

  // Panel DY
  c->cd(2); gPad->SetLeftMargin(0.14); gPad->SetTopMargin(0.08);
  ymx=TMath::Max(hDY_raw->GetMaximum(),hDY_dig->GetMaximum())*1.50;
  hDY_raw->SetTitle(";#DeltaY = Y^{L1}_{hit} - Y_{extrap} (cm);Entries / 0.2 cm");
  hDY_raw->GetYaxis()->SetRangeUser(0,ymx);
  hDY_raw->Draw("hist"); hDY_dig->Draw("hist same");
  makeLine(-5,ymx)->Draw(); makeLine(5,ymx)->Draw();
  tx.DrawLatex(0.16,0.88,Form("Mean: %.3f / %.3f cm",sumY_raw/n,sumY_dig/n));
  tx.DrawLatex(0.16,0.82,Form("RMS:  %.2f / %.2f cm",hDY_raw->GetRMS(),hDY_dig->GetRMS()));
  tx.DrawLatex(0.16,0.74,"995 #mu^{-}, 5-50 GeV, |#eta|<1.2");
  tx.DrawLatex(0.16,0.68,"Bar pitch = 5.2 cm");
  tx.DrawLatex(0.16,0.62,Form("Eff: %.1f%% / %.1f%%",100.*nWindow_raw/n,100.*nWindow_dig/n));
  leg->Draw();

  // Panel DZ
  c->cd(3); gPad->SetLeftMargin(0.14); gPad->SetTopMargin(0.08);
  ymx=TMath::Max(hDZ_raw->GetMaximum(),hDZ_dig->GetMaximum())*1.50;
  hDZ_raw->SetTitle(";#DeltaZ = Z^{L1}_{bar} - Z_{extrap} (cm);Entries / 0.2 cm");
  hDZ_raw->GetYaxis()->SetRangeUser(0,ymx);
  hDZ_raw->Draw("hist"); hDZ_dig->Draw("hist same");
  makeLine(-5,ymx)->Draw(); makeLine(5,ymx)->Draw();
  tx.DrawLatex(0.16,0.88,Form("Mean: %.3f / %.3f cm",sumZ_raw/n,sumZ_dig/n));
  tx.DrawLatex(0.16,0.82,Form("RMS:  %.2f / %.2f cm",hDZ_raw->GetRMS(),hDZ_dig->GetRMS()));
  tx.DrawLatex(0.16,0.74,"995 #mu^{-}, 5-50 GeV, |#eta|<1.2");
  tx.DrawLatex(0.16,0.68,"Bar pitch = 5.2 cm");
  tx.DrawLatex(0.16,0.62,Form("Eff: %.1f%% / %.1f%%",100.*nWindow_raw/n,100.*nWindow_dig/n));
  leg->Draw();

  c->SaveAs("/tmp/residuos_final_v4.png");
  printf("Guardado: /tmp/residuos_final_v4.png\n");
}
