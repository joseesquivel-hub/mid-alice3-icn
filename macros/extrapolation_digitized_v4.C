void extrapolation_digitized_v4() {

  TFile* fmap = TFile::Open("/tmp/HitMap_O2.root");
  TGraph* gx = (TGraph*)fmap->Get("gx");
  TGraph* gy = (TGraph*)fmap->Get("gy");
  TGraph* gz = (TGraph*)fmap->Get("gz");
  if(!gx || !gy || !gz){ printf("ERROR: HitMap no cargado\n"); return; }

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

  // Residuos sin digitizacion
  TH1F* hDX_raw = new TH1F("hDX_raw","#DeltaX sin digit.;#DeltaX (cm);Entries",100,-10,10);
  TH1F* hDY_raw = new TH1F("hDY_raw","#DeltaY sin digit.;#DeltaY (cm);Entries",100,-10,10);
  TH1F* hDZ_raw = new TH1F("hDZ_raw","#DeltaZ sin digit.;#DeltaZ (cm);Entries",100,-10,10);
  // Residuos con digitizacion
  // DX,DY: barra L0 vs extrapolacion en L1
  // DZ: barra L1 vs extrapolacion en L1
  TH1F* hDX_dig = new TH1F("hDX_dig","#DeltaX con digit. (L0);#DeltaX (cm);Entries",100,-10,10);
  TH1F* hDY_dig = new TH1F("hDY_dig","#DeltaY con digit. (L0);#DeltaY (cm);Entries",100,-10,10);
  TH1F* hDZ_dig = new TH1F("hDZ_dig","#DeltaZ con digit. (L1);#DeltaZ (cm);Entries",100,-10,10);

  const double massparticle   = 0.10566;
  const double massparticleKG = 1.883531627e-28;
  const double clight         = 2.99792458e8;
  const double chargeparticle = -1.602176634e-19;
  const double Bz             = -1.619;

  int nMatched=0, nWindow_raw=0, nWindow_dig=0;
  int nNoBarL0=0, nNoBarL1=0;

  for(int iev=0; iev<thits->GetEntries(); iev++){
    thits->GetEntry(iev);
    tmc->GetEntry(iev);

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

    // Posicion MC en L0
    double x0_mc=h0.GetStartX(), y0_mc=h0.GetStartY(), z0_mc=h0.GetStartZ();
    double phi_hit0=TMath::ATan2(y0_mc,x0_mc);

    // Digitizacion L0: phi -> rango Z -> distancia XY
    double minXY0=1e9;
    int bestBar0=-1;
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
    if(bestBar0<0){ nNoBarL0++; continue; }

    // Centro de barra L0 digitalizada — X,Y de la barra, Z del hit MC
    double x0_dig = gx->GetY()[bestBar0];
    double y0_dig = gy->GetY()[bestBar0];
    double z0_dig = z0_mc; // Z no se digitaliza en L0 — viene del hit

    // Momento en L0
    double px_ini=h0.GetPx(), py_ini=h0.GetPy(), pz_ini=h0.GetPz();
    double p_tot=sqrt(px_ini*px_ini+py_ini*py_ini+pz_ini*pz_ini);
    double E=sqrt(p_tot*p_tot+massparticle*massparticle);
    double gamma=E/massparticle;
    double vx_0=px_ini*clight/(massparticle*gamma);
    double vy_0=py_ini*clight/(massparticle*gamma);
    double vz_0=pz_ini*clight/(massparticle*gamma);
    double omega_c=chargeparticle*Bz/(massparticleKG*gamma);

    // Extrapolacion desde posicion MC
    double x=x0_mc/100.,y=y0_mc/100.,z=z0_mc/100.,dt=0.;
    int step=0;
    while(sqrt(x*x+y*y)<=3.11 && step<100000){
      x=x0_mc/100.+(1./omega_c)*(vx_0*sin(omega_c*dt)+vy_0*(1.-cos(omega_c*dt)));
      y=y0_mc/100.+(1./omega_c)*(vx_0*(cos(omega_c*dt)-1.)+vy_0*sin(omega_c*dt));
      z=z0_mc/100.+vz_0*dt;
      dt+=(2.*TMath::Pi()/abs(omega_c))/7000.;
      step++; if(abs(z)>6.) break;
    }
    double xex_raw=x*100., yex_raw=y*100., zex_raw=z*100.;

    // Extrapolacion desde posicion digitalizada
    // X,Y del centro de barra L0; Z del hit MC
    x=x0_dig/100.; y=y0_dig/100.; z=z0_dig/100.; dt=0.; step=0;
    while(sqrt(x*x+y*y)<=3.11 && step<100000){
      x=x0_dig/100.+(1./omega_c)*(vx_0*sin(omega_c*dt)+vy_0*(1.-cos(omega_c*dt)));
      y=y0_dig/100.+(1./omega_c)*(vx_0*(cos(omega_c*dt)-1.)+vy_0*sin(omega_c*dt));
      z=z0_dig/100.+vz_0*dt;
      dt+=(2.*TMath::Pi()/abs(omega_c))/7000.;
      step++; if(abs(z)>6.) break;
    }
    double xex_dig=x*100., yex_dig=y*100., zex_dig=z*100.;

    // Hit mas cercano en L1 MC
    double minDist1=1e9, bx1=0, by1=0, bz1=0;
    for(int idx : hitsL1){
      auto& h1=(*hits)[idx];
      double x2=h1.GetStartX(), y2=h1.GetStartY(), z2=h1.GetStartZ();
      double d=sqrt((x2-xex_raw)*(x2-xex_raw)+(y2-yex_raw)*(y2-yex_raw)+(z2-zex_raw)*(z2-zex_raw));
      if(d<minDist1){ minDist1=d; bx1=x2; by1=y2; bz1=z2; }
    }
    if(minDist1>50.) continue;

    // Digitizacion L1: stave mas cercano en phi, luego barra mas cercana en Z
    double phi1=TMath::ATan2(by1,bx1);
    double minDphi=1e9;
    int bestStave=-1;
    for(int istave=0; istave<16; istave++){
      int idx=nBarsL0+istave*barsPerStaveL1;
      double bx=gx->GetY()[idx], by=gy->GetY()[idx];
      double phi_bar=TMath::ATan2(by,bx);
      double dphi=TMath::Abs(phi1-phi_bar);
      if(dphi>TMath::Pi()) dphi=TMath::TwoPi()-dphi;
      if(dphi<minDphi){ minDphi=dphi; bestStave=istave; }
    }
    if(bestStave<0){ nNoBarL1++; continue; }

    double minDZ=1e9;
    int bestBar1=-1;
    int staveStart=nBarsL0+bestStave*barsPerStaveL1;
    for(int ib=staveStart; ib<staveStart+barsPerStaveL1; ib++){
      double dZ=TMath::Abs(gz->GetY()[ib]-bz1);
      if(dZ<minDZ){ minDZ=dZ; bestBar1=ib; }
    }
    if(bestBar1<0){ nNoBarL1++; continue; }

    // Z digitalizada de L1
    double bz1_dig = gz->GetY()[bestBar1];

    // Residuos sin digitizacion
    double dX_raw=bx1-xex_raw, dY_raw=by1-yex_raw, dZ_raw=bz1-zex_raw;

    // Residuos con digitizacion:
    // DX,DY: efecto de digitalizar L0 en la extrapolacion a L1
    // DZ: efecto de digitalizar L1
    double dX_dig=bx1-xex_dig;   // X del hit MC vs extrapolacion desde barra L0
    double dY_dig=by1-yex_dig;   // Y del hit MC vs extrapolacion desde barra L0
    double dZ_dig=bz1_dig-zex_raw; // Z barra L1 vs extrapolacion desde MC

    hDX_raw->Fill(dX_raw); hDY_raw->Fill(dY_raw); hDZ_raw->Fill(dZ_raw);
    hDX_dig->Fill(dX_dig); hDY_dig->Fill(dY_dig); hDZ_dig->Fill(dZ_dig);

    if(abs(dX_raw)<5. && abs(dY_raw)<5. && abs(dZ_raw)<5.) nWindow_raw++;
    if(abs(dX_dig)<5. && abs(dY_dig)<5. && abs(dZ_dig)<5.) nWindow_dig++;
    nMatched++;
  }

  printf("\nEventos con coincidencia: %d\n", nMatched);
  printf("Sin barra en L0: %d\n", nNoBarL0);
  printf("Sin barra en L1: %d\n", nNoBarL1);
  printf("Dentro ventana +-5cm (sin digit.): %d (%.1f%%)\n", nWindow_raw, 100.*nWindow_raw/nMatched);
  printf("Dentro ventana +-5cm (con digit.): %d (%.1f%%)\n", nWindow_dig, 100.*nWindow_dig/nMatched);

  gStyle->SetOptStat(1);
  TCanvas* c = new TCanvas("c","Residuos digitizados v4",1800,600);
  c->Divide(3,2);
  c->cd(1); hDX_raw->SetLineColor(kBlue); hDX_raw->Draw();
  c->cd(2); hDY_raw->SetLineColor(kBlue); hDY_raw->Draw();
  c->cd(3); hDZ_raw->SetLineColor(kBlue); hDZ_raw->Draw();
  c->cd(4); hDX_dig->SetLineColor(kRed);  hDX_dig->Draw();
  c->cd(5); hDY_dig->SetLineColor(kRed);  hDY_dig->Draw();
  c->cd(6); hDZ_dig->SetLineColor(kRed);  hDZ_dig->Draw();
  c->SaveAs("/tmp/residuos_digitizados_v4.png");
  printf("Guardado: residuos_digitizados_v4.png\n");
}
