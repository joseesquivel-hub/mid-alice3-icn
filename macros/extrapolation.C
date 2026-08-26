void extrapolation() {

  TFile* fhits = TFile::Open("/tmp/SimALICE3/o2sim_HitsMI3.root");
  TFile* fmc   = TFile::Open("/tmp/SimALICE3/o2sim_Kine.root");
  TTree* thits = (TTree*)fhits->Get("o2sim");
  TTree* tmc   = (TTree*)fmc->Get("o2sim");

  std::vector<o2::MCTrack>* tracks = nullptr;
  std::vector<o2::itsmft::Hit>* hits = nullptr;
  tmc->SetBranchAddress("MCTrack", &tracks);
  thits->SetBranchAddress("MI3Hit", &hits);

  TH1F* hDX = new TH1F("hDX","#DeltaX;#DeltaX (cm);Entries",100,-10,10);
  TH1F* hDY = new TH1F("hDY","#DeltaY;#DeltaY (cm);Entries",100,-10,10);
  TH1F* hDZ = new TH1F("hDZ","#DeltaZ;#DeltaZ (cm);Entries",100,-10,10);
  TH1F* hDR = new TH1F("hDR","#DeltaR 3D;#DeltaR (cm);Entries",100,0,15);

  const double massparticle    = 0.10566;         // GeV/c2
  const double massparticleKG  = 1.883531627e-28; // kg
  const double clight          = 2.99792458e8;    // m/s
  const double chargeparticle  = -1.602176634e-19;// C (mu-)
  const double Bz              = -1.619;          // T (B2 en absorbedor)

  int nMatched = 0;

  for(int iev=0; iev<thits->GetEntries(); iev++){
    thits->GetEntry(iev);
    tmc->GetEntry(iev);

    std::vector<int> hitsL0, hitsL1;
    for(int i=0; i<(int)hits->size(); i++){
      int layer = (*hits)[i].GetDetectorID();
      if(layer==0) hitsL0.push_back(i);
      else if(layer==1) hitsL1.push_back(i);
    }
    if(hitsL0.empty() || hitsL1.empty()) continue;

    auto& h0 = (*hits)[hitsL0[0]];
    int tid = h0.GetTrackID();
    if(tid<0 || tid>=(int)tracks->size()) continue;
    auto& tr = (*tracks)[tid];
    if(abs(tr.GetPdgCode())!=13) continue;

    // Posicion en L0 — metros
    double xini = h0.GetStartX() / 100.;
    double yini = h0.GetStartY() / 100.;
    double zini = h0.GetStartZ() / 100.;

    // Momento en el hit L0 — GeV/c (punto de desconfianza: ver nota)
    double px_ini = h0.GetPx();
    double py_ini = h0.GetPy();
    double pz_ini = h0.GetPz();
    double p_tot  = sqrt(px_ini*px_ini + py_ini*py_ini + pz_ini*pz_ini);
    double E      = sqrt(p_tot*p_tot + massparticle*massparticle);
    double gamma  = E / massparticle;

    double vx_0   = px_ini * clight / (massparticle * gamma);
    double vy_0   = py_ini * clight / (massparticle * gamma);
    double vz_0   = pz_ini * clight / (massparticle * gamma);
    double omega_c = chargeparticle * Bz / (massparticleKG * gamma);

    // Integrar hasta R = 3.11 m (L1 = R=311 cm)
    double x=xini, y=yini, z=zini, dt=0.;
    int step=0;
    while(sqrt(x*x+y*y) <= 3.11 && step < 100000){
      x = xini + (1./omega_c)*(vx_0*sin(omega_c*dt) + vy_0*(1.-cos(omega_c*dt)));
      y = yini + (1./omega_c)*(vx_0*(cos(omega_c*dt)-1.) + vy_0*sin(omega_c*dt));
      z = zini + vz_0*dt;
      dt += (2.*TMath::Pi()/abs(omega_c))/7000.;
      step++;
      if(abs(z)>6.) break;
    }

    double xex=x*100., yex=y*100., zex=z*100.;

    // Hit mas cercano en L1
    double minDist=1e9, bx=0, by=0, bz=0;
    for(int idx : hitsL1){
      auto& h1 = (*hits)[idx];
      double x2=h1.GetStartX(), y2=h1.GetStartY(), z2=h1.GetStartZ();
      double d=sqrt((x2-xex)*(x2-xex)+(y2-yex)*(y2-yex)+(z2-zex)*(z2-zex));
      if(d<minDist){ minDist=d; bx=x2; by=y2; bz=z2; }
    }
    if(minDist>50.) continue;

    hDX->Fill(bx-xex);
    hDY->Fill(by-yex);
    hDZ->Fill(bz-zex);
    hDR->Fill(minDist);
    nMatched++;
  }

  printf("Coincidencias L0-L1: %d / %d eventos\n", nMatched, (int)thits->GetEntries());

  gStyle->SetOptStat(1);
  TCanvas* c = new TCanvas("c","Residuos MID",1200,400);
  c->Divide(4,1);
  c->cd(1); hDX->Draw();
  c->cd(2); hDY->Draw();
  c->cd(3); hDZ->Draw();
  c->cd(4); hDR->Draw();
  c->SaveAs("/tmp/residuos_extrapolacion.png");
  printf("Guardado: /tmp/residuos_extrapolacion.png\n");
}
