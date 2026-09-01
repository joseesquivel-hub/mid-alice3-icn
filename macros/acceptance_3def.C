void acceptance_3def() {

  TFile* fhits = TFile::Open("/tmp/SimALICE3/o2sim_HitsMI3.root");
  TFile* fmc   = TFile::Open("/tmp/SimALICE3/o2sim_Kine.root");
  TTree* thits = (TTree*)fhits->Get("o2sim");
  TTree* tmc   = (TTree*)fmc->Get("o2sim");

  std::vector<o2::MCTrack>* tracks = nullptr;
  std::vector<o2::itsmft::Hit>* hits = nullptr;
  tmc->SetBranchAddress("MCTrack", &tracks);
  thits->SetBranchAddress("MI3Hit", &hits);

  // Histogramas — tres definiciones
  TH1F* hGen  = new TH1F("hGen", "",30,-1.5,1.5);  // generados
  TH1F* hDef1 = new TH1F("hDef1","",30,-1.5,1.5);  // hit en L0 o L1
  TH1F* hDef2 = new TH1F("hDef2","",30,-1.5,1.5);  // hit en L0 Y L1
  TH1F* hDef3 = new TH1F("hDef3","",30,-1.5,1.5);  // hit en L0 Y L1, mismo stave

  int nEv = thits->GetEntries();

  for(int iev=0; iev<nEv; iev++){
    thits->GetEntry(iev);
    tmc->GetEntry(iev);

    // Muon primario generado
    double eta_gen = -999.;
    int muonTID = -1;
    for(int it=0; it<(int)tracks->size(); it++){
      auto& tr = (*tracks)[it];
      if(tr.getMotherTrackId()>=0) continue;
      if(abs(tr.GetPdgCode())!=13) continue;
      double p = tr.GetP();
      if(p<1e-9) continue;
      eta_gen = TMath::ATanH(tr.Pz()/p);
      muonTID = it;
      break;
    }
    if(muonTID<0) continue;
    hGen->Fill(eta_gen);

    // Clasificar hits por capa
    bool hitL0=false, hitL1=false;
    int staveL0=-1, staveL1=-1;

    for(auto& h:*hits){
      if(h.GetTrackID()!=muonTID) continue;
      int lay=h.GetDetectorID();
      if(lay==0){
        hitL0=true;
        // Stave aproximado por phi
        double phi=TMath::ATan2(h.GetStartY(),h.GetStartX());
        if(phi<0) phi+=2.*TMath::Pi();
        staveL0=(int)(phi/(2.*TMath::Pi()/16.));
      }
      if(lay==1){
        hitL1=true;
        double phi=TMath::ATan2(h.GetStartY(),h.GetStartX());
        if(phi<0) phi+=2.*TMath::Pi();
        staveL1=(int)(phi/(2.*TMath::Pi()/16.));
      }
    }

    // Definicion 1: hit en L0 O L1
    if(hitL0 || hitL1) hDef1->Fill(eta_gen);

    // Definicion 2: hit en L0 Y L1
    if(hitL0 && hitL1) hDef2->Fill(eta_gen);

    // Definicion 3: hit en L0 Y L1, mismo stave
    if(hitL0 && hitL1 && staveL0==staveL1) hDef3->Fill(eta_gen);
  }

  // Calcular eficiencias
  printf("\n--- Aceptancia preliminar (1000 eventos, g=1mm) ---\n");
  printf("Def 1 (L0 OR L1):           %.1f%%\n", 100.*hDef1->Integral()/hGen->Integral());
  printf("Def 2 (L0 AND L1):          %.1f%%\n", 100.*hDef2->Integral()/hGen->Integral());
  printf("Def 3 (L0 AND L1, same stave): %.1f%%\n", 100.*hDef3->Integral()/hGen->Integral());

  // Figura
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c","MID Acceptance — 3 definitions",900,600);
  gPad->SetLeftMargin(0.12);

  TH1F* hEff1 = (TH1F*)hDef1->Clone("hEff1");
  TH1F* hEff2 = (TH1F*)hDef2->Clone("hEff2");
  TH1F* hEff3 = (TH1F*)hDef3->Clone("hEff3");
  hEff1->Divide(hGen); hEff2->Divide(hGen); hEff3->Divide(hGen);

  hEff1->SetLineColor(kGreen+2);  hEff1->SetLineWidth(2);
  hEff2->SetLineColor(kAzure+2);  hEff2->SetLineWidth(2);
  hEff3->SetLineColor(kRed+1);    hEff3->SetLineWidth(2);

  hEff1->SetTitle("MID Geometric Acceptance (preliminary);#eta;Acceptance");
  hEff1->GetYaxis()->SetRangeUser(0,1.2);
  hEff1->Draw("hist");
  hEff2->Draw("hist same");
  hEff3->Draw("hist same");

  TLegend* leg = new TLegend(0.13,0.20,0.55,0.38);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.038);
  leg->AddEntry(hEff1,Form("L0 OR L1  (%.1f%%)",100.*hDef1->Integral()/hGen->Integral()),"l");
  leg->AddEntry(hEff2,Form("L0 AND L1 (%.1f%%)",100.*hDef2->Integral()/hGen->Integral()),"l");
  leg->AddEntry(hEff3,Form("L0 AND L1, same stave (%.1f%%)",100.*hDef3->Integral()/hGen->Integral()),"l");
  leg->Draw();

  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.036);
  tx.DrawLatex(0.13,0.88,"1000 #mu^{-}, 1.5-5 GeV, |#eta|<1.2 (preliminary)");

  c->SaveAs("/tmp/acceptance_3def.png");
  printf("Guardado: /tmp/acceptance_3def.png\n");
}
