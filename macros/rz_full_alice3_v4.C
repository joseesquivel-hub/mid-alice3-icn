void rz_full_alice3_v4() {

  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetGridColor(kGray);

  TCanvas* c = new TCanvas("c","ALICE 3 RZ",1400,700);
  c->SetFillColor(kWhite);
  gPad->SetFillColor(kWhite);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.10);
  gPad->SetTopMargin(0.09);
  gPad->SetBottomMargin(0.20);
  gPad->SetGridx(); gPad->SetGridy();

  TH2F* frame = new TH2F("frame",
    "ALICE 3 MID - R-Z Layout (current implementation);z (m);R (m)",
    120,-6,6,70,0,3.5);
  frame->SetStats(0);
  frame->GetXaxis()->SetNdivisions(12);
  frame->GetYaxis()->SetNdivisions(7);
  frame->GetXaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetXaxis()->SetTitleSize(0.040);
  frame->GetYaxis()->SetTitleSize(0.040);
  frame->Draw();

  auto box = [](double zmin, double zmax, double rmin, double rmax,
                int fillcol, int linecol, double alpha=0.5) {
    TBox* b = new TBox(zmin,rmin,zmax,rmax);
    b->SetFillColorAlpha(fillcol,alpha);
    b->SetLineColor(linecol);
    b->SetLineWidth(1);
    b->Draw();
    return b;
  };

  // ============================================================
  // LINEAS DE ETA
  // ============================================================
  double etaVals[5]={0.5,1.0,1.5,2.0,2.5};
  int etaColors[5]={kGray+1,kGray+1,kGray+2,kGray+1,kGray+1};
  int etaStyles[5]={3,3,2,3,3}; // eta=1.5 punteado distinto
  TLatex etaLab; etaLab.SetTextColor(kGray+1); etaLab.SetTextSize(0.022);

  for(int ie=0;ie<5;ie++){
    double eta=etaVals[ie];
    double theta=2*TMath::ATan(TMath::Exp(-eta));
    double zend=6.0, rend=zend*TMath::Tan(theta);
    if(rend>3.5){ rend=3.5; zend=rend/TMath::Tan(theta); }
    TLine* lp=new TLine(0,0,zend,rend);
    lp->SetLineColor(etaColors[ie]); lp->SetLineStyle(etaStyles[ie]);
    lp->SetLineWidth(ie==2?2:1); lp->Draw();
    TLine* lm=new TLine(0,0,-zend,rend);
    lm->SetLineColor(etaColors[ie]); lm->SetLineStyle(etaStyles[ie]);
    lm->SetLineWidth(ie==2?2:1); lm->Draw();
    etaLab.DrawLatex(zend+0.05,rend,Form("#eta=%.1f",eta));
    etaLab.DrawLatex(-zend-0.40,rend,Form("#eta=-%.1f",eta));
  }

  // Lineas especiales eta=1.28 y 1.30 (borde MID)
  for(double eta : {1.28, 1.30}){
    double theta=2*TMath::ATan(TMath::Exp(-eta));
    double zend=4.99, rend=zend*TMath::Tan(theta);
    TLine* lp=new TLine(0,0,zend,rend);
    lp->SetLineColor(kOrange+1); lp->SetLineStyle(2); lp->SetLineWidth(1); lp->Draw();
    TLine* lm=new TLine(0,0,-zend,rend);
    lm->SetLineColor(kOrange+1); lm->SetLineStyle(2); lm->SetLineWidth(1); lm->Draw();
  }
  TLatex etaMID; etaMID.SetTextColor(kOrange+1); etaMID.SetTextSize(0.020);
  etaMID.DrawLatex(4.8,1.18,"#eta=1.28/1.30");

  // ============================================================
  // TRK
  // ============================================================
  int cTRK=kCyan+2;
  box(-0.64,0.64,0.068,0.081,cTRK,cTRK,0.7);
  box(-0.64,0.64,0.088,0.100,cTRK,cTRK,0.7);
  box(-0.64,0.64,0.117,0.129,cTRK,cTRK,0.7);
  box(-0.64,0.64,0.187,0.202,cTRK,cTRK,0.7);
  box(-0.64,0.64,0.299,0.311,cTRK,cTRK,0.7);
  box(-1.28,1.28,0.449,0.461,cTRK,cTRK,0.7);
  box(-1.28,1.28,0.599,0.611,cTRK,cTRK,0.7);
  box(-1.28,1.28,0.799,0.811,cTRK,cTRK,0.7);
  TLatex lTRK; lTRK.SetTextColor(cTRK); lTRK.SetTextSize(0.026); lTRK.SetTextFont(62);
  lTRK.DrawLatex(0.70,0.83,"TRK");

  // ============================================================
  // TOF
  // ============================================================
  int cTOF=kOrange+7;
  box(-0.645,0.645,0.210,0.222,cTOF,cTOF,0.9);
  box(-3.40,3.40,0.920,0.932,cTOF,cTOF,0.9);
  TLatex lTOF; lTOF.SetTextColor(cTOF); lTOF.SetTextSize(0.024); lTOF.SetTextFont(62);
  lTOF.DrawLatex(0.70,0.222,"iTOF");
  lTOF.DrawLatex(3.45,0.935,"oTOF");

  // ============================================================
  // RICH
  // ============================================================
  int cRICH=kGreen+2;
  box(-2.00,2.00,1.050,1.309,cRICH,cRICH,0.4);
  TLatex lRICH; lRICH.SetTextColor(cRICH); lRICH.SetTextSize(0.028); lRICH.SetTextFont(62);
  lRICH.DrawLatex(-0.5,1.17,"RICH");

  // ============================================================
  // IMAN
  // ============================================================
  box(-4.00,4.00,1.400,1.410,kGray+1,kGray+2,0.9); // InnerWrap Al
  box(-4.00,4.00,1.970,2.000,kGray+1,kGray+2,0.9); // OuterWrap Al
  box(-4.00,4.00,1.410,1.970,kPink-9,kPink-9,0.3);  // Criostato
  box(-4.00,4.00,1.600,1.700,kOrange-3,kOrange,0.8);// Soporte Al
  TLatex lMAG; lMAG.SetTextColor(kOrange+1); lMAG.SetTextSize(0.026); lMAG.SetTextFont(62);
  lMAG.DrawLatex(-1.0,1.85,"MAGNET (2T)");

  // ============================================================
  // ABSORBEDOR
  // ============================================================
  int cABS=kAzure-4, cABSL=kAzure+2;
  box(-3.00,3.00,2.20,2.90,cABS,cABSL,0.5);
  box(3.00,5.00,2.45,2.90,cABS,cABSL,0.5);
  box(-5.00,-3.00,2.45,2.90,cABS,cABSL,0.5);
  TLatex lABS; lABS.SetTextColor(kAzure+3); lABS.SetTextSize(0.028); lABS.SetTextFont(62);
  lABS.DrawLatex(-0.8,2.52,"ABSORBER (Fe)");

  // ============================================================
  // MID — una sola banda declarada como L0+L1
  // ============================================================
  TBox* bmid0=new TBox(-4.99,3.005,4.99,3.025);
  bmid0->SetFillColorAlpha(kMagenta-4,0.9);
  bmid0->SetLineColor(kMagenta); bmid0->SetLineWidth(2); bmid0->Draw();
  TBox* bmid1=new TBox(-5.25,3.030,5.25,3.050);
  bmid1->SetFillColorAlpha(kViolet+2,0.9);
  bmid1->SetLineColor(kViolet+2); bmid1->SetLineWidth(2); bmid1->Draw();
  TLatex lMID; lMID.SetTextSize(0.024);

  // ============================================================
  // LEYENDA
  // ============================================================
  TLegend* leg=new TLegend(0.08,0.01,0.90,0.10);
  leg->SetNColumns(6);
  leg->SetBorderSize(1); leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite); leg->SetTextSize(0.024);
  TBox* bAl=new TBox(0,0,1,1); bAl->SetFillColor(kGray+1);
  TBox* bSup=new TBox(0,0,1,1); bSup->SetFillColor(kOrange-3);
  TBox* bFe=new TBox(0,0,1,1); bFe->SetFillColorAlpha(kAzure-4,0.5);
  TBox* bL0=new TBox(0,0,1,1); bL0->SetFillColorAlpha(kMagenta-4,0.9);
  TBox* bL1=new TBox(0,0,1,1); bL1->SetFillColorAlpha(kViolet+2,0.9);
  TBox* bTRK=new TBox(0,0,1,1); bTRK->SetFillColorAlpha(kCyan+2,0.7);
  TBox* bTOF=new TBox(0,0,1,1); bTOF->SetFillColorAlpha(kOrange+7,0.9);
  TBox* bRCH=new TBox(0,0,1,1); bRCH->SetFillColorAlpha(kGreen+2,0.4);
  leg->AddEntry(bTRK,"TRK","f");
  leg->AddEntry(bTOF,"TOF (inner+outer)","f");
  leg->AddEntry(bRCH,"RICH","f");
  leg->AddEntry(bFe,"Fe absorber (Rext=290 const.)","f");
  leg->AddEntry(bL0,"MID L0 (|z|<499 cm)","f");
  leg->AddEntry(bL1,"MID L1 (|z|<525 cm)","f");
  leg->Draw();

  // Nota grosor MID
  TLatex note; note.SetNDC(); note.SetTextColor(kGray+1); note.SetTextSize(0.020);
  note.DrawLatex(0.08,0.115,"MID layer thickness exaggerated for visibility");

  c->SaveAs("/tmp/rz_full_alice3_v4.png");
  printf("Guardado: rz_full_alice3_v4.png\n");
}
