// generate_hitmap.C
// Genera HitMap_O2.root con posiciones exactas de barras del MID
// Geometria validada por A. Ortiz (ICN-UNAM/ALICE-CERN), agosto 2026
// Formula de posicionamiento basada en Ian Perez Garcia (DetectorConstruction.cc)
// Nota: barStep = 5.2 cm = 2*barWidth + barSpacing = width + gap (Ian)

void generate_hitmap() {

  const int nLayers   = 2;
  const int nStaves   = 16;
  const int nModules  = 10;
  const float barStep    = 5.2f;   // paso entre centros de barras
  const float barWidth   = 2.5f;   // ancho de barra
  const float barSpacing = 0.2f;   // gap entre barras

  const float R[2]        = {301.f, 311.f};
  const int   nBars[2]    = {23, 20};
  const float barLen[2]   = {49.9f, 61.75f};
  const float staveLen[2] = {500.f, 525.f};
  const float sumWidth[2] = {
    (barStep * nBars[0]) / 2.f,  // L0: 59.8 cm
    (barStep * nBars[1]) / 2.f   // L1: 52.0 cm
  };

  TGraph* gx = new TGraph();
  TGraph* gy = new TGraph();
  TGraph* gz = new TGraph();
  gx->SetName("gx"); gy->SetName("gy"); gz->SetName("gz");

  int globalIdx = 0;

  for(int ilay=0; ilay<nLayers; ilay++){
    float r    = R[ilay];
    int   nb   = nBars[ilay];
    float bl   = barLen[ilay];
    float sw   = sumWidth[ilay];
    float zOff = -staveLen[ilay];

    for(int istave=0; istave<nStaves; istave++){
      float phi = TMath::TwoPi() * istave / nStaves;

      for(int imod=0; imod<nModules; imod++){
        // Posicion Z del centro del modulo — MIDLayer.cxx lineas 194-196
        float zMod;
        if(ilay==0){
          zMod = zOff + imod * 2 * bl + bl;
        } else {
          zMod = zOff + imod * 2 * sw + sw;
        }

        int counter  = 1;
        int halfBars = nb / 2;

        for(int ibar=0; ibar<nb; ibar++){
          float x_pos, y_pos, z_pos;

          if(ilay==0){
            // L0: barras en phi — formula polar de Ian
            // paso = barStep = 5.2 cm (= width+gap en Ian)
            float omega;
            if(ibar <= halfBars){
              omega = TMath::ATan((float)ibar * barStep / r);
              x_pos = TMath::Sqrt((float)(ibar*barStep*ibar*barStep) + r*r) * TMath::Cos(phi - omega);
              y_pos = TMath::Sqrt((float)(ibar*barStep*ibar*barStep) + r*r) * TMath::Sin(phi - omega);
            } else {
              omega = TMath::ATan((float)counter * barStep / r);
              x_pos = TMath::Sqrt((float)(counter*barStep*counter*barStep) + r*r) * TMath::Cos(phi + omega);
              y_pos = TMath::Sqrt((float)(counter*barStep*counter*barStep) + r*r) * TMath::Sin(phi + omega);
              counter++;
            }
            z_pos = zMod;

          } else {
            // L1: barras en Z — posicion en Z dentro del modulo
            x_pos = r * TMath::Cos(phi);
            y_pos = r * TMath::Sin(phi);
            float barCenter = -nb * barStep / 2.f + (2*ibar+1) * barStep/2.f;
            z_pos = zMod + barCenter;
          }

          gx->SetPoint(globalIdx, globalIdx, x_pos);
          gy->SetPoint(globalIdx, globalIdx, y_pos);
          gz->SetPoint(globalIdx, globalIdx, z_pos);
          globalIdx++;
        }
      }
    }
  }

  int nL0 = nStaves * nModules * nBars[0];
  int nL1 = nStaves * nModules * nBars[1];
  printf("Total barras: %d\n", globalIdx);
  printf("  L0: %d (%d staves x %d modulos x %d barras)\n", nL0, nStaves, nModules, nBars[0]);
  printf("  L1: %d (%d staves x %d modulos x %d barras)\n", nL1, nStaves, nModules, nBars[1]);
  printf("  Separacion L0/L1 en indice: %d\n", nL0);

  // Verificacion rapida — primeras y ultimas barras
  printf("\nPrimeras 3 barras L0:\n");
  for(int i=0; i<3; i++){
    printf("  idx=%d  X=%.2f  Y=%.2f  Z=%.2f\n",
      i, gx->GetY()[i], gy->GetY()[i], gz->GetY()[i]);
  }
  printf("Primeras 3 barras L1:\n");
  for(int i=nL0; i<nL0+3; i++){
    printf("  idx=%d  X=%.2f  Y=%.2f  Z=%.2f\n",
      i, gx->GetY()[i], gy->GetY()[i], gz->GetY()[i]);
  }

  TFile* fout = TFile::Open("/home/franciscoe/MID-Geometry-Studies/HitMap_O2.root","RECREATE");
  gx->Write(); gy->Write(); gz->Write();
  fout->Close();
  printf("\nGuardado: HitMap_O2.root\n");
}
