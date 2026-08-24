void color_mid_geometry() {
  TGeoManager::Import("o2sim_geometry.root");
  if(!gGeoManager){ printf("ERROR: no se pudo cargar o2sim_geometry.root\n"); return; }

  TGeoVolume* v;
  v = gGeoManager->GetVolume("cave");        if(v) v->SetVisibility(0);
  v = gGeoManager->GetVolume("innerVacuum"); if(v) v->SetVisibility(0);
  v = gGeoManager->GetVolume("outerVacuum"); if(v) v->SetVisibility(0);
  v = gGeoManager->GetVolume("innerWrap");
  if(v){ v->SetLineColor(kGray+1); v->SetFillColor(kGray+1); v->SetTransparency(20); }
  v = gGeoManager->GetVolume("coils");
  if(v){ v->SetLineColor(kRed+1); v->SetFillColor(kRed+1); v->SetTransparency(10); }
  v = gGeoManager->GetVolume("restMaterial");
  if(v){ v->SetLineColor(kOrange+1); v->SetFillColor(kOrange+1); v->SetTransparency(20); }
  v = gGeoManager->GetVolume("outerWrap");
  if(v){ v->SetLineColor(kGray+2); v->SetFillColor(kGray+2); v->SetTransparency(20); }
  v = gGeoManager->GetVolume("Absorber");
  if(v){ v->SetLineColor(kAzure-4); v->SetFillColor(kAzure-4); v->SetTransparency(40); }

  int nL0=0, nL1=0;
  TObjArray* vols = gGeoManager->GetListOfVolumes();
  for(int i=0; i<vols->GetEntries(); i++){
    TGeoVolume* vol = (TGeoVolume*)vols->At(i);
    TString name = vol->GetName();
    if(name.Contains("MIDSensor")){
      if(name.Contains("_L0_")){ vol->SetLineColor(kViolet-4); vol->SetFillColor(kViolet-4); nL0++; }
      else { vol->SetLineColor(kViolet+2); vol->SetFillColor(kViolet+2); nL1++; }
    }
  }
  printf("MIDSensor L0: %d   L1: %d\n", nL0, nL1);
  gGeoManager->Export("geom_mid_colored.root");
  printf("Exportado: geom_mid_colored.root\n");
  printf("Abrir en: https://root.cern/js/latest/\n");
}
