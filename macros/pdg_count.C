void pdg_count() {
  TFile* fhits = TFile::Open("o2sim_HitsMI3.root");
  TFile* fmc   = TFile::Open("o2sim_Kine.root");
  TTree* thits = (TTree*)fhits->Get("o2sim");
  TTree* tmc   = (TTree*)fmc->Get("o2sim");

  std::vector<o2::MCTrack>* tracks = nullptr;
  std::vector<o2::itsmft::Hit>* hits = nullptr;
  tmc->SetBranchAddress("MCTrack", &tracks);
  thits->SetBranchAddress("MI3Hit", &hits);

  std::map<int,int> pdgCount;
  int total = 0;

  for(int i=0; i<thits->GetEntries(); i++){
    thits->GetEntry(i);
    tmc->GetEntry(i);
    for(int j=0; j<(int)hits->size(); j++){
      int tid = (*hits)[j].GetTrackID();
      if(tid < (int)tracks->size()){
        int pdg = abs((*tracks)[tid].GetPdgCode());
        pdgCount[pdg]++;
        total++;
      }
    }
  }

  printf("Total hits: %d\n", total);
  printf("PDG        Hits\n");
  std::map<int,int>::iterator it;
  for(it=pdgCount.begin(); it!=pdgCount.end(); ++it)
    printf("%-10d %d\n", it->first, it->second);
}
