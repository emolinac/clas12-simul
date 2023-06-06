// Script that takes a .dat file and transforms it into a TTree in a ROOT file

// Req: .dat file has to have the following format
//      <event_index> <particle_id> <parent_id> <px> <py> <pz> <E> <x> <y> <z>

// author : Esteban Molina (May 2022)

void MACRO_dat2root(){
  // Open dat file
  const char* file_name = "lepto_out.dat";
  std::ifstream file(file_name);

  // Create Tree
  TTree* t = new TTree("ntuple_thrown_raw","");

  // Make the tree read the .dat file
  t->ReadFile(file_name,"event_index/D:PID:parent_PID:px:py:pz:E:x:y:z");

  // Create target root file
  TFile* f = new TFile("dat2root.root","RECREATE");
  f->cd();

  t->Write();
}
