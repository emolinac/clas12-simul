// Script that takes a .dat file and transforms it into a TTree in a ROOT file

// Req: .dat file has to have the following format
//      <event_index> <particle_id> <parent_id> <px> <py> <pz> <E> <x> <y> <z>

// author : Esteban Molina (May 2022)

// Constants
Float_t kEbeam = 11.;
Float_t kMprt  = 0.938272;
Float_t kMpi   = 0.13960;
Float_t kMntr  = 0.939565;

void setBranchesAddresses(TTree* t, Float_t* event_index, Float_t* PID, Float_t* parent_PID, Float_t* Px, Float_t* Py, Float_t* Pz, Float_t* E, Float_t* x, Float_t* y, Float_t* z){
  t->SetBranchAddress("event_index", event_index);
  t->SetBranchAddress("PID",         PID);
  t->SetBranchAddress("parent_PID",  parent_PID);
  t->SetBranchAddress("Px",          Px);
  t->SetBranchAddress("Py",          Py);
  t->SetBranchAddress("Pz",          Pz);
  t->SetBranchAddress("E",           E);
  t->SetBranchAddress("x",           x);
  t->SetBranchAddress("y",           y);
  t->SetBranchAddress("z",           z);
}

class DISkinematics{
}

// Fundamental variables
Float_t ThetaLab(Float_t Px, Float_t Py, Float_t Pz){
  TVector3 v(Px, Py, Pz);
  return v.Theta()*TMath::RadToDeg();
}
Float_t PhiLab(Float_t Px, Float_t Py, Float_t Pz){
  TVector3 v(Px, Py, Pz);
  return v.Phi()*TMath::RadToDeg();
}
Float_t P(Float_t Px, Float_t Py, Float_t Pz){
  return Px*Px+Py*Py+Pz*Pz;
}

// Leptonic variables
Float_t Q2(Float_t Momentum, Float_t ThetaLab){
  return 4.*kEbeam*TMath::Power(TMath::Sin(ThetaLab*TMath::DegToRad()/2.),2);
}
Float_t Xb(Float_t Q2, Float_t Nu){
  return Q2/2./kMprt/Nu;
}
Float_t Nu(Float_t Eel){
  return kEbeam - Eel;
}

// Hadronic variables
Float_t Pt(){}
Float_t Pl(){}
Float_t Zh(){}
Float_t PhiPQ(){}
Float_t ThetaPQ(){}

void MACRO_dat2thrown(){
  // Open dat file
  const char* file_name = "lepto_out.dat";
  std::ifstream file(file_name);

  // Create Tree
  TTree* t				= new TTree("ntuple_thrown_raw","");
  TNtuple* ntuple_thrown_electrons	= new TNtuple("ntuple_thrown_electrons",
						      "",
						      "Q2:Xb:Nu:ThetaLab:PhiLab:P:Px:Py:Pz:E");
  TNtuple* ntuple_thrown		= new TNtuple("ntuple_thrown",
						      "",
						      "Q2:Xb:Nu:Zh:Pt:Pt2:Pl:Pl2:PhiPQ:Xf:ThetaLab:PhiLab:P:Px:Py:Pz:E:ThetaLab_el:PhiLab_el:P_el:Px_el:Py_el:Pz_el:E_el");
  // Make the tree read the .dat file
  t->ReadFile(file_name,"event_index/F:PID:parent_PID:Px:Py:Pz:E:x:y:z");

  //Process the tree
  Float_t event_index, PID, parent_PID, Px, Py, Pz, E, x, y, z;
  setBranchesAddresses(t, &event_index, &PID, &parent_PID, &Px, &Py, &Pz, &E, &x, &y, &z);

  Int_t Nentries = t->GetEntries();
  for(Int_t entry = 0 ; entry < Nentries ; entry++){
    t->GetEntry(entry);
    std::cout<<"event_index="<<event_index<<"   Pz="<<Pz<<std::endl;
    if(PID==11){
      // Calculate leptonic variables
    }
    else{
      // Calculate hadronic variables
    }
  }



  // // Create target root file
  // TFile* f = new TFile("dat2root.root","RECREATE");
  // f->cd();

  // t->Write();
}
