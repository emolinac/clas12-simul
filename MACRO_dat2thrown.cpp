// Script that takes a .dat file and transforms it into a TTree in a ROOT file

// Req: .dat file has to have the following format
//      <event_index> <particle_id> <parent_id> <px> <py> <pz> <E> <x> <y> <z>

// author : Esteban Molina (May 2022)

void setBranchesAddresses(TTree* t, Float_t* event_index, Float_t* PID, Float_t* parent_PID, Float_t* Px, Float_t* Py, Float_t* Pz, Float_t* E, Float_t* x, Float_t* y, Float_t* z);

// Constants
const Float_t kEbeam = 11.;
const Float_t kMassGamma = 0.000000;
const Float_t kMassPositron = 5.109998e-4;
const Float_t kMassElectron = 5.109998e-4;
const Float_t kMassNeutrino = 0.000000;
const Float_t kMassMuonPlus = 0.105658;
const Float_t kMassMuonMinus = 0.105658;
const Float_t kMassPiZero = 0.134977;
const Float_t kMassPiPlus = 0.139570;
const Float_t kMassPiMinus = 0.139570;
const Float_t kMassKaonZeroLong = 0.497614;
const Float_t kMassKaonPlus = 0.493677;
const Float_t kMassKaonMinus = 0.493677;
const Float_t kMassNeutron = 0.939565;
const Float_t kMassProton = 0.938272;
const Float_t kMassKaonZeroShort = 0.497614;
const Float_t kMassEta = 0.547853;
const Float_t kMassOmega = 0.782650;
const Float_t kMassKaonZero = 0.497614;

// LEPTONIC CLASS

class LeptonicKinematics{
  Float_t Q2, Xb, Nu, Px_el, Py_el, Pz_el, P_el, ThetaLab_el, PhiLab_el;

public:
  LeptonicKinematics(Float_t, Float_t, Float_t);
  ~LeptonicKinematics();

private:
  Float_t getQ2()		{return Q2;}
  Float_t getXb()		{return Xb;}
  Float_t getNu()		{return Nu;}
  Float_t getP_el()		{return P_el;}
  Float_t getThetaLab_el()	{return ThetaLab_el;}
  Float_t getPhiLab_el()	{return PhiLab_el;}
  Float_t getPx_el()		{return Px_el;}
  Float_t getPy_el()		{return Py_el;}
  Float_t getPz_el()		{return Pz_el;}  
};

LeptonicKinematics::LeptonicKinematics(Float_t Px, Float_t Py, Float_t Pz){
  // Class constructor
  TVector3 v(Px, Py, Pz);

  // momentum
  P_el		= TMath::Sqrt(v.Mag());
  Px_el		= Px;
  Py_el		= Py;
  Pz_el		= Pz;
  // direction
  ThetaLab_el	= v.Theta()*TMath::RadToDeg();
  PhiLab_el	= v.Phi()*TMath::RadToDeg();
  // leptonic
  Q2		= 4.*kEbeam*TMath::Power(TMath::Sin(ThetaLab_el*TMath::DegToRad()/2.),2);
  Nu		= kEbeam - P_el;
  Xb		= Q2/2./kMprt/Nu;
}

LeptonicKinematics::~LeptonicKinematics(){}

// HADRONIC CLASS

class HadronicKinematics{
  Float_t Px_h, Py_h, Pz_h, P_h, ThetaLab_h, PhiLab_h, Zh, Pl2, Pt2, PhiPQ, ThetaPQ, PID_h;

public:
  HadronicKinematics(Float_t, Float_t, Float_t, Float_t);
  ~HadronicKinematics();
  friend class LeptonicKinematics;

private:
  Float_t getMass_h(Float_t PID_h);
  
  Float_t getThetaLab_h()	{return ThetaLab_h;}
  Float_t getPhiLab_h()	        {return PhiLab_h;}
  Float_t getP_h()		{return P_h;}
  Float_t getPx_h()		{return Px_h;}
  Float_t getPy_h()		{return Py_h;}
  Float_t getPz_h()		{return Pz_h;}

  Float_t getThetaPQ(LeptonicKinematics* lk, Float_t Px_h, Float_t Py_h, Float_t Pz_h);
  Float_t getPhiPQ(LeptonicKinematics* lk, Float_t Px_h, Float_t Py_h, Float_t Pz_h);
  Float_t getCosThetaPQ(LeptonicKinematics* lk, Float_t Px_h, Float_t Py_h, Float_t Pz_h);
  Float_t getZh();
  Float_t getPl2();
  Float_t getPt2();
};

HadronicKinematics::HadronicKinematics(Float_t Px, Float_t Py, Float_t Pz, Float_t PID){
  // Class constructor
  TVector3 v(Px, Py, Pz);
  PID_h = PID;
  
  // momentum
  P_h		= TMath::Sqrt(v.Mag());
  Px_h		= Px;
  Py_h		= Py;
  Pz_h		= Pz;
  // direction
  ThetaLab_h	= v.Theta()*TMath::RadToDeg();
  PhiLab_h	= v.Phi()*TMath::RadToDeg();
}

HadronicKinematics::~HadronicKinematics(){}

Float_t HadronicKinetics::getMass_h(Float_t PID_h){
const Float_t kMassGamma = 0.000000;
const Float_t kMassPositron = 5.109998e-4;
const Float_t kMassElectron = 5.109998e-4;
const Float_t kMassNeutrino = 0.000000;
const Float_t kMassMuonPlus = 0.105658;
const Float_t kMassMuonMinus = 0.105658;
const Float_t kMassPiZero = 0.134977;
const Float_t kMassPiPlus = 0.139570;
const Float_t kMassPiMinus = 0.139570;
const Float_t kMassKaonZeroLong = 0.497614;
const Float_t kMassKaonPlus = 0.493677;
const Float_t kMassKaonMinus = 0.493677;
const Float_t kMassNeutron = 0.939565;
const Float_t kMassProton = 0.938272;
const Float_t kMassKaonZeroShort = 0.497614;
const Float_t kMassEta = 0.547853;
const Float_t kMassOmega = 0.782650;
const Float_t kMassKaonZero = 0.497614;

if 
}

Float_t HadronicKinematics::getPhiPQ(LeptonicKinematics* lk, Float_t Px, Float_t Py, Float_t Pz) {
  // Returns the azimuthal angle of the particle w.r.t. the virtual photon direction
  // First, it Z-rotates the virtual photon momentum to have Y-component=0
  // Second, it Z-rotates the particle momentum by the same amount
  // Third, it Y-rotates the virtual photon to have X-component=0
  // Lastly, it Y-rotates the particle momentum by the same amount
  // In the end, the values of the particle momentum components will be w.r.t to the virtual photon momentum
  TVector3 Vpi(Px, Py, Pz);
  TVector3 Vvirt(-lk.Px(), -lk.Py(), kEbeam - lk.Pz());
  Double_t phi_z = TMath::Pi() - Vvirt.Phi();
  Vvirt.RotateZ(phi_z);
  Vpi.RotateZ(phi_z);
  TVector3 Vhelp(0., 0., 1.);
  Double_t phi_y = Vvirt.Angle(Vhelp);
  Vvirt.RotateY(phi_y);
  Vpi.RotateY(phi_y);
  return Vpi.Phi()*TMath::RadToDeg();
}

Float_t HadronicKinematics::getThetaPQ(LeptonicKinematics* lk) {
  // Return the polar angle of the particle w.r.t. the virtual photon direction
  // It's defined as the angle between both particle's momentum
  TVector3 Vpi(this.Px_h, this.Py_h, this.Pz_h);
  TVector3 Vvirt(-lk.Px(), -lk.Py(), kEbeam - lk.Pz());
  return Vvirt.Angle(Vpi)*TMath::RadToDeg();
}

Float_t HadronicKinematics::getCosThetaPQ(LeptonicKinematics* lk) {
  // Returns the cosine of ThetaPQ for the particle
  return (this.Pz_h*(kEbeam - lk.Pz()) - this.Px*lk.Px() - this.Py*lk.Py())/(TMath::Sqrt(lk.Nu()*lk.Nu() + lk.Q2())*this.P_h());
}

Float_t HadronicKinematics::Zh(LeptonicKinematics* lk) {
  // Returns the energy fraction of the particle
  return TMath::Sqrt(mass*mass + this.P_h()*this.P_h())/lk.Nu();
}
Float_t HadronicKinematics::Pt2(LeptonicKinematics* lk) {
  // Returns the square of the transverse momentum component w.r.t. the virtual photon direction
  return Momentum(k, kind)*Momentum(k, kind)*(1 - CosThetaPQ(k, kind)*CosThetaPQ(k, kind));
}
Float_t HadronicKinematics::Pl2(LeptonicKinematics* lk) {
  // Returns the square of the longitudinal momentum component w.r.t. the virtual photon direction
  return Momentum(k, kind)*Momentum(k, kind)*CosThetaPQ(k, kind)*CosThetaPQ(k, kind);
}
 
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
    LeptonicKinematics* lk = new LeptonicKinematics(Px,Py,Pz);
    if(PID==11){
      // Calculate leptonic variables
      std::cout<<"Q2   = "<<lk->getQ2()<<std::endl;
      std::cout<<"Xb   = "<<lk->getXb()<<std::endl;
      std::cout<<"Nu   = "<<lk->getNu()<<std::endl;
      std::cout<<"P_el = "<<lk->getP_el()<<std::endl;
      std::cout<<"Theta= "<<lk->getThetaLab_el()<<std::endl;
      std::cout<<"Phi  = "<<lk->getPhiLab_el()<<std::endl;
    }
    else{
      // Calculate hadronic variables
    }
    delete lk;
  }



  // // Create target root file
  // TFile* f = new TFile("dat2root.root","RECREATE");
  // f->cd();

  // t->Write();
}


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
