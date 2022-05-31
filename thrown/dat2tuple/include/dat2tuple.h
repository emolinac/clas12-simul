#include "TTree.h"
#include "TNtuple.h"
#include "TVector3.h"
#include "TMath.h"

#include "constants.h"
//####################################################################################################################//
//########################################       COMMON FUNCTIONS        #############################################//
//####################################################################################################################//

void setBranchesAddresses(TTree* t, float* event_index, float* PID, float* parent_PID, float* Px, float* Py, float* Pz, float* E, float* x, float* y, float* z){
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

//####################################################################################################################//
//########################################     CLASSES DECLARATIONS      #############################################//
//####################################################################################################################//

// LEPTONIC CLASS

class LeptonicKinematics{
  float Q2, Xb, Nu, Px_el, Py_el, Pz_el, P_el, ThetaLab_el, PhiLab_el;

public:
  LeptonicKinematics(float, float, float);
  ~LeptonicKinematics();

  //protected:
  float getQ2()	        	{return Q2;}
  float getXb()	         	{return Xb;}
  float getNu()		        {return Nu;}
  float getP_el()		{return P_el;}
  float getThetaLab_el()	{return ThetaLab_el;}
  float getPhiLab_el()      	{return PhiLab_el;}
  float getPx_el()		{return Px_el;}
  float getPy_el()		{return Py_el;}
  float getPz_el()		{return Pz_el;}  
};

// HADRONIC CLASS

class HadronicKinematics{
  float Px_h, Py_h, Pz_h, P_h, ThetaLab_h, PhiLab_h, PID_h;

public:
  friend class LeptonicKinematics;
  
  HadronicKinematics(float, float, float, float);
  ~HadronicKinematics();

  float getMass_h(float);
  
  float getThetaLab_h()	        {return ThetaLab_h;}
  float getPhiLab_h()	        {return PhiLab_h;}
  float getP_h()		{return P_h;}
  float getPx_h()		{return Px_h;}
  float getPy_h()		{return Py_h;}
  float getPz_h()		{return Pz_h;}

  float getThetaPQ(LeptonicKinematics* lk);
  float getPhiPQ(LeptonicKinematics* lk);
  float getCosThetaPQ(LeptonicKinematics* lk);
  float getZh(LeptonicKinematics* lk);
  float getPl2(LeptonicKinematics* lk);
  float getPt2(LeptonicKinematics* lk);
};

//####################################################################################################################//
//########################################       METHOD DEFINITIONS      #############################################//
//####################################################################################################################//

// Leptonic class

LeptonicKinematics::LeptonicKinematics(float Px, float Py, float Pz){
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
  Xb		= Q2/2./kMassProton/Nu;
}

LeptonicKinematics::~LeptonicKinematics(){}

// Hadronic class

HadronicKinematics::HadronicKinematics(float Px, float Py, float Pz, float PID){
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

float HadronicKinematics::getMass_h(float PID_h){
  if(PID_h == 211){
    return kMassPiPlus;
  } else if(PID_h == -211){
    return kMassPiMinus;
  } else if(PID_h == 111){
    return kMassPiZero;
  } else if(PID_h == 221){
    return kMassEta;
  } else if(PID_h == 223){
    return kMassOmega;
  } else if(PID_h == 2212){
    return kMassProton;
  } else if(PID_h == 2112){
    return kMassNeutron;
  } else if(PID_h == 311){
    return kMassKaonZero;
  } else if(PID_h == 321){
    return kMassKaonPlus;    
  } else if(PID_h == -321){
    return kMassKaonMinus;    
  } else if(PID_h == 22){
    return kMassGamma;
  } else{
    return 0;
  }
}

float HadronicKinematics::getPhiPQ(LeptonicKinematics* lk) {
  // Returns the azimuthal angle of the particle w.r.t. the virtual photon direction
  // First, it Z-rotates the virtual photon momentum to have Y-component=0
  // Second, it Z-rotates the particle momentum by the same amount
  // Third, it Y-rotates the virtual photon to have X-component=0
  // Lastly, it Y-rotates the particle momentum by the same amount
  // In the end, the values of the particle momentum components will be w.r.t to the virtual photon momentum
  TVector3 Vpi(this->Px_h, this->Py_h, this->Pz_h);
  TVector3 Vvirt(-lk->getPx_el(), -lk->getPy_el(), kEbeam - lk->getPz_el());
  Double_t phi_z = TMath::Pi() - Vvirt.Phi();
  Vvirt.RotateZ(phi_z);
  Vpi.RotateZ(phi_z);
  TVector3 Vhelp(0., 0., 1.);
  Double_t phi_y = Vvirt.Angle(Vhelp);
  Vvirt.RotateY(phi_y);
  Vpi.RotateY(phi_y);
  return Vpi.Phi()*TMath::RadToDeg();
}

float HadronicKinematics::getThetaPQ(LeptonicKinematics* lk) {
  // Return the polar angle of the particle w.r.t. the virtual photon direction
  // It's defined as the angle between both particle's momentum
  TVector3 Vpi(this->Px_h, this->Py_h, this->Pz_h);
  TVector3 Vvirt(-lk->getPx_el(), -lk->getPy_el(), kEbeam - lk->getPz_el());
  return Vvirt.Angle(Vpi)*TMath::RadToDeg();
}

float HadronicKinematics::getCosThetaPQ(LeptonicKinematics* lk) {
  // Returns the cosine of ThetaPQ for the particle
  return (this->Pz_h*(kEbeam - lk->getPz_el()) - this->Px_h*lk->getPx_el() - this->Py_h*lk->getPy_el())/(TMath::Sqrt(lk->getNu()*lk->getNu() + lk->getQ2())*this->getP_h());

  return (Pz(k, kind)*(kEbeam - Pz(0, kind)) - Px(k, kind)*Px(0, kind) - Py(k, kind)*Py(0, kind))/(TMath::Sqrt(Nu(kind)*Nu(kind) + Q2(kind))*Momentum(k, kind));
}

float HadronicKinematics::getZh(LeptonicKinematics* lk) {
  // Returns the energy fraction of the particle
  // TEMP: PI+ MASS
  return TMath::Sqrt(this->getMass_h(this->PID_h)*this->getMass_h(this->PID_h) + this->getP_h()*this->getP_h())/lk->getNu();
}
float HadronicKinematics::getPt2(LeptonicKinematics* lk) {
  // Returns the square of the transverse momentum component w.r.t. the virtual photon direction
  return this->getP_h()*this->getP_h()*(1 - TMath::Power(this->getCosThetaPQ(lk),2));
}
float HadronicKinematics::getPl2(LeptonicKinematics* lk) {
  // Returns the square of the longitudinal momentum component w.r.t. the virtual photon direction
  return this->getP_h()*this->getP_h()*TMath::Power(this->getCosThetaPQ(lk),2);
}
