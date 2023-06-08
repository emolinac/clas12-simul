#include "TTree.h"
#include "TNtuple.h"
#include "TVector3.h"
#include "TMath.h"
#include "constants.h"

#include <iostream>
//####################################################################################################################//
//########################################       COMMON FUNCTIONS        #############################################//
//####################################################################################################################//

void setBranchesAddresses(TTree* t, double* event_index, double* PID, double* parent_PID, double* Px, double* Py, double* Pz, double* E, double* x, double* y, double* z){
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
  double Q2, Xb, Nu, W, y, Px_el, Py_el, Pz_el, P_el, ThetaLab_el, PhiLab_el;

public:
  LeptonicKinematics(double, double, double);
  ~LeptonicKinematics();

  //protected:
  double getQ2()	        {return Q2;}
  double getXb()	        {return Xb;}
  double getNu()		{return Nu;}
  double getW()		        {return W;}
  double gety()		        {return y;}
  double getP_el()		{return P_el;}
  double getThetaLab_el()	{return ThetaLab_el;}
  double getPhiLab_el()      	{return PhiLab_el;}
  double getPx_el()		{return Px_el;}
  double getPy_el()		{return Py_el;}
  double getPz_el()		{return Pz_el;}  
};

// HADRONIC CLASS

class HadronicKinematics{
  double Px_h, Py_h, Pz_h, P_h, ThetaLab_h, PhiLab_h, PID_h;

public:
  friend class LeptonicKinematics;
  
  HadronicKinematics(double, double, double, double);
  ~HadronicKinematics();

  double getMass_h(double);
  
  double getThetaLab_h()	        {return ThetaLab_h;}
  double getPhiLab_h()	        {return PhiLab_h;}
  double getP_h()		{return P_h;}
  double getPx_h()		{return Px_h;}
  double getPy_h()		{return Py_h;}
  double getPz_h()		{return Pz_h;}

  double getThetaPQ(LeptonicKinematics* lk);
  double getPhiPQ(LeptonicKinematics* lk);
  double getCosThetaPQ(LeptonicKinematics* lk);
  double getZh(LeptonicKinematics* lk);
  double getPl2(LeptonicKinematics* lk);
  double getPt2(LeptonicKinematics* lk);
};

//####################################################################################################################//
//########################################       METHOD DEFINITIONS      #############################################//
//####################################################################################################################//

// Leptonic class

LeptonicKinematics::LeptonicKinematics(double Px, double Py, double Pz){
  // Class constructor
  TVector3 v(Px, Py, Pz);

  // momentum
  P_el		= TMath::Sqrt(Px*Px + Py*Py + Pz*Pz);
  Px_el		= Px;
  Py_el		= Py;
  Pz_el		= Pz;

  // direction
  ThetaLab_el	= v.Theta()*TMath::RadToDeg();
  PhiLab_el	= v.Phi()*TMath::RadToDeg();

  // leptonic
  Q2		= 4.*kEbeam*P_el*TMath::Sin(v.Theta()/2.)*TMath::Sin(v.Theta()/2.);
  Nu		= kEbeam - P_el;
  Xb		= Q2/2./kMassProton/Nu;
  W     = TMath::Sqrt(kMassProton*kMassProton + 2.*kMassProton*Nu - Q2);
  y     = Nu/kEbeam;
}

LeptonicKinematics::~LeptonicKinematics(){}

// Hadronic class

HadronicKinematics::HadronicKinematics(double Px, double Py, double Pz, double PID){
  // Class constructor
  TVector3 v(Px, Py, Pz);
  PID_h = PID;
  
  // momentum
  //  P_h		= v.Mag();
  P_h           = TMath::Sqrt(Px*Px + Py*Py + Pz*Pz);
  Px_h		= Px;
  Py_h		= Py;
  Pz_h		= Pz;
  // direction
  ThetaLab_h	= v.Theta()*TMath::RadToDeg();
  PhiLab_h	= v.Phi()*TMath::RadToDeg();
}

HadronicKinematics::~HadronicKinematics(){}

double HadronicKinematics::getMass_h(double PID_h){
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

double HadronicKinematics::getPhiPQ(LeptonicKinematics* lk) {
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

double HadronicKinematics::getThetaPQ(LeptonicKinematics* lk) {
  // Return the polar angle of the particle w.r.t. the virtual photon direction
  // It's defined as the angle between both particle's momentum
  TVector3 Vpi(this->Px_h, this->Py_h, this->Pz_h);
  TVector3 Vvirt(-lk->getPx_el(), -lk->getPy_el(), kEbeam - lk->getPz_el());
  return Vvirt.Angle(Vpi)*TMath::RadToDeg();
}

double HadronicKinematics::getCosThetaPQ(LeptonicKinematics* lk) {
  // Returns the cosine of ThetaPQ for the particle
  double Px_h = this->Px_h;
  double Py_h = this->Py_h;
  double Pz_h = this->Pz_h;
  //double Ph_mag = this->getP_h();
  double Ph_mag = sqrt(Px_h*Px_h + Py_h*Py_h + Pz_h*Pz_h);
  
  double Px_q = - lk->getPx_el();
  double Py_q = - lk->getPy_el();
  double Pz_q = (kEbeam - lk->getPz_el());
  double Pq_mag = sqrt(lk->getNu()*lk->getNu() + lk->getQ2());
  
  double result = (Pz_h*Pz_q + Px_h*Px_q + Py_h*Py_q)/(Pq_mag*Ph_mag);

  if(result > 1){
    std::cout<<" Numerator   = "<<(this->Pz_h * (kEbeam - lk->getPz_el()) - this->Px_h * lk->getPx_el() - this->Py_h * lk->getPy_el())<<std::endl;
    std::cout<<" Denominator = "<<(TMath::Sqrt(lk->getNu()*lk->getNu() + lk->getQ2())*this->getP_h())<<std::endl;;
  }
  return result;
}

double HadronicKinematics::getZh(LeptonicKinematics* lk) {
  // Returns the energy fraction of the particle
  double mass = this->getMass_h(this->PID_h);
  double P_h  = this->getP_h();
  double Nu   = lk->getNu();

  return TMath::Sqrt(mass*mass + P_h*P_h)/Nu;
}
double HadronicKinematics::getPt2(LeptonicKinematics* lk) {
  // Returns the square of the transverse momentum component w.r.t. the virtual photon direction
  double P_h        = this->getP_h();
  double CosThetaPQ = this->getCosThetaPQ(lk);

  return P_h*P_h*(1. - TMath::Power(CosThetaPQ,2));
}
double HadronicKinematics::getPl2(LeptonicKinematics* lk) {
  // Returns the square of the longitudinal momentum component w.r.t. the virtual photon direction
  double P_h        = this->getP_h();
  double CosThetaPQ = this->getCosThetaPQ(lk);
  
  return P_h*P_h*CosThetaPQ*CosThetaPQ;
}
