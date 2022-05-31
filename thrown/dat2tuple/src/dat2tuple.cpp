// Script that takes a .dat file and transforms it into a TNtupe in a ROOT file

// Req: .dat file has to have the following format
//      <event_index> <particle_id> <parent_id> <px> <py> <pz> <E> <x> <y> <z>
// Pro Tip: Use lepto2dat.pl (wink wink)

// author : Esteban Molina (May 2022)

#include "dat2tuple.h"
#include "TFile.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>

int main(int argc, char** argv){

  if(argc < 3){
    std::cout<<"Number of arguments is not correct!"<<std::endl;
    return 0;
  }
  
  // Input variables
  const char* file_in  = argv[1];
  const char* file_out = argv[2];

  // Open dat file
  std::ifstream file(file_in);

  // Create Tree that reads file
  TTree* t = new TTree("ntuple_thrown_raw","");
  // Make the tree read the .dat file
  t->ReadFile(file_in,"event_index/F:PID:parent_PID:Px:Py:Pz:E:x:y:z");

  // Create final ntuples
  TNtuple* ntuple_thrown_electrons	= new TNtuple("ntuple_thrown_electrons","","Q2:Xb:Nu:ThetaLab:PhiLab:P:Px:Py:Pz");
  TNtuple* ntuple_thrown		= new TNtuple("ntuple_thrown"          ,"","Q2:Xb:Nu:Zh:Pt2:Pl2:PhiPQ:ThetaLab:PhiLab:P:Px:Py:Pz:ThetaLab_el:PhiLab_el:P_el:Px_el:Py_el:Pz_el");
  
  //Process the tree
  Float_t event_index, PID, parent_PID, Px, Py, Pz, E, x, y, z;
  Float_t elP[3];
  setBranchesAddresses(t, &event_index, &PID, &parent_PID, &Px, &Py, &Pz, &E, &x, &y, &z);

  Int_t Nentries = t->GetEntries();
  for(Int_t entry = 0 ; entry < Nentries ; entry++){
    t->GetEntry(entry);
    if(PID==11){
      // Calculate leptonic variables
      LeptonicKinematics lk(Px,Py,Pz);
      elP[0] = Px;
      elP[1] = Py;
      elP[2] = Pz;

      ntuple_thrown_electrons->Fill(lk.getQ2(), lk.getXb(), lk.getNu(), lk.getThetaLab_el(), lk.getPhiLab_el(), lk.getP_el(), Px, Py, Pz);
    }
    else{
      // Calculate hadronic variables
      LeptonicKinematics lk(elP[0],elP[1],elP[2]);
      HadronicKinematics hk(Px,Py,Pz,PID);

      Float_t vars_h[19] = {lk.getQ2(), lk.getXb(), lk.getNu(), hk.getZh(&lk), hk.getPt2(&lk), hk.getPl2(&lk), hk.getPhiPQ(&lk), hk.getThetaLab_h(), hk.getPhiLab_h(), hk.getP_h(), hk.getPx_h(), hk.getPy_h(), hk.getPz_h(), lk.getThetaLab_el(), lk.getPhiLab_el(),
  			    lk.getP_el(), lk.getPx_el(), lk.getPy_el(), lk.getPz_el()};
      
      ntuple_thrown->Fill(vars_h);
    }
  }

  // Create target root file
  TFile* f = new TFile(file_out,"RECREATE");

  f->cd();
  t->Write();
  ntuple_thrown->Write();
  ntuple_thrown_electrons->Write();
  f->Close();
  
  gROOT->cd();
  delete f;
  delete t;
  delete ntuple_thrown;
  delete ntuple_thrown_electrons;

  return 0;
}
