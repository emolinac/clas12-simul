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
  t->ReadFile(file_in,"event_index/D:PID:parent_PID:Px:Py:Pz:E:x:y:z");

  // Create final ntuples
  TNtuple* ntuple_thrown_electrons	= new TNtuple("ntuple_thrown_electrons","","Q2:x_{bjorken}:#nu:W:y:#theta:#phi:p:p_{x}:p_{y}:p_{z}:vz");
  TNtuple* ntuple_thrown	         	= new TNtuple("ntuple_thrown"          ,"","Q2:x_{bjorken}:#nu:W:y:z_{h}:Pt2:Pl2:#theta_{PQ}:#phi_{PQ}:#theta:#phi:p:p_{x}:p_{y}:p_{z}:#theta_{el}:#phi_{el}:p_{el}:p_{xel}:p_{yel}:p_{zel}:pid");
  
  //Process the tree
  Double_t event_index, PID, parent_PID, Px, Py, Pz, E, x, y, z;
  Double_t elP[3];
  setBranchesAddresses(t, &event_index, &PID, &parent_PID, &Px, &Py, &Pz, &E, &x, &y, &z);

  Int_t Nentries = t->GetEntries();
  for(Int_t entry1 = 0 ; entry1 < Nentries ; entry1++){
    t->GetEntry(entry1);
    if(PID==11 && parent_PID==0){
      // Calculate leptonic variables
      LeptonicKinematics lk(Px,Py,Pz);
      elP[0] = Px;
      elP[1] = Py;
      elP[2] = Pz;

      ntuple_thrown_electrons->Fill(lk.getQ2(), lk.getXb(), lk.getNu(), lk.getW(), lk.gety(), lk.getThetaLab_el(), lk.getPhiLab_el(), lk.getP_el(), Px, Py, Pz, z);
    }
    else if(PID != 11 && PID != 22 && PID !=-11){
      // Calculate hadronic variables
      LeptonicKinematics lk(elP[0],elP[1],elP[2]);
      HadronicKinematics hk(Px,Py,Pz,PID);

      float vars_h[23] = {(float) lk.getQ2(), (float) lk.getXb(), (float) lk.getNu(), (float) lk.getW(), (float) lk.gety(), (float) hk.getZh(&lk), (float) hk.getPt2(&lk),
			  (float) hk.getPl2(&lk), (float) hk.getThetaPQ(&lk), (float) hk.getPhiPQ(&lk), (float) hk.getThetaLab_h(),
			  (float) hk.getPhiLab_h(), (float) hk.getP_h(), (float) hk.getPx_h(), (float) hk.getPy_h(), (float) hk.getPz_h(),
			  (float) lk.getThetaLab_el(), (float) lk.getPhiLab_el(), (float) lk.getP_el(), (float) lk.getPx_el(), (float) lk.getPy_el(), (float) lk.getPz_el(), (float) PID};

      ntuple_thrown->Fill(vars_h);
    }
  }

  // Create target root file
  TFile* f = new TFile(file_out,"RECREATE");

  f->cd();
  //  t->Write();
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
