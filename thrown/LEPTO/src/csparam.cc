#include <stdio.h>
#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TStyle.h>

using namespace std;

void csparam(double*, double*, double*, double*, double*, double*, double*);
Double_t FermiCS(Double_t*, Double_t*);

extern "C" // csparam_() is to be called externally in F
{
  void csparam_(double *choice, double *res, double *cost, double *phi,
double *kx, double *ky, double *kz){
    csparam(choice, res, cost, phi, kx, ky, kz);
 } 
}

extern "C" int qp1_();
int main() {
  return qp1_();
}

void  csparam(double *choice, double *res, double *cost, double *phi, 
	      double *kx, double *ky, double *kz ){  

  TF1 *fun = new TF1("fun", FermiCS,0,5,1);
  if(choice[0] == 2.)fun->SetParameter(0,0);    //D
  if(choice[0] == 12.)fun->SetParameter(0,1);  //C 
  if(choice[0] == 56.)fun->SetParameter(0,2);  //Fe 
  if(choice[0] == 207.)fun->SetParameter(0,3); //Pb
  
  fun->SetNpx(1000);
  
  gRandom = new TRandom3(0); //! magic line, takes time as a random seed
  //gRandom = new TRandom3(8968);
  res[0] = (fun->GetRandom(0,0.5));
  
  cost[0] = gRandom->Uniform(-1,1);
  phi[0] =  gRandom->Uniform(0, 2*TMath::Pi());
  
  kx[0] = res[0]*TMath::Sin(TMath::ACos(cost[0]))*TMath::Sin(phi[0]);
  ky[0] = res[0]*TMath::Sin(TMath::ACos(cost[0]))*TMath::Cos(phi[0]);
  kz[0] = res[0]*cost[0];
  
}

const float p[4][11] = {
  // parameters numbers        1       2      3      4       5      6        7       8        9        10          11 
  {157.4 ,  1.24 , 18.3 , 0.234 , 1.27 ,  0.0 , 0.00623 , 0.22 ,    0.0 , 27.7387  ,      0.0  },
  { 2.61 ,  2.66 , 3.54 ,  0.0  ,  0.0 ,  0.0 ,  0.426  , 1.60 , 0.0237 ,    0.22  ,  2.70513  },
  { 3.57 ,  4.97 ,  0.0 , 19.8  , 15.0 ,  0.0 ,  0.230  , 1.20 , 0.0286 ,    0.22  ,  2.83756  },
  {  1.8 ,  4.77 ,  0.0 , 25.5  ,  0.0 , 40.3 ,  0.275  , 1.01 , 0.0304 ,    0.22  ,  2.01735  }
};
// p[0] -> D; p[1] -> C; p[2] -> Fe; p[3] -> Pb; 

Double_t FermiCS(Double_t *x, Double_t *c){ // c is param choice[0]
  Double_t result =0.;
  Double_t gev2fm=5.0677;
  Double_t k=x[0]*gev2fm;
  // Double_t k=x[0];
  int i = TMath::FloorNint(c[0]);
  // std::cout<<__LINE__<<"\t"<<std::endl;
  if(i!=0){
    Double_t n0=p[i][0]*TMath::Exp(-p[i][1]*k*k)*(1.0+p[i][2]*k*k+p[i][3]*TMath::Power(k,4)+p[i][4]*TMath::Power(k,6)+p[i][5]*TMath::Power(k,8));
    Double_t n1=p[i][6]*TMath::Exp(-p[i][7]*k*k)+p[i][8]*TMath::Exp(-p[i][9]*k*k);
    result = k*k*(n1 +n0)*gev2fm;// /p[i][10];
  }
  if(i==0){
    Double_t n0 = p[i][0]*TMath::Exp(-p[i][1]*k*k)/TMath::Power(1+p[i][2]*k*k,2)+p[i][3]*TMath::Exp(-p[i][4]*k*k)/TMath::Power(1+p[i][5]*k*k,2)+p[i][6]*TMath::Exp(-p[i][7]*k*k)/TMath::Power(1+p[i][8]*k*k,2);
    result = k*k*n0*gev2fm;// /p[i][9];
  }
  //std::cout<<__LINE__<<"\t"<<std::endl;
  // Double_t gev2fm=5.0677;
  // Double_t k=x[0]*gev2fm;
  return result;
}
