#ifndef FUNCTIONSTRINGCOLLECTOR_H
#define FUNCTIONSTRINGCOLLECTOR_H
#include "FunctionClassBase.h"
#include "TFormula.h"

using namespace std;
 
class FunctionStringCollector {

 public:
  //Data members
  TString funcstring;  
  TString type;

  //Constructor and destructor
  FunctionStringCollector();
  ~FunctionStringCollector() {}

  //Member functions
  void Add(FunctionBase*);
  void DrawFormula1D(TString Var,int=10,int=10,int=10,double=1.0);
  void DrawFormula2D(TString Var1, TString Var2,int=10,int=10,int=10,double=1.0);
  void SetPositionXY(double x, double y, double &posx, double &posy);

};

FunctionStringCollector::FunctionStringCollector() {

  funcstring = "1.";
  type       = "Constant";

  //cout << "The function is now: " << funcstring << endl;  

}

void FunctionStringCollector::Add(FunctionBase* F) {

  stringstream scoef;
  TString coef;
  scoef << F->meanC;
  scoef >> coef;

  funcstring = funcstring + " + " + coef + "*" + F->formula;

  //cout << "The function is now: " << funcstring << endl;  

}

void FunctionStringCollector::DrawFormula1D(TString Var, int Nt, int Np, int Nm, double Amount){
  //Var = t1/t2/ph1/ph2/dp/sp
  double t1, t2, ph1, ph2, dp, sp;   
  Double_t pos(0.);
  Double_t angv[4];
  Double_t momv[2];
  Double_t pi = TMath::Pi();
  Double_t BinFactor = 1.0;

  Double_t func = 0.;

  //Creating function object
  Function F("F",this->funcstring,1.,1.,"Default");
 
  TProfile* prof;
  if(Var == "t1" || Var == "t2"){  
    prof = new TProfile("prof","",Nt,0,pi);
    BinFactor = Amount/Nt; 
  }
  else if(Var == "ph1" || Var == "ph2"){  
    prof = new TProfile("prof","",Np,-pi,pi); 
    BinFactor = Amount/Np;
  }
  else if(Var == "dp" || Var == "sp"){  
    prof = new TProfile("prof","",Nm,-pi,pi); 
    BinFactor = Amount/Nm;
  }
  else 
  {
    cout << "No valid argument supplied. Should be: t1/t2/ph1/ph2/dp/sp" << endl;
    return;
  }
  
  //Normalization factor
  if(Amount == 1.0) BinFactor = 1.0;

  //These loop represents a uniform sampling for the angular and momentum distributions
  for(int t1i=0; t1i<Nt; t1i++){
   for(int t2i=0; t2i<Nt; t2i++){
    for(int p1i=0; p1i<Np; p1i++){
     for(int p2i=0; p2i<Np; p2i++){
      for(int m1i=0; m1i<Nm; m1i++){
       for(int m2i=0; m2i<Nm; m2i++){

         //Make sure to always pick the center of the "bin"
         t1  = t1i*(pi/Nt) + 0.5*(pi/Nt);
         t2  = t2i*(pi/Nt) + 0.5*(pi/Nt);
         ph1 = (2.*p1i*(pi/Np) - pi) + 0.5*(2.*pi/Np) ;  
         ph2 = (2.*p2i*(pi/Np) - pi) + 0.5*(2.*pi/Np) ;  
         dp  = (2.*m1i*(pi/Nm) - pi) + 0.5*(2.*pi/Nm) ;  
         sp  = (2.*m2i*(pi/Nm) - pi) + 0.5*(2.*pi/Nm) ;  

         angv[0] = t1;
         angv[1] = t2;
         angv[2] = ph1;
         angv[3] = ph2;
         momv[0] = dp;
         momv[1] = sp;

         func = (F.f->EvalPar(angv,momv))*BinFactor;

         if(Var == "t1")  pos = t1;
         if(Var == "t2")  pos = t2;
         if(Var == "ph1") pos = ph1;
         if(Var == "ph2") pos = ph2;
         if(Var == "dp")  pos = dp;
         if(Var == "sp")  pos = sp;

         prof->Fill(pos,func); 
       }
      }
     }
    }
   }
  }

  prof->GetXaxis()->SetTitle(Var);
  prof->Draw();

}

void FunctionStringCollector::SetPositionXY(double x, double y, double &posx, double &posy){
  posx = x;
  posy = y;
}

void FunctionStringCollector::DrawFormula2D(TString Var1, TString Var2, int Nt, int Np, int Nm, double Amount){
  //Var = t1/t2/ph1/ph2/dp/sp
  double t1, t2, ph1, ph2, dp, sp;
  Double_t posx(0.);
  Double_t posy(0.);
  Double_t angv[4];
  Double_t momv[2];
  Double_t pi = TMath::Pi();
  Double_t BinFactor = 1.0; 

  Double_t func = 0.;
 
  //Creating function object
  Function F("F",this->funcstring,1.,1.,"Default");

  double xdo(-pi), ydo(-pi); 
  if(Var1 == "t1" || Var1 == "t2") xdo = 0;  
  if(Var2 == "t1" || Var2 == "t2") ydo = 0;  

  TProfile2D* prof;
  if(Var1 == "t1" && Var2 == "t2"){
    prof = new TProfile2D("prof","",Nt,xdo,pi,Nt,ydo,pi);
    BinFactor = Amount/(Nt*Nt); 
  }
  else if((Var1 == "t1" || Var1 == "t2") && (Var2 == "ph1" || Var2 == "ph2")){
    prof = new TProfile2D("prof","",Nt,xdo,pi,Np,ydo,pi);
    BinFactor = Amount/(Nt*Np); 
  }
  else if((Var1 == "ph1" || Var1 == "ph2") && (Var2 == "t1" || Var2 == "t2")){
    prof = new TProfile2D("prof","",Np,xdo,pi,Nt,ydo,pi); 
    BinFactor = Amount/(Np*Nt);
  }
  else if((Var1 == "t1" || Var1 == "t2") && (Var2 == "dp" || Var2 == "sp")){
    prof = new TProfile2D("prof","",Nt,xdo,pi,Nm,ydo,pi); 
    BinFactor = Amount/(Nt*Nm);
  }
  else if((Var1 == "dp" || Var1 == "sp") && (Var2 == "t1" || Var2 == "t2")){
    prof = new TProfile2D("prof","",Nm,xdo,pi,Nt,ydo,pi); 
    BinFactor = Amount/(Nm*Nt);
  }
  else if(Var1 == "ph1" && Var2 == "ph2"){
    prof = new TProfile2D("prof","",Np,xdo,pi,Np,ydo,pi); 
    BinFactor = Amount/(Np*Np);
  }
  else if((Var1 == "ph1" || Var1 == "ph2") && (Var2 == "dp" || Var2 == "sp")){
    prof = new TProfile2D("prof","",Np,xdo,pi,Nm,ydo,pi); 
    BinFactor = Amount/(Np*Nm);
  }
  else if((Var1 == "dp" || Var1 == "sp") && (Var2 == "ph1" || Var2 == "ph2")){
    prof = new TProfile2D("prof","",Nm,xdo,pi,Np,ydo,pi); 
    BinFactor = Amount/(Nm*Np);
  }
  else if(Var1 == "dp" && Var2 == "sp"){
    prof = new TProfile2D("prof","",Nm,xdo,pi,Nm,ydo,pi); 
    BinFactor = Amount/(Nm*Nm);
  }
  else{
    cout << "No valid arguments supplied. Should be two from these: t1/t2/ph1/ph2/dp/sp" << endl;
    return;
  }

  if(Amount == 1.0) BinFactor = 1.0;

  //These loop represents a uniform sampling for the angular and momentum distributions
  for(int t1i=0; t1i<Nt; t1i++){
   for(int t2i=0; t2i<Nt; t2i++){
    for(int p1i=0; p1i<Np; p1i++){
     for(int p2i=0; p2i<Np; p2i++){
      for(int m1i=0; m1i<Nm; m1i++){
       for(int m2i=0; m2i<Nm; m2i++){

         //Make sure to always pick the center of the "bin"
         t1  = t1i*(pi/Nt) + 0.5*(pi/Nt);
         t2  = t2i*(pi/Nt) + 0.5*(pi/Nt);
         ph1 = (2.*p1i*(pi/Np) - pi) + 0.5*(2.*pi/Np) ;  
         ph2 = (2.*p2i*(pi/Np) - pi) + 0.5*(2.*pi/Np) ;  
         dp  = (2.*m1i*(pi/Nm) - pi) + 0.5*(2.*pi/Nm) ;  
         sp  = (2.*m2i*(pi/Nm) - pi) + 0.5*(2.*pi/Nm) ;  

         angv[0] = t1;
         angv[1] = t2;
         angv[2] = ph1;
         angv[3] = ph2;
         momv[0] = dp;
         momv[1] = sp;

         if(Var1 == "t1"  && Var2 == "t2")  SetPositionXY(t1,t2,posx,posy);  
         if(Var1 == "t1"  && Var2 == "ph1") SetPositionXY(t1,ph1,posx,posy);
         if(Var1 == "t1"  && Var2 == "ph2") SetPositionXY(t1,ph2,posx,posy); 
         if(Var1 == "t1"  && Var2 == "dp")  SetPositionXY(t1,dp,posx,posy);  
         if(Var1 == "t1"  && Var2 == "sp")  SetPositionXY(t1,sp,posx,posy);  
         if(Var1 == "t2"  && Var2 == "t1")  SetPositionXY(t2,t1,posx,posy);  
         if(Var1 == "t2"  && Var2 == "ph1") SetPositionXY(t2,ph1,posx,posy);
         if(Var1 == "t2"  && Var2 == "ph2") SetPositionXY(t2,ph2,posx,posy); 
         if(Var1 == "t2"  && Var2 == "dp")  SetPositionXY(t2,dp,posx,posy);  
         if(Var1 == "t2"  && Var2 == "sp")  SetPositionXY(t2,sp,posx,posy);  
         if(Var1 == "ph1" && Var2 == "t1")  SetPositionXY(ph1,t1,posx,posy);  
         if(Var1 == "ph1" && Var2 == "t2")  SetPositionXY(ph1,t2,posx,posy);
         if(Var1 == "ph1" && Var2 == "ph2") SetPositionXY(ph1,ph2,posx,posy); 
         if(Var1 == "ph1" && Var2 == "dp")  SetPositionXY(ph1,dp,posx,posy);  
         if(Var1 == "ph1" && Var2 == "sp")  SetPositionXY(ph1,sp,posx,posy);  
         if(Var1 == "ph2" && Var2 == "t1")  SetPositionXY(ph2,t1,posx,posy);  
         if(Var1 == "ph2" && Var2 == "t2")  SetPositionXY(ph2,t2,posx,posy);
         if(Var1 == "ph2" && Var2 == "ph1") SetPositionXY(ph2,ph1,posx,posy); 
         if(Var1 == "ph2" && Var2 == "dp")  SetPositionXY(ph2,dp,posx,posy);  
         if(Var1 == "ph2" && Var2 == "sp")  SetPositionXY(ph2,sp,posx,posy);  
         if(Var1 == "dp"  && Var2 == "t1")  SetPositionXY(dp,t1,posx,posy);  
         if(Var1 == "dp"  && Var2 == "t2")  SetPositionXY(dp,t2,posx,posy);
         if(Var1 == "dp"  && Var2 == "ph1") SetPositionXY(dp,ph1,posx,posy); 
         if(Var1 == "dp"  && Var2 == "ph2") SetPositionXY(dp,ph2,posx,posy);  
         if(Var1 == "dp"  && Var2 == "sp")  SetPositionXY(dp,sp,posx,posy);  
         if(Var1 == "sp"  && Var2 == "t1")  SetPositionXY(sp,t1,posx,posy);  
         if(Var1 == "sp"  && Var2 == "t2")  SetPositionXY(sp,t2,posx,posy);
         if(Var1 == "sp"  && Var2 == "ph1") SetPositionXY(sp,ph1,posx,posy); 
         if(Var1 == "sp"  && Var2 == "ph2") SetPositionXY(sp,ph2,posx,posy);  
         if(Var1 == "sp"  && Var2 == "dp")  SetPositionXY(sp,dp,posx,posy);  

         func = (F.f->EvalPar(angv,momv))*BinFactor;
         prof->Fill(posx,posy,func);
       }
      }
     }
    }
   }
  }

  prof->GetXaxis()->SetTitle(Var1);
  prof->GetYaxis()->SetTitle(Var2);
  prof->Draw("colz");

}

#endif

