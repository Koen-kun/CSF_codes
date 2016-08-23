#ifndef FUNCTIONCOLLECTOR_H
#define FUNCTIONCOLLECTOR_H
#include "FunctionClassBase.h"
#include "TFormula.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"

using namespace std;
 
class FunctionCollector {

 public:
  //Data members
  vector<FunctionBase*> FuncVec;
  vector<double> dphvec;
  vector<double> sphvec;
  vector<double> dtvec;
  vector<double> stvec;
  vector<double> dpvec;
  vector<double> spvec;
  vector<double> wvec;

  //Constructor and destructor
  FunctionCollector();
  FunctionCollector(vector<FunctionBase*>);
  ~FunctionCollector() {}

  //Member functions
  void Add(FunctionBase*);
  void FillVarVectors(double var,TString svar);

  //void DrawFormula1D(TString Var,int=10,int=10,int=10,double=1.0);
  //void DrawFormula2D(TString Var1, TString Var2,int=10,int=10,int=10,double=1.0);
  TProfile*   DrawFormula1D(TString Var,int=10,int=10,int=10);
  TProfile2D* DrawFormula2D(TString Var1, TString Var2,int=10,int=10,int=10);
  void SetPositionXY(double x, double y, double &posx, double &posy);

  //Comparing Data with Funtions
  void DataCompare1D(TString Var,int=10);
  void DataCompare2D(TString Var1, TString Var2,int=10);

};


FunctionCollector::FunctionCollector() {
}

FunctionCollector::FunctionCollector(vector<FunctionBase*> Functions) {
  FuncVec = Functions;  
}

void FunctionCollector::Add(FunctionBase* F) {
  FuncVec.push_back(F);
}

void FunctionCollector::FillVarVectors(double var, TString svar){

  if(svar == "dph") dphvec.push_back(var);
  else if(svar == "sph") sphvec.push_back(var);
  else if(svar == "dt")  dtvec.push_back(var);
  else if(svar == "st")  stvec.push_back(var);
  else if(svar == "dp")  dpvec.push_back(var);
  else if(svar == "sp")  spvec.push_back(var);
  else if(svar == "w")   wvec.push_back(var);
  else {
    cout << "Error: Var " << svar << " can't be stored in VarVectors!" << endl;
    return;
  }
}

TProfile* FunctionCollector::DrawFormula1D(TString Var, int Nt, int Np, int Nm){

  double dph, sph, dt, st, dp, sp;   
  Double_t pos(0.);
  Double_t angv[4];
  Double_t momv[2];
  Double_t pi = TMath::Pi();

  Double_t func = 0.;
 
  TProfile* prof;
  TString profname = "prof_" + Var;
  if(Var == "dph" || Var == "sph"){  
    prof = new TProfile(profname,profname,Np,-pi,pi); 
  }
  else if(Var == "dt" || Var == "st"){  
    prof = new TProfile(profname,profname,Nt,-pi,pi); 
  }
  else if(Var == "dp" || Var == "sp"){  
    prof = new TProfile(profname,profname,Nm,-pi,pi); 
  }
  else 
  {
    cout << "No valid argument supplied. Should be: dph/sph/dt/st/dp/sp" << endl;
    prof = new TProfile(profname,profname,1,0,1);
    return prof;
  }
  
  //These loop represents a uniform sampling for the angular and momentum distributions
  for(int dti=0; dti<Nt; dti++){
   for(int sti=0; sti<Nt; sti++){
    for(int dphi=0; dphi<Np; dphi++){
     for(int sphi=0; sphi<Np; sphi++){
      for(int dpi=0;  dpi<Nm;  dpi++){
       for(int spi=0;  spi<Nm;  spi++){

         //Make sure to always pick the center of the "bin"
         dt  = (dti  * (2.*pi/Nt) - pi) + 0.5*(2.*pi/Nt) ;  
         st  = (sti  * (2.*pi/Nt) - pi) + 0.5*(2.*pi/Nt) ;  
         dph = (dphi * (2.*pi/Np) - pi) + 0.5*(2.*pi/Np) ;  
         sph = (sphi * (2.*pi/Np) - pi) + 0.5*(2.*pi/Np) ;  
         dp  = (dpi  * (2.*pi/Nm) - pi) + 0.5*(2.*pi/Nm) ;  
         sp  = (spi  *( 2.*pi/Nm) - pi) + 0.5*(2.*pi/Nm) ;  

         angv[0] = dt;
         angv[1] = st;
         angv[2] = dph;
         angv[3] = sph;
         momv[0] = dp;
         momv[1] = sp;

         func = 1.0;
         for(vector<FunctionBase*>::iterator it = FuncVec.begin(); it != FuncVec.end(); it++){
           func += (*it)->meanC * ((*it)->f->EvalPar(angv,momv));
         }

         if(Var == "dt")  pos = dt;
         if(Var == "st")  pos = st;
         if(Var == "dph") pos = dph;
         if(Var == "sph") pos = sph;
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

  return prof;

}

void FunctionCollector::SetPositionXY(double x, double y, double &posx, double &posy){
  posx = x;
  posy = y;
}

TProfile2D* FunctionCollector::DrawFormula2D(TString Var1, TString Var2, int Nt, int Np, int Nm){

  double dt, st, dph, sph, dp, sp;
  Double_t posx(0.);
  Double_t posy(0.);
  Double_t angv[4];
  Double_t momv[2];
  Double_t pi = TMath::Pi();

  Double_t func = 0.;
 
  double xdo(-pi), ydo(-pi); 

  TProfile2D* prof;
  TString profname = "2dprof_" + Var1 + "_" + Var2;
  if(Var1 == "dt" && Var2 == "st"){
    prof = new TProfile2D(profname,profname,Nt,xdo,pi,Nt,ydo,pi);
  }
  else if((Var1 == "dt" || Var1 == "st") && (Var2 == "dph" || Var2 == "sph")){
    prof = new TProfile2D(profname,profname,Nt,xdo,pi,Np,ydo,pi);
  }
  else if((Var1 == "dph" || Var1 == "sph") && (Var2 == "dt" || Var2 == "st")){
    prof = new TProfile2D(profname,profname,Np,xdo,pi,Nt,ydo,pi); 
  }
  else if((Var1 == "dt" || Var1 == "st") && (Var2 == "dp" || Var2 == "sp")){
    prof = new TProfile2D(profname,profname,Nt,xdo,pi,Nm,ydo,pi); 
  }
  else if((Var1 == "dp" || Var1 == "sp") && (Var2 == "dt" || Var2 == "st")){
    prof = new TProfile2D(profname,profname,Nm,xdo,pi,Nt,ydo,pi); 
  }
  else if(Var1 == "dph" && Var2 == "sph"){
    prof = new TProfile2D(profname,profname,Np,xdo,pi,Np,ydo,pi); 
  }
  else if((Var1 == "dph" || Var1 == "sph") && (Var2 == "dp" || Var2 == "sp")){
    prof = new TProfile2D(profname,profname,Np,xdo,pi,Nm,ydo,pi); 
  }
  else if((Var1 == "dp" || Var1 == "sp") && (Var2 == "dph" || Var2 == "sph")){
    prof = new TProfile2D(profname,profname,Nm,xdo,pi,Np,ydo,pi); 
  }
  else if(Var1 == "dp" && Var2 == "sp"){
    prof = new TProfile2D(profname,profname,Nm,xdo,pi,Nm,ydo,pi); 
  }
  else{
    cout << "No valid arguments supplied. Should be two from these: dt/st/dph/sph/dp/sp" << endl;
    prof = new TProfile2D(profname,profname,1,0,1,1,0,1);
    return prof;
  }

  //These loop represents a uniform sampling for the angular and momentum distributions
  for(int dti=0; dti<Nt; dti++){
   for(int sti=0; sti<Nt; sti++){
    for(int dphi=0; dphi<Np; dphi++){
     for(int sphi=0; sphi<Np; sphi++){
      for(int dpi=0;  dpi<Nm;  dpi++){
       for(int spi=0;  spi<Nm;  spi++){

         //Make sure to always pick the center of the "bin"
         dt  = (dti  * (2.*pi/Nt) - pi) + 0.5*(2.*pi/Nt) ;  
         st  = (sti  * (2.*pi/Nt) - pi) + 0.5*(2.*pi/Nt) ;  
         dph = (dphi * (2.*pi/Np) - pi) + 0.5*(2.*pi/Np) ;  
         sph = (sphi * (2.*pi/Np) - pi) + 0.5*(2.*pi/Np) ;  
         dp  = (dpi  * (2.*pi/Nm) - pi) + 0.5*(2.*pi/Nm) ;  
         sp  = (spi  *( 2.*pi/Nm) - pi) + 0.5*(2.*pi/Nm) ;  

         angv[0] = dt;
         angv[1] = st;
         angv[2] = dph;
         angv[3] = sph;
         momv[0] = dp;
         momv[1] = sp;

         if(Var1 == "dt"  && Var2 == "st")  SetPositionXY(dt,st,posx,posy);  
         if(Var1 == "dt"  && Var2 == "dph") SetPositionXY(dt,dph,posx,posy);
         if(Var1 == "dt"  && Var2 == "sph") SetPositionXY(dt,sph,posx,posy); 
         if(Var1 == "dt"  && Var2 == "dp")  SetPositionXY(dt,dp,posx,posy);  
         if(Var1 == "dt"  && Var2 == "sp")  SetPositionXY(dt,sp,posx,posy);  
         if(Var1 == "st"  && Var2 == "dt")  SetPositionXY(st,dt,posx,posy);  
         if(Var1 == "st"  && Var2 == "dph") SetPositionXY(st,dph,posx,posy);
         if(Var1 == "st"  && Var2 == "sph") SetPositionXY(st,sph,posx,posy); 
         if(Var1 == "st"  && Var2 == "dp")  SetPositionXY(st,dp,posx,posy);  
         if(Var1 == "st"  && Var2 == "sp")  SetPositionXY(st,sp,posx,posy);  
         if(Var1 == "dph" && Var2 == "dt")  SetPositionXY(dph,dt,posx,posy);  
         if(Var1 == "dph" && Var2 == "st")  SetPositionXY(dph,st,posx,posy);
         if(Var1 == "dph" && Var2 == "sph") SetPositionXY(dph,sph,posx,posy); 
         if(Var1 == "dph" && Var2 == "dp")  SetPositionXY(dph,dp,posx,posy);  
         if(Var1 == "dph" && Var2 == "sp")  SetPositionXY(dph,sp,posx,posy);  
         if(Var1 == "sph" && Var2 == "dt")  SetPositionXY(sph,dt,posx,posy);  
         if(Var1 == "sph" && Var2 == "st")  SetPositionXY(sph,st,posx,posy);
         if(Var1 == "sph" && Var2 == "dph") SetPositionXY(sph,dph,posx,posy); 
         if(Var1 == "sph" && Var2 == "dp")  SetPositionXY(sph,dp,posx,posy);  
         if(Var1 == "sph" && Var2 == "sp")  SetPositionXY(sph,sp,posx,posy);  
         if(Var1 == "dp"  && Var2 == "dt")  SetPositionXY(dp,dt,posx,posy);  
         if(Var1 == "dp"  && Var2 == "st")  SetPositionXY(dp,st,posx,posy);
         if(Var1 == "dp"  && Var2 == "dph") SetPositionXY(dp,dph,posx,posy); 
         if(Var1 == "dp"  && Var2 == "sph") SetPositionXY(dp,sph,posx,posy);  
         if(Var1 == "dp"  && Var2 == "sp")  SetPositionXY(dp,sp,posx,posy);  
         if(Var1 == "sp"  && Var2 == "dt")  SetPositionXY(sp,dt,posx,posy);  
         if(Var1 == "sp"  && Var2 == "st")  SetPositionXY(sp,st,posx,posy);
         if(Var1 == "sp"  && Var2 == "dph") SetPositionXY(sp,dph,posx,posy); 
         if(Var1 == "sp"  && Var2 == "sph") SetPositionXY(sp,sph,posx,posy);  
         if(Var1 == "sp"  && Var2 == "dp")  SetPositionXY(sp,dp,posx,posy);  

         func = 1.0;
         for(vector<FunctionBase*>::iterator it = FuncVec.begin(); it != FuncVec.end(); it++){
           func += (*it)->meanC * ((*it)->f->EvalPar(angv,momv));
         }

         prof->Fill(posx,posy,func);
       }
      }
     }
    }
   }
  }

  prof->GetXaxis()->SetTitle(Var1);
  prof->GetYaxis()->SetTitle(Var2);

  return prof;

}

void FunctionCollector::DataCompare1D(TString Var, int N){

  double Nt = 10;
  double Np = 10;
  double Nm = 10;
  Double_t pi = TMath::Pi();
  double xdo = -pi;

  if(Var == "dt" || Var == "st"){
    Nt = N;
  } 
  else if(Var == "dph" || Var == "sph"){
    Np = N;
  } 
  else if(Var == "dp" || Var == "sp"){
    Nm = N;
  } 

  TH1D* hist;
  TString histname = "hist_" + Var;
  hist = new TH1D(histname,histname,N,xdo,pi);

  vector<double> varvec;
  if(Var == "dt")       varvec = dtvec; 
  else if(Var == "st")  varvec = stvec; 
  else if(Var == "dph") varvec = dphvec; 
  else if(Var == "sph") varvec = sphvec; 
  else if(Var == "dp")  varvec = dpvec; 
  else if(Var == "sp")  varvec = spvec; 
  else{
    cout << "Error: Var " << Var << " can't be drawn!" << endl;
    return;
  }

  vector<double>::iterator wit = wvec.begin();
  for(vector<double>::iterator it = varvec.begin(); it != varvec.end(); it++){
    hist->Fill(*it,*wit); 
    wit++;
  }

  TProfile* Funcprof;
  Funcprof = DrawFormula1D(Var,Nt,Np,Nm);

  double Histint;
  double Funcint;

  Histint = hist->Integral();
  Funcint = Funcprof->Integral();
  Funcprof->Scale(Histint/Funcint);

  TString cname = "can_" + Var;
  TCanvas* W = new TCanvas(cname,cname,1100,1000);
  W->Divide(1); 
  W->cd(1); 

  Funcprof->SetLineColor(2);
  hist->Draw();
  Funcprof->Draw("same");

}

void FunctionCollector::DataCompare2D(TString Var1, TString Var2, int N){

  double Nt = 10;
  double Np = 10;
  double Nm = 10;
  Double_t pi = TMath::Pi();
  double xdo = -pi;
  double ydo = -pi;

  if(Var1 == "dt" || Var1 == "st"){
    Nt = N;
  } 
  else if(Var1 == "dph" || Var1 == "sph"){
    Np = N;
  } 
  else if(Var1 == "dp" || Var1 == "sp"){
    Nm = N;
  } 

  TH2D* hist;
  TString histname = "2dhist_" + Var1 + "_" + Var2;
  hist = new TH2D(histname,histname,N,xdo,pi,N,ydo,pi);

  vector<double> varvec1;
  vector<double> varvec2;

  if(Var1 == "dt")       varvec1 = dtvec; 
  else if(Var1 == "st")  varvec1 = stvec; 
  else if(Var1 == "dph") varvec1 = dphvec; 
  else if(Var1 == "sph") varvec1 = sphvec; 
  else if(Var1 == "dp")  varvec1 = dpvec; 
  else if(Var1 == "sp")  varvec1 = spvec; 
  else{
    cout << "Error: Var1 " << Var1 << " can't be drawn!" << endl;
    return;
  }
  if(Var2 == "dt")       varvec2 = dtvec; 
  else if(Var2 == "st")  varvec2 = stvec; 
  else if(Var2 == "dph") varvec2 = dphvec; 
  else if(Var2 == "sph") varvec2 = sphvec; 
  else if(Var2 == "dp")  varvec2 = dpvec; 
  else if(Var2 == "sp")  varvec2 = spvec; 
  else{
    cout << "Error: Var2 " << Var2 << " can't be drawn!" << endl;
    return;
  }

  vector<double>::iterator wit  = wvec.begin();
  vector<double>::iterator var2 = varvec2.begin();
  for(vector<double>::iterator var1 = varvec1.begin(); var1 != varvec1.end(); var1++){
    hist->Fill(*var1,*var2,*wit);
    var2++; 
    wit++;
  }

  TProfile2D* Funcprof;
  Funcprof = DrawFormula2D(Var1,Var2,Nt,Np,Nm);

  double Histint;
  double Funcint;
  Histint = hist->Integral();
  Funcint = Funcprof->Integral();
  Funcprof->Scale(Histint/Funcint);

  TH2D* Rat = (TH2D*) hist->Clone();
  Rat->Divide(Funcprof);

  TString cname = "2Dcan_" + Var1 + "_" + Var2;
  TCanvas* W = new TCanvas(cname,cname,1100,1000);
  W->Divide(3,1); 
  W->cd(1); 
  hist->Draw("colz");

  W->cd(2);
  Funcprof->Draw("colz");

  W->cd(3);
  Rat->SetMaximum(2.0);
  Rat->SetMinimum(0.0);
  Rat->Draw("colz");

}

#endif

