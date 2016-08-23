#ifndef FUNCTIONTYPECOLLECTOR_H
#define FUNCTIONTYPECOLLECTOR_H
#include "FunctionClassBase.h"
#include "TFormula.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH1.h"
#include "TH2.h"
#include <map>

using namespace std;
 
class FunctionTypeCollector {

 public:
  //Data members
  map<TString,TString> FuncStringMap;
  map<TString,vector<FunctionBase*> > FuncVectorMap;

  //Constructor and destructor
  FunctionTypeCollector();
  FunctionTypeCollector(vector<FunctionBase*>);
  ~FunctionTypeCollector() {}

  //Member functions
  void Add(FunctionBase*);
  void DrawFormula1D(TString Var,int=10,double=1.0);
  void DrawFormula2D(TString Var1, TString Var2,int=10,double=1.0);

};


FunctionTypeCollector::FunctionTypeCollector() {
}

FunctionTypeCollector::FunctionTypeCollector(vector<FunctionBase*> Functions) {

  for(vector<FunctionBase*>::iterator it = Functions.begin(); it != Functions.end(); it++){
    //Adding the function when it's not created yet
    if(FuncStringMap.find((*it)->type) != FuncStringMap.end()){
      TString oldfunc, newfunc;
      oldfunc = FuncStringMap[(*it)->type];   
      stringstream scoef;
      TString coef;
      scoef << (*it)->meanC;
      scoef >> coef;
      newfunc = oldfunc + " + " + coef + "*" + (*it)->formula;
      FuncStringMap[(*it)->type] = newfunc;
      cout << "New Function of type: " << (*it)->type << ": " << FuncStringMap[(*it)->type] << endl;
      //Adding Function to Vector for corresponding type
      FuncVectorMap[(*it)->type].push_back(*it);
    }
    else if(FuncStringMap.find((*it)->type) == FuncStringMap.end()){
      stringstream scoef;
      TString coef;
      scoef << (*it)->meanC;
      scoef >> coef;
      FuncStringMap[(*it)->type] = "1 + " + coef + "*" + (*it)->formula;
      cout << "New Function of type: " << (*it)->type << ": " << FuncStringMap[(*it)->type] << endl;
      //Adding Function to Vector for corresponding type
      FuncVectorMap[(*it)->type].push_back(*it);
    }
  }

}

void FunctionTypeCollector::Add(FunctionBase* F) {

  if(FuncStringMap.find(F->type) != FuncStringMap.end()){
    TString oldfunc, newfunc;
    oldfunc = FuncStringMap[F->type];   
    stringstream scoef;
    TString coef;
    scoef << F->meanC;
    scoef >> coef;
    newfunc = oldfunc + " + " + coef + "*" + F->formula;
    FuncStringMap[F->type] = newfunc;
    cout << "New Function of type: " << F->type << ": " << FuncStringMap[F->type] << endl;
    //Adding Function to Vector for corresponding type
    FuncVectorMap[F->type].push_back(F);
  }
  else if(FuncStringMap.find(F->type) == FuncStringMap.end()){
    stringstream scoef;
    TString coef;
    scoef << F->meanC;
    scoef >> coef;
    FuncStringMap[F->type] = "1 + " + coef + "*" + F->formula;
    cout << "New Function of type: " << F->type << ": " << FuncStringMap[F->type] << endl;
    //Adding Function to Vector for corresponding type
    FuncVectorMap[F->type].push_back(F);
  }

}

void FunctionTypeCollector::DrawFormula1D(TString Var, int N, double Amount){
  //Var = t1/t2/ph1/ph2/dp/sp
  double t1, t2, ph1, ph2, dp, sp;   
  Double_t pos(0.);
  Double_t angv[4];
  Double_t momv[2];
  Double_t pi = TMath::Pi();
  Double_t BinFactor = 1.0;
  TString TypeSelect = "";

  Double_t func = 0.;

  TProfile* prof;
  if(Var == "t1" || Var == "t2"){  
    prof = new TProfile("prof","",N,0,pi);
    BinFactor  = Amount/N;
    TypeSelect = "Theta"; 
    ph1 = 0.0;
    ph2 = 0.0;
    dp  = 0.0;
    sp  = 0.0;
  }
  else if(Var == "ph1" || Var == "ph2"){  
    prof = new TProfile("prof","",N,-pi,pi); 
    BinFactor  = Amount/N;
    TypeSelect = "Phi";
    t1  = 0.5*pi;
    t2  = 0.5*pi;
    dp  = 0.0;
    sp  = 0.0;
  }
  else if(Var == "dp" || Var == "sp"){  
    prof = new TProfile("prof","",N,-pi,pi); 
    BinFactor  = Amount/N;
    TypeSelect = "Mom";
    t1  = 0.5*pi;
    t2  = 0.5*pi;
    ph1 = 0.0;
    ph2 = 0.0; 
  }
  else 
  {
    cout << "No valid argument supplied. Should be: t1/t2/ph1/ph2/dp/sp" << endl;
    return;
  }

  //Selecting Function Vector type
  vector<FunctionBase*> FuncSelectVec;

  //Checking if Function of Type TypeSelect exists 
  if(FuncVectorMap.find(TypeSelect) != FuncVectorMap.end()){
    FuncSelectVec = FuncVectorMap[TypeSelect];
  }
  else if(FuncVectorMap.find(TypeSelect) == FuncVectorMap.end()){
     cout << "Error! Function of type " << TypeSelect << " doesn't exist in the FunctionVectorMap!" << endl;
     return;
  }
 
  //Normalization factor
  if(Amount == 1.0) BinFactor = 1.0;

  //These loop represents a uniform sampling for the angular and momentum distributions
  for(int var1i=0; var1i<N; var1i++){
   for(int var2i=0; var2i<N; var2i++){

     //Make sure to always pick the center of the "bin"
     if(TypeSelect == "Theta"){
       t1  = var1i*(pi/N) + 0.5*(pi/N);
       t2  = var2i*(pi/N) + 0.5*(pi/N);
     }
     if(TypeSelect == "Phi"){
       ph1 = (2.*var1i*(pi/N) - pi) + 0.5*(2.*pi/N) ;  
       ph2 = (2.*var2i*(pi/N) - pi) + 0.5*(2.*pi/N) ;  
     }
     if(TypeSelect == "Mom"){
       dp  = (2.*var1i*(pi/N) - pi) + 0.5*(2.*pi/N) ;  
       sp  = (2.*var2i*(pi/N) - pi) + 0.5*(2.*pi/N) ;  
     }

     angv[0] = t1;
     angv[1] = t2;
     angv[2] = ph1;
     angv[3] = ph2;
     momv[0] = dp;
     momv[1] = sp;

     func = 1.0;
     for(vector<FunctionBase*>::iterator it = FuncSelectVec.begin(); it != FuncSelectVec.end(); it++){
       func += (*it)->meanC * ((*it)->f->EvalPar(angv,momv))*BinFactor;
     }

     if(Var == "t1")  pos = t1;
     if(Var == "t2")  pos = t2;
     if(Var == "ph1") pos = ph1;
     if(Var == "ph2") pos = ph2;
     if(Var == "dp")  pos = dp;
     if(Var == "sp")  pos = sp;

     prof->Fill(pos,func); 
   }
  }

  prof->GetXaxis()->SetTitle(Var);
  prof->Draw();

}

void FunctionTypeCollector::DrawFormula2D(TString Var1, TString Var2, int N, double Amount){
  //Var = t1/t2/ph1/ph2/dp/sp
  double t1, t2, ph1, ph2, dp, sp;
  Double_t angv[4];
  Double_t momv[2];
  Double_t pi = TMath::Pi();
  TString TypeSelect = "";

  Double_t func = 0.;
 
  double xdo(-pi), ydo(-pi); 
  if(Var1 == "t1" || Var1 == "t2") xdo = 0;  
  if(Var2 == "t1" || Var2 == "t2") ydo = 0;  

  TH2D* hist;
  if(Var1 == "t1" && Var2 == "t2"){
    hist = new TH2D("hist","",N,xdo,pi,N,ydo,pi);
    TypeSelect = "Theta"; 
    ph1 = 0.0;
    ph2 = 0.0;
    dp  = 0.0;
    sp  = 0.0;
  }
  else if(Var1 == "ph1" && Var2 == "ph2"){
    hist = new TH2D("hist","",N,xdo,pi,N,ydo,pi); 
    TypeSelect = "Phi";
    t1  = 0.5*pi;
    t2  = 0.5*pi;
    dp  = 0.0;
    sp  = 0.0;
  }
  else if(Var1 == "dp" && Var2 == "sp"){
    hist = new TH2D("hist","",N,xdo,pi,N,ydo,pi); 
    TypeSelect = "Mom";
    t1  = 0.5*pi;
    t2  = 0.5*pi;
    ph1 = 0.0;
    ph2 = 0.0;
  }
  else{
    cout << "No valid arguments supplied. Should be two from these: t1/t2/ph1/ph2/dp/sp" << endl;
    cout << "This class only accepts: t1 vs t2 / ph1 vs ph2 / dp vs sp" << endl;
    return;
  }

  //Selecting Function Vector type
  vector<FunctionBase*> FuncSelectVec;

  //Checking if Function of Type TypeSelect exists 
  if(FuncVectorMap.find(TypeSelect) != FuncVectorMap.end()){
    FuncSelectVec = FuncVectorMap[TypeSelect];
  }
  else if(FuncVectorMap.find(TypeSelect) == FuncVectorMap.end()){
     cout << "Error! Function of type " << TypeSelect << " doesn't exist in the FunctionVectorMap!" << endl;
     return;
  }

  //These loop represents a uniform sampling for the angular and momentum distributions
  for(int var1i=0; var1i<N; var1i++){
   for(int var2i=0; var2i<N; var2i++){

    //Make sure to always pick the center of the "bin"
    if(TypeSelect == "Theta"){
       t1  = var1i*(pi/N) + 0.5*(pi/N);
       t2  = var2i*(pi/N) + 0.5*(pi/N);
     }
     if(TypeSelect == "Phi"){
       ph1 = (2.*var1i*(pi/N) - pi) + 0.5*(2.*pi/N) ;  
       ph2 = (2.*var2i*(pi/N) - pi) + 0.5*(2.*pi/N) ;  
     }
     if(TypeSelect == "Mom"){
       dp  = (2.*var1i*(pi/N) - pi) + 0.5*(2.*pi/N) ;  
       sp  = (2.*var2i*(pi/N) - pi) + 0.5*(2.*pi/N) ;  
     }

     angv[0] = t1;
     angv[1] = t2;
     angv[2] = ph1;
     angv[3] = ph2;
     momv[0] = dp;
     momv[1] = sp;

     func = 1.0;
     for(vector<FunctionBase*>::iterator it = FuncSelectVec.begin(); it != FuncSelectVec.end(); it++){
       func += (*it)->meanC * ((*it)->f->EvalPar(angv,momv));
     }
     hist->SetBinContent(var1i+1,var2i+1,func);

   }
  }

  hist->GetXaxis()->SetTitle(Var1);
  hist->GetYaxis()->SetTitle(Var2);
  hist->Scale(Amount);
  hist->Draw("colz");

}

#endif

