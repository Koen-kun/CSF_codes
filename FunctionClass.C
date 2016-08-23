#ifndef FUNCTIONCLASS_H
#define FUNCTIONCLASS_H

#include <iostream>
#include "TFormula.h"
#include "FunctionClassBase.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TProfile2D.h"

using namespace std;

class  Function : public FunctionBase {

  public:
    Function(TString name, TString form, Double_t pref, TString type = "Constant");
    Function(TString name, TString form, Double_t pref, Double_t lim, TString type = "Constant");
    ~Function();

    void Project_and_Fill(Double_t [],Double_t [], Double_t=1.0);
    void Normalize();
    void Print();

    Double_t calcIntegral(Function*,int=10,int=10,int=10,bool=false);
    Double_t calcMeanError(TH1D*);

    bool OrthoCheck(Function*,bool=false);
    bool PrefacCheck(bool=false);

    void DrawFormula1D(TString Var,int=10,int=10,int=10,double=1.0);
    void DrawFormula2D(TString Var1, TString Var2,int=10,int=10,int=10,double=1.0);
    void SetPositionXY(double x, double y, double &posx, double &posy);
};

Function::Function(TString name, TString form, Double_t pref, TString _type) : FunctionBase(name,form,pref,_type){
  hC = new TH1D("hC_"+formName,formName,100,-limit,limit);
  hC->Sumw2();
}

Function::Function(TString name, TString form, Double_t pref, Double_t lim, TString _type) : FunctionBase(name,form,pref,lim,_type) {
  hC = new TH1D("hC_"+formName,formName,100,-limit,limit);
  hC->Sumw2();
}

Function::~Function(){
if(!isNormalized) delete hC;
}

void Function::Project_and_Fill(Double_t angv[], Double_t momv[], Double_t weight){
  Double_t app(0.);
  
  //Projecting out term
  app = f->EvalPar(angv,momv) * preFactor;

  //Summing
  events += weight;
  sumC   += app*weight;
  sumsqC += app*app*weight*weight; 

  //Filling histogram
  hC->Fill(app,weight);
}

void Function::Normalize(){
  if(events == 0 || isNormalized) return;

  Double_t N = events;

  if(formName == "fR") meanC = sumC;
  else meanC = sumC/N;

  if(formName == "fR") sigC = sqrt(sumC);
  else sigC = calcMeanError(hC);

  delete hC;
  isNormalized = true;
}

void Function::Print(){

   if(formName == "fR"){
     printf("\n====================================================================================================================================\n");
     printf("|                            Rate Information\n");
     printf("|                            R: % 6.6f +- %f",sumC,sqrt(sumC));
     printf("\n====================================================================================================================================\n");
   }
   else{
     printf("\n====================================================================================================================================\n");
     printf("|                            Formula: %s \n",formName.Data());
     printf("|                            Expresion: %2.1f * %s \n",preFactor,formula.Data());
     printf("|                            Legend: x = theta1, y = theta2, z = phi1, t = phi2, [0] = dp, [1] = sp\n");
     printf("|                            Coefficient: % 6.6f +- %f, sig: % 6.1f",meanC,sigC,meanC/sigC);
     printf("\n====================================================================================================================================\n");
   }
}

Double_t Function::calcIntegral(Function* F2, int Nt, int Nph, int Nm, bool notify){

   Double_t t1, t2, ph1, ph2, dp, sp;
   Double_t angv[4];
   Double_t momv[2];
   Double_t pi = TMath::Pi();

   Double_t func = 0.;
   Double_t fact = 0.;
   Double_t intg = 0.;
   Double_t intc = 0.;

   //These loops represent a uniform sampling for the angular and momentum distributions
   for(int t1i=0; t1i<Nt; t1i++){
    for(int t2i=0;  t2i<Nt; t2i++) {
     for(int ph1i=0; ph1i<Nph; ph1i++){
      for(int ph2i=0; ph2i<Nph; ph2i++){
       for(int  dpi=0;  dpi<Nm; dpi++) {
        for(int  spi=0;  spi<Nm; spi++) { 

          //Make sure to always pick the center of the "bin"
          t1  = t1i*(pi/Nt) + 0.5*(pi/Nt);
          t2  = t2i*(pi/Nt) + 0.5*(pi/Nt);
          ph1 = (2*ph1i*(pi/Nph) - pi) + 0.5*(2*pi/Nph);
          ph2 = (2*ph2i*(pi/Nph) - pi) + 0.5*(2*pi/Nph);
          dp  = (2*dpi*(pi/Nm) - pi) + 0.5*(2*pi/Nm);
          sp  = (2*spi*(pi/Nm) - pi) + 0.5*(2*pi/Nm);

          angv[0] = t1;
          angv[1] = t2;
          angv[2] = ph1;
          angv[3] = ph2;
          momv[0] = dp; 
          momv[1] = sp; 

          func  = 1.0;
          func += F2->f->EvalPar(angv,momv);
          fact  = func * this->f->EvalPar(angv,momv) * this->preFactor;
          intg += fact;
          intc += func;

        }
       }
      }
     }
    }
   }

   intg = intg/intc;

   //Usefull for Debugging
   if(notify){
     printf("\n====================================================================================================================================\n");
     printf("|                            Orhtogonality check for the Function1 with Function2:\n");
     printf("|                            Function1: %2.1f * %s \n",this->preFactor,this->formula.Data());
     printf("|                            Function2: %s \n",F2->formula.Data());
     printf("|                            Integral: % 6.6f",intg);
     printf("\n====================================================================================================================================\n");
   }

   return intg;

}

bool Function::OrthoCheck(Function* F2, bool notify){
  bool passed = false;
  Double_t check = 0.;
  check = calcIntegral(F2,10,10,10,notify);
  if(check < 0.1) passed = true;

  if(notify){
    if(!passed) cout << "Functions are NOT orthogonal!\n" << endl;
    if(passed) cout  << "Functions are orthogonal!\n" << endl;
  }

  return passed;
}

bool Function::PrefacCheck(bool notify){
  bool passed = false;
  Double_t check = 0.;
  check = calcIntegral(this,10,10,10,notify);
  if(check > 0.9 && check < 1.1) passed = true;

  if(notify){
    if(!passed) cout << "Prefactor should be multiplied with: " << 1.0/check << "\n" << endl;
    if(passed) cout  << "Prefactor is OK!\n" << endl;
  }

  return passed;
}

Double_t Function::calcMeanError(TH1D* hist){

  Double_t sig(0.);
  Double_t mean(0.);
  TH1D* hpull = (TH1D*) hist->Clone();
  hpull->Reset();
  TH1D* hmeans =  new TH1D("hmeans","hmeans",1000,-4,4);

  double binCon(0.);
  double binErr(0.);
  double newBinCon(0.);
  TRandom3 ran;
  int bins = hist->GetNbinsX();
  for(int k=0; k<10000; k++){
    for(int i=0; i<bins; i++){
      binCon    = 0.;
      binErr    = 0.;
      newBinCon = 0.;
      if(hist->GetBinContent(i+1) != 0){
        binCon = hist->GetBinContent(i+1);
        binErr = hist->GetBinError(i+1);

        newBinCon = ran.Gaus(binCon,binErr);
      }
      hpull->SetBinContent((i+1),newBinCon);
      hpull->SetBinError((i+1),binErr);
    }
    mean = hpull->GetMean();
    hmeans->Fill(mean);
    hpull->Reset();
  }

  sig = hmeans->GetRMS();

  delete hmeans; //You can write it away if commenting out
  delete hpull;

  return sig;
}

void Function::DrawFormula1D(TString Var, int Nt, int Np, int Nm, double Amount){
  //Var = t1/t2/ph1/ph2/dp/sp
  double t1, t2, ph1, ph2, dp, sp;   
  Double_t pos(0.);
  Double_t angv[4];
  Double_t momv[2];
  Double_t pi = TMath::Pi();
  Double_t BinFactor = 1.0;

  Double_t func = 0.;
 
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

         func = (this->f->EvalPar(angv,momv))*BinFactor;

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

void Function::SetPositionXY(double x, double y, double &posx, double &posy){
  posx = x;
  posy = y;
}

void Function::DrawFormula2D(TString Var1, TString Var2, int Nt, int Np, int Nm, double Amount){
  //Var = t1/t2/ph1/ph2/dp/sp
  double t1, t2, ph1, ph2, dp, sp;
  Double_t posx(0.);
  Double_t posy(0.);
  Double_t angv[4];
  Double_t momv[2];
  Double_t pi = TMath::Pi();
  Double_t BinFactor = 1.0; 

  Double_t func = 0.;
 
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

         func = (this->f->EvalPar(angv,momv))*BinFactor;
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

