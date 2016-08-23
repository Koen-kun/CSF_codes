#ifndef FUNCTIONBASE_H
#define FUNCTIONBASE_H
#include <iostream>
#include "TFormula.h"
#include "TH1.h"
	
class  FunctionBase {

  public:
    //Formula id, Formula name, Formula expression and Formula prefactor
    TFormula* f;
    TString   formName;
    TString   formula;
    TString   type; //"Phi", "Theta", "DeltaP", "SumP" and "DeltaPxSumP" //In case you want to classify it
    Double_t  preFactor;
    Double_t  limit;

    //Coefficient and error
    Double_t sumC; 
    Double_t sumsqC;
    Double_t meanC; 
    Double_t sigC;

    //Amount of events
    Double_t events;

    //Histogram
    TH1D* hC;
   
    //Bool for normalization check
    bool isNormalized;

    //Constructors and destructors
    FunctionBase();
    FunctionBase(TString,TString,Double_t,TString="Constant");
    FunctionBase(TString,TString,Double_t,Double_t,TString="Constant");
    ~FunctionBase() {delete f;}

    //Member functions
    virtual void Project_and_Fill(Double_t [],Double_t [], Double_t=1.0) {}
    virtual void Normalize() {}
    virtual void Print() {}
 
};

FunctionBase::FunctionBase(){
  f = new TFormula("f1","1");
  formName  = "dummy";
  formula   = "1";
  type      = "Constant";
  preFactor = 1.;
  limit     = 1.;
  sumC      = 0.;
  sumsqC    = 0.;
  meanC     = 0.;
  sigC      = 0.;
  events    = 0.;
}

FunctionBase::FunctionBase(TString name, TString form, Double_t pref, TString _type) : sumC(0.), sumsqC(0.), meanC(0.), sigC(0.), events(0.), isNormalized(false){
  f = new TFormula(name,form);
  formName  = name;
  formula   = form;
  type      = _type;
  preFactor = pref;
  limit     = pref;//Assuming that the prefactor sets the max and min value for the coefficient
}

FunctionBase::FunctionBase(TString name, TString form, Double_t pref, Double_t lim, TString _type) : sumC(0.), sumsqC(0.), meanC(0.), sigC(0.), events(0.), isNormalized(false){
  f = new TFormula(name,form);
  formName  = name;
  formula   = form;
  type      = _type;
  preFactor = pref;
  limit     = lim;
}

#endif


