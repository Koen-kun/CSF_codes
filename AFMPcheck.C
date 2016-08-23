//////////////////////////////////////////////////////////////////////////////////////////////
////             									  ////
////                   A.F.M.P - AWESOME FOURIER MODE PROJECTOR                           ////
//// 											  ////
//// ------------------------------------------------------------------------------------ ////         
////											  ////
////                              Author: Koen Oussoren					  ////
////                         Group: ATLAS experiment / NIKHEF				  ////
////     Mail: koeno@nikhef.nl / koen.pieter.oussoren@cern.ch / kp.oussoren@gmail.com     ////
////                            Date: February 3th, 2015                                  ////
////											  //// 
//////////////////////////////////////////////////////////////////////////////////////////////

#define AFMP_cxx
#include "AFMPcheck.h"
#include "TMath.h"
#include <TFormula.h>
#include <sstream>
#include <fstream>
#include <iostream>

#include "FunctionClassBase.h"
#include "FunctionClass.C"
//#include "FunctionStringCollector.C"
//#include "FunctionTypeCollector.C"
#include "FunctionCollector.C"

//#include "Sig_terms.C"
#include "Sig_terms_check.C"

using namespace std;

double pi = TMath::Pi();

void AFMP::ModesLoop()
{
   if (fChain == 0) return;

   //Defining basic construction functions
   //Convention: x = dt, y = st, z = dph, t = sph, [0] = dp, [1] = sp

   //FunctionCollector container
   FunctionCollector Fvec4;

   vector<FunctionBase*> functions;
   vector<TString> coefs;
   TString histname = prefix; //"Modes"
   TString title;

   do_Terms(functions);

   TString coef;
   for(vector<FunctionBase*>::iterator it = functions.begin(); it != functions.end(); ++it){
     coef = (*it)->formName;
     coefs.push_back(coef);
   }

   cout << "\nCheck functions size:          " << functions.size() << endl;
   cout << "Check coefs size:              " << coefs.size() << endl;

   Long64_t nentries = fChain->GetEntriesFast();

   Double_t p1, p2, t1, t2, ph1, ph2, weightfac;
   Double_t dt, st;
   Double_t dph, sph;
   Double_t dp, sp, spmean;

   Double_t angv[4];
   Double_t momv[2];

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<10000;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      weightfac = 1;
      if(weight<0) weightfac = -1;

      p1  = csl1P;
      p2  = csl2P;
      t1  = thetal1;
      t2  = thetal2;
      ph1 = phil1;
      ph2 = phil2;
      spmean = 62.5;

      dp  = (p1-p2)*pi/(62.5);
      sp  = (p1+p2-spmean)*pi/(42); 
      dt  = t1 - t2;
      st  = t1 + t2 -pi;
      dph = ph1 - ph2;
      sph = ph1 + ph2;

      if(dph < -pi) dph += 2*pi;
      else if(dph > pi) dph -= 2*pi;
      if(sph < -pi) sph += 2*pi;
      else if(sph > pi) sph -= 2*pi;

      angv[0] = dt;
      angv[1] = st;
      angv[2] = dph;
      angv[3] = sph;
      momv[0] = dp;
      momv[1] = sp;

      for (vector<FunctionBase*>::iterator it = functions.begin(); it != functions.end(); ++it){
        (*it)->Project_and_Fill(angv,momv,weightfac);
      }

      Fvec4.FillVarVectors(dt,"dt");
      Fvec4.FillVarVectors(st,"st");
      Fvec4.FillVarVectors(dph,"dph");
      Fvec4.FillVarVectors(sph,"sph");
      Fvec4.FillVarVectors(dp,"dp");
      Fvec4.FillVarVectors(sp,"sp");
      Fvec4.FillVarVectors(weightfac,"w");

   } //End EventLoop

   for(vector<FunctionBase*>::iterator it = functions.begin(); it != functions.end(); ++it){

     (*it)->Normalize();
     (*it)->Print();

     Fvec4.Add(*it);

   }

  functions.erase(functions.begin(),functions.end());

  cout << "Check functions size: " << functions.size() << endl;

  Fvec4.DataCompare1D("sp",25);
  Fvec4.DataCompare2D("dp","sp",25);
  Fvec4.DataCompare2D("dph","sph",25);
  Fvec4.DataCompare2D("dt","st",25);

}


