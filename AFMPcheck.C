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

//  This class projects out predefined coefficients by looping over H->WW->lvlv event  
//  Created with the MadGraph or aMC@NLO generators. There are 6 sum and difference
//  variables to chose from based on the angular and momentum of the charged leptons
//
//  -  The FunctionClass contains all the functionality to project out certain Fourier
//  mode and do a proper normalization when looping over the event. FunctionClassBass
//  is a virtual base class. I used it in some older version to let other Function
//  classes inherited from it, but it in this version only FunctionClass inherits 
//  from it. 
//
//  - The FunctionCollector class is a convenient container class to store multiple
//  Function objects. It is also able to store the original sum and difference
//  variables as data members so they can later on be used to compare the a
//  theorectical function with the original data. In the code below it stores the
//  functions from Sig_terms.C and plots a 1D and 2D data/theory comparison at the end 
//
//  - f(x) = a0 + a1*cos(x) + a2*cos(2*x) + a3*cos(3*x) + ... (Assuming a symmetric function)
//  The code projects out the coefficients a0, a1, a2, ..., by means of integration
//  After the integration step the coefficients have to be normalized and you'll end
//  up with something like:
//  f(x) = 1 + 3*cos(x) + .5*cos(3x), which can be compared to the original distribution
//  of x to see if all went alright

//Some basic headers from the ROOT libraries
#define AFMP_cxx
#include "AFMPcheck.h"
#include "TMath.h"
#include <TFormula.h>
#include <sstream>
#include <fstream>
#include <iostream>

//The files for the Function classes
#include "FunctionClassBase.h"
#include "FunctionClass.C"
#include "FunctionCollector.C"

//File containing a collection of Functions that can be used
#include "Sig_terms.C"

//Some regular stuff
using namespace std;
double pi = TMath::Pi();

//This is the main body of the AFMP class
void AFMP::ModesLoop()
{
   if (fChain == 0) return;

   //Defining basic construction functions
   //Convention: x = DeltaTheta, y = SumTheta, z = DeltaPhi, 
   //            t = SumPhi, [0] = DeltaP, [1] = SumP

   //This is an example function that works with DeltaPhi (z) and we will project out
   //the relevant coefficient in the data sample for the "cos(z)" term (argument 2). The Function
   //class needs a couple of things to be initialized: the term to be projected out
   //(of course), the proper normalization factor called the prefactor (argument 3) and some other stuff
   //which is not relevant for now. The right prefactor is very important, because
   //it takes care of an additinal factor when integrating out the coefficient while looping
   //over the events
   Function *f1 = new Function("f1","cos(1*z)",2.,2.,"Phi");

   //We also introduce 3 other functions, for which we will project out the coefficients
   //However, we don't want to introduce coefficients that are correlated with the ones we
   //had so we have to make sure that the new functions are orthogonal with respect to f1.
   //The 3rd new funtion looks a bit complicated so we might want to check if we have the
   //right prefactor for this one.
   Function *test1 = new Function("test1","4*cos(1*z)",2.,2.,"Phi"); //This one obviously correlated with f1
   Function *test2 = new Function("test2","4*cos(2*z)",2.,2.,"Phi"); //This one should be orthogonal w.r.t. f1
   Function *test3 = new Function("test3","cos(2*z) + cos(3*z)",2.,2.,"Phi"); //Do we have the right prefactor here?

   //FunctionCollector container
   FunctionCollector Fvec4;

   //We declare a vector object that contains Function objects
   vector<FunctionBase*> functions;
   vector<TString> coefs;
   TString histname = prefix; //"Modes"
   TString title;

   //The vector functions is being filled with Function objects defined in
   //Sig_terms.C (6 in total)
   do_Terms(functions);

   //Some initialization for a vector containing the names of the Function objects
   TString coef;
   for(vector<FunctionBase*>::iterator it = functions.begin(); it != functions.end(); ++it){
     coef = (*it)->formName;
     coefs.push_back(coef);
   }

   //Check output
   cout << "\nCheck functions size:          " << functions.size() << endl;
   cout << "Check coefs size:              " << coefs.size() << endl;

   //Load the entries of the Ntuple
   Long64_t nentries = fChain->GetEntriesFast();

   //We need to define some variables, before we can calculate the sum and difference
   Double_t p1, p2, t1, t2, ph1, ph2, weightfac;
   Double_t dt, st;
   Double_t dph, sph;
   Double_t dp, sp, spmean;

   //These two arrays are necessary for ROOT's TFormula magic 
   Double_t angv[4];
   Double_t momv[2];

   //The actual event loop starts here
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<10000;jentry++) {

      //Loading in (next) event
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //Setting all relevant variables to correspond with this event
      weightfac = 1;
      if(weight<0) weightfac = -1;

      p1  = csl1P;
      p2  = csl2P;
      t1  = thetal1;
      t2  = thetal2;
      ph1 = phil1;
      ph2 = phil2;
      spmean = 62.5;

      //Calculate sum and differences
      dp  = (p1-p2)*pi/(62.5);
      sp  = (p1+p2-spmean)*pi/(42); 
      dt  = t1 - t2;
      st  = t1 + t2 -pi;
      dph = ph1 - ph2;
      sph = ph1 + ph2;

      //Remap to range between -pi to pi
      if(dph < -pi) dph += 2*pi;
      else if(dph > pi) dph -= 2*pi;
      if(sph < -pi) sph += 2*pi;
      else if(sph > pi) sph -= 2*pi;

      //Put the sum and difference variables in the two arrays
      angv[0] = dt;
      angv[1] = st;
      angv[2] = dph;
      angv[3] = sph;
      momv[0] = dp;
      momv[1] = sp;

      //Here we do the actual projection of a coefficient for a certain term (function) per event
      //We have our functions f1, test1, test2 and test3
      f1->Project_and_Fill(angv,momv,weightfac);
      test1->Project_and_Fill(angv,momv,weightfac);
      test2->Project_and_Fill(angv,momv,weightfac);
      test3->Project_and_Fill(angv,momv,weightfac);

      //We also do the projection for all the Function objects inside the Function vector
      for (vector<FunctionBase*>::iterator it = functions.begin(); it != functions.end(); ++it){
        (*it)->Project_and_Fill(angv,momv,weightfac);
      }

      //Fill the Function Collector class with the original data
      Fvec4.FillVarVectors(dph,"dph");
      Fvec4.FillVarVectors(sph,"sph");
      Fvec4.FillVarVectors(weightfac,"w");

      //Stuff which is not used for now
      //Fvec4.FillVarVectors(dt,"dt");
      //Fvec4.FillVarVectors(st,"st");
      //Fvec4.FillVarVectors(dp,"dp");
      //Fvec4.FillVarVectors(sp,"sp");

   } //End EventLoop

   //Normalizing the coefficients, that is w.r.t. a0 set to 1
   f1->Normalize();
   test1->Normalize();
   test2->Normalize();
   test3->Normalize();

   //Printing out the results of the projection
   cout << "\nPrinting out the results of the projection:" << endl;
   f1->Print();
   test1->Print();
   test2->Print();
   test3->Print();

   //Doing the normalization for the Functions in the vector and adding them to the Function Collector
   for(vector<FunctionBase*>::iterator it = functions.begin(); it != functions.end(); ++it){
     (*it)->Normalize();
     //(*it)->Print();
     Fvec4.Add(*it);
   }

  //We don't need the functions vector anymore, so we can clean it up and free memory
  functions.erase(functions.begin(),functions.end());

  //Checking if the function is really empty after cleanup
  cout << "Check functions size: " << functions.size() << endl;

  //Now we've actually done the projection we can check if f1 is orthogonal to
  //test1, test2 or test3
  cout << "\nDoing the orthogonality check:" << endl;
  f1->OrthoCheck(test1,true);
  f1->OrthoCheck(test2,true);
  f1->OrthoCheck(test3,true);

  //And also let's check if the prefactor for test3 was the right guess
  cout << "\nTesting the prefactor:" << endl;
  test3->PrefacCheck(true);

  //We can use the Function Collector to create a 1D or 2D (theoretical)
  //distribution based on the projected coeffients and compare it to the
  //original 1D or 2D distribution to see if we found the most significant
  //coefficients and if the projection went alright
  cout << "\nDrawing a 1D and 2D distribution:" << endl;
  Fvec4.DataCompare1D("dph",25);
  Fvec4.DataCompare2D("dph","sph",25);

  //Not needed for now
  //Fvec4.DataCompare2D("dp","sp",25);
  //Fvec4.DataCompare2D("dph","sph",25);
  //Fvec4.DataCompare2D("dt","st",25);

}

//The End of the Example



