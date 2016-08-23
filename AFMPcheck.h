////////////////////////////////////////////////////////////////////////////////////////////////
//////                                                                                      ////
//////                   A.F.M.P - AWESOME FOURIER MODE PROJECTOR                           ////
//////                                 Header file                                          ////
////// ------------------------------------------------------------------------------------ ////         
//////                                                                                      ////
//////                              Author: Koen Oussoren                                   ////
//////                         Group: ATLAS experiment / NIKHEF                             ////
//////     Mail: koeno@nikhef.nl / koen.pieter.oussoren@cern.ch / kp.oussoren@gmail.com     ////
//////                            Date: October 30th, 2014                                  ////
//////                                                                                      //// 
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef AFMP_h
#define AFMP_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TString.h"

TString prefix = "Modes"; //Default root output file prefix

class AFMP {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        weight;
   //Double_t        metrew_0j;
   //Double_t        metrew_1j;
   //Double_t        metrew_0n1j;
   Int_t           passed;

   //Kinematics
   //leptons
   Double_t        l1Pt;
   Double_t        l2Pt;
   Double_t        l1Phi;
   Double_t        l2Phi;
   Double_t        l1Eta;
   Double_t        l2Eta;
   Double_t        l1E;
   Double_t        l2E;

   //Higgs
   Double_t        HPt;
   Double_t        HM;

   //jets
   Double_t        nJets;
   Double_t        j1Pt;
   Double_t        j2Pt;
   Double_t        j1Phi;
   Double_t        j2Phi;
   Double_t        j1Eta;
   Double_t        j2Eta;
   Double_t        j1E;
   Double_t        j2E;

   //CS frame variables
   //Double_t        rthetal1;
   //Double_t        rthetal2;
   //Double_t        rphil1;
   //Double_t        rphil2;
   //Double_t        rl1P;
   //Double_t        rl2P;

   Double_t        thetal1;
   Double_t        thetal2;
   Double_t        phil1;
   Double_t        phil2;
   Double_t        csl1P;
   Double_t        csl2P;

   //Double_t        rthetaj1;
   //Double_t        rthetaj2;
   //Double_t        rphij1;
   //Double_t        rphij2;
   //Double_t        rj1P;
   //Double_t        rj2P;

   //Other
   Double_t        Mll;
   Double_t        DPhill;
   Double_t        sMET;
   Double_t        Ptll;
   Double_t        MT;
   Double_t        METtr;
   Double_t        MWT1;
   Double_t        MWT2;

   // List of branches
   TBranch        *b_weight;   //!
   //TBranch        *b_metrew_0j;   //!
   //TBranch        *b_metrew_1j;   //!
   //TBranch        *b_metrew_0n1j;   //!
   TBranch        *b_passed;   //!

   TBranch        *b_l1Pt;   //!
   TBranch        *b_l2Pt;   //!
   TBranch        *b_l1Phi;   //!
   TBranch        *b_l2Phi;   //!
   TBranch        *b_l1Eta;   //!
   TBranch        *b_l2Eta;   //!
   TBranch        *b_l1E;   //!
   TBranch        *b_l2E;   //!

   TBranch        *b_HPt;   //!
   TBranch        *b_HM;   //!

   TBranch        *b_nJets;   //!
   TBranch        *b_j1Pt;   //!
   TBranch        *b_j2Pt;   //!
   TBranch        *b_j1Phi;   //!
   TBranch        *b_j2Phi;   //!
   TBranch        *b_j1Eta;   //!
   TBranch        *b_j2Eta;   //!
   TBranch        *b_j1E;   //!
   TBranch        *b_j2E;   //!

   //TBranch        *b_rthetal1;   //!
   //TBranch        *b_rthetal2;   //!
   //TBranch        *b_rphil1;   //!
   //TBranch        *b_rphil2;   //!
   //TBranch        *b_rl1P;   //!
   //TBranch        *b_rl2P;   //!

   TBranch        *b_thetal1;   //!
   TBranch        *b_thetal2;   //!
   TBranch        *b_phil1;   //!
   TBranch        *b_phil2;   //!
   TBranch        *b_csl1P;   //!
   TBranch        *b_csl2P;   //!

   //TBranch        *b_rthetaj1;   //!
   //TBranch        *b_rthetaj2;   //!
   //TBranch        *b_rphij1;   //!
   //TBranch        *b_rphij2;   //!
   //TBranch        *b_rj1P;   //!
   //TBranch        *b_rj2P;   //!

   TBranch        *b_Mll;   //!
   TBranch        *b_DPhill;   //!
   TBranch        *b_sMET;   //!
   TBranch        *b_Ptll;   //!
   TBranch        *b_MT;   //!
   TBranch        *b_METtr;   //!
   TBranch        *b_MWT1;   //!
   TBranch        *b_MWT2;   //!

   AFMP();
   //AFMP(TTree *tree=0);
   //AFMP(TString s);
   AFMP(TTree *tree);
   virtual ~AFMP();
   //virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     ModesLoop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AFMP_cxx
AFMP::AFMP(){

  TChain* chain = new TChain("Tree");

  //chain->Add("/data/atlas/users/koeno/Spin_studies/shower/output_trees/h0pww_aMCatNLO_corrected_200k_tree.root");
  chain->Add("/data/atlas/users/koeno/Spin_studies/shower/output_trees/h0mww_aMCatNLO_corrected_200k_tree.root");
  std::cout << "NOTE: using default tree: ...../h0pww_aMCatNLO_corrected_200k_tree.root" << std::endl;

  Init(chain);
}

AFMP::AFMP(TTree *tree) : fChain(0) 
{
   if (tree == 0) {
        std::cout << "ERROR: NO tree selected or tree is NULL pointer" << std::endl;
        return;
   }
   std::cout << "Initializing TTree" << std::endl;
   Init(tree);
}

AFMP::~AFMP()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AFMP::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AFMP::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AFMP::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("weight", &weight, &b_weight);
   //fChain->SetBranchAddress("metrew_0j", &metrew_0j, &b_metrew_0j);
   //fChain->SetBranchAddress("metrew_1j", &metrew_1j, &b_metrew_1j);
   //fChain->SetBranchAddress("metrew_0n1j", &metrew_0n1j, &b_metrew_0n1j);
   fChain->SetBranchAddress("passed", &passed, &b_passed);

   fChain->SetBranchAddress("l1Pt", &l1Pt, &b_l1Pt);
   fChain->SetBranchAddress("l2Pt", &l2Pt, &b_l2Pt);
   fChain->SetBranchAddress("l1Phi", &l1Phi, &b_l1Phi);
   fChain->SetBranchAddress("l2Phi", &l2Phi, &b_l2Phi);
   fChain->SetBranchAddress("l1Eta", &l1Eta, &b_l1Eta);
   fChain->SetBranchAddress("l2Eta", &l2Eta, &b_l2Eta);
   fChain->SetBranchAddress("l1E", &l1E, &b_l1E);
   fChain->SetBranchAddress("l2E", &l2E, &b_l2E);

   fChain->SetBranchAddress("HPt", &HPt, &b_HPt);
   fChain->SetBranchAddress("HM", &HM, &b_HM);

   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("j1Pt", &j1Pt, &b_j1Pt);
   fChain->SetBranchAddress("j2Pt", &j2Pt, &b_j2Pt);
   fChain->SetBranchAddress("j1Phi", &j1Phi, &b_j1Phi);
   fChain->SetBranchAddress("j2Phi", &j2Phi, &b_j2Phi);
   fChain->SetBranchAddress("j1Eta", &j1Eta, &b_j1Eta);
   fChain->SetBranchAddress("j2Eta", &j2Eta, &b_j2Eta);
   fChain->SetBranchAddress("j1E", &j1E, &b_j1E);
   fChain->SetBranchAddress("j2E", &j2E, &b_j2E);

   //fChain->SetBranchAddress("rthetal1", &rthetal1, &b_rthetal1);
   //fChain->SetBranchAddress("rthetal2", &rthetal2, &b_rthetal2);
   //fChain->SetBranchAddress("rphil1", &rphil1, &b_rphil1);
   //fChain->SetBranchAddress("rphil2", &rphil2, &b_rphil2);
   //fChain->SetBranchAddress("rl1P", &rl1P, &b_rl1P);
   //fChain->SetBranchAddress("rl2P", &rl2P, &b_rl2P);

   fChain->SetBranchAddress("thetal1", &thetal1, &b_thetal1);
   fChain->SetBranchAddress("thetal2", &thetal2, &b_thetal2);
   fChain->SetBranchAddress("phil1",   &phil1,   &b_phil1);
   fChain->SetBranchAddress("phil2",   &phil2,   &b_phil2);
   fChain->SetBranchAddress("csl1P",   &csl1P,   &b_csl1P);
   fChain->SetBranchAddress("csl2P",   &csl2P,   &b_csl2P);

   //fChain->SetBranchAddress("rthetaj1", &rthetaj1, &b_rthetaj1);
   //fChain->SetBranchAddress("rthetaj2", &rthetaj2, &b_rthetaj2);
   //fChain->SetBranchAddress("rphij1", &rphij1, &b_rphij1);
   //fChain->SetBranchAddress("rphij2", &rphij2, &b_rphij2);
   //fChain->SetBranchAddress("rj1P", &rj1P, &b_rj1P);
   //fChain->SetBranchAddress("rj2P", &rj2P, &b_rj2P);

   fChain->SetBranchAddress("Mll", &Mll, &b_Mll);
   fChain->SetBranchAddress("DPhill", &DPhill, &b_DPhill);
   fChain->SetBranchAddress("sMET", &sMET, &b_sMET);
   fChain->SetBranchAddress("Ptll", &Ptll, &b_Ptll);
   fChain->SetBranchAddress("MT", &MT, &b_MT);


   Notify();
}

Bool_t AFMP::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AFMP::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
// Int_t AFMP::Cut(Long64_t entry)
//{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
//   return 1;
//}
#endif // #ifdef AFMP_cxx
