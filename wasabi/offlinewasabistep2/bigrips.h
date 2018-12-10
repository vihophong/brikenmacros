//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec  3 18:19:50 2018 by ROOT version 6.08/00
// from TTree tree/tree
// found on file: bigrips/run2152.root
//////////////////////////////////////////////////////////

#ifndef bigrips_h
#define bigrips_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class bigrips {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
 //TreeData        *bigrips;
   ULong64_t       ts;
   ULong64_t       sts;
   Double_t        tof;
   Double_t        zet;
   Double_t        aoq;
   Double_t        f5x;
   Double_t        f11x;
   Double_t        f11y;
   Double_t        f11dt;
   Double_t        beta;

   // List of branches
   TBranch        *b_bigrips_ts;   //!
   TBranch        *b_bigrips_sts;   //!
   TBranch        *b_bigrips_tof;   //!
   TBranch        *b_bigrips_zet;   //!
   TBranch        *b_bigrips_aoq;   //!
   TBranch        *b_bigrips_f5x;   //!
   TBranch        *b_bigrips_f11x;   //!
   TBranch        *b_bigrips_f11y;   //!
   TBranch        *b_bigrips_f11dt;   //!
   TBranch        *b_bigrips_beta;   //!

   bigrips(char* infile, TTree *tree=0);
   virtual ~bigrips();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   void GetTsEntry(std::multimap<Long64_t,Long64_t> &mts);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

