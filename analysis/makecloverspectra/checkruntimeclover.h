//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Dec 23 14:38:56 2018 by ROOT version 6.08/00
// from TTree gamma/tree gamma
// found on file: rootfiles/eucal_prep_run1314.root
//////////////////////////////////////////////////////////

#ifndef checkruntimeclover_h
#define checkruntimeclover_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TVectorD.h>

// Header file for the classes stored in the TTree if any.
#include "TObject.h"

class checkruntimeclover {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
 //CloverHit       *gamma;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   UShort_t        fhitsadded;
   TVector3        fpos;
   Short_t         fid;
   Short_t         fdaqid;
   Short_t         fclover;
   Short_t         fcloverleaf;
   Int_t           fadc;
   Double_t        fen;
   ULong64_t       fts;

   // List of branches
   TBranch        *b_gamma_fUniqueID;   //!
   TBranch        *b_gamma_fBits;   //!
   TBranch        *b_gamma_fhitsadded;   //!
   TBranch        *b_gamma_fpos;   //!
   TBranch        *b_gamma_fid;   //!
   TBranch        *b_gamma_fdaqid;   //!
   TBranch        *b_gamma_fclover;   //!
   TBranch        *b_gamma_fcloverleaf;   //!
   TBranch        *b_gamma_fadc;   //!
   TBranch        *b_gamma_fen;   //!
   TBranch        *b_gamma_fts;   //!

   checkruntimeclover(char* infile,TTree *tree=0);
   virtual ~checkruntimeclover();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef checkruntimeclover_cxx
checkruntimeclover::checkruntimeclover(char* infile,TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(infile);
      if (!f || !f->IsOpen()) {
         f = new TFile(infile);
      }
      f->GetObject("gamma",tree);

   }
   Init(tree);
}

checkruntimeclover::~checkruntimeclover()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t checkruntimeclover::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t checkruntimeclover::LoadTree(Long64_t entry)
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

void checkruntimeclover::Init(TTree *tree)
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

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_gamma_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_gamma_fBits);
   fChain->SetBranchAddress("fhitsadded", &fhitsadded, &b_gamma_fhitsadded);
   fChain->SetBranchAddress("fpos", &fpos, &b_gamma_fpos);
   fChain->SetBranchAddress("fid", &fid, &b_gamma_fid);
   fChain->SetBranchAddress("fdaqid", &fdaqid, &b_gamma_fdaqid);
   fChain->SetBranchAddress("fclover", &fclover, &b_gamma_fclover);
   fChain->SetBranchAddress("fcloverleaf", &fcloverleaf, &b_gamma_fcloverleaf);
   fChain->SetBranchAddress("fadc", &fadc, &b_gamma_fadc);
   fChain->SetBranchAddress("fen", &fen, &b_gamma_fen);
   fChain->SetBranchAddress("fts", &fts, &b_gamma_fts);
   Notify();
}

Bool_t checkruntimeclover::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void checkruntimeclover::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t checkruntimeclover::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef checkruntimeclover_cxx
