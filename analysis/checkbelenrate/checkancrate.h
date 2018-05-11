//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 24 16:30:42 2018 by ROOT version 6.08/00
// from TTree anc/tree anc
// found on file: belen82.root
//////////////////////////////////////////////////////////

#ifndef checkancrate_h
#define checkancrate_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TObject.h"

class checkancrate {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
 //BELENHit        *anc;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   UShort_t        fhitsadded;
   TVector3        fpos;
   TVector3        frndpos;
   Short_t         fdaqid;
   Short_t         fid;
   Short_t         fring;
   Short_t         ftype;
   Int_t           fadc;
   Double_t        fen;
   ULong64_t       fts;
   Double_t        ff11time;
   Double_t        fdvetotime;
   Double_t        fvetotime;

   // List of branches
   TBranch        *b_anc_fUniqueID;   //!
   TBranch        *b_anc_fBits;   //!
   TBranch        *b_anc_fhitsadded;   //!
   TBranch        *b_anc_fpos;   //!
   TBranch        *b_anc_frndpos;   //!
   TBranch        *b_anc_fdaqid;   //!
   TBranch        *b_anc_fid;   //!
   TBranch        *b_anc_fring;   //!
   TBranch        *b_anc_ftype;   //!
   TBranch        *b_anc_fadc;   //!
   TBranch        *b_anc_fen;   //!
   TBranch        *b_anc_fts;   //!
   TBranch        *b_anc_ff11time;   //!
   TBranch        *b_anc_fdvetotime;   //!
   TBranch        *b_anc_fvetotime;   //!

   checkancrate(char* infile,TTree *tree=0);
   virtual ~checkancrate();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef checkancrate_cxx
checkancrate::checkancrate(char* infile,TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(infile);
      if (!f || !f->IsOpen()) {
         f = new TFile(infile);
      }
      f->GetObject("anc",tree);
   }
   Init(tree);
}

checkancrate::~checkancrate()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t checkancrate::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t checkancrate::LoadTree(Long64_t entry)
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

void checkancrate::Init(TTree *tree)
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

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_anc_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_anc_fBits);
   fChain->SetBranchAddress("fhitsadded", &fhitsadded, &b_anc_fhitsadded);
   fChain->SetBranchAddress("fpos", &fpos, &b_anc_fpos);
   fChain->SetBranchAddress("frndpos", &frndpos, &b_anc_frndpos);
   fChain->SetBranchAddress("fdaqid", &fdaqid, &b_anc_fdaqid);
   fChain->SetBranchAddress("fid", &fid, &b_anc_fid);
   fChain->SetBranchAddress("fring", &fring, &b_anc_fring);
   fChain->SetBranchAddress("ftype", &ftype, &b_anc_ftype);
   fChain->SetBranchAddress("fadc", &fadc, &b_anc_fadc);
   fChain->SetBranchAddress("fen", &fen, &b_anc_fen);
   fChain->SetBranchAddress("fts", &fts, &b_anc_fts);
   fChain->SetBranchAddress("ff11time", &ff11time, &b_anc_ff11time);
   fChain->SetBranchAddress("fdvetotime", &fdvetotime, &b_anc_fdvetotime);
   fChain->SetBranchAddress("fvetotime", &fvetotime, &b_anc_fvetotime);
   Notify();
}

Bool_t checkancrate::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void checkancrate::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t checkancrate::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef checkancrate_cxx
