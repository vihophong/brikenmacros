//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar  5 16:32:13 2018 by ROOT version 6.08/00
// from TTree treeimpSn136/treeimpSn136
// found on file: /home/phong/rootfiles/decaysep_0.5cps/final_decay_brips3106.root
//////////////////////////////////////////////////////////

#ifndef checkimprate_h
#define checkimprate_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TObject.h"

class checkimprate {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
 //IonBetaMult     *implant;
   UInt_t          fUniqueID;
   UInt_t          fBits;
 //vector<BELENHit*> neuf;
 //vector<BELENHit*> neub;
 //vector<CloverHit*> clover;
 //vector<BELENHit*> anc;
   UShort_t        nneuf;
   UShort_t        nneub;
   UShort_t        nclover;
   UShort_t        nanc;
   UShort_t        nbeam;
   ULong64_t       fts;
   UInt_t          fevt;
   UChar_t         fid;
   UShort_t        fmult;
   Double_t        fx;
   Double_t        fy;
   Short_t         fz;
   Short_t         fminx;
   Short_t         fmaxx;
   Short_t         fminy;
   Short_t         fmaxy;
   Double_t        fex;
   Double_t        fey;
   ULong64_t       fxts;
   ULong64_t       fyts;
   Int_t           ftdiffmin;
   Int_t           ftdiffmax;
   UShort_t        fnx;
   UShort_t        fny;
   UShort_t        fncx;
   UShort_t        fncy;
   UShort_t        fnz;
   UShort_t        fniz;
   UChar_t         frflag;
   UChar_t         fdrflag;
   UShort_t        fsumexyrank;
   Double_t        fmindx;
   Double_t        fmindy;
   Double_t        ftw;
   Double_t        fdtion;
   UShort_t        fdz;
   Double_t        fdeltax;
   Double_t        fdeltay;

   // List of branches
   TBranch        *b_implant_fUniqueID;   //!
   TBranch        *b_implant_fBits;   //!
   TBranch        *b_implant_nneuf;   //!
   TBranch        *b_implant_nneub;   //!
   TBranch        *b_implant_nclover;   //!
   TBranch        *b_implant_nanc;   //!
   TBranch        *b_implant_nbeam;   //!
   TBranch        *b_implant_fts;   //!
   TBranch        *b_implant_fevt;   //!
   TBranch        *b_implant_fid;   //!
   TBranch        *b_implant_fmult;   //!
   TBranch        *b_implant_fx;   //!
   TBranch        *b_implant_fy;   //!
   TBranch        *b_implant_fz;   //!
   TBranch        *b_implant_fminx;   //!
   TBranch        *b_implant_fmaxx;   //!
   TBranch        *b_implant_fminy;   //!
   TBranch        *b_implant_fmaxy;   //!
   TBranch        *b_implant_fex;   //!
   TBranch        *b_implant_fey;   //!
   TBranch        *b_implant_fxts;   //!
   TBranch        *b_implant_fyts;   //!
   TBranch        *b_implant_ftdiffmin;   //!
   TBranch        *b_implant_ftdiffmax;   //!
   TBranch        *b_implant_fnx;   //!
   TBranch        *b_implant_fny;   //!
   TBranch        *b_implant_fncx;   //!
   TBranch        *b_implant_fncy;   //!
   TBranch        *b_implant_fnz;   //!
   TBranch        *b_implant_fniz;   //!
   TBranch        *b_implant_frflag;   //!
   TBranch        *b_implant_fdrflag;   //!
   TBranch        *b_implant_fsumexyrank;   //!
   TBranch        *b_implant_fmindx;   //!
   TBranch        *b_implant_fmindy;   //!
   TBranch        *b_implant_ftw;   //!
   TBranch        *b_implant_fdtion;   //!
   TBranch        *b_implant_fdz;   //!
   TBranch        *b_implant_fdeltax;   //!
   TBranch        *b_implant_fdeltay;   //!

   checkimprate(char* infile,char* ri, TTree *tree=0);
   char identifier[1000];
   virtual ~checkimprate();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Int_t selectz);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef checkimprate_cxx
checkimprate::checkimprate(char* infile,char* ri, TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(infile);
      if (!f || !f->IsOpen()) {
         f = new TFile(infile);
      }
      char tempchar[500];
      sprintf(tempchar,"treeimp%s",ri);
      f->GetObject(tempchar,tree);

      sprintf(identifier,"%s\t%s",infile,ri);

   }
   Init(tree);
}

checkimprate::~checkimprate()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t checkimprate::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t checkimprate::LoadTree(Long64_t entry)
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

void checkimprate::Init(TTree *tree)
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

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_implant_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_implant_fBits);
   fChain->SetBranchAddress("nneuf", &nneuf, &b_implant_nneuf);
   fChain->SetBranchAddress("nneub", &nneub, &b_implant_nneub);
   fChain->SetBranchAddress("nclover", &nclover, &b_implant_nclover);
   fChain->SetBranchAddress("nanc", &nanc, &b_implant_nanc);
   fChain->SetBranchAddress("nbeam", &nbeam, &b_implant_nbeam);
   fChain->SetBranchAddress("fts", &fts, &b_implant_fts);
   fChain->SetBranchAddress("fevt", &fevt, &b_implant_fevt);
   fChain->SetBranchAddress("fid", &fid, &b_implant_fid);
   fChain->SetBranchAddress("fmult", &fmult, &b_implant_fmult);
   fChain->SetBranchAddress("fx", &fx, &b_implant_fx);
   fChain->SetBranchAddress("fy", &fy, &b_implant_fy);
   fChain->SetBranchAddress("fz", &fz, &b_implant_fz);
   fChain->SetBranchAddress("fminx", &fminx, &b_implant_fminx);
   fChain->SetBranchAddress("fmaxx", &fmaxx, &b_implant_fmaxx);
   fChain->SetBranchAddress("fminy", &fminy, &b_implant_fminy);
   fChain->SetBranchAddress("fmaxy", &fmaxy, &b_implant_fmaxy);
   fChain->SetBranchAddress("fex", &fex, &b_implant_fex);
   fChain->SetBranchAddress("fey", &fey, &b_implant_fey);
   fChain->SetBranchAddress("fxts", &fxts, &b_implant_fxts);
   fChain->SetBranchAddress("fyts", &fyts, &b_implant_fyts);
   fChain->SetBranchAddress("ftdiffmin", &ftdiffmin, &b_implant_ftdiffmin);
   fChain->SetBranchAddress("ftdiffmax", &ftdiffmax, &b_implant_ftdiffmax);
   fChain->SetBranchAddress("fnx", &fnx, &b_implant_fnx);
   fChain->SetBranchAddress("fny", &fny, &b_implant_fny);
   fChain->SetBranchAddress("fncx", &fncx, &b_implant_fncx);
   fChain->SetBranchAddress("fncy", &fncy, &b_implant_fncy);
   fChain->SetBranchAddress("fnz", &fnz, &b_implant_fnz);
   fChain->SetBranchAddress("fniz", &fniz, &b_implant_fniz);
   fChain->SetBranchAddress("frflag", &frflag, &b_implant_frflag);
   fChain->SetBranchAddress("fdrflag", &fdrflag, &b_implant_fdrflag);
   fChain->SetBranchAddress("fsumexyrank", &fsumexyrank, &b_implant_fsumexyrank);
   fChain->SetBranchAddress("fmindx", &fmindx, &b_implant_fmindx);
   fChain->SetBranchAddress("fmindy", &fmindy, &b_implant_fmindy);
   fChain->SetBranchAddress("ftw", &ftw, &b_implant_ftw);
   fChain->SetBranchAddress("fdtion", &fdtion, &b_implant_fdtion);
   fChain->SetBranchAddress("fdz", &fdz, &b_implant_fdz);
   fChain->SetBranchAddress("fdeltax", &fdeltax, &b_implant_fdeltax);
   fChain->SetBranchAddress("fdeltay", &fdeltay, &b_implant_fdeltay);
   Notify();
}

Bool_t checkimprate::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void checkimprate::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t checkimprate::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef checkimprate_cxx
