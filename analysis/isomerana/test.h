//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr  3 14:37:19 2019 by ROOT version 6.08/00
// from TTree treeIn134/treeIn134
// found on file: /home/phong/briken17/makerootfiles/rootfiles/decay_isomer_wab_slewcorr_gc/decay_brips3107.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class test {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       isomer_evt;
   ULong64_t       isomer_ts;
   Double_t        isomer_ion_t;
   Double_t        isomer_ion_x;
   Double_t        isomer_ion_y;
   Double_t        isomer_ion_ex;
   Double_t        isomer_ion_ey;
   Double_t        isomer_zet;
   Double_t        isomer_aoq;
   Double_t        isomer_beta;
   Double_t        isomer_F11L_T;
   Double_t        isomer_F11L_E;
   Double_t        isomer_F11R_T;
   Double_t        isomer_F11R_E;
   Double_t        isomer_F7_T;
   Double_t        isomer_F7_E;
   Double_t        isomer_veto_T;
   Double_t        isomer_veto_E;
   Double_t        isomer_de_T;
   Double_t        isomer_de_E;
   Short_t         isomer_ion_z;
   Short_t         isomer_multx;
   Short_t         isomer_multy;
   Short_t         isomer_multz;
   Int_t           gc_hit;
   Double_t        gc_E[8];   //[gc_hit]
   Double_t        gc_T[8];   //[gc_hit]
   Double_t        gc_Tslew[8];   //[gc_hit]
   Int_t           gc_ch[8];   //[gc_hit]
   Int_t           gc1_hit;
   Double_t        gc1_E[5];   //[gc1_hit]
   Double_t        gc1_T[5];   //[gc1_hit]
   Double_t        gc1_Tslew[5];   //[gc1_hit]
   Int_t           gc1_ch[5];   //[gc1_hit]
   Int_t           gc2_hit;
   Double_t        gc2_E[5];   //[gc2_hit]
   Double_t        gc2_T[5];   //[gc2_hit]
   Double_t        gc2_Tslew[5];   //[gc2_hit]
   Int_t           gc2_ch[5];   //[gc2_hit]
   Int_t           ab1_hit;
   Double_t        ab1_E[2];   //[ab1_hit]
   Double_t        ab1_T[2];   //[ab1_hit]
   Double_t        ab1_Tslew[2];   //[ab1_hit]
   Int_t           ab1_ch[2];   //[ab1_hit]
   Short_t         ab1_mult[2];   //[ab1_hit]
   Int_t           ab2_hit;
   Double_t        ab2_E[3];   //[ab2_hit]
   Double_t        ab2_T[3];   //[ab2_hit]
   Double_t        ab2_Tslew[3];   //[ab2_hit]
   Int_t           ab2_ch[3];   //[ab2_hit]
   Short_t         ab2_mult[3];   //[ab2_hit]
   Int_t           neu_hit;
   Double_t        neu_E[1];   //[neu_hit]
   Double_t        neu_T[1];   //[neu_hit]
   Int_t           neu_ch[1];   //[neu_hit]

   // List of branches
   TBranch        *b_isomer;   //!
   TBranch        *b_gc_hit;   //!
   TBranch        *b_gc_E;   //!
   TBranch        *b_gc_T;   //!
   TBranch        *b_gc_Tslew;   //!
   TBranch        *b_gc_ch;   //!
   TBranch        *b_gc1_hit;   //!
   TBranch        *b_gc1_E;   //!
   TBranch        *b_gc1_T;   //!
   TBranch        *b_gc1_Tslew;   //!
   TBranch        *b_gc1_ch;   //!
   TBranch        *b_gc2_hit;   //!
   TBranch        *b_gc2_E;   //!
   TBranch        *b_gc2_T;   //!
   TBranch        *b_gc2_Tslew;   //!
   TBranch        *b_gc2_ch;   //!
   TBranch        *b_ab1_hit;   //!
   TBranch        *b_ab1_E;   //!
   TBranch        *b_ab1_T;   //!
   TBranch        *b_ab1_Tslew;   //!
   TBranch        *b_ab1_ch;   //!
   TBranch        *b_ab1_mult;   //!
   TBranch        *b_ab2_hit;   //!
   TBranch        *b_ab2_E;   //!
   TBranch        *b_ab2_T;   //!
   TBranch        *b_ab2_Tslew;   //!
   TBranch        *b_ab2_ch;   //!
   TBranch        *b_ab2_mult;   //!
   TBranch        *b_neu_hit;   //!
   TBranch        *b_neu_E;   //!
   TBranch        *b_neu_T;   //!
   TBranch        *b_neu_ch;   //!

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/phong/briken17/makerootfiles/rootfiles/decay_isomer_wab_slewcorr_gc/decay_brips3107.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/home/phong/briken17/makerootfiles/rootfiles/decay_isomer_wab_slewcorr_gc/decay_brips3107.root");
      }
      f->GetObject("treeIn134",tree);

   }
   Init(tree);
}

test::~test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
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

void test::Init(TTree *tree)
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

   fChain->SetBranchAddress("isomer", &isomer_evt, &b_isomer);
   fChain->SetBranchAddress("gc_hit", &gc_hit, &b_gc_hit);
   fChain->SetBranchAddress("gc_E", gc_E, &b_gc_E);
   fChain->SetBranchAddress("gc_T", gc_T, &b_gc_T);
   fChain->SetBranchAddress("gc_Tslew", gc_Tslew, &b_gc_Tslew);
   fChain->SetBranchAddress("gc_ch", gc_ch, &b_gc_ch);
   fChain->SetBranchAddress("gc1_hit", &gc1_hit, &b_gc1_hit);
   fChain->SetBranchAddress("gc1_E", gc1_E, &b_gc1_E);
   fChain->SetBranchAddress("gc1_T", gc1_T, &b_gc1_T);
   fChain->SetBranchAddress("gc1_Tslew", gc1_Tslew, &b_gc1_Tslew);
   fChain->SetBranchAddress("gc1_ch", gc1_ch, &b_gc1_ch);
   fChain->SetBranchAddress("gc2_hit", &gc2_hit, &b_gc2_hit);
   fChain->SetBranchAddress("gc2_E", gc2_E, &b_gc2_E);
   fChain->SetBranchAddress("gc2_T", gc2_T, &b_gc2_T);
   fChain->SetBranchAddress("gc2_Tslew", gc2_Tslew, &b_gc2_Tslew);
   fChain->SetBranchAddress("gc2_ch", gc2_ch, &b_gc2_ch);
   fChain->SetBranchAddress("ab1_hit", &ab1_hit, &b_ab1_hit);
   fChain->SetBranchAddress("ab1_E", ab1_E, &b_ab1_E);
   fChain->SetBranchAddress("ab1_T", ab1_T, &b_ab1_T);
   fChain->SetBranchAddress("ab1_Tslew", ab1_Tslew, &b_ab1_Tslew);
   fChain->SetBranchAddress("ab1_ch", ab1_ch, &b_ab1_ch);
   fChain->SetBranchAddress("ab1_mult", ab1_mult, &b_ab1_mult);
   fChain->SetBranchAddress("ab2_hit", &ab2_hit, &b_ab2_hit);
   fChain->SetBranchAddress("ab2_E", ab2_E, &b_ab2_E);
   fChain->SetBranchAddress("ab2_T", ab2_T, &b_ab2_T);
   fChain->SetBranchAddress("ab2_Tslew", ab2_Tslew, &b_ab2_Tslew);
   fChain->SetBranchAddress("ab2_ch", ab2_ch, &b_ab2_ch);
   fChain->SetBranchAddress("ab2_mult", ab2_mult, &b_ab2_mult);
   fChain->SetBranchAddress("neu_hit", &neu_hit, &b_neu_hit);
   fChain->SetBranchAddress("neu_E", &neu_E, &b_neu_E);
   fChain->SetBranchAddress("neu_T", &neu_T, &b_neu_T);
   fChain->SetBranchAddress("neu_ch", &neu_ch, &b_neu_ch);
   Notify();
}

Bool_t test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef test_cxx
