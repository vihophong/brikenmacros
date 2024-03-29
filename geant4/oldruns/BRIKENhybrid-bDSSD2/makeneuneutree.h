//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec 18 16:39:02 2017 by ROOT version 6.08/00
// from TTree BRIKEN/BRIKEN
// found on file: tempHist100k.root
//////////////////////////////////////////////////////////

#ifndef makeneuneutree_h
#define makeneuneutree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class makeneuneutree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        id;
   Double_t        detGroup;
   Double_t        tubeGroup;
   Double_t        x;
   Double_t        y;
   Double_t        z;
   Double_t        E;
   Double_t        T;
   Double_t        partEnergy;
   Double_t        id2;
   Double_t        detGroup2;
   Double_t        tubeGroup2;
   Double_t        x2;
   Double_t        y2;
   Double_t        z2;
   Double_t        E2;
   Double_t        T2;
   Double_t        partEnergy2;

   // List of branches
   TBranch        *b_id;   //!
   TBranch        *b_detGroup;   //!
   TBranch        *b_tubeGroup;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_E;   //!
   TBranch        *b_T;   //!
   TBranch        *b_partEnergy;   //!
   TBranch        *b_id2;   //!
   TBranch        *b_detGroup2;   //!
   TBranch        *b_tubeGroup2;   //!
   TBranch        *b_x2;   //!
   TBranch        *b_y2;   //!
   TBranch        *b_z2;   //!
   TBranch        *b_E2;   //!
   TBranch        *b_T2;   //!
   TBranch        *b_partEnergy2;   //!

   makeneuneutree(TTree *tree=0);
   virtual ~makeneuneutree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(char *outfile, Double_t neutronrate);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef makeneuneutree_cxx
makeneuneutree::makeneuneutree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("tempHist1000k.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("tempHist1000k.root");
      }
      f->GetObject("BRIKEN",tree);

   }
   Init(tree);
}

makeneuneutree::~makeneuneutree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t makeneuneutree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t makeneuneutree::LoadTree(Long64_t entry)
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

void makeneuneutree::Init(TTree *tree)
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

   fChain->SetBranchAddress("id", &id, &b_id);
   fChain->SetBranchAddress("detGroup", &detGroup, &b_detGroup);
   fChain->SetBranchAddress("tubeGroup", &tubeGroup, &b_tubeGroup);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("E", &E, &b_E);
   fChain->SetBranchAddress("T", &T, &b_T);
   fChain->SetBranchAddress("partEnergy", &partEnergy, &b_partEnergy);
   fChain->SetBranchAddress("id2", &id2, &b_id2);
   fChain->SetBranchAddress("detGroup2", &detGroup2, &b_detGroup2);
   fChain->SetBranchAddress("tubeGroup2", &tubeGroup2, &b_tubeGroup2);
   fChain->SetBranchAddress("x2", &x2, &b_x2);
   fChain->SetBranchAddress("y2", &y2, &b_y2);
   fChain->SetBranchAddress("z2", &z2, &b_z2);
   fChain->SetBranchAddress("E2", &E2, &b_E2);
   fChain->SetBranchAddress("T2", &T2, &b_T2);
   fChain->SetBranchAddress("partEnergy2", &partEnergy2, &b_partEnergy2);
   Notify();
}

Bool_t makeneuneutree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void makeneuneutree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t makeneuneutree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef makeneuneutree_cxx
