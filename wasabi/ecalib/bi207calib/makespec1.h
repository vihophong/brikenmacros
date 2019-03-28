//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 18 12:26:39 2018 by ROOT version 6.08/00
// from TChain group1/
//////////////////////////////////////////////////////////

#ifndef makespec1_h
#define makespec1_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class makespec1 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           evt[128];
   Double_t        dgtz_e[128];
   Double_t        dgtz_bl[128];
   Int_t           dgtz_ch[128];
   UShort_t        dgtz_nsample[128];
   ULong64_t       dgtz_ts[128];
   UShort_t        dgtz_waveform[128][300];
   UShort_t        dgtz_sample[128][300];

   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_dgtz_e;   //!
   TBranch        *b_dgtz_bl;   //!
   TBranch        *b_dgtz_ch;   //!
   TBranch        *b_dgtz_nsample;   //!
   TBranch        *b_dgtz_ts;   //!
   TBranch        *b_dgtz_waveform;   //!
   TBranch        *b_dgtz_sample;   //!

   makespec1(TTree *tree=0);
   virtual ~makespec1();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef makespec1_cxx
makespec1::makespec1(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("group1",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("group1","");
      chain->Add("preparation/prerun00042.root/group1");
      chain->Add("preparation/prerun00043.root/group1");
      chain->Add("preparation/prerun00044.root/group1");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

makespec1::~makespec1()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t makespec1::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t makespec1::LoadTree(Long64_t entry)
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

void makespec1::Init(TTree *tree)
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

   fChain->SetBranchAddress("evt", evt, &b_evt);
   fChain->SetBranchAddress("dgtz_e", dgtz_e, &b_dgtz_e);
   fChain->SetBranchAddress("dgtz_bl", dgtz_bl, &b_dgtz_bl);
   fChain->SetBranchAddress("dgtz_ch", dgtz_ch, &b_dgtz_ch);
   fChain->SetBranchAddress("dgtz_nsample", dgtz_nsample, &b_dgtz_nsample);
   fChain->SetBranchAddress("dgtz_ts", dgtz_ts, &b_dgtz_ts);
   fChain->SetBranchAddress("dgtz_waveform", dgtz_waveform, &b_dgtz_waveform);
   fChain->SetBranchAddress("dgtz_sample", dgtz_sample, &b_dgtz_sample);
   Notify();
}

Bool_t makespec1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void makespec1::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t makespec1::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef makespec1_cxx
