#include "bigrips.h"
bigrips::bigrips(char* infile, TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(infile);
      if (!f || !f->IsOpen()) {
         f = new TFile(infile);
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

bigrips::~bigrips()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t bigrips::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

void bigrips::GetTsEntry(std::multimap<Long64_t,Long64_t> &mts)
{
   cout<<"Reading BigRIPS timestamp table"<<endl;
   Long64_t nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb =0;
   for(Long64_t jentry = 0;jentry < nentries ;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if(ientry<0) break;
      nb = fChain->GetEvent(jentry); nbytes += nb;
      Double_t percent_complete=(Double_t)jentry/(Double_t)nentries*100;
      if (jentry%50000==0) cout<<percent_complete<<" % completed"<<endl;
      if (zet>0) mts.insert(std::pair<Long64_t,Long64_t> ((Long64_t)ts,jentry));
   }
}

Long64_t bigrips::LoadTree(Long64_t entry)
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

void bigrips::Init(TTree *tree)
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

   fChain->SetBranchAddress("ts", &ts, &b_bigrips_ts);
   fChain->SetBranchAddress("sts", &sts, &b_bigrips_sts);
   fChain->SetBranchAddress("tof", &tof, &b_bigrips_tof);
   fChain->SetBranchAddress("zet", &zet, &b_bigrips_zet);
   fChain->SetBranchAddress("aoq", &aoq, &b_bigrips_aoq);
   fChain->SetBranchAddress("f5x", &f5x, &b_bigrips_f5x);
   fChain->SetBranchAddress("f11x", &f11x, &b_bigrips_f11x);
   fChain->SetBranchAddress("f11y", &f11y, &b_bigrips_f11y);
   fChain->SetBranchAddress("f11dt", &f11dt, &b_bigrips_f11dt);
   fChain->SetBranchAddress("beta", &beta, &b_bigrips_beta);
   Notify();
}

Bool_t bigrips::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void bigrips::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t bigrips::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
