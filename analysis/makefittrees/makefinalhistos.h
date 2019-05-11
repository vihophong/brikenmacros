//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Apr 28 21:06:36 2019 by ROOT version 6.08/00
// from TChain tree/
//////////////////////////////////////////////////////////

#ifndef makefinalhistos_h
#define makefinalhistos_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.


const Double_t lowlimit=-10;
const Double_t uplimit=20;
const Int_t binning=30000;
const Int_t kMaxGamma = 500;
const Int_t kMaxNeutron = 500;
const Double_t outerinnersep = 14.5;//cm

#define MaxID 1000
#define MaxIndex1 30
#define MaxIndex2 10
#define belenClockResolution 20 //in ns


class makefinalhistos {
public :
   char nameri[200];
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TTree          *fChainImp;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       decay_evt;
   ULong64_t       decay_ts;
   Double_t        decay_t;
   Double_t        decay_x;
   Double_t        decay_y;
   Double_t        decay_ex;
   Double_t        decay_ey;
   Double_t        decay_ion_x;
   Double_t        decay_ion_y;
   Double_t        decay_ion_ex;
   Double_t        decay_ion_ey;
   Double_t        decay_zet;
   Double_t        decay_aoq;
   Double_t        decay_beta;
   Double_t        decay_deltaxy;
   Short_t         decay_z;
   Short_t         decay_ion_z;
   Short_t         decay_multx;
   Short_t         decay_multy;
   Short_t         decay_multz;
   Short_t         decay_ndecay;
   Short_t         decay_isbump;
   Int_t           gc_hit;
   Double_t        gc_E[kMaxGamma];   //[gc_hit]
   Double_t        gc_T[kMaxGamma];   //[gc_hit]
   Int_t           gc_ch[kMaxGamma];   //[gc_hit]
   Int_t           neu_hit;
   Double_t        neu_E[kMaxNeutron];   //[neu_hit]
   Double_t        neu_T[kMaxNeutron];   //[neu_hit]
   Int_t           neu_ch[kMaxNeutron];   //[neu_hit]
   Double_t        neu_x[kMaxNeutron];   //[neu_hit]
   Double_t        neu_y[kMaxNeutron];   //[neu_hit]
   Int_t           neub_hit;
   Double_t        neub_E[kMaxNeutron];   //[neub_hit]
   Double_t        neub_T[kMaxNeutron];   //[neub_hit]
   Int_t           neub_ch[kMaxNeutron];   //[neub_hit]
   Double_t        neub_x[kMaxNeutron];   //[neub_hit]
   Double_t        neub_y[kMaxNeutron];   //[neub_hit]

   // List of branches
   TBranch        *b_decay;   //!
   TBranch        *b_gc_hit;   //!
   TBranch        *b_gc_E;   //!
   TBranch        *b_gc_T;   //!
   TBranch        *b_gc_ch;   //!
   TBranch        *b_neu_hit;   //!
   TBranch        *b_neu_E;   //!
   TBranch        *b_neu_T;   //!
   TBranch        *b_neu_ch;   //!
   TBranch        *b_neu_x;   //!
   TBranch        *b_neu_y;   //!
   TBranch        *b_neub_hit;   //!
   TBranch        *b_neub_E;   //!
   TBranch        *b_neub_T;   //!
   TBranch        *b_neub_ch;   //!
   TBranch        *b_neub_x;   //!
   TBranch        *b_neub_y;   //!

   makefinalhistos(char* listfile,char* pid,TTree *tree=0);
   virtual ~makefinalhistos();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     MakeFinalHisto(char* outfile, Double_t decaytmin=0.05,Double_t decaytmax=1,Int_t layer=-1);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef makefinalhistos_cxx
makefinalhistos::makefinalhistos(char* listfile,char* pid,TTree *tree) : fChain(0)
{
    sprintf(nameri,"%s",pid);

    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
       fChainImp = 0;
       if (tree == 0) {
          // The following code should be used if you want this class to access a chain
          // of trees.
          char tempchar1[1000];
          sprintf(tempchar1,"tree%s",pid);
          TChain * chain = new TChain(tempchar1,"");
          sprintf(tempchar1,"treeimp%s",pid);
          TChain* chainimp=new TChain(tempchar1,"");
          std::ifstream ifs(listfile);
          string filelist[1000];

          Int_t nfiles=0;
          while (!ifs.eof()){
              ifs>>filelist[nfiles];
              cout<<filelist[nfiles]<<endl;
              nfiles++;
          }
          nfiles=nfiles-1;
          cout<<"There are "<<nfiles<<" files in total!"<<endl;

          for (Int_t i=0;i<nfiles;i++){
              char tempchar2[1000];
              sprintf(tempchar2,"%s/tree%s",filelist[i].c_str(),pid);
              chain->Add(tempchar2);
              sprintf(tempchar2,"%s/treeimp%s",filelist[i].c_str(),pid);
              chainimp->Add(tempchar2);
          }

          tree = chain;
          fChainImp=chainimp;
       }
       Init(tree);
}

makefinalhistos::~makefinalhistos()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t makefinalhistos::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t makefinalhistos::LoadTree(Long64_t entry)
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

void makefinalhistos::Init(TTree *tree)
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

   fChain->SetBranchAddress("decay", &decay_evt, &b_decay);
   fChain->SetBranchAddress("gc_hit", &gc_hit, &b_gc_hit);
   fChain->SetBranchAddress("gc_E", &gc_E, &b_gc_E);
   fChain->SetBranchAddress("gc_T", &gc_T, &b_gc_T);
   fChain->SetBranchAddress("gc_ch", &gc_ch, &b_gc_ch);
   fChain->SetBranchAddress("neu_hit", &neu_hit, &b_neu_hit);
   fChain->SetBranchAddress("neu_E", &neu_E, &b_neu_E);
   fChain->SetBranchAddress("neu_T", &neu_T, &b_neu_T);
   fChain->SetBranchAddress("neu_ch", &neu_ch, &b_neu_ch);
   fChain->SetBranchAddress("neu_x", &neu_x, &b_neu_x);
   fChain->SetBranchAddress("neu_y", &neu_y, &b_neu_y);
   fChain->SetBranchAddress("neub_hit", &neub_hit, &b_neub_hit);
   fChain->SetBranchAddress("neub_E", &neub_E, &b_neub_E);
   fChain->SetBranchAddress("neub_T", &neub_T, &b_neub_T);
   fChain->SetBranchAddress("neub_ch", &neub_ch, &b_neub_ch);
   fChain->SetBranchAddress("neub_x", &neub_x, &b_neub_x);
   fChain->SetBranchAddress("neub_y", &neub_y, &b_neub_y);
   Notify();
}

Bool_t makefinalhistos::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void makefinalhistos::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t makefinalhistos::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef makefinalhistos_cxx
