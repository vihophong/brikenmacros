//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec  6 16:01:56 2017 by ROOT version 6.08/00
// from TChain treeIn134/
//////////////////////////////////////////////////////////

#ifndef makemlhtree_h
#define makemlhtree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TObject.h"


#define MaxID 1000
#define MaxIndex1 30
#define MaxIndex2 10
#define belenClockResolution 20 //in ns

class makemlhtree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

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
   Double_t        gc_E[9];   //[gc_hit]
   Double_t        gc_T[9];   //[gc_hit]
   Int_t           gc_ch[9];   //[gc_hit]
   Int_t           neu_hit;
   Double_t        neu_E[13];   //[neu_hit]
   Double_t        neu_T[13];   //[neu_hit]
   Int_t           neu_ch[13];   //[neu_hit]
   Double_t        neu_x[13];   //[neu_hit]
   Double_t        neu_y[13];   //[neu_hit]
   Int_t           neub_hit;
   Double_t        neub_E[16];   //[neub_hit]
   Double_t        neub_T[16];   //[neub_hit]
   Int_t           neub_ch[16];   //[neub_hit]
   Double_t        neub_x[16];   //[neub_hit]
   Double_t        neub_y[16];   //[neub_hit]

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


   TTree          *fChainImp;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrentImp; //!current Tree number in a TChain

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
   Double_t        fex;
   Double_t        fey;
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
   TBranch        *b_implant_fex;   //!
   TBranch        *b_implant_fey;   //!
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

   makemlhtree(char* listfile,char* pid,TTree *tree=0);
   virtual ~makemlhtree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(char* outfile, Int_t binning=1200,Int_t layer=-1);
   virtual void     PlotNeuHitPattern(char* outfile,Double_t decaytmin=0.05,Double_t decaytmax=1);

   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef makemlhtree_cxx
makemlhtree::makemlhtree(char* listfile,char* pid,TTree *tree) : fChain(0)
{
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

makemlhtree::~makemlhtree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t makemlhtree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t makemlhtree::LoadTree(Long64_t entry)
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

void makemlhtree::Init(TTree *tree)
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
   fChain->SetBranchAddress("gc_E", gc_E, &b_gc_E);
   fChain->SetBranchAddress("gc_T", gc_T, &b_gc_T);
   fChain->SetBranchAddress("gc_ch", gc_ch, &b_gc_ch);
   fChain->SetBranchAddress("neu_hit", &neu_hit, &b_neu_hit);
   fChain->SetBranchAddress("neu_E", neu_E, &b_neu_E);
   fChain->SetBranchAddress("neu_T", neu_T, &b_neu_T);
   fChain->SetBranchAddress("neu_ch", neu_ch, &b_neu_ch);
   fChain->SetBranchAddress("neu_x", neu_x, &b_neu_x);
   fChain->SetBranchAddress("neu_y", neu_y, &b_neu_y);
   fChain->SetBranchAddress("neub_hit", &neub_hit, &b_neub_hit);

   fCurrentImp=-1;
   fChainImp->SetMakeClass(1);
   fChainImp->SetBranchAddress("fUniqueID", &fUniqueID, &b_implant_fUniqueID);
   fChainImp->SetBranchAddress("fBits", &fBits, &b_implant_fBits);
   fChainImp->SetBranchAddress("nneuf", &nneuf, &b_implant_nneuf);
   fChainImp->SetBranchAddress("nneub", &nneub, &b_implant_nneub);
   fChainImp->SetBranchAddress("nclover", &nclover, &b_implant_nclover);
   fChainImp->SetBranchAddress("nanc", &nanc, &b_implant_nanc);
   fChainImp->SetBranchAddress("nbeam", &nbeam, &b_implant_nbeam);
   fChainImp->SetBranchAddress("fts", &fts, &b_implant_fts);
   fChainImp->SetBranchAddress("fevt", &fevt, &b_implant_fevt);
   fChainImp->SetBranchAddress("fid", &fid, &b_implant_fid);
   fChainImp->SetBranchAddress("fmult", &fmult, &b_implant_fmult);
   fChainImp->SetBranchAddress("fx", &fx, &b_implant_fx);
   fChainImp->SetBranchAddress("fy", &fy, &b_implant_fy);
   fChainImp->SetBranchAddress("fz", &fz, &b_implant_fz);
   fChainImp->SetBranchAddress("fex", &fex, &b_implant_fex);
   fChainImp->SetBranchAddress("fey", &fey, &b_implant_fey);
   fChainImp->SetBranchAddress("fnx", &fnx, &b_implant_fnx);
   fChainImp->SetBranchAddress("fny", &fny, &b_implant_fny);
   fChainImp->SetBranchAddress("fncx", &fncx, &b_implant_fncx);
   fChainImp->SetBranchAddress("fncy", &fncy, &b_implant_fncy);
   fChainImp->SetBranchAddress("fnz", &fnz, &b_implant_fnz);
   fChainImp->SetBranchAddress("fniz", &fniz, &b_implant_fniz);
   fChainImp->SetBranchAddress("frflag", &frflag, &b_implant_frflag);
   fChainImp->SetBranchAddress("fdrflag", &fdrflag, &b_implant_fdrflag);
   fChainImp->SetBranchAddress("fsumexyrank", &fsumexyrank, &b_implant_fsumexyrank);
   fChainImp->SetBranchAddress("fmindx", &fmindx, &b_implant_fmindx);
   fChainImp->SetBranchAddress("fmindy", &fmindy, &b_implant_fmindy);
   fChainImp->SetBranchAddress("ftw", &ftw, &b_implant_ftw);
   fChainImp->SetBranchAddress("fdtion", &fdtion, &b_implant_fdtion);
   fChainImp->SetBranchAddress("fdz", &fdz, &b_implant_fdz);

   Notify();
}

Bool_t makemlhtree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void makemlhtree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t makemlhtree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef makemlhtree_cxx
