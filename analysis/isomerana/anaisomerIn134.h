//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Feb 17 17:24:46 2019 by ROOT version 6.08/00
// from TChain treeIn134/
//////////////////////////////////////////////////////////

#ifndef anaisomer_h
#define anaisomer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.



const Int_t kMaxGamma = 200;
const Int_t kMaxNeutron = 200;

class anaisomer {
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
   Double_t        gc_E[kMaxGamma];   //[gc_hit]
   Double_t        gc_T[kMaxGamma];   //[gc_hit]
   Int_t           gc_ch[kMaxGamma];   //[gc_hit]
   Int_t           gc1_hit;
   Double_t        gc1_E[kMaxGamma];   //[gc1_hit]
   Double_t        gc1_T[kMaxGamma];   //[gc1_hit]
   Int_t           gc1_ch[kMaxGamma];   //[gc1_hit]
   Int_t           gc2_hit;
   Double_t        gc2_E[kMaxGamma];   //[gc2_hit]
   Double_t        gc2_T[kMaxGamma];   //[gc2_hit]
   Int_t           gc2_ch[kMaxGamma];   //[gc2_hit]
   Int_t           ab1_hit;
   Double_t        ab1_E[kMaxGamma];   //[ab1_hit]
   Double_t        ab1_T[kMaxGamma];   //[ab1_hit]
   Double_t        ab1_Tslew[kMaxGamma];   //[ab1_hit]
   Int_t           ab1_ch[kMaxGamma];   //[ab1_hit]
   Short_t         ab1_mult[kMaxGamma];   //[ab1_hit]
   Int_t           ab2_hit;
   Double_t        ab2_E[kMaxGamma];   //[ab2_hit]
   Double_t        ab2_T[kMaxGamma];   //[ab2_hit]
   Double_t        ab2_Tslew[kMaxGamma];   //[ab2_hit]
   Int_t           ab2_ch[kMaxGamma];   //[ab2_hit]
   Short_t         ab2_mult[kMaxGamma];   //[ab2_hit]
   Int_t           neu_hit;
   Double_t        neu_E[kMaxNeutron];   //[neu_hit]
   Double_t        neu_T[kMaxNeutron];   //[neu_hit]
   Int_t           neu_ch[kMaxNeutron];   //[neu_hit]

   // List of branches
   TBranch        *b_isomer;   //!
   TBranch        *b_gc_hit;   //!
   TBranch        *b_gc_E;   //!
   TBranch        *b_gc_T;   //!
   TBranch        *b_gc_ch;   //!
   TBranch        *b_gc1_hit;   //!
   TBranch        *b_gc1_E;   //!
   TBranch        *b_gc1_T;   //!
   TBranch        *b_gc1_ch;   //!
   TBranch        *b_gc2_hit;   //!
   TBranch        *b_gc2_E;   //!
   TBranch        *b_gc2_T;   //!
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

   anaisomer(TTree *tree=0);
   virtual ~anaisomer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     MakeTree(char *outfile);
   virtual void     MakeGGMatrix(char *outfile);
   virtual void     CheckTimeSlewCorr();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef anaisomer_cxx
anaisomer::anaisomer(TTree *tree) : fChain(0) 
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
      f->GetObject("treeIn134",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("treeIn134","");
      chain->Add("rootfilesexp/decay_brips3004.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3005.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3006.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3007.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3008.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3009.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3010.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3011.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3012.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3013.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3014.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3016.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3017.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3018.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3019.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3020.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3021.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3022.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3023.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3024.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3025.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3026.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3027.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3028.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3029.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3030.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3031.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3032.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3033.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3034.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3035.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3036.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3038.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3039.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3040.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3041.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3043.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3044.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3046.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3051.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3052.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3053.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3054.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3055.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3056.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3057.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3058.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3059.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3060.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3061.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3062.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3063.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3066.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3067.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3068.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3069.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3070.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3071.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3072.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3073.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3074.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3075.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3076.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3078.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3079.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3080.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3081.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3083.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3084.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3085.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3087.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3088.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3089.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3090.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3091.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3093.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3094.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3095.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3096.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3097.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3098.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3099.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3100.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3101.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3102.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3103.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3104.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3105.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3106.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3107.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3108.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3109.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3110.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3111.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3112.root/treeIn134");
      chain->Add("rootfilesexp/decay_brips3113.root/treeIn134");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

anaisomer::~anaisomer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t anaisomer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t anaisomer::LoadTree(Long64_t entry)
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

void anaisomer::Init(TTree *tree)
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
   fChain->SetBranchAddress("gc_ch", gc_ch, &b_gc_ch);
   fChain->SetBranchAddress("gc1_hit", &gc1_hit, &b_gc1_hit);
   fChain->SetBranchAddress("gc1_E", gc1_E, &b_gc1_E);
   fChain->SetBranchAddress("gc1_T", gc1_T, &b_gc1_T);
   fChain->SetBranchAddress("gc1_ch", gc1_ch, &b_gc1_ch);
   fChain->SetBranchAddress("gc2_hit", &gc2_hit, &b_gc2_hit);
   fChain->SetBranchAddress("gc2_E", gc2_E, &b_gc2_E);
   fChain->SetBranchAddress("gc2_T", gc2_T, &b_gc2_T);
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
   fChain->SetBranchAddress("neu_E", neu_E, &b_neu_E);
   fChain->SetBranchAddress("neu_T", neu_T, &b_neu_T);
   fChain->SetBranchAddress("neu_ch", neu_ch, &b_neu_ch);
   Notify();
}

Bool_t anaisomer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void anaisomer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t anaisomer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef anaisomer_cxx
