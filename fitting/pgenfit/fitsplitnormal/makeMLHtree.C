#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TROOT.h"

void makeMLHtree(char* inputfile, char* outputfile, Int_t findex=18,Int_t f1nindex=18,Int_t f2nindex=36 ){
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Thu Feb 20 15:28:52 2020 by ROOT version6.14/00)
//   from TTree treeout/treeout
//   found on file: outFit/In135out.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*) gROOT->GetListOfFiles()->FindObject(inputfile);
   if (!f) {
      f = new TFile(inputfile);
   }
   TTree* tree;
   f->GetObject("treeout",tree);

//Declaration of leaves types
   Double_t        pVal[150];
   Double_t        pCentralVal[150];
   Double_t        pValError[150];
   Int_t           ispVary[150];
   Int_t           ipVal[150];
   Double_t        nsigCentralVal;
   Double_t        nsigVal;
   Double_t        nbkgVal;
   Double_t        nbkgError;
   Double_t        bkg1nratioVal;
   Double_t        bkg1nratioError;
   Double_t        bkg2nratioVal;
   Double_t        bkg2nratioError;
   Double_t        slope1pos;
   Double_t        slope2pos;
   Double_t        slope3pos;
   Double_t        slope1posError;
   Double_t        slope2posError;
   Double_t        slope3posError;
   Double_t        fFitTime;
   Int_t           fitStatus;
   Int_t           fitCovQual;
   Int_t           fitNumInvalidNLL;
   Double_t        fitEdm;
   Double_t        fitMinNll;

   // Set branch addresses.
   tree->SetBranchAddress("pVal",pVal);
   tree->SetBranchAddress("pCentralVal",pCentralVal);
   tree->SetBranchAddress("pValError",pValError);
   tree->SetBranchAddress("ispVary",ispVary);
   tree->SetBranchAddress("ipVal",ipVal);
   tree->SetBranchAddress("nsigCentralVal",&nsigCentralVal);
   tree->SetBranchAddress("nsigVal",&nsigVal);
   tree->SetBranchAddress("nbkgVal",&nbkgVal);
   tree->SetBranchAddress("nbkgError",&nbkgError);
   tree->SetBranchAddress("bkg1nratioVal",&bkg1nratioVal);
   tree->SetBranchAddress("bkg1nratioError",&bkg1nratioError);
   tree->SetBranchAddress("bkg2nratioVal",&bkg2nratioVal);
   tree->SetBranchAddress("bkg2nratioError",&bkg2nratioError);
   tree->SetBranchAddress("slope1pos",&slope1pos);
   tree->SetBranchAddress("slope2pos",&slope2pos);
   tree->SetBranchAddress("slope3pos",&slope3pos);
   tree->SetBranchAddress("slope1posError",&slope1posError);
   tree->SetBranchAddress("slope2posError",&slope2posError);
   tree->SetBranchAddress("slope3posError",&slope3posError);
   tree->SetBranchAddress("fFitTime",&fFitTime);
   tree->SetBranchAddress("fitStatus",&fitStatus);
   tree->SetBranchAddress("fitCovQual",&fitCovQual);
   tree->SetBranchAddress("fitNumInvalidNLL",&fitNumInvalidNLL);
   tree->SetBranchAddress("fitEdm",&fitEdm);
   tree->SetBranchAddress("fitMinNll",&fitMinNll);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// tree->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   TFile* fileout=new TFile(outputfile,"recreate");

   TTree* outtree=new TTree("tree","tree");
   Double_t x=0;
   outtree->Branch("x",&x,"x/D");

   Long64_t nentries = tree->GetEntries();

   Long64_t nbytes = 0;
   for (Long64_t i=0; i<nentries;i++) {
      nbytes += tree->GetEntry(i);
      x=pVal[findex];
      if (pVal[f1nindex]+pVal[f2nindex]<1&&fitStatus==0)
        outtree->Fill();
   }
   outtree->Write();
   fileout->Close();
}
