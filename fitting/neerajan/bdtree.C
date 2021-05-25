#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"

void bdtree(){
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Thu Mar  4 19:46:59 2021 by ROOT version6.14/00)
//   from TTree bdtree/decay data tree
//   found on file: Rb102_decaydata.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TTree* tree;
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Rb102_decaydata.root");
   if (!f) {
      f = new TFile("Rb102_decaydata.root");
   }
    f->GetObject("bdtree",tree);
   TH1F* hdecayall = (TH1F*) f->Get("decayall");
   TH1F* hdecay1nbkg = (TH1F*) f->Get("decay1n_bkg");
   TH1F* hdecay2nbkg = (TH1F*) f->Get("decay2n_bkg");

   TH1F* hdecayallc = (TH1F*) hdecayall->Clone();
   hdecayallc->SetName("hdecay");
   TH1F* hdecay1nbkgc = (TH1F*) hdecay1nbkg->Clone();
   hdecay1nbkgc->SetName("hdecay1nbwd");
   TH1F* hdecay2nbkgc = (TH1F*) hdecay2nbkg->Clone();
   hdecay2nbkgc->SetName("hdecay2nbwd");
   TH1F* hdecaygt0nbwd = (TH1F*) hdecay1nbkg->Clone();
   hdecaygt0nbwd->Add(hdecay2nbkg);
   hdecaygt0nbwd->SetName("hdecaygt0nbwd");

//Declaration of leaves types
   Double_t        decay_time;
   Int_t           neu_mult;
   Int_t           bkg_neu_mult;

   // Set branch addresses.
   tree->SetBranchAddress("decay_time",&decay_time);
   tree->SetBranchAddress("neu_mult",&neu_mult);
   tree->SetBranchAddress("bkg_neu_mult",&bkg_neu_mult);

   TFile* outfile=new TFile("Rb102_decaydata_rmk.root","RECREATE");
   TTree* treereco = new TTree("tree","tree");
   TTree* treerecob = new TTree("treebw","treebw");
   Double_t x;
   Int_t y;
   Int_t y_2;
   treereco->Branch("x",&x,"x/D");
   treereco->Branch("y",&y,"y/I");
   treerecob->Branch("x",&x,"x/D");
   treerecob->Branch("y",&y_2,"y/I");

   Long64_t nentries = tree->GetEntries();

   Long64_t nbytes = 0;
   for (Long64_t i=0; i<nentries;i++) {
      nbytes += tree->GetEntry(i);
      x = decay_time/1000.;
      y = neu_mult;
      y_2 = bkg_neu_mult;
      treereco->Fill();
      treerecob->Fill();
   }
   hdecayallc->Write();
   hdecay1nbkgc->Write();
   hdecay2nbkgc->Write();
   hdecaygt0nbwd->Write();
   treereco->Write();
   treerecob->Write();

   outfile->Close();
}
