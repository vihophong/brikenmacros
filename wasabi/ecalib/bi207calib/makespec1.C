#define makespec1_cxx
#include "makespec1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void makespec1::Loop()
{
//   In a ROOT session, you can do:
//      root> .L makespec1.C
//      root> makespec1 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Double_t xych[]={48,64,112,128};
   Long64_t nentries = fChain->GetEntries();


   TFile* outfile=new TFile("spec.root","recreate");
   TH1F *h[32];
   TH1F *hh[32];
   TH1F *horigin[32];
   for (Int_t i=0;i<32;i++){
       h[i]=new TH1F(Form("h%d",i),Form("h%d",i),500,10,2010);
       hh[i]=new TH1F(Form("hh%d",i),Form("hh%d",i),500,10,2010);
       horigin[i]=new TH1F(Form("horigin%d",i),Form("horigin%d",i),500,10,2010);

   }

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (jentry%10000==0) cout<<jentry<<"/"<<nentries<<endl;
       Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      Int_t multx=0;
      Int_t multy=0;

      Int_t multx1=0;
      Int_t multy1=0;

      Int_t k=0;

      for (Int_t i=0;i<64;i++){
          if(dgtz_e[i]>50){
              multx1++;
          }
          if (i>=xych[0]&&i<xych[1]){
              if(dgtz_e[i]>50){
                  multx++;
              }
              horigin[k]->Fill(dgtz_e[i]);
              k++;
          }

      }



      for (Int_t i=64;i<128;i++){
          if(dgtz_e[i]>50){
              multy1++;
          }
          if (i>=xych[2]&&i<xych[3]){
              if(dgtz_e[i]>50){
                  multy++;
              }
              horigin[k]->Fill(dgtz_e[i]);
              k++;
          }
      }




      k=0;
      for (Int_t i=xych[0];i<xych[1];i++){
          if (multx==1) h[k]->Fill(dgtz_e[i]);
          if (multx1==1) hh[k]->Fill(dgtz_e[i]);
          k++;
      }

      for (Int_t i=xych[2];i<xych[3];i++){
          if (multy==1) h[k]->Fill(dgtz_e[i]);
          if (multy1==1) hh[k]->Fill(dgtz_e[i]);
          k++;
      }


      //if (jentry>50000) break;
      // if (Cut(ientry) < 0) continue;
   }


   for (Int_t i=0;i<32;i++){
       h[i]->Write();
       hh[i]->Write();
       horigin[i]->Write();
   }
   outfile->Close();
}
