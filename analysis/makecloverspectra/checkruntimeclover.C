#define checkruntimeclover_cxx
#include "checkruntimeclover.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void checkruntimeclover::Loop()
{
//   In a ROOT session, you can do:
//      root> .L checkruntimeclover.C
//      root> checkruntimeclover t
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

   Long64_t nentries = fChain->GetEntriesFast();


   Long64_t ts_prev=0;

   Long64_t entrystart=0;

   Long64_t entrysep=0;



   Long64_t tsbeg=1e18;

   Long64_t tsend=0;


   Long64_t runtime=0;

   Long64_t nbytes = 0, nb = 0;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (ts_prev>fts) {
          cout<<"Time jump at entry number "<< jentry<<endl;
          entrysep=jentry;
      }
      if ((Long64_t)fts-ts_prev>0.05e12&&jentry<10000&&jentry!=0) {
          cout<<"Start entry at "<< jentry<<endl;
          entrystart=jentry;
      }
      ts_prev=fts;

      if (fts>tsend) tsend=fts;
      //tsend=fts;
      if (fts<tsbeg) tsbeg=fts;
      // if (Cut(ientry) < 0) continue;
   }

   TH1F* hts=new TH1F("hts","hts",20000,tsbeg-1e9,tsend+1e9);

   TH1F* hclover[8];

   Double_t lowlim[]={1,2,3,4,6205,6600,7,5050};
   Double_t uplim[]={1,2,3,4,6240,6640,7,5150};
   Double_t crate[]={0,0,0,0,0,0,0,0};
   for (Int_t i=0;i<8;i++){
       hclover[i]=new TH1F(Form("h%d",i+1),Form("h%d",i+1),8000,0,8000);
   }

   tsbeg=0;
   tsend=0;

   ts_prev=0;

   Int_t ncounts=0;

   Int_t k=0;

   if (entrysep==0) entrysep=entrystart;
   for (Long64_t jentry=entrystart; jentry<entrysep;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (ts_prev>fts) cout<<"1. Time jump at entry number "<< jentry<<endl;
      ts_prev=fts;
      hts->Fill(fts);
      // if (Cut(ientry) < 0) continue;

      for (Int_t i=0;i<8;i++){
          if ((fclover==i/4+1)&&(fcloverleaf==i%4+1)) {
              hclover[i]->Fill(fen);
              crate[i]+=1.;
          }
      }
      if (k==0) tsbeg=fts;
      tsend=fts;
      k++;

      ncounts++;
   }
   runtime=tsend-tsbeg;

   cout<<"runtime 1 in second = "<<(Double_t)runtime/1e9<<endl;


   tsbeg=0;
   tsend=0;

   ts_prev=0;
   k=0;
   for (Long64_t jentry=entrysep; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (ts_prev>fts) cout<<"2. Time jump at entry number "<< jentry<<endl;
      ts_prev=fts;

      hts->Fill(fts);
      for (Int_t i=0;i<8;i++){
          if ((fclover==i/4+1)&&(fcloverleaf==i%4+1)) {
              hclover[i]->Fill(fen);
              crate[i]+=1.;
          }
      }

      if (k==0) tsbeg=fts;
      tsend=fts;
      k++;
      ncounts++;
   }

   runtime+=tsend-tsbeg;
   cout<<"runtime 2 in second = "<<(Double_t)(tsend-tsbeg)/1e9<<endl;

   cout<<"runtime total in second = "<<(Double_t)runtime/1e9<<endl;

   cout<<"ncounts= "<<ncounts<<endl;

   TFile* fout=new TFile("out.root","recreate");
   fout->cd();

   TVectorD runtimes(1);
   runtimes[0] = runtime;

   TCanvas * c1=new TCanvas("c1","c1",900,700);
   c1->Divide(2,4);
   for (Int_t i=0;i<8;i++){
       c1->cd(i+1);
       c1->cd(i+1)->SetLogy();
       hclover[i]->Draw();
       hclover[i]->Write();

       cout<<"live time "<<i<<"\t"<<(Double_t)hclover[i]->Integral(hclover[i]->FindBin(lowlim[i]),hclover[i]->FindBin(uplim[i]))/10/runtime*1e9*100<<"%, rate = "<<crate[i]/runtime*1e9<<" cps"<<endl;
   }
   hts->Write();
   runtimes.Write("runtimes");
   fout->Close();

}
