#define checkancrate_cxx
#include "checkancrate.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void checkancrate::Loop()
{
//   In a ROOT session, you can do:
//      root> .L checkancrate.C
//      root> checkancrate t
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

   Int_t nf11r=0;
   Int_t nfdtveto=0;
   Long64_t tsstartf11r=0;
   Long64_t tsstopf11r=0;

   Int_t nf11l=0;
   Long64_t tsstartf11l=0;
   Long64_t tsstopf11l=0;


   Long64_t tsstartdtveto=0;
   Long64_t tsstopdtveto=0;

   Int_t ntsjump=0;
   Long64_t tsprev=0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (fring==1&&fid==1){
          if (tsstartf11r==0) tsstartf11r=fts;
          tsstopf11r=fts;
          nf11r++;
      }
      if (fring==1&&fid==2){
          if (tsstartf11l==0) tsstartf11l=fts;
          tsstopf11l=fts;
          nf11l++;
      }
      if (fring==4){
          if (tsstartdtveto==0) tsstartdtveto=fts;
          tsstopdtveto=fts;
          nfdtveto++;
      }
      if (fts<tsprev) ntsjump++;
      tsprev = fts;
      // if (Cut(ientry) < 0) continue;
   }
   cout<<"f11r rate = "<<(Double_t)nf11r/(Double_t)(tsstopf11r-tsstartf11r)*1e9<<endl;
   cout<<"f11l rate = "<<(Double_t)nf11l/(Double_t)(tsstopf11l-tsstartf11l)*1e9<<endl;
   cout<<"aida veto rate = "<<(Double_t)nfdtveto/(Double_t)(tsstopdtveto-tsstartdtveto)*1e9<<endl;
   cout<<"number of time jump = "<<ntsjump<<endl;
}
