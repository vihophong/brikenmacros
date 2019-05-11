#define countn_cxx
#include "countn.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void countn::Loop()
{
//   In a ROOT session, you can do:
//      root> .L countn.C
//      root> countn t
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

   Long64_t nbytes = 0, nb = 0;

   Double_t arr[1000];
   for (Int_t i=0;i<1000;i++) arr[i]=0;
   Int_t pt=1;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (id>139) continue;
      Double_t r=sqrt(x2*x2+y2*y2);
      if (jentry==0) arr[0]=r;
      Int_t isfound=0;
      for (Int_t i=0;i<1000;i++){
          if (r==arr[i]){
              isfound++;
          }
      }
      if (isfound==0){
          arr[pt]=r;
          pt++;
      }
      //cout<<r<<endl;
      // if (Cut(ientry) < 0) continue;
   }
   for (Int_t i=0;i<pt;i++) {
       cout<<arr[i]<<endl;
   }
}
