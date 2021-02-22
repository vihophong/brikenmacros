#define makeneuneutree_cxx
#include "makeneuneutree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

typedef struct {
    Double_t ts;
    Double_t id;
    Double_t e;
    Double_t x,y,z;
    Double_t rx,ry,rz;
    Double_t vtall;
    Double_t vtf11;//f11
    Double_t vtap;//aida plastic
}neu_neu;

std::multimap < unsigned long long, unsigned int> fneumap1;
std::multimap < unsigned long long, unsigned int> fneumap2;

std::multimap < unsigned long long, unsigned int>::iterator fneumap1_it;
std::multimap < unsigned long long, unsigned int>::iterator fneumap2_it;

void makeneuneutree::Loop(char* outfile, Double_t neutronrate)
{
//   In a ROOT session, you can do:
//      root> .L makeneuneutree.C
//      root> makeneuneutree t
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

   Long64_t twlow=5e6;
   Long64_t twup=5e6;

   Long64_t nentries = fChain->GetEntries();
   cout<<nentries<<endl;
   //nentries=1000000*20;
   //nentries=1000000;

   std::multimap < unsigned long long, unsigned int> fneumap1;
   std::multimap < unsigned long long, unsigned int> fneumap2;

   std::multimap < unsigned long long, unsigned int>::iterator fneumap1_it;
   std::multimap < unsigned long long, unsigned int>::iterator fneumap2_it;

   TFile* ofile = new TFile(outfile,"recreate");
   ofile->cd();
   TTree* tree = new TTree("tree","tree");

   neu_neu neu1;
   neu_neu neu2;
   tree->Branch("neu1",&neu1,"ts/D:id/D:e/D:x/D:y/D:z/D:rx/D:ry/D:rz/D:vtall/D:vtf11/D:vtap/D");
   tree->Branch("neu2",&neu2,"ts/D:id/D:e/D:x/D:y/D:z/D:rx/D:ry/D:rz/D:vtall/D:vtf11/D:vtap/D");

   TRandom3* rseed=new TRandom3;

   //! accumulate data
   Long64_t nbytes = 0, nb = 0;
   unsigned long long tabs=1e8;//absolute time in nano seconds
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;      
      unsigned long long dtevt=rseed->Exp(1/neutronrate)*1e9; //ns unit
      tabs=tabs+dtevt;
      if (T>0){
          fneumap1.insert(make_pair(tabs+(unsigned long long)T,(unsigned int)id));
      }
      if (T2>0){
          fneumap2.insert(make_pair(tabs+(unsigned long long)T2,(unsigned int)id2));
      }
   }
   cout<<fneumap1.size()<<endl;
   cout<<fneumap2.size()<<endl;

   //! make correlation
   for (fneumap1_it=fneumap1.begin();fneumap1_it!=fneumap1.end();fneumap1_it++){
       Long64_t ts=(Long64_t)fneumap1_it->first;
       unsigned int entry=(unsigned int)fneumap1_it->second;
       neu1.ts=(Double_t)ts;
       neu1.id=entry;

       Long64_t ts1 = ts - twlow;
       Long64_t ts2 = ts + twup;
       Long64_t corrts = 0;
       Long64_t check_time = 0;
       fneumap2_it = fneumap2.lower_bound(ts1);
       while(fneumap2_it!=fneumap2.end()&&fneumap2_it->first<ts2){
           corrts = (Long64_t) fneumap2_it->first;
           unsigned int correntry=(unsigned int)fneumap2_it->second;
           if (corrts!=check_time&&corrts!=ts){
               neu2.ts=(Double_t)corrts;
               neu2.id=correntry;
               tree->Fill();
           }
           fneumap2_it++;
       }
   }
   tree->Write();
   ofile->Close();


}
