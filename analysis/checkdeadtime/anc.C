#define anc_cxx
#include "anc.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void anc::Loop()
{
//   In a ROOT session, you can do:
//      root> .L anc.C
//      root> anc t
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


   std::multimap < unsigned long long, unsigned int> fF11RMap;
   std::multimap < unsigned long long, unsigned int>::iterator fF11MapR_it;


   std::multimap < unsigned long long, unsigned int> fAllveto;
   std::multimap < unsigned long long, unsigned int>::iterator fAllveto_it;

   std::multimap < unsigned long long, unsigned int> fAIDAVeto;
   std::multimap < unsigned long long, unsigned int>::iterator fAIDAVeto_it;

   std::multimap < unsigned long long, unsigned int> fdtpulserMap;
   std::multimap < unsigned long long, unsigned int>::iterator fdtpulserMap_it;

   TH1F* h1=new TH1F("h1","h1",200,0,200000);
   TH1F* h2=new TH1F("h2","h2",200,0,200000000);
   TH1F* h3=new TH1F("h3","h3",200,-1000,1000);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   std::cout<<"collect map"<<std::endl;

   Long64_t tprev=0;
   Long64_t tbegin=0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (fring==5) {
          fdtpulserMap.insert(std::make_pair(fts,jentry));
      }
      if ((fring==1&&fid==1)||(fring==1&&fid==2)||(fring==4)){
          if (tbegin==0) tbegin=fts;
          fAllveto.insert(std::make_pair(fts,jentry));
          h2->Fill(fts-tprev);
          tprev=fts;
      }

      if (fring==1&&fid==1){
          if (tbegin==0) tbegin=fts;
          fF11RMap.insert(std::make_pair(fts,jentry));
          h2->Fill(fts-tprev);
          tprev=fts;
      }

      if (fring==4){
          if (tbegin==0) tbegin=fts;
          fAIDAVeto.insert(std::make_pair(fts,jentry));
          h2->Fill(fts-tprev);
          tprev=fts;
      }

   }
   std::cout<<fdtpulserMap.size()<<"^^"<<fAllveto.size()<<"~~"<<fF11RMap.size()<<"~~"<<fAIDAVeto.size()<<std::endl;
   std::cout<<"time="<<(Double_t)(tprev-tbegin)/1e9<<std::endl;
   Int_t ncorr=0;
   for (fdtpulserMap_it=fdtpulserMap.begin();fdtpulserMap_it!=fdtpulserMap.end();fdtpulserMap_it++){
       unsigned long long ts=fdtpulserMap_it->first;
       //! Correlate imp with bigrips
       Long64_t ts1 = (Long64_t)ts - 400000;
       Long64_t ts2 = (Long64_t)ts + 0;
       Long64_t corrts = 0;
       Long64_t check_time = 0;
       fAllveto_it = fAllveto.lower_bound(ts1);
       while(fAllveto_it!=fAllveto.end()&&fAllveto_it->first<ts2){
           corrts = (Long64_t) fAllveto_it->first;
           if (corrts!=check_time){
               check_time=corrts;
               h1->Fill(ts-corrts);
               ncorr++;
               break;
           }
           fAllveto_it++;
       }
   }


   Int_t ncorr2=0;
   for (fF11MapR_it=fF11RMap.begin();fF11MapR_it!=fF11RMap.end();fF11MapR_it++){
       unsigned long long ts=fF11MapR_it->first;
       //! Correlate imp with bigrips
       Long64_t ts1 = (Long64_t)ts - 400;
       Long64_t ts2 = (Long64_t)ts + 400;
       Long64_t corrts = 0;
       Long64_t check_time = 0;
       fAIDAVeto_it = fAIDAVeto.lower_bound(ts1);
       while(fAIDAVeto_it!=fAIDAVeto.end()&&fAIDAVeto_it->first<ts2){
           corrts = (Long64_t) fAIDAVeto_it->first;
           if (corrts!=check_time){
               check_time=corrts;
               h3->Fill(ts-corrts);
               ncorr2++;
               break;
           }
           fAIDAVeto_it++;
       }
   }
   std::cout<<ncorr<<"-"<<ncorr2<<std::endl;
   h3->Draw();
}
