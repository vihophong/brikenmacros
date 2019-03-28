#define anaisomer_cxx
#include "anaisomer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

Double_t fcn_gen(Double_t *x, Double_t *par) {
    Double_t returnval=par[3]+(par[0]-par[3])/(1+pow(x[0]/par[2],par[1]));
    return returnval;
}
void anaisomer::MakeGGMatrix()
{
    //   In a ROOT session, you can do:
    //      root> .L anaisomer.C
    //      root> anaisomer t
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

       Long64_t nentries = fChain->GetEntries();
       TTree *ggtree = new TTree("ggtree","ggtree");
       ggtree->Branch("t",&t,"t/D");

       ggtree->Branch("gg_matrix",&ab_matrix,"gg_matrix/I");
       ggtree->Branch("gg1_E",ab1_E,"gg1_E[ab_matrix]/D");
       ggtree->Branch("gg2_E",ab2_E,"gg2_E[ab_matrix]/D");
       ggtree->Branch("gg1_T",ab1_T,"gg1_T[ab_matrix]/D");
       ggtree->Branch("gg2_T",ab2_T,"gg2_T[ab_matrix]/D");

       ggtree->Branch("ab_hit",&ab_hit,"ab_hit/I");
       ggtree->Branch("ab_E",ab_E,"ab_E[ab_hit]/D");
       ggtree->Branch("ab_T",ab_T,"ab_T[ab_hit]/D");


       Long64_t nbytes = 0, nb = 0;
       for (Long64_t jentry=0; jentry<nentries;jentry++) {
          Long64_t ientry = LoadTree(jentry);
          if (ientry < 0) break;
          nb = fChain->GetEntry(jentry);   nbytes += nb;
          // if (Cut(ientry) < 0) continue;
          for (Int_t i=0;i<ab1_hit;i++){
              if (isomer_F11L_T<800&&isomer_F11L_T>600&&isomer_F11R_T<800&&isomer_F11R_T>600&&isomer_ion_z>=0&&isomer_ion_z<5){

              }
          }

       }

       TCanvas*c1=new TCanvas("c1","c1",10,10,700,900);
       c1->Divide(1,2);
       c1->cd(1);
       gPad->SetBottomMargin(0.001);
       gPad->SetTopMargin(0.1);
       gPad->SetRightMargin(0.01);
       h1->Draw("colz");
       c1->cd(2);
       gPad->SetBottomMargin(0.1);
       gPad->SetTopMargin(0.001);
       gPad->SetRightMargin(0.01);
       h2->Draw("colz");
}

void anaisomer::CheckTimeSlewCorr()
{
//   In a ROOT session, you can do:
//      root> .L anaisomer.C
//      root> anaisomer t
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

   Long64_t nentries = fChain->GetEntries();


   // time slew function
   Double_t a[8];
   Double_t b[8];
   Double_t c[8];
   Double_t d[8];
   std::ifstream ifs("slewcorrparms.txt");
   for (Int_t i=0;i<8;i++){
       Int_t temp;
       ifs>>temp>>a[i]>>b[i]>>c[i]>>d[i];
   }

   //TH1F* h1=new TH1F("h1","",20000,0.45,0.51);
   TH2F* h1=new TH2F("h1","",8000,0,4000,300,-1000,5000);
   TH2F* h2=new TH2F("h2","",8000,0,4000,300,-1000,5000);
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      for (Int_t i=0;i<gc1_hit;i++){
          if (isomer_F11L_T<800&&isomer_F11L_T>600&&isomer_F11R_T<800&&isomer_F11R_T>600&&isomer_ion_z>=0&&isomer_ion_z<5&&gc1_ch[i]==1){
              h1->Fill(gc1_E[i],gc1_T[i]);
              h2->Fill(gc1_E[i],gc1_T[i]-(d[gc1_ch[i]-1]+(a[gc1_ch[i]-1]-d[gc1_ch[i]-1])/(1+pow(gc1_E[i]/c[gc1_ch[i]-1],b[gc1_ch[i]-1]))));
          }
      }
   }

   TCanvas*c1=new TCanvas("c1","c1",10,10,700,900);
   c1->Divide(1,2);
   c1->cd(1);
   gPad->SetBottomMargin(0.001);
   gPad->SetTopMargin(0.1);
   gPad->SetRightMargin(0.01);
   h1->Draw("colz");
   c1->cd(2);
   gPad->SetBottomMargin(0.1);
   gPad->SetTopMargin(0.001);
   gPad->SetRightMargin(0.01);
   h2->Draw("colz");
}
