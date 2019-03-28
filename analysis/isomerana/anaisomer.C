#define anaisomer_cxx
#include "anaisomer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVectorD.h>

typedef struct{
    Double_t ab_E;
    Double_t ab_T;//gamma time in ns
    Int_t ab_ch;//first hit channel
} gammaab;
void anaisomer::MakeTree(char* outfile)
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

       Int_t ab_hit;
       Double_t ab_E[2000];   //[ab_hit]
       Double_t ab_T[2000];   //[ab_hit]
       Int_t ab_ch[2000];   //[ab_hit]

       Long64_t nentries = fChain->GetEntries();

       TFile *rootfile = new TFile(outfile,"RECREATE");
       TTree *ggtree = new TTree("tree","tree");

       ggtree->Branch("ab_hit",&ab_hit,"ab_hit/I");
       ggtree->Branch("ab_E",ab_E,"ab_E[ab_hit]/D");
       ggtree->Branch("ab_T",ab_T,"ab_T[ab_hit]/D");
       ggtree->Branch("ab_ch",ab_ch,"ab_ch[ab_hit]/I");


       Long64_t nbytes = 0, nb = 0;

       Int_t ioncnt=0;

       TVectorD ioncountcontainer(2);
       ioncountcontainer[0]=0;
       ioncountcontainer[1]=0;

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
          Long64_t ientry = LoadTree(jentry);
          if (ientry < 0) break;
          nb = fChain->GetEntry(jentry);   nbytes += nb;
          // if (Cut(ientry) < 0) continue;
          //!sort the timestamp
          std::multimap<Double_t, gammaab*> abmap;
          std::multimap<Double_t, gammaab*>::iterator abmap_it;
          abmap.clear();

          //! count how many implanted ions
          if (isomer_F11L_T<800&&isomer_F11L_T>600&&isomer_F11R_T<800&&isomer_F11R_T>600&&isomer_ion_z>=0&&isomer_ion_z<5)
          {
              ioncnt++;
          }

          //! add ab1 to map
          for (Int_t i=0;i<ab1_hit;i++){
              if (isomer_F11L_T<800&&isomer_F11L_T>600&&isomer_F11R_T<800&&isomer_F11R_T>600&&isomer_ion_z>=0&&isomer_ion_z<5&&ab1_mult[i]<4){
                  gammaab* ab=new gammaab;
                  ab->ab_E=ab1_E[i];
                  ab->ab_ch=ab1_ch[i];
                  ab->ab_T=ab1_T[i]-ab1_Tslew[i];
                  abmap.insert(make_pair(ab->ab_T,ab));
              }
          }
          for (Int_t i=0;i<ab2_hit;i++){
              if (isomer_F11L_T<800&&isomer_F11L_T>600&&isomer_F11R_T<800&&isomer_F11R_T>600&&isomer_ion_z>=0&&isomer_ion_z<5&&ab2_mult[i]<4){
                  gammaab* ab=new gammaab;
                  ab->ab_E=ab2_E[i];
                  ab->ab_ch=ab2_ch[i]+4;
                  ab->ab_T=ab2_T[i]-ab2_Tslew[i];
                  abmap.insert(make_pair(ab->ab_T,ab));
              }
          }

          //! make gamma-gamma matrix
          Int_t k=0;
          Int_t k1=0;
          Int_t k2=0;

          for (abmap_it=abmap.begin();abmap_it!=abmap.end();abmap_it++){
              gammaab* ab1=(gammaab*) abmap_it->second;
              ab_E[k1]=ab1->ab_E;
              ab_T[k1]=ab1->ab_T;
              ab_ch[k1]=ab1->ab_ch;
              k1++;
              k2=0;
          }
          ab_hit=k1;
          ggtree->Fill();
       }
       ioncountcontainer[0]=ioncnt;
       ioncountcontainer[1]=nentries;

       ioncountcontainer.Write("nions");
       ggtree->Write();
       rootfile->Close();
}

void anaisomer::MakeGGMatrix(char* outfile)
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

       Int_t gg_matrix;
       Double_t gg1_E[10000];
       Double_t gg2_E[10000];
       Double_t gg1_T[10000];
       Double_t gg2_T[10000];
       Int_t gg1_ch[10000];
       Int_t gg2_ch[10000];

       Int_t ab_hit;
       Double_t ab_E[2000];   //[ab_hit]
       Double_t ab_T[2000];   //[ab_hit]
       Int_t ab_ch[2000];   //[ab_hit]

       Long64_t nentries = fChain->GetEntries();

       TFile *rootfile = new TFile(outfile,"RECREATE");
       TTree *ggtree = new TTree("ggtree","ggtree");
       ggtree->Branch("gg_matrix",&gg_matrix,"gg_matrix/I");
       ggtree->Branch("gg1_E",gg1_E,"gg1_E[gg_matrix]/D");
       ggtree->Branch("gg2_E",gg2_E,"gg2_E[gg_matrix]/D");
       ggtree->Branch("gg1_T",gg1_T,"gg1_T[gg_matrix]/D");
       ggtree->Branch("gg2_T",gg2_T,"gg2_T[gg_matrix]/D");
       ggtree->Branch("gg1_ch",gg1_ch,"gg1_ch[gg_matrix]/I");
       ggtree->Branch("gg2_ch",gg2_ch,"gg2_ch[gg_matrix]/I");

       ggtree->Branch("ab_hit",&ab_hit,"ab_hit/I");
       ggtree->Branch("ab_E",ab_E,"ab_E[ab_hit]/D");
       ggtree->Branch("ab_T",ab_T,"ab_T[ab_hit]/D");
       ggtree->Branch("ab_ch",ab_ch,"ab_ch[ab_hit]/I");



       Long64_t nbytes = 0, nb = 0;
       for (Long64_t jentry=0; jentry<nentries;jentry++) {
          Long64_t ientry = LoadTree(jentry);
          if (ientry < 0) break;
          nb = fChain->GetEntry(jentry);   nbytes += nb;
          // if (Cut(ientry) < 0) continue;
          //!sort the timestamp
          std::multimap<Double_t, gammaab*> abmap;
          std::multimap<Double_t, gammaab*> abmap2;
          std::multimap<Double_t, gammaab*>::iterator abmap_it;
          std::multimap<Double_t, gammaab*>::iterator abmap2_it;
          abmap.clear();
          abmap2.clear();
          //! add ab1 to map
          for (Int_t i=0;i<ab1_hit;i++){
              if (isomer_F11L_T<800&&isomer_F11L_T>600&&isomer_F11R_T<800&&isomer_F11R_T>600&&isomer_ion_z>=0&&isomer_ion_z<5&&ab1_mult[i]<4){
                  gammaab* ab=new gammaab;
                  ab->ab_E=ab1_E[i];
                  ab->ab_ch=ab1_ch[i];
                  ab->ab_T=ab1_T[i]-ab1_Tslew[i];

                  gammaab* abb=new gammaab;
                  abb->ab_E=ab1_E[i];
                  abb->ab_ch=ab1_ch[i];
                  abb->ab_T=ab1_T[i]-ab1_Tslew[i];
                  abmap.insert(make_pair(ab->ab_T,ab));
                  abmap2.insert(make_pair(abb->ab_T,abb));
              }
          }
          for (Int_t i=0;i<ab2_hit;i++){
              if (isomer_F11L_T<800&&isomer_F11L_T>600&&isomer_F11R_T<800&&isomer_F11R_T>600&&isomer_ion_z>=0&&isomer_ion_z<5&&ab2_mult[i]<4){
                  gammaab* ab=new gammaab;
                  ab->ab_E=ab2_E[i];
                  ab->ab_ch=ab2_ch[i]+4;
                  ab->ab_T=ab2_T[i]-ab2_Tslew[i];

                  gammaab* abb=new gammaab;
                  abb->ab_E=ab2_E[i];
                  abb->ab_ch=ab2_ch[i];
                  abb->ab_T=ab2_T[i]-ab2_Tslew[i];
                  abmap.insert(make_pair(ab->ab_T,ab));
                  abmap2.insert(make_pair(abb->ab_T,abb));
              }
          }

          //! make gamma-gamma matrix
          Int_t k=0;
          Int_t k1=0;
          Int_t k2=0;

          for (abmap_it=abmap.begin();abmap_it!=abmap.end();abmap_it++){
              gammaab* ab1=(gammaab*) abmap_it->second;
              ab_E[k1]=ab1->ab_E;
              ab_T[k1]=ab1->ab_T;
              ab_ch[k1]=ab1->ab_ch;

              k1++;
              k2=0;
              for (abmap2_it=abmap2.begin();abmap2_it!=abmap2.end();abmap2_it++){
                  gammaab* ab2=(gammaab*) abmap2_it->second;
                  if (k1==k2) continue;
                  gg1_ch[k]=ab1->ab_ch;
                  gg1_T[k]=ab1->ab_T;
                  gg1_E[k]=ab1->ab_E;
                  gg2_ch[k]=ab2->ab_ch;
                  gg2_T[k]=ab2->ab_T;
                  gg2_E[k]=ab2->ab_E;
                  k2++;
                  k++;
              }
          }
          ab_hit=k1;
          gg_matrix=k;
          ggtree->Fill();
       }

       ggtree->Write();
       rootfile->Close();


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
