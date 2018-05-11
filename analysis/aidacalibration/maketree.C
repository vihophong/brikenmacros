#define maketree_cxx
#include "maketree.h"

#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <unistd.h>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "TH1.h"
#include "TH2.h"

#include "TCutG.h"
#include "TLatex.h"
#include "TMath.h"
#include "TCanvas.h"


typedef struct {
    Int_t x,y,z;
    Double_t ex,ey,adcx,adcy;
} datatype;

void maketree::ReadPulserCalibTable(char *infile)
{

    //clean up
    for (Int_t i=0;i<6;i++){
        for (Int_t j=0;j<256;j++){
            dssd_cal[i][j][0]=0.;
            dssd_cal[i][j][1]=1.;
        }
    }

    ifstream inpf(infile);
    if (inpf.fail()){
        cout<<"No Calibration table is given"<<endl;
        return;
    }

    cout<<"Start reading calibration table: "<<infile<<endl;
    Int_t dssd_index,strip_index;
    Double_t cal1,cal2;
    Int_t mm=0;

    while (inpf.good()){
    //for (Int_t i=0;i<100;i++){
        inpf>>dssd_index>>strip_index>>cal1>>cal2;
        dssd_cal[dssd_index][strip_index][0]=cal1;
        dssd_cal[dssd_index][strip_index][1]=cal2;
        mm++;
    }
    cout<<"Read "<<mm<<" line"<<endl;
    inpf.close();
}

void maketree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L maketree.C
//      root> maketree t
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

    ReadPulserCalibTable("dummy.txt");
    //! ouput files, trees and histograms
    TFile* outfile=new TFile("outfile.root","recreate");
    datatype data;
    TTree* tree=new TTree("tree","tree");
    tree->Branch("x",&data.x,"x/I");
    tree->Branch("y",&data.y,"y/I");
    tree->Branch("z",&data.z,"z/I");
    tree->Branch("ex",&data.ex,"ex/D");
    tree->Branch("ey",&data.ey,"ey/D");
    tree->Branch("adcx",&data.adcx,"adcx/D");
    tree->Branch("adcy",&data.adcy,"adcy/D");

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<10000000;jentry++) {
      if (jentry%1000000==0) cout<<jentry<<" / "<<nentries<<" complete , "<<(Double_t)jentry/(Double_t)nentries*100.<< "% \r"<<flush;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      // EX>1500&&EY>1500&&ID==5&&nx==1&&ny==1&&EY-EX<400&&EY-EX>-400
      if (aida_ID==4&&aida_nx==1&&aida_ny==1&&aida_x>1&&aida_x<126&&aida_y>1&&aida_y<126&&aida_EY-aida_EX<1500&&aida_EY-aida_EX>-1500&&aida_EX>800){
          data.ex=aida_EX;
          data.ey=aida_EY;
          data.x=aida_x;
          data.y=aida_y;
          data.z=aida_z;
          data.adcx=(aida_EX-dssd_cal[data.z][data.x][0])/dssd_cal[data.z][data.x][1];
          data.adcy=(aida_EY-dssd_cal[data.z][data.y+128][0])/dssd_cal[data.z][data.y+128][1];
          tree->Fill();// optional cut
      }
   }
   tree->Write();
   outfile->Close();
}


void maketree::CheckCalib(char *calfile, char *outname)
{
//   In a ROOT session, you can do:
//      root> .L maketree.C
//      root> maketree t
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

    ReadPulserCalibTable(calfile);
    //! ouput files, trees and histograms
    TFile* outfile=new TFile(outname,"recreate");
    datatype data;
    TTree* tree=new TTree("tree","tree");
    tree->Branch("x",&data.x,"x/I");
    tree->Branch("y",&data.y,"y/I");
    tree->Branch("z",&data.z,"z/I");
    tree->Branch("ex",&data.ex,"ex/D");
    tree->Branch("ey",&data.ey,"ey/D");
    tree->Branch("adcx",&data.adcx,"adcx/D");
    tree->Branch("adcy",&data.adcy,"adcy/D");

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   Int_t ntemp=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<10000000;jentry++) {
      if (jentry%1000000==0) cout<<jentry<<" / "<<nentries<<" complete , "<<(Double_t)jentry/(Double_t)nentries*100.<< "% \r"<<flush;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      // EX>1500&&EY>1500&&ID==5&&nx==1&&ny==1&&EY-EX<400&&EY-EX>-400
      if (aida_ID==4&&aida_nx==1&&aida_ny==1){
          data.x=aida_x;
          data.y=aida_y;
          data.z=aida_z;
          data.adcx=aida_EX;
          data.adcy=aida_EY;
          data.ex=aida_EX*dssd_cal[data.z][data.x][1]+dssd_cal[data.z][data.x][0];
          data.ey=aida_EY*dssd_cal[data.z][data.y+128][1]+dssd_cal[data.z][data.y+128][0];
          tree->Fill();// optional cut
          ntemp++;
      }
   }
   cout<<ntemp<<" entries written!"<<endl;
   tree->Write();
   outfile->Close();
}
