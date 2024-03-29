#define fitchannel_cxx
#include "fitchannel.h"

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

#include "TGraph.h"

typedef struct {
    Int_t x,y,z;
    Double_t ex,ey,adcx,adcy;
} datatype;


void fitchannel::Loop(Int_t dssd,Double_t ecut,Int_t npoint,char* outname,Int_t startch=0,Int_t stopch=16384)
{
//   In a ROOT session, you can do:
//      Root > .L fitchannel.C
//      Root > fitchannel t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
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

   Long64_t nbytes = 0, nb = 0;

   //! ouput files, trees and histograms
   char tempout[500];
   sprintf(tempout,"%s.root",outname);
   TFile* outfile=new TFile(tempout,"recreate");

   TGraph* fitGraph[16384];//from 0

   sprintf(tempout,"%s.txt",outname);
   std::ofstream ofs(tempout);

   TH2F* hfit=new TH2F("hfit","hfit",2000,0,32000,2000,0,32000);
   TH1F* hdiff=new TH1F("hdiff","hdiff",2000,-1000,1000);

   Int_t current_ii=0;
   Int_t np=0;

   Int_t nptotal = 0;
   for (Int_t ii=0;ii<16384;ii++){
       //! reject outer channels
       Int_t x_c=ii/128;
       Int_t y_c=ii%128;
       bool foolflag=false;
       if (x_c>1&&x_c<126&&y_c>1&&y_c<126){
           foolflag=true;
       }
       if (!foolflag) continue;
       nptotal++;
   }

   std::multimap < double,std::pair< int, double> > emap;
   std::multimap < double,std::pair< int, double> >::iterator emap_it;
   emap.clear();
   Int_t ik=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Int_t index=x*128+y;
      if (z==dssd&&ex>ecut&&ey>ecut)
        emap.insert(std::make_pair(adcx,std::make_pair(index,adcy)));
   }
   cout<<"emap size"<<emap.size()<<endl;

   //for (int ii=0;ii<16384;ii++){
   for (int ii=startch;ii<stopch;ii++){

       //! reject outer channels
       Int_t x_c=ii/128;
       Int_t y_c=ii%128;
       bool foolflag=false;
       if (x_c>1&&x_c<126&&y_c>1&&y_c<126){
           foolflag=true;
       }
       if (!foolflag) continue;
       np++;
       cout<<"Working on channel "<<ii<<" - number of calculated pixels = "<<np<<" / "<<nptotal<<"..."<<endl;

       Double_t exx[5000];
       Double_t eyy[5000];
       Int_t k=0;
       Double_t exprev=-1000;
       for(emap_it = emap.begin(); emap_it != emap.end(); emap_it++){
           if (emap_it->second.first==ii){
               if (exprev==-1000||emap_it->first!=exprev){
                   if (k>4999) break;
                   exx[k]=emap_it->first;
                   eyy[k]=emap_it->second.second;
                   exprev=emap_it->first;
                   hfit->Fill(exx[k],eyy[k]);
                   hdiff->Fill(eyy[k]-exx[k]);
                   k++;
               }
           }
       }
       if (k<npoint){
           ofs<<ii<<"\t"<<-99999<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<emap.size()<<"\t"<<0<<endl;
           cout<<ii<<"\t"<<-99999<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<emap.size()<<"\t"<<0<<endl;
           fitGraph[ii]=new TGraph(1);
       }else{
           //fit
           fitGraph[ii]=new TGraph(k,exx,eyy);
           fitGraph[ii]->Fit("pol1","Q");
           TF1* f=(TF1*) fitGraph[ii]->GetFunction("pol1");
           ofs<<ii<<"\t"<<f->GetParameter(0)<<"\t"<<f->GetParameter(1)<<"\t"<<f->GetParError(0)<<"\t"<<f->GetParError(1)<<"\t"<<k<<"\t"<<(Double_t)f->GetChisquare()/(Double_t)f->GetNDF()<<endl;
           cout<<ii<<"\t"<<f->GetParameter(0)<<"\t"<<f->GetParameter(1)<<"\t"<<f->GetParError(0)<<"\t"<<f->GetParError(1)<<"\t"<<k<<"\t"<<(Double_t)f->GetChisquare()/(Double_t)f->GetNDF()<<endl;
       }
       current_ii=ii;
       fitGraph[ii]->SetName(Form("g%d",ii));
       fitGraph[ii]->Write();
   }

   hfit->Write();
   hdiff->Write();
   outfile->Close();
}
void fitchannel::CheckCalib(char *calibfile, char *outfilename)
{

    //! read calib table

    Double_t dssd_cal[6][256][2];
    //clean up
    for (Int_t i=0;i<6;i++){
        for (Int_t j=0;j<256;j++){
            dssd_cal[i][j][0]=0.;
            dssd_cal[i][j][1]=1.;
        }
    }

    ifstream inpf(calibfile);
    if (inpf.fail()){
        cout<<"No Calibration table is given"<<endl;
        return;
    }

    cout<<"Start reading calibration table: "<<calibfile<<endl;
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

    //! ouput files, trees and histograms
    TFile* outfile=new TFile(outfilename,"recreate");
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
        if (jentry%1000000==0) cout<<jentry<<" / "<<nentries<<" complete , "<<(Double_t)jentry/(Double_t)nentries*100.<< "% \r"<<flush;
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        data.x=x;
        data.y=y;
        data.z=z;
        data.ex=dssd_cal[z][x][0]+adcx*dssd_cal[z][x][1];
        data.ey=dssd_cal[z][y+128][0]+adcy*dssd_cal[z][y+128][1];
        data.adcx=adcx;
        data.adcy=adcy;
        tree->Fill();
    }
    tree->Write();
    outfile->Close();

}
