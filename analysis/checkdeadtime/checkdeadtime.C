#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TString.h"

#include <fstream>

TTree* chain;
TH1F* hdeadtimepulser;
TH1F* hdeadtimepulser_dtpulser;

void checkdeadtime(char* filename)
{
    TFile* fin=TFile::Open(filename);
    fin->GetObject("treeneuveto",chain);
    fin->GetObject("h1deadtime",hdeadtimepulser);
    fin->GetObject("h1deadtime3",hdeadtimepulser_dtpulser);

    //Declaration of leaves types
    Int_t           id;
    Long64_t        ts;

    // Set branch addresses.
    chain->SetBranchAddress("id",&id);
    chain->SetBranchAddress("ts",&ts);

    Long64_t mints=0;
    Long64_t maxts=0;

    Long64_t nentries = chain->GetEntries();
    cout<<nentries<<endl;
    Long64_t ts_prev=0;
    Int_t ntsreset=0;
    Long64_t nbytes = 0;
    for (Long64_t i=0; i<nentries;i++) {
       nbytes += chain->GetEntry(i);
       if (id==0){
           if (ts<mints||mints==0) mints=ts;
           if (ts>maxts||maxts==0) maxts=ts;
           if (ts_prev>ts) ntsreset++;
           ts_prev=ts;
       }
    }
    cout<<mints<<endl;
    cout<<maxts<<endl;

    std::ofstream str("out_deadtime.txt",std::ios::app);

    Bool_t ts_reset_flag=false;
    if (maxts<mints||ntsreset>0) {
        str<<filename<<"\ttime_jump_detected\t"<<ntsreset<<endl;
        ts_reset_flag=true;
    }
    TH2I* h2=new TH2I("h2","h2",20000,0,nentries,20000,mints-10000,maxts+10000);


    Long64_t beg_dt=0;
    Long64_t end_dt=0;
    Long64_t time_window=400000;//400 us
    Int_t ntout=0;
    Long64_t totaldeadtime=0;
    for (Long64_t i=0; i<nentries;i++) {
       nbytes += chain->GetEntry(i);
       if (id==0){
          //h2->Fill(i,ts);
          if (beg_dt==0) {
              beg_dt=ts;
              end_dt=ts+time_window;
          }else{
              if (ts<end_dt) {//within windows
                  end_dt=ts+time_window;
              }else{//timeout
                  Long64_t deadtime=end_dt-beg_dt;
                  if (deadtime<0) cout<<"something wrong!"<<endl;
                  totaldeadtime+=deadtime;
                  ntout++;
                  beg_dt=ts;
                  end_dt=ts+time_window;
              }
          }

       }//if id=0
    }//end entries loop

    if(!ts_reset_flag) str<<filename<<"\t"<<nentries<<"\t"<<ntout<<"\t"<<totaldeadtime<<"\t"<<maxts-mints<<"\t"<<(Double_t)totaldeadtime/(Double_t)(maxts-mints)*100<<"\t"<<hdeadtimepulser->GetMean()<<"\t"<<hdeadtimepulser_dtpulser->GetMean()<<endl;

    //h2->SaveAs("h2ts.root");

}
