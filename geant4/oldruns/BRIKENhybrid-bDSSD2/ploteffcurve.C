#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TCanvas.h>
#include <iostream>
#include <TGraph.h>
void ploteffcurve(char* path,Double_t Ntotal=1000000){


  Int_t nRun=28;

  Double_t  ke[]={5,15,25,27,35,45,52.6739,56.7,65,79.7871,85,95,105,115,125,135,145,200,300,400,500,538,600,700,800,900,1000,1332};
  Double_t  effab[100];
  Double_t  effgc[100];


  for (Int_t i=1;i<nRun;i++){
      TFile* h=new TFile(Form("%stempHist%d.root",path,i));
      TTree* BRIKEN=(TTree*) h->Get("BRIKEN");
      //GetTotal fire on All

      Double_t NcountsAB=BRIKEN->Draw("",Form("E>%f&&id==140",ke[i]-1),"goff");
      Double_t NcountsGC=0;

      for (Int_t j=0;j<8;j++){
          NcountsGC+=BRIKEN->Draw("",Form("E>%f&&id==%d",ke[i]-1,141+j),"goff");
         //cout<<BRIKEN->Draw("",Form("E>%f&&id==%d",ke[i]-1,141+j),"goff")<<endl;
      }
      effab[i]=NcountsAB/Ntotal*100;
      effgc[i]=NcountsGC/Ntotal*100;


      cout<<ke[i]<<"\t"<<effgc[i]<<"\t"<<effab[i]<<endl;
      h->Close("R");
   }
   TGraph* grAB=new TGraph(nRun,ke,effab);
   grAB->SetName("grAB");
   TGraph* grGC=new TGraph(nRun,ke,effgc);
   grGC->SetName("grGC");
   grAB->SetMarkerStyle(20);
   grGC->SetMarkerStyle(20);
   grAB->SetLineColor(2);
   grAB->Draw("APL");
   grGC->Draw("PL SAME");
}
