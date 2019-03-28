#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TCutG.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"

//! Global function
Double_t fcn_gen(Double_t *x, Double_t *par) {
    Double_t returnval=par[3]+(par[0]-par[3])/(1+pow(x[0]/par[2],par[1]));
    return returnval;
}

void maketimeslewfunctions()
{
    TFile* f1=new TFile("hist/timeslewhist.root");
    TH2F* hslew[8];
    hslew[0]=(TH2F*)f1->Get("hslewgc1_1");
    hslew[1]=(TH2F*)f1->Get("hslewgc1_2");
    hslew[2]=(TH2F*)f1->Get("hslewgc1_3");
    hslew[3]=(TH2F*)f1->Get("hslewgc1_4");
    hslew[4]=(TH2F*)f1->Get("hslewgc2_1");
    hslew[5]=(TH2F*)f1->Get("hslewgc2_2");
    hslew[6]=(TH2F*)f1->Get("hslewgc2_3");
    hslew[7]=(TH2F*)f1->Get("hslewgc2_4");
    TCanvas*c1=new TCanvas("c1","c1",900,700);
    c1->Divide(4,2);
    for (Int_t i=0;i<8;i++){
        c1->cd(i+1)->SetLogz();
        hslew[i]->GetYaxis()->SetRangeUser(200,1200);
        hslew[i]->GetXaxis()->SetRangeUser(0,4000);
        hslew[i]->Draw("colz");
    }

    // Projection on 8000 bin histograms
    Double_t emin[]={20,40,20,20,20,20,20,40};
    Double_t emax[]={4000,4000,4000,4000,4000,4000,4000,4000};
    Double_t ewidth=40;
    Double_t estep=1;

    Double_t fcnarrx[8][8000];
    Double_t fcnarry[8][8000];
    Double_t npntfcnarr[8];


    for (Int_t i=0;i<8;i++){
        Int_t nproj=(Int_t)((emax[i]-emin[i])/estep);
        Double_t eminj=emin[i];

        //while(emaxj>emax[i]){
        Int_t k=0;
        for (Int_t j=0;j<nproj;j++){
            Double_t emaxj=eminj+ewidth;
            if (emaxj>emax[i]) break;
            Double_t slewval=hslew[i]->ProjectionY("p1",hslew[i]->GetXaxis()->FindBin(eminj),hslew[i]->GetXaxis()->FindBin(emaxj))->GetBinCenter(hslew[i]->ProjectionY("p1",hslew[i]->GetXaxis()->FindBin(eminj),hslew[i]->GetXaxis()->FindBin(emaxj))->GetMaximumBin());
            //cout<<j<<"\t"<<eminj<<"\t"<<emaxj<<"\t"<<(emaxj+eminj)/2<<"\t";
            //cout<<slewval<<endl;

            fcnarrx[i][k]=(emaxj+eminj)/2;
            fcnarry[i][k]=slewval;
            eminj+=estep;
            k++;
        }
        npntfcnarr[i]=k;
    }

    TGraph* gr[8];
    for (Int_t i=0;i<8;i++){
        gr[i]=new TGraph(npntfcnarr[i],fcnarrx[i],fcnarry[i]);
    }


    TCanvas*c2=new TCanvas("c2","c2",900,700);
    c2->Divide(4,2);
    TF1* fB[8];


    std::ofstream ofs("slewcorrparms.txt");



    for (Int_t i=0;i<8;i++){
        fB[i]=new TF1(Form("fB%d",i),fcn_gen,0,3980,4);
        fB[i]->SetParameter(0,800);//first y
        fB[i]->SetParameter(1,0.5);//fix
        fB[i]->SetParameter(2,20.);//xmin
        fB[i]->SetParameter(3,350);//last y

        /*
        fB[i]->SetParameter(0,400);//first y
        fB[i]->SetParameter(1,0.5);//fix
        fB[i]->SetParameter(2,20);//xmin
        fB[i]->SetParameter(3,-300);//last y
        */

        fB[i]->SetNpx(500);
        fB[i]->SetLineWidth(2);
        fB[i]->SetLineColor(8);

        c2->cd(i+1);
        gr[i]->SetTitle("");
        gr[i]->GetYaxis()->SetRangeUser(200,1200);
        gr[i]->GetXaxis()->SetRangeUser(0,4000);
        gr[i]->SetMarkerStyle(20);
        gr[i]->SetMarkerColor(2);
        gr[i]->SetMarkerSize(.5);
        gr[i]->Draw("AP");
        fB[i]->SetLineColor(1);
        fB[i]->Draw("same");
        gr[i]->Fit(fB[i],"LRE0");
        ofs<<i<<"\t"<<fB[i]->GetParameter(0)<<"\t"<<fB[i]->GetParameter(1)<<"\t"<<fB[i]->GetParameter(2)<<"\t"<<fB[i]->GetParameter(3)<<endl;
    }
    ofs.close();
}
