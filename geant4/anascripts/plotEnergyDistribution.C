#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <string>
#include <TLatex.h>
#include <TStyle.h>
void plotEnergyDistribution(){

    Int_t nQRPA=12;
    Int_t nGT=7;
    TString  nameriQRPA[]={"Ag129","Ag130","Cd131","Cd132","In131","In132","In133","In134","In135","In136","Sn136","Sn137"};
    TString  nameriQRPAlatex[]={"^{129}Ag","^{130}Ag","^{131}Cd","^{132}Cd","^{131}In","^{132}In","^{133}In","^{134}In","^{135}In","^{136}In","^{136}Sn","^{137}Sn"};
    TString  nameriGT[]={"Ag129","Ag130","Cd131","Cd132","In133","In134","In135"};

    Int_t xGT[]={1,2,3,4,7,8,9};

    TH1F* hQRPA[100];
    for (Int_t i=0;i<nQRPA;i++){
        TFile* h=new TFile(Form("%sQRPA.root",(char*)nameriQRPA[i].Data()));
        TTree* BRIKEN=(TTree*) h->Get("BRIKEN");
        BRIKEN->Draw(Form("partEnergy>>hQRPA%d(100,0,10)",i),"id<140","goff");
        hQRPA[i]=(TH1F*)gDirectory->Get(Form("hQRPA%d",i));
        hQRPA[i]->SetTitle(nameriQRPAlatex[i]);
        hQRPA[i]->Scale(1/(hQRPA[i]->Integral()));//normanilzation
        hQRPA[i]->SetLineColor(1);
        hQRPA[i]->SetLineWidth(2);

     }

    TH1F* hGT[100];
    for (Int_t i=0;i<nGT;i++){
        TFile* h=new TFile(Form("%sGT.root",(char*)nameriGT[i].Data()));
        TTree* BRIKEN=(TTree*) h->Get("BRIKEN");
        BRIKEN->Draw(Form("partEnergy>>hGT%d(100,0,10)",i),"id<140","goff");
        hGT[i]=(TH1F*)gDirectory->Get(Form("hGT%d",i));
        hGT[i]->SetTitle(nameriQRPAlatex[xGT[i]-1]);
        hGT[i]->Scale(1/(hGT[i]->Integral()));//normanilzation
        hGT[i]->SetLineColor(2);
        hGT[i]->SetLineWidth(2);
     }

    TH1F* hQRPAoriginal[100];
    for (Int_t i=0;i<nQRPA;i++){
        TFile* h=new TFile(Form("neuspecs/qrpahfspec/%s.mac.root",(char*)nameriQRPA[i].Data()));
        hQRPAoriginal[i]=(TH1F*) h->Get("hSpec");
        hQRPAoriginal[i]->SetName(Form("%soriginal",(char*)nameriQRPA[i].Data()));
        hQRPAoriginal[i]->SetTitle(nameriQRPAlatex[i]);
        hQRPAoriginal[i]->Scale(1/(hQRPAoriginal[i]->Integral()));//normanilzation
        hQRPAoriginal[i]->SetLineColor(4);
        hQRPAoriginal[i]->SetLineWidth(2);
    }

    TH1F* hGToriginal[100];
    for (Int_t i=0;i<nGT;i++){
        TFile* h=new TFile(Form("neuspecs/minatospec/%s.mac.root",(char*)nameriGT[i].Data()));
        hGToriginal[i]=(TH1F*) h->Get("hSpec");
        hGToriginal[i]->SetName(Form("%soriginal",(char*)nameriGT[i].Data()));
        hGToriginal[i]->SetTitle(nameriQRPAlatex[xGT[i]-1]);
        hGToriginal[i]->Scale(1/(hQRPAoriginal[i]->Integral()));//normanilzation
        hGToriginal[i]->SetLineColor(5);
        hGToriginal[i]->SetLineWidth(2);
    }

    TCanvas * c1=new TCanvas("c1","c1",900,700);
    c1->Divide(4,3);

    for (Int_t i=0;i<nQRPA;i++){
        c1->cd(i+1);
        gStyle->SetOptStat(0);
        gStyle->SetTitleSize(0.1,"t");
        hQRPA[i]->GetYaxis()->SetTitle("Probability (normalized)");
        hQRPA[i]->GetYaxis()->SetTitleOffset(0.9);
        hQRPA[i]->GetYaxis()->SetTitleSize(0.05);
        hQRPA[i]->GetXaxis()->SetTitle("Energy (MeV)");
        hQRPA[i]->GetXaxis()->SetTitleSize(0.05);

        hQRPA[i]->Draw("hist");
        for (Int_t j=0;j<nGT;j++){
            if (xGT[j]-1==i){
                hGT[j]->Draw("hist same");
            }
        }
    }
    c1->cd(1);
    TLegend* lg=new TLegend;
    lg->AddEntry(hQRPA[0],"QRPA-HF");
    lg->AddEntry(hGT[0],"pnQRPA");
    lg->SetLineColor(0);
    lg->SetTextSize(0.08);
    lg->Draw();

    TCanvas * c2=new TCanvas("c2","c2",900,700);
    c2->cd();
    hQRPA[6]->Draw("hist");
    hQRPAoriginal[6]->Draw("hist same");
    TLegend* lg2=new TLegend;
    lg2->AddEntry(hQRPA[6],"GEANT4 spectrum");
    lg2->AddEntry(hQRPAoriginal[6],"Original spectrum");
    lg2->SetLineColor(0);
    lg2->Draw();



}
