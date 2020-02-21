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
void ploteff(Double_t Ntotal=1000000){

    TFile* halle=new TFile("effcurvealle.root");

    TGraph* grAllE=(TGraph*) halle->Get("grAllE");

    TCanvas * c1=new TCanvas("c1","c1",900,700);
    c1->cd();
    c1->SetLogx();
    grAllE->Draw("APL");

    Int_t nQRPA=12;
    Int_t nGT=7;

    TString  nameriQRPA[]={"Ag129","Ag130","Cd131","Cd132","In131","In132","In133","In134","In135","In136","Sn136","Sn137"};
    TString  nameriQRPAlatex[]={"^{129}Ag","^{130}Ag","^{131}Cd","^{132}Cd","^{131}In","^{132}In","^{133}In","^{134}In","^{135}In","^{136}In","^{136}Sn","^{137}Sn"};
    TString  nameriGT[]={"Ag129","Ag130","Cd131","Cd132","In133","In134","In135"};

    Double_t xQRPA[100];
    Double_t effQRPA[100];
    for (Int_t i=0;i<nQRPA;i++){
        TFile* h=new TFile(Form("%sQRPA.root",(char*)nameriQRPA[i].Data()));
        TTree* BRIKEN=(TTree*) h->Get("BRIKEN");
        //GetTotal fire on All
        Double_t Ncounts=BRIKEN->Draw("","id<140","goff");
        effQRPA[i]=Ncounts/Ntotal*100;
        xQRPA[i]=i+1;
        cout<<nameriQRPA[i]<<"\t"<<effQRPA[i]<<endl;
        h->Close("R");
     }

    Double_t xGT[]={1,2,3,4,7,8,9};
    Double_t effGT[100];
    for (Int_t i=0;i<nGT;i++){
        TFile* h=new TFile(Form("%sGT.root",(char*)nameriGT[i].Data()));
        TTree* BRIKEN=(TTree*) h->Get("BRIKEN");
        //GetTotal fire on All
        Double_t Ncounts=BRIKEN->Draw("","id<140","goff");
        effGT[i]=Ncounts/Ntotal*100;

        cout<<nameriGT[i]<<"GT\t"<<effGT[i]<<endl;
        h->Close("R");
     }

    TCanvas * c2=new TCanvas("c2","c2",900,700);
    c2->cd();

    TGraph* grISOQRPA=new TGraph(nQRPA,xQRPA,effQRPA);
    grISOQRPA->GetYaxis()->SetRangeUser(50,72);
    grISOQRPA->SetLineWidth(0);
    grISOQRPA->SetMarkerSize(0);
    grISOQRPA->SetName("grISOQRPA");
    grISOQRPA->SetFillColor(40);
    grISOQRPA->Draw("AB");


    TGraph* grISOGT=new TGraph(nGT,xGT,effGT);
    grISOGT->SetName("grISOGT");
    grISOGT->SetMarkerStyle(29);
    grISOGT->SetMarkerColor(2);
    grISOGT->SetLineWidth(0);
    grISOGT->SetMarkerSize(3);
    grISOGT->Draw("PSAME");


    Int_t nHitDistr=11;
    Double_t xHitDistr[]={1,2,3,4,5,6,7,8,9,11,12};
    Double_t effHitDistr[100];
    Double_t effHitDistrErr[100];

    std::ifstream ifs("effhitdistributions.txt");
    for (Int_t i=0;i<nHitDistr;i++){
        std::string tempstr;
        ifs>>tempstr>>effHitDistr[i]>>effHitDistrErr[i];
        cout<<nameriQRPA[(Int_t)xHitDistr[i]-1]<<"Distr\t"<<effHitDistr[i]<<"\t"<<effHitDistrErr[i]<<endl;
    }
    TGraphErrors* grHitDistr=new TGraphErrors(nHitDistr,xHitDistr,effHitDistr,0,effHitDistrErr);
    grHitDistr->SetName("grHitDistr");
    grHitDistr->SetMarkerStyle(20);
    grHitDistr->SetMarkerColor(6);
    grHitDistr->SetMarkerSize(2);
    grHitDistr->Draw("PSAME");

    TLatex latexdraw;
    latexdraw.SetTextSize(0.03);
    for (Int_t i=0;i<nQRPA;i++){
        latexdraw.DrawLatex(xQRPA[i]-0.4,50.2,nameriQRPAlatex[i]);
    }

    TLegend* lg=new TLegend;
    lg->SetLineColor(0);
    lg->AddEntry(grISOQRPA,"QRPA-HF (M.Mumpower)");
    lg->AddEntry(grISOGT,"pnQRPA (F.Minato)");
    lg->AddEntry(grHitDistr,"Hit Distribution");
    lg->Draw();

}
