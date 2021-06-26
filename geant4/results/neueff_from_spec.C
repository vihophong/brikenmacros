#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TCanvas.h>
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMarker.h>

#include <fstream>
Double_t evaleff(Double_t E, char* infileeff)
{
    std::ifstream inpf1(infileeff);

    Double_t enearr1[1000];
    Double_t effarr1[1000];
    Double_t efferrarr1[1000];

    Double_t tempnum=0;
    Int_t nlines=0;
    while (inpf1.good()){
        inpf1>>tempnum>>tempnum>>tempnum>>enearr1[nlines]>>effarr1[nlines]>>efferrarr1[nlines];
        nlines++;
    }
    nlines--;
    //cout<<nlines<<endl;
    for (Int_t i=0;i<nlines;i++){
        enearr1[i]=enearr1[i];
        effarr1[i]=effarr1[i];
        efferrarr1[i]=efferrarr1[i];
    }

    TGraphErrors* gr1=new TGraphErrors(nlines,enearr1,effarr1,0,efferrarr1);
    //gr1->SetMarkerStyle(21);
    //gr1->SetMarkerColor(3);
    gr1->SetLineColor(2);
    gr1->SetFillColor(0);
    //gr1->SetMarkerSize(1.2);
    gr1->SetLineWidth(2);

    //TCanvas* c2 = new TCanvas("c2","c2",900,700);
    //c2->cd()->SetLogx();
    //c2->cd()->SetGrid();

    //gr1->Draw("APL");
    Double_t eff_eval=gr1->Eval(E);
    //TMarker* mk=new TMarker(E,eff_eval,20);
    //mk->SetMarkerSize(2);
    //mk->Draw();
    return eff_eval;
}
void neueff_from_spec(char* infilespec,char* specname=(char*)"hSpecRebin", char* infileeff=(char*)"mysim_briken_wClover_noAIDA.txt",Double_t Qbn=0.,Double_t Qb=0.)
{
    std::ofstream outfile((char*)"caleff.txt",std::ofstream::out | std::ofstream::app);
    TFile* file1=TFile::Open(infilespec);
    TH1F* h1=(TH1F*)file1->Get(specname);
    Double_t eff1=0;
    for (Int_t i=0;i<h1->GetNbinsX();i++){
        eff1+=h1->GetBinContent(i+1)*evaleff(h1->GetBinCenter(i+1),infileeff);
    }
    Double_t eff2 = evaleff(h1->GetMean(),infileeff);
    outfile<<infilespec<<"\t"<<infileeff<<"\t"<<Qbn<<"\t"<<Qb<<"\t"<<eff1<<"\t"<<eff2<<std::endl;
    cout<<infilespec<<"\t"<<infileeff<<"\t"<<Qbn<<"\t"<<Qb<<"\t"<<eff1<<"\t"<<eff2<<std::endl;
}

