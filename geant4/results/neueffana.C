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
void compare_eff(char* infile1,char* infile2)
{
    std::ifstream inpf1(infile1);

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
    cout<<nlines<<endl;
    for (Int_t i=0;i<nlines;i++){
        enearr1[i]=enearr1[i]*1000;
        effarr1[i]=effarr1[i]*100;
        efferrarr1[i]=efferrarr1[i]*100;
    }


    TGraphErrors* gr1=new TGraphErrors(nlines,enearr1,effarr1,0,efferrarr1);
    //gr1->SetMarkerStyle(21);
    //gr1->SetMarkerColor(3);
    gr1->SetLineColor(2);
    gr1->SetFillColor(0);
    //gr1->SetMarkerSize(1.2);
    gr1->SetLineWidth(2);


    std::ifstream inpf2(infile2);

    Double_t enearr2[1000];
    Double_t effarr2[1000];
    Double_t efferrarr2[1000];

    tempnum=0;
    nlines=0;
    while (inpf2.good()){
        inpf2>>tempnum>>tempnum>>tempnum>>enearr2[nlines]>>effarr2[nlines]>>efferrarr2[nlines];
        nlines++;
    }
    nlines--;
    cout<<nlines<<endl;
    for (Int_t i=0;i<nlines;i++){
        enearr2[i]=enearr2[i]*1000;
        effarr2[i]=effarr2[i]*100;
        efferrarr2[i]=efferrarr2[i]*100;
    }


    TGraphErrors* gr2=new TGraphErrors(nlines,enearr2,effarr2,0,efferrarr2);
    //gr2->SetMarkerStyle(21);
    //gr2->SetMarkerColor(3);
    gr2->SetLineColor(3);
    gr2->SetFillColor(0);
    //gr2->SetMarkerSize(1.2);
    gr2->SetLineWidth(2);



    Double_t effcomparr[500];
    for(Int_t i=0;i<nlines;i++){
        for(Int_t j=0;j<nlines;j++){
            if (enearr1[i]==enearr2[j]){
                effcomparr[i]=effarr2[j]-effarr1[i];
            }
        }
    }

    TGraph* gr3=new TGraph(nlines,enearr1,effcomparr);
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerColor(1);
    gr3->SetLineColor(1);
    gr3->SetMarkerSize(1.2);
    gr3->SetLineWidth(2);





    TCanvas* c2 = new TCanvas("c2","c2",900,700);
    c2->Divide(1,2);
    c2->cd(1);
    c2->cd(1)->SetLogx();
    c2->cd(1)->SetGrid();
    gr1->Draw("APL");
    gr2->Draw("PL SAME");
    auto leg = new TLegend(0.1,0.7,0.5,0.9);
    leg->AddEntry(gr1,infile1);
    leg->AddEntry(gr2,infile2);
    leg->Draw();

    c2->cd(2);
    c2->cd(2)->SetLogx();
    c2->cd(2)->SetGrid();
    gr3->GetYaxis()->SetRangeUser(-20,20);
    char tmp[500];
    sprintf(tmp,"%s - %s",infile1, infile2);
    gr3->SetTitle(tmp);
    gr3->Draw("APL");
    cout<<endl;

}

void compare_eff_check(char* infile1)
{
    std::ifstream inpf1(infile1);

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
    cout<<nlines<<endl;
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

    TFile* file1=TFile::Open("mysim_briken_wClover_wAIDA_4hitdistr.root");
    TGraph* greff=(TGraph*)file1->Get("effvsr");
    greff->SetFillColor(0);

    TCanvas* c2 = new TCanvas("c2","c2",900,700);
    c2->cd()->SetLogx();
    c2->cd()->SetGrid();
    gr1->Draw("APL");
    greff->Draw("P SAME");

    auto leg = new TLegend(0.1,0.7,0.5,0.9);
    leg->AddEntry(gr1,infile1);
    leg->AddEntry(greff,"mysim_briken_wClover_wAIDA_4hitdistr.root");
    leg->Draw();
}

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

void caleff(char* infilespec , char* infileeff=(char*)"mysim_briken_wClover_noAIDA.txt",char* specname=(char*)"hSpecRebin")
{
    std::ofstream outfile((char*)"caleff.txt",std::ofstream::out | std::ofstream::app);

    TFile* file1=TFile::Open(infilespec);
    TH1F* h1=(TH1F*)file1->Get(specname);
    Double_t eff1=0;
    for (Int_t i=0;i<h1->GetNbinsX();i++){
        eff1+=h1->GetBinContent(i+1)*evaleff(h1->GetBinCenter(i+1),infileeff);
    }
    Double_t eff2 = evaleff(h1->GetMean(),infileeff);
    outfile<<infilespec<<"\t"<<infileeff<<"\t"<<eff1<<"\t"<<eff2<<std::endl;
}
void caleffall()
{
    caleff((char*)"neuspecs/qrpahfspec/Ag129.mac.root");
    caleff((char*)"neuspecs/qrpahfspec/Ag130.mac.root");
    caleff((char*)"neuspecs/qrpahfspec/Cd131.mac.root");
    caleff((char*)"neuspecs/qrpahfspec/Cd132.mac.root");
    caleff((char*)"neuspecs/qrpahfspec/In131.mac.root");
    caleff((char*)"neuspecs/qrpahfspec/In132.mac.root");
    caleff((char*)"neuspecs/qrpahfspec/In133.mac.root");
    caleff((char*)"neuspecs/qrpahfspec/In134.mac.root");
    caleff((char*)"neuspecs/qrpahfspec/In135.mac.root");
    caleff((char*)"neuspecs/qrpahfspec/In136.mac.root");
    caleff((char*)"neuspecs/qrpahfspec/Sn136.mac.root");
    caleff((char*)"neuspecs/qrpahfspec/Sn137.mac.root");
    caleff((char*)"neuspecs/minatospec/Ag129.mac.root");
    caleff((char*)"neuspecs/minatospec/Ag130.mac.root");
    caleff((char*)"neuspecs/minatospec/Cd131.mac.root");
    caleff((char*)"neuspecs/minatospec/Cd132.mac.root");
    caleff((char*)"neuspecs/minatospec/In133.mac.root");
    caleff((char*)"neuspecs/minatospec/In134.mac.root");
    caleff((char*)"neuspecs/minatospec/In135.mac.root");
}

