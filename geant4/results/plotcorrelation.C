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
void plotcorrelation(char* infileeff)
{
    std::ifstream inpf1(infileeff);

    Double_t enearr1[1000];
    Double_t effarr1[1000];

    Double_t tempnum=0;
    std::string tmpstring1,tmpstring2;
    Int_t nlines=0;
    while (inpf1.good()){
        inpf1>>tmpstring1>>tmpstring2>>enearr1[nlines]>>effarr1[nlines]>>tempnum;
        if (effarr1[nlines]<0.6) cout<<tmpstring1<<endl;
        nlines++;
    }
    nlines--;
    //cout<<nlines<<endl;


    TGraph* gr1=new TGraph(nlines,enearr1,effarr1);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerColor(2);
    gr1->SetMarkerSize(1.2);
    gr1->SetLineColor(2);
    gr1->SetFillColor(0);
    gr1->SetLineWidth(2);

    TCanvas* c2 = new TCanvas("c2","c2",900,700);
    c2->cd()->SetGrid();
    gr1->Draw("AP");
}
