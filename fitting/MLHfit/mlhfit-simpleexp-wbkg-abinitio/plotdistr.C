#include <TROOT.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TLine.h>
#include <TGraphErrors.h>
using namespace std;
#define kmaxsamples 10000
void plotdistr(){

    Int_t nsampleslimit=400;
    Double_t inputvalhl=0.14;
    Double_t inputvalp1n=0.5;
    Double_t inputvalp2n=0.25;

    Double_t sampleno[kmaxsamples];
    Double_t hl[kmaxsamples];
    Double_t p1n[kmaxsamples];
    Double_t p2n[kmaxsamples];

    Double_t hlerr[kmaxsamples];
    Double_t p1nerr[kmaxsamples];
    Double_t p2nerr[kmaxsamples];

    Int_t nsamples=0;


    //generate histo with background
    TH1F* hhl=new TH1F("hhl","hhl",5000,0.1,0.2);
    TH1F* hp1n=new TH1F("hp1n","hp1n",5000,0,1);
    TH1F* hp2n=new TH1F("hp2n","hp2n",5000,0,1);

    ifstream ifs("out.txt");

    for (Int_t i=0;i<8500;i++){
        ifs>>hl[nsamples]>>p1n[nsamples]>>p2n[nsamples]>>hlerr[nsamples]>>p1nerr[nsamples]>>p2nerr[nsamples];
        hlerr[nsamples]=TMath::Log(2)/hl[nsamples]/hl[nsamples]*hlerr[nsamples];
        hl[nsamples]=TMath::Log(2)/hl[nsamples];

        hhl->Fill(hl[nsamples]);
        hp1n->Fill(p1n[nsamples]);
        hp2n->Fill(p2n[nsamples]);

        sampleno[nsamples]=nsamples;
        nsamples++;
    }

    if (nsamples>nsampleslimit) nsamples=nsampleslimit;

    TCanvas* c1=new TCanvas("c1","c1",900,700);
    c1->Divide(3,2);
    c1->cd(1);
    hhl->Draw();
    c1->cd(2);
    hp1n->Draw();
    c1->cd(3);
    hp2n->Draw();

    TGraphErrors* grhl=new TGraphErrors(nsamples,sampleno,hl,0,hlerr);
    TGraphErrors* grp1n=new TGraphErrors(nsamples,sampleno,p1n,0,p1nerr);
    TGraphErrors* grp2n=new TGraphErrors(nsamples,sampleno,p2n,0,p2nerr);

    grhl->SetMarkerStyle(20);
    grhl->SetMarkerSize(1);
    grp1n->SetMarkerStyle(20);
    grp1n->SetMarkerSize(1);
    grp2n->SetMarkerStyle(20);
    grp2n->SetMarkerSize(1);

    TLine* lhl=new TLine(sampleno[0],inputvalhl,sampleno[nsamples-1],inputvalhl);
    TLine* lp1n=new TLine(sampleno[0],inputvalp1n,sampleno[nsamples-1],inputvalp1n);
    TLine* lp2n=new TLine(sampleno[0],inputvalp2n,sampleno[nsamples-1],inputvalp2n);

    lhl->SetLineWidth(2);
    lp1n->SetLineWidth(2);
    lp2n->SetLineWidth(2);
    lhl->SetLineColor(2);
    lp1n->SetLineColor(2);
    lp2n->SetLineColor(2);

    c1->cd(4);
    grhl->Draw("AP");
    lhl->Draw("same");
    c1->cd(5);
    grp1n->Draw("AP");
    lp1n->Draw("same");
    c1->cd(6);
    grp2n->Draw("AP");
    lp2n->Draw("same");
}

