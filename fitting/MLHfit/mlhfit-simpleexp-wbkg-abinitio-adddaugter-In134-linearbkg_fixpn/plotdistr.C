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
  gStyle->SetOptStat(0);
    Int_t nsampleslimit=1000;
    Double_t inputvalhl=140;
    Double_t inputvalp1n=50;
    Double_t inputvalp2n=25;

    Double_t sampleno[kmaxsamples];
    Double_t hl[kmaxsamples];
    Double_t p1n[kmaxsamples];
    Double_t p2n[kmaxsamples];

    Double_t hlerr[kmaxsamples];
    Double_t p1nerr[kmaxsamples];
    Double_t p2nerr[kmaxsamples];

    Int_t nsamples=0;


    //generate histo with background
    TH1F* hhl=new TH1F("hhl","hhl",100,100,220);
    TH1F* hp1n=new TH1F("hp1n","hp1n",100,30,100);
    TH1F* hp2n=new TH1F("hp2n","hp2n",100,10,60);

    hhl->GetXaxis()->SetTitle("T_{1/2} (ms)");
    hhl->GetYaxis()->SetTitle("Probability (a.u)");
    hp1n->GetXaxis()->SetTitle("P_{1n} (%)");
    hp1n->GetYaxis()->SetTitle("Probability (a.u)");
    hp2n->GetXaxis()->SetTitle("P_{2n} (%)");
    hp2n->GetYaxis()->SetTitle("Probability (a.u)");


    ifstream ifs("out.txt");

    for (Int_t i=0;i<223;i++){
        ifs>>hl[nsamples]>>p1n[nsamples]>>p2n[nsamples]>>hlerr[nsamples]>>p1nerr[nsamples]>>p2nerr[nsamples];
        hlerr[nsamples]=TMath::Log(2)/hl[nsamples]/hl[nsamples]*hlerr[nsamples]*1000;
        hl[nsamples]=TMath::Log(2)/hl[nsamples]*1000;
        p1n[nsamples]=p1n[nsamples]*100;
        p2n[nsamples]=p2n[nsamples]*100;
        p1nerr[nsamples]=p1nerr[nsamples]*100;
        p2nerr[nsamples]=p2nerr[nsamples]*100;

        hhl->Fill(hl[nsamples]);
        hp1n->Fill(p1n[nsamples]);
        hp2n->Fill(p2n[nsamples]);

        sampleno[nsamples]=nsamples;
        nsamples++;
    }

    if (nsamples>nsampleslimit) nsamples=nsampleslimit;

    TCanvas* c1=new TCanvas("c1","c1",900,700);

    TLine* lhl1=new TLine(inputvalhl,0,inputvalhl,hhl->GetMaximum());
    TLine* lp1n1=new TLine(inputvalp1n,0,inputvalp1n,hp1n->GetMaximum());
    TLine* lp2n1=new TLine(inputvalp2n,0,inputvalp2n,hp2n->GetMaximum());

    c1->Divide(3,2);
    c1->cd(1);
    hhl->SetLineWidth(3);
    hhl->SetFillColor(kGray);
    hhl->Draw();
    hhl->Fit("gaus","LE+");
    lhl1->SetLineColor(2);
    lhl1->SetLineWidth(2);
    lhl1->SetLineStyle(2);
    lhl1->Draw("sane");
    c1->cd(2);
    hp1n->SetLineWidth(3);
    hp1n->SetFillColor(kGray);
    hp1n->Draw();
    hp1n->Fit("gaus","LE+");
    lp1n1->SetLineColor(2);
    lp1n1->SetLineWidth(2);
    lp1n1->SetLineStyle(2);
    lp1n1->Draw("sane");
    c1->cd(3);
    hp2n->SetLineWidth(3);
    hp2n->SetFillColor(kGray);
    hp2n->Draw();
    hp2n->Fit("gaus","LE+");
    lp2n1->SetLineColor(2);
    lp2n1->SetLineWidth(2);
    lp2n1->SetLineStyle(2);
    lp2n1->Draw("sane");

    TGraphErrors* grhl=new TGraphErrors(nsamples,sampleno,hl,0,hlerr);
    TGraphErrors* grp1n=new TGraphErrors(nsamples,sampleno,p1n,0,p1nerr);
    TGraphErrors* grp2n=new TGraphErrors(nsamples,sampleno,p2n,0,p2nerr);

    grhl->GetYaxis()->SetTitle("T_{1/2} (ms)");
    grhl->GetXaxis()->SetTitle("Random sample");
    grp1n->GetYaxis()->SetTitle("P_{1n} (%)");
    grp1n->GetXaxis()->SetTitle("Random sample");
    grp2n->GetYaxis()->SetTitle("P_{2n} (%)");
    grp2n->GetXaxis()->SetTitle("Random sample");

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

