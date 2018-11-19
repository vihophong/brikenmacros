#include <vector>
#include "TGraph.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TH1.h"
#include "TH2.h"
#include "TVirtualPad.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace  std;

void plainPlot(TCanvas* c1,Double_t xrange[],Double_t yrange[])
{

  Double_t minhalflife=0.0001;//100 ns

  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(12);
  gStyle->SetOptStat(0);

  Int_t NumRI;

  Int_t nprot;
  Int_t nneut;
  Int_t nmass;
  Double_t hlval;

  //TCanvas* c1=new TCanvas("c1","",900,700) ;

  
  NumRI = 5346;
  //  NumRI = 24;

  //  ofstream fout2;
  //  fout2.open("zzz.dat");

  ifstream fdat;
  fdat.open("FRDM-QRPA12-halflife.txt");
  
  cout << "Get data" << endl;

  TH2F *hchart = new TH2F("hist","",185,-0.5,184.5,127,-0.5,126.5);


  for (Int_t i=0; i<NumRI; i++) {
    fdat >> nprot >> nneut >> hlval;
    nmass = nneut+nprot;
    hchart->Fill(nneut,nprot,hlval);
  }  
  
  c1->SetLogz(0);
  hchart->SetTitleSize(0.04);
  hchart->GetXaxis()->SetTitleOffset(1.0);
  hchart->GetYaxis()->SetTitleOffset(1.2);
  hchart->GetYaxis()->CenterTitle();
  hchart->GetXaxis()->SetLabelSize(0.03);
  hchart->GetYaxis()->SetLabelSize(0.03);
  //  hchart->GetXaxis()->SetTitleSize(1.1);
  hchart->GetYaxis()->SetTitle("N_{Proton}");
  hchart->GetXaxis()->SetTitle("N_{Neutron}");

  hchart->GetXaxis()->SetRangeUser(xrange[0],xrange[1]);
  hchart->GetYaxis()->SetRangeUser(yrange[0],yrange[1]);

  hchart->SetMinimum(minhalflife);


  c1->SetLogz();
  hchart->SetLineWidth(10);
  hchart->SetLineColor(1);
  hchart->Draw("COLZ");



  //! draw border of isotopes
  TBox b2;
  b2.SetFillStyle(0);
  b2.SetLineColor(2);
  b2.SetLineWidth(1);
  fdat.seekg(0, ios::beg);
  for (Int_t i=0; i<NumRI; i++) {
    fdat >> nprot >> nneut >> hlval;
    if(nprot>=yrange[0]&&nprot<=yrange[1]&&nneut>=xrange[0]&&nneut<=xrange[1]) b2.DrawBox(nneut-0.5,nprot-0.5,nneut+0.5,nprot+0.5);
  }
  fdat.close();



  //! Drawing magic number
  Double_t dd = 0.5;
  TLine a1;
  //  a1.SetLineWidth(1.5);
  a1.SetLineWidth(3.0);
  a1.SetLineColor(7);

  Int_t magicn[]={8,20,28,50,82,126};

  for (Int_t i=0;i<6;i++){
      a1.DrawLine(magicn[i]-dd,yrange[0]-dd,magicn[i]-dd,yrange[1]+dd); a1.DrawLine(magicn[i]+1-dd,yrange[0]-dd,magicn[i]+1-dd,yrange[1]+dd);
      a1.DrawLine(xrange[0]-dd,magicn[i]-dd,xrange[1]+dd,magicn[i]-dd); a1.DrawLine(xrange[0]-dd,magicn[i]+1-dd,xrange[1]+dd,magicn[i]+1-dd);
  }


  //! Drawing r-process path
  TLine a0;
  a0.SetLineWidth(4);
  a0.SetLineStyle(1);
  a0.SetLineColor(3);
  Double_t nn;
  Double_t pp1;
  Double_t pp2;
  ifstream rpathfile("r-process_path.txt");
   while (rpathfile.good()){
       rpathfile>>nn>>pp1>>pp2;
       if (nn>=xrange[0]&&nn<=xrange[1]){
           Bool_t isplot=true;
           if (pp1<yrange[0]) pp1=yrange[0];
           if (pp1>yrange[1]) isplot=false;
           if (pp2<yrange[0]) isplot=false;
           if (pp2>yrange[1]) pp2=yrange[1];
           if (isplot){
            a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5);a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5);
            a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5);a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);
           }
       }
   }

  ifstream fdat7;
  fdat7.open("stable.csv");
  cout << "Stable data " << endl;
  Int_t n7=287;
  Double_t xx7[300];
  Double_t yy7[300];
  
  for (Int_t i=0; i<n7; i++) {
    fdat7 >> yy7[i] >> xx7[i];    
  }  
  TGraph *gr7 = new TGraph(n7,xx7,yy7);
  gr7->SetMarkerStyle(21);
  gr7->SetMarkerColor(1);
  gr7->SetMarkerSize(1.8);
  gr7->Draw("PS");
}


