#include <vector>
#include "TGraph.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TVirtualPad.h"
#include <stdio.h>
#include <stdlib.h>


void pn(Int_t choice) 
{

  c1 = new TCanvas("c1","test",10,10,1000,700);
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(12);
  gStyle->SetOptStat(0);

  Int_t NumRI;

  Int_t nprot;
  Int_t nneut;
  Int_t nmass;
  Double_t qval;
  Double_t pn0;
  Double_t pn1; 
  Double_t pn2;
  Double_t pn3;

  Double_t pn;
  
  
  NumRI = 5346;
  //  NumRI = 24;

  //  ofstream fout2;
  //  fout2.open("zzz.dat");

  ifstream fdat;
  fdat.open("Moller_FRDM_tpnff.dat");
  
  cout << "Get data" << endl;

  TH2F *hchart = new TH2F("hist","",185,-0.5,184.5,127,-0.5,126.5);
  TH2F *hchartP1 = new TH2F("hist1","",185,-0.5,184.5,127,-0.5,126.5);
  TH2F *hchartP2 = new TH2F("hist2","",185,-0.5,184.5,127,-0.5,126.5);
  TH2F *hchartP3 = new TH2F("hist3","",185,-0.5,184.5,127,-0.5,126.5);
  

  for (Int_t i=0; i<NumRI; i++) {
    fdat >> nprot >> nneut >> qval >> pn0 >> pn1 >> pn2 >> pn3;
    nmass = nneut+nprot;
/*
    cout << " Np=" << nprot 
	 << " Nn=" << nneut 
	 << " Pn0=" << pn0
	 << " Pn1=" << pn1
	 << " Pn2=" << pn2 
	 << endl;
	 */
    pn = pn0*0. + pn1*1. + pn2*2. + pn3*3.;
    //    pn = pn0*45.+pn1*63.4*pn2*71.56+pn3*75.96;
    //hchart->Fill(nneut,nprot,(pn1+pn2+pn3)*100.);
    hchart->Fill(nneut,nprot,pn0);
    hchartP1->Fill(nneut,nprot,pn1);
    hchartP2->Fill(nneut,nprot,pn2);
    hchartP3->Fill(nneut,nprot,pn3);
    //    fout2 << nprot << " " << nneut << " " << pn << endl;
    //    TArrow ar(nneut,nprot,nneut-pn,nprot-1);
    //    ar.DrawClone("same");
    //    ar.Delete();
    //    ar->Delete[];

    //    TArrow ar1(60,40,70,80,0.01);
    //    ar1.DrawClone();



  }

  //  fout2->close();

  fdat.close();

  
  
  
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


  Double_t xrange[2]={56,84};
  Double_t yrange[2]={25,55};
  hchart->GetXaxis()->SetRangeUser(xrange[0],xrange[1]);
  hchart->GetYaxis()->SetRangeUser(yrange[0],yrange[1]);

  //hchart->SetMinimum(3);
  //  hchart->SetMinimum(45.);
  //hchart->SetMaximum(50);

  hchart->SetFillColor(16);
  hchart->SetLineWidth(1);
  //hchart->SetMarkerColor(5);
  hchart->Draw("box1");
  hchartP1->SetFillColor(46);
  hchartP1->Draw("box1same");
  hchartP2->SetFillColor(6);
  hchartP2->Draw("box1same");
  hchartP3->SetFillColor(36);
  hchartP3->Draw("box1same");
  //hchart->DrawClone("boxsame");
  



  Double_t dd = 0.5;

  Double_t nn;
  Double_t pp1;
  Double_t pp2;
  //  nn=82;
  //  pp=44; TBox *b1 = new TBox(nn,pp,nn+1.,pp+1.);
  TBox b1;
  b1.SetFillColor(11);
  //  b1.SetFillColor(8);
  b1.SetLineColor(2);
  b1.SetLineWidth(1);
  //  b1.SetLineStyle(2);

  TLine a0;
  a0.SetLineWidth(4);
  //  a0.SetLineStyle(3);
  a0.SetLineStyle(1);
  a0.SetLineColor(3);

  nn=50; pp1=28; pp2=29; 


  //  nn=51; pp1=28; pp2=30; b1.DrawBox(nn-0.03-dd,pp1-dd,nn+1.03-dd,pp2+1-dd);
  //  nn=53; pp1=29; pp2=32; b1.DrawBox(nn-0.03-dd,pp1-dd,nn+1.03-dd,pp2+1-dd);
  //  nn=55; pp1=30; pp2=33; b1.DrawBox(nn-0.03-dd,pp1-dd,nn+1.03-dd,pp2+1-dd);

  
  nn=50; pp1=28; pp2=29; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=52; pp1=28; pp2=32; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=54; pp1=30; pp2=33; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=56; pp1=32; pp2=34; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=58; pp1=32; pp2=34; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=60; pp1=33; pp2=35; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=62; pp1=33; pp2=36; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=64; pp1=34; pp2=37; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=66; pp1=35; pp2=38; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=68; pp1=36; pp2=40; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=70; pp1=38; pp2=40; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=72; pp1=39; pp2=41; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=74; pp1=39; pp2=41; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=76; pp1=39; pp2=42; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=78; pp1=40; pp2=42; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=80; pp1=41; pp2=44; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=82; pp1=42; pp2=49; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=84; pp1=47; pp2=50; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=86; pp1=47; pp2=50; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=88; pp1=48; pp2=52; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=90; pp1=50; pp2=52; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=92; pp1=51; pp2=52; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=94; pp1=51; pp2=53; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=96; pp1=51; pp2=54; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=98; pp1=52; pp2=54; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=100; pp1=52; pp2=56;
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=102; pp1=53; pp2=56; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=104; pp1=54; pp2=58; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=106; pp1=55; pp2=59; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=108; pp1=55; pp2=60; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=110; pp1=56; pp2=61; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=112; pp1=57; pp2=61; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=114; pp1=57; pp2=62; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=116; pp1=58; pp2=61; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=118; pp1=58; pp2=60; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=120; pp1=59; pp2=61; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=122; pp1=59; pp2=61; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=124; pp1=60; pp2=62; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=126; pp1=61; pp2=72; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=128; pp1=69; pp2=72; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=130; pp1=70; pp2=72; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=132; pp1=70; pp2=72; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=134; pp1=70; pp2=72; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=136; pp1=70; pp2=72; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=138; pp1=70; pp2=73; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=140; pp1=71; pp2=73; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=142; pp1=71; pp2=74; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=144; pp1=72; pp2=75; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=146; pp1=73; pp2=76; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=148; pp1=74; pp2=76; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=150; pp1=74; pp2=78; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=152; pp1=76; pp2=79; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=154; pp1=77; pp2=80; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=156; pp1=76; pp2=81; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=158; pp1=78; pp2=82; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=160; pp1=79; pp2=82; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=162; pp1=80; pp2=82; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=164; pp1=80; pp2=83; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=166; pp1=80; pp2=83; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=168; pp1=81; pp2=84; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=170; pp1=81; pp2=84; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=172; pp1=81; pp2=85; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=174; pp1=82; pp2=86; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=176; pp1=83; pp2=86; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=178; pp1=83; pp2=86; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=180; pp1=84; pp2=86; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  nn=182; pp1=85; pp2=85; 
  a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5); 
  a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);

  
 

  TLine a1;
  //  a1.SetLineWidth(1.5);
  a1.SetLineWidth(3.0);
  a1.SetLineColor(9);

  // N=8
  //  a1.DrawLine(8-dd,2-dd,8-dd,16-dd); a1.DrawLine(9-dd,2-dd,9-dd,16-dd);
  // N=20
  //  a1.DrawLine(20-dd,7-dd,20-dd,30-dd); a1.DrawLine(21-dd,7-dd,21-dd,30-dd);

  // N=28
  //  a1.DrawLine(28-dd,12-dd,28-dd,35-dd); a1.DrawLine(29-dd,12-dd,29-dd,35-dd);

  // N=50
  a1.DrawLine(50-dd,25-dd,50-dd,54-dd); a1.DrawLine(51-dd,25-dd,51-dd,54-dd);

  // N=82 
  a1.DrawLine(82-dd,40-dd,82-dd,73-dd); a1.DrawLine(83-dd,40-dd,83-dd,73-dd);

  // N=126
  a1.DrawLine(126-dd,70-dd,126-dd,95-dd); a1.DrawLine(127-dd,70-dd,127-dd,95-dd);

  // Z=8
  //  a1.DrawLine(3-dd,8-dd,20-dd,8-dd); a1.DrawLine(3-dd,9-dd,20-dd,9-dd);

  // Z=20
  //  a1.DrawLine(15-dd,20-dd,40-dd,20-dd); a1.DrawLine(15-dd,21-dd,40-dd,21-dd);

  // Z=28
  a1.DrawLine(45-dd,28-dd,60-dd,28-dd); a1.DrawLine(45-dd,29-dd,60-dd,29-dd);

  // Z=50
  a1.DrawLine(48-dd,50-dd,100-dd,50-dd); a1.DrawLine(48-dd,51-dd,100-dd,51-dd);

  // Z=82
  a1.DrawLine(95-dd,82-dd,145-dd,82-dd); a1.DrawLine(95-dd,83-dd,145-dd,83-dd);
    

  ifstream fdat7;

  fdat7.open("stable.csv");
  cout << "Stable data " << endl;
  Int_t n7=287;
  Double_t xx7[300];
  Double_t yy7[300];
  
  for (Int_t i=0; i<n7; i++) {
    fdat7 >> yy7[i] >> xx7[i];
    /*
	cout << " i=" << i
	 << " np=" << yy7[i]
	 << " nn=" << xx7[i]
	 << endl;
	 */
  }
  

  TGraph *gr7 = new TGraph(n7,xx7,yy7);
  gr7->SetMarkerStyle(21);
  gr7->SetMarkerColor(1);
  gr7->SetMarkerSize(1.7);
  gr7->Draw("PS");

	  




  ifstream fdat3;
  fdat3.open("Moller_FRDM_tpnff.dat");

  for (Int_t i=0; i<NumRI; i++) {
    fdat3 >> nprot >> nneut >> qval >> pn0 >> pn1 >> pn2 >> pn3;
    nmass = nneut+nprot;
    //    if (nprot>40 && nprot<52 && nneut>60 && nneut<110) {

    pn = pn0*0. + pn1*1. + pn2*2. + pn3*3.;
    //cout << i << " " << pn << endl;

    if (pn>0) {
      TArrow ar(nneut,nprot,nneut-pn-1,nprot+1,0.008,"|>");
      if (nneut-pn-1<xrange[0]||nneut-pn-1>xrange[1]||nprot+1<yrange[0]||nprot+1>yrange[1]||nneut<xrange[0]||nneut>xrange[1]||nprot<yrange[0]||nprot>yrange[1]) continue;
      //ar.DrawClone("");
    }
    //    }
  }
  fdat3.close();

  //KTUY05 deformation parameters
  ifstream fdat4;
  fdat4.open("KTUY05_m246.dat");
  TString temp;
  fdat4>>temp>>temp>>temp>>temp>>temp>>temp>>temp;

  Double_t Mcal,Esh,alpha2,alpha4,alpha6;


  TH2F *deformChartKTUY = new TH2F("deformKTUY","",185,-0.5,184.5,127,-0.5,126.5);
  for (Int_t i=0; i<NumRI; i++) {
    fdat4 >> nprot >> nneut >> Mcal >> Esh >> alpha2 >> alpha4 >> alpha6;
    nmass = nneut+nprot;
    //    if (nprot>40 && nprot<52 && nneut>60 && nneut<110) {
    deformChartKTUY->Fill(nneut,nprot,alpha2*sqrt(TMath::Pi()*20));
    //pn = pn0*0. + pn1*1. + pn2*2. + pn3*3.;
    //cout << i << " " << alpha2 << endl;

    //    }
  }
  deformChartKTUY->SetMinimum(-0.4);
  deformChartKTUY->SetMaximum(0.4);
  if (choice==0) deformChartKTUY->Draw("cont1zsame");
  fdat4.close();

  //FRDM deformation parameters
  FILE * pFile;
  pFile = fopen ("FRDM-beta.txt","r");
  TH2F *deformChartFRDM = new TH2F("deformFRDM","",185,-0.5,184.5,127,-0.5,126.5);
  float beta2;

  char tmp[100];
  for (Int_t i=0; i<NumRI; i++) {
    fscanf (pFile, "%5d", &nprot);
    fscanf (pFile, "%5d", &nneut);
    fscanf (pFile, "%5d", &nmass);
    for (Int_t j=0;j<15;j++) {
      fscanf (pFile,"%10c",tmp);
      if (j==5) beta2=atof(tmp);
    }

    deformChartFRDM->Fill(nneut,nprot,beta2);
  }
  deformChartFRDM->SetMinimum(-0.4);
  deformChartFRDM->SetMaximum(0.4);
  if (choice==1) deformChartFRDM->Draw("cont1zsame");
  fclose(pFile);

  gStyle->SetHatchesLineWidth(5);
  ifstream fdat2;
  fdat2.open("newPn_setting1.txt");
    cout << "New Pn " << endl;
   TH2F *newPn1 = new TH2F("newPn","",185,-0.5,184.5,127,-0.5,126.5);
  for (Int_t i=0; i<98; i++) {
          fdat2 >> nprot >> nneut;
          newPn1->Fill(nneut,nprot,1);
  }
 newPn1->SetFillColor(1);
 //newPn1->SetLineWidth(3.0);
 newPn1->SetFillStyle(3021);
newPn1->Draw("boxsame");
fdat2.close();

ifstream fdat5;
  fdat5.open("newPn_setting2_plus.txt");
    cout << "New Pn " << endl;
   TH2F *newPn2 = new TH2F("newPn","",185,-0.5,184.5,127,-0.5,126.5);
  for (Int_t i=0; i<98; i++) {
          fdat5 >> nprot >> nneut;
          newPn2->Fill(nneut,nprot,1);
  }
 newPn2->SetFillColor(4);
 //newPn2->SetLineWidth(3.0);
 newPn2->SetFillStyle(3005);
newPn2->Draw("boxsame");
fdat5.close();

  ifstream fdat8;



  //fdat8.open("MeasuredPn.txt");
  fdat8.open("overlap2.txt");
  cout << "measured Pn " << endl;

  Double_t xx8[50];
  Double_t yy8[50];

  for (Int_t i=0; i<15; i++) {
    fdat8 >> yy8[i] >> xx8[i];

  }
    TGraph *gr8 = new TGraph(15,xx8,yy8);
  gr8->SetMarkerStyle(29);
  gr8->SetMarkerColor(5);
  gr8->SetLineWidth(3);
  gr8->SetMarkerSize(1.7);
  gr8->Draw("PS");


	  
  
}


