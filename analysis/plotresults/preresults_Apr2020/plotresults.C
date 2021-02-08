#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TAxis.h>
void plotresults(char* element)
{
  TMultiGraph* mg1=new TMultiGraph();
  TGraphAsymmErrors* grP1nchain=new TGraphAsymmErrors(Form("p1nchain%s.txt",element),"%*lg %lg %lg %lg %lg");
  grP1nchain->SetTitle("P1n stat.");

  TGraphErrors* grP1nchain_stat=new TGraphErrors(Form("p1nchain_stat%s.txt",element),"%*lg %lg %lg %lg");
  grP1nchain_stat->SetTitle("P1n sys.");

  TGraph* grQRPAHF1n=new TGraph(Form("p1nQRPAHF%s.txt",element),"%*lg %lg %lg");
  grQRPAHF1n->SetTitle("QRPA-HF P1n");

  grP1nchain->SetLineColor(2);
  grP1nchain->SetMarkerStyle(20);
  grP1nchain->SetMarkerColor(2);
  grP1nchain->SetFillColor(0);

  grP1nchain_stat->SetLineColor(4);
  grP1nchain_stat->SetMarkerSize(0);
  grP1nchain_stat->SetMarkerColor(2);
  grP1nchain_stat->SetFillColor(0);

  grQRPAHF1n->SetLineColor(3);
  grQRPAHF1n->SetFillColor(0);

  grQRPAHF1n->SetMarkerStyle(21);
  grQRPAHF1n->SetMarkerColor(3);



  mg1->Add(grP1nchain);
  mg1->Add(grP1nchain_stat);
  //mg1->Add(grQRPAHF1n);


  TMultiGraph* mg2=new TMultiGraph();
  TGraphAsymmErrors* grp2nchain=new TGraphAsymmErrors(Form("p2nchain%s.txt",element),"%*lg %lg %lg %lg %lg");
  grp2nchain->SetTitle("P2n stat.");

  TGraphErrors* grp2nchain_stat=new TGraphErrors(Form("p2nchain_stat%s.txt",element),"%*lg %lg %lg %lg");
  grp2nchain_stat->SetTitle("P2n sys.");

  TGraph* grQRPAHF2n=new TGraph(Form("p2nQRPAHF%s.txt",element),"%*lg %lg %lg");
  grQRPAHF2n->SetTitle("QRPA-HF P2n");

  grp2nchain->SetLineColor(2);
  grp2nchain->SetMarkerStyle(20);
  grp2nchain->SetMarkerColor(2);
  grp2nchain->SetFillColor(0);

  grp2nchain_stat->SetLineColor(4);
  grp2nchain_stat->SetMarkerSize(0);
  grp2nchain_stat->SetMarkerColor(2);
  grp2nchain_stat->SetFillColor(0);

  grQRPAHF2n->SetLineColor(3);
  grQRPAHF2n->SetFillColor(0);

  grQRPAHF2n->SetMarkerStyle(21);
  grQRPAHF2n->SetMarkerColor(3);



  mg2->Add(grp2nchain);
  mg2->Add(grp2nchain_stat);
  //mg2->Add(grQRPAHF2n);


  TCanvas* c1=new TCanvas("c1","c1",900,900);
  //c1->cd();
  c1->Divide(1,2);
  c1->cd(1);
  mg1->Draw("ap");
  mg1->GetXaxis()->SetTitle("Mass number");
  mg1->GetYaxis()->SetTitle("P1n (%)");
  mg1->GetYaxis()->SetRangeUser(0,120);
  TLegend* leg1=c1->cd(1)->BuildLegend();
  leg1->SetLineColor(0);
  grQRPAHF1n->Draw("l same");

  c1->cd(2);
  mg2->Draw("ap");
  mg2->GetXaxis()->SetTitle("Mass number");
  mg2->GetYaxis()->SetTitle("P1n (%)");
  mg2->GetYaxis()->SetRangeUser(0,120);
  TLegend* leg2=c1->cd(2)->BuildLegend();
  leg2->SetLineColor(0);
  grQRPAHF2n->Draw("l same");
}
