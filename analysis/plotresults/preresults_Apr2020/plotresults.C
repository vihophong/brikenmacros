#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
void plotresults(char* element)
{
  TMultiGraph* mg1=new TMultiGraph();
  TGraphAsymmErrors* grP1nchain=new TGraphAsymmErrors(Form("p1nchain.txt",element),"%*lg %lg %lg %lg %lg");
  grP1nchain->SetTitle("P1n Current value");

  TGraphErrors* grP1nchain_stat=new TGraphErrors(Form("p1nchain_stat.txt",element),"%*lg %lg %lg %lg");
  grP1nchain_stat->SetTitle("P1n Prev. value");

  TGraph* grQRPAHF=new TGraph(Form("p1nQRPAHF.txt",element),"%*lg %lg %lg");
  grQRPAHF->SetTitle("QRPA-HF P1n");

  grP1nchain->SetLineColor(2);
  grP1nchain->SetMarkerStyle(20);
  grP1nchain->SetMarkerColor(2);

  grP1nchain_stat->SetLineColor(4);
  grP1nchain_stat->SetMarkerSize(0);
  grP1nchain_stat->SetMarkerColor(2);

  grQRPAHF->SetLineColor(3);
  grQRPAHF->SetMarkerStyle(21);
  grQRPAHF->SetMarkerColor(3);


  mg1->Add(grP1nchain);
  mg1->Add(grP1nchain_stat);
  //mg1->Add(grQRPAHF);


  /*
  TMultiGraph* mg2=new TMultiGraph();
  TGraphAsymmErrors* grP2nchain=new TGraphAsymmErrors(Form("p2nchain.txt",element),"%*lg %lg %lg %lg");
  grP2nchain->SetTitle("P2n Current value");
  TGraphErrors* grP2nchain_stat=new TGraphErrors(Form("p2nchain_stat.txt",element),"%*lg %lg %lg %lg");
  grP2nchain_stat->SetTitle("P2n Prev. value");
  TGraph* grQRPAHFP2n=new TGraph(Form("p2nQRPAHF.txt",element),"%*lg %lg %lg");
  grQRPAHFP2n->SetTitle("QRPA-HF P2n");

  grP2nchain->SetLineColor(2);
  grP2nchain->SetMarkerStyle(20);
  grP2nchain->SetMarkerColor(2);

  grP2nchain_stat->SetLineColor(4);
  grP2nchain_stat->SetMarkerSize(0);
  grP2nchain_stat->SetMarkerColor(2);

  grQRPAHFP2n->SetLineColor(3);
  grQRPAHFP2n->SetMarkerStyle(21);
  grQRPAHFP2n->SetMarkerColor(3);


  mg2->Add(grP2nchain);
  mg2->Add(grP2nchain_stat);
  mg2->Add(grQRPAHFP2n);

  */

  TCanvas* c1=new TCanvas("c1","c1",900,900);
  c1->Divide(1,2);
  c1->cd(1);
  mg1->Draw("ap");
  c1->cd(1)->BuildLegend();
  grQRPAHF->Draw("pl same");
  c1->cd(2);
  //mg2->Draw("ap");
  //c1->cd(2)->BuildLegend();
  //grQRPAHFP2n->Draw("pl same");
}
