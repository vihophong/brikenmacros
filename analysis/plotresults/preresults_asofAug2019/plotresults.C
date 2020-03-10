#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
void plotresults()
{
    TMultiGraph* mg1=new TMultiGraph();
  TGraphErrors* grAg=new TGraphErrors("agchain.txt","%*lg %lg %lg %lg");
  grAg->SetTitle("Ag Isotopes");
  TGraphErrors* grCd=new TGraphErrors("cdchain.txt","%*lg %lg %lg %lg");
  grCd->SetTitle("Cd Isotopes");
  TGraphErrors* grIn=new TGraphErrors("inchain.txt","%*lg %lg %lg %lg");
  grIn->SetTitle("In Isotopes");
  TGraphErrors* grSn=new TGraphErrors("snchain.txt","%*lg %lg %lg %lg");
  grSn->SetTitle("Sn Isotopes");

  grAg->SetLineColor(2);
  grAg->SetMarkerStyle(20);
  grAg->SetMarkerColor(2);

  grCd->SetLineColor(3);
  grCd->SetMarkerStyle(21);
  grCd->SetMarkerColor(3);

  grIn->SetLineColor(4);
  grIn->SetMarkerStyle(22);
  grIn->SetMarkerColor(4);

  grSn->SetLineColor(5);
  grSn->SetMarkerStyle(23);
  grSn->SetMarkerColor(5);


  mg1->Add(grAg);
  mg1->Add(grCd);
  mg1->Add(grIn);
  mg1->Add(grSn);



  TMultiGraph* mg2=new TMultiGraph();
  TGraphErrors* grAgP2n=new TGraphErrors("agchainP2n.txt","%*lg %lg %lg %lg");
  TGraphErrors* grCdP2n=new TGraphErrors("cdchainP2n.txt","%*lg %lg %lg %lg");
  TGraphErrors* grInP2n=new TGraphErrors("inchainP2n.txt","%*lg %lg %lg %lg");
  TGraphErrors* grSnP2n=new TGraphErrors("snchainP2n.txt","%*lg %lg %lg %lg");

  grAgP2n->SetLineColor(2);
  grAgP2n->SetMarkerStyle(20);
  grAgP2n->SetMarkerColor(2);

  grCdP2n->SetLineColor(3);
  grCdP2n->SetMarkerStyle(21);
  grCdP2n->SetMarkerColor(3);

  grInP2n->SetLineColor(4);
  grInP2n->SetMarkerStyle(22);
  grInP2n->SetMarkerColor(4);

  grSnP2n->SetLineColor(5);
  grSnP2n->SetMarkerStyle(23);
  grSnP2n->SetMarkerColor(5);

  mg2->Add(grAgP2n);
  mg2->Add(grCdP2n);
  mg2->Add(grInP2n);
  mg2->Add(grSnP2n);


  TCanvas* c1=new TCanvas("c1","c1",900,900);
  c1->Divide(1,2);
  c1->cd(1);
  mg1->Draw("apl");
  c1->cd(1)->BuildLegend();
  c1->cd(2);
  mg2->Draw("apl");
}
