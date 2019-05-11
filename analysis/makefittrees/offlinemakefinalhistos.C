#include "stdio.h"
#include "string.h"
#include "TROOT.h"
#include "TString.h"
#include "TH1.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
void offlinemakefinalhistos(char* input,char* riname,char* output,Double_t tmin, Double_t tmax,Int_t layer=-1)
{
  gROOT->ProcessLine(".L makefinalhistos.C");
  gROOT->ProcessLine(Form("makefinalhistos o(\"%s\",\"%s\");",input,riname));
  gROOT->ProcessLine(Form("o.MakeFinalHisto(\"%s\",%f,%f,%d);",output,tmin,tmax,layer));
}
