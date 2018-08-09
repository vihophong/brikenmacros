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
void offlinemakehitdistr(char* input,char* riname,char* output,Double_t tmin, Double_t tmax)
{
  gROOT->ProcessLine(".L makemlhtree.C");
  gROOT->ProcessLine(Form("makemlhtree o(\"%s\",\"%s\");",input,riname));
  gROOT->ProcessLine(Form("o.PlotNeuHitPattern(\"%s\",%f,%f);",output,tmin,tmax));
}
