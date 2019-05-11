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
void offlinefit(char* fitname,char* infile, char* parms,char* outfile)
{
  gROOT->ProcessLine(".L /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/mlhfitv2/mlhfitv2.C");
  gROOT->ProcessLine(Form("mlhfitv2(\"%s\",\"%s\",\"%s\",\"%s\")",fitname,infile,parms,outfile));
}
