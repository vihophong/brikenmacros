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
void offlinefitall(char* thr, char* ri)
{
  gROOT->ProcessLine(".L fitDecayLL3hists_syserr.C");
  gROOT->ProcessLine(Form("fitDecayLL3hists_wrndcoin(\"%s_%s\",\"../../../brikenmacros-rootfiles/fitting/fithistslowin/histlowin%s_%s.root\",\"parms/parms%s.txt\",1000,\"fitresults/%s_%skeV.root\")",ri,thr,ri,thr,ri,ri,thr));
  //cout<<Form("fitDecayLL3hists_wrndcoin(\"%s_%s\",\"../../../brikenmacros-rootfiles/fitting/fithistslowin/histlowin%s_%s.root\",\"parms/parms%slowin.txt\",0,\"fitresults/%s_%skeV.root\")",ri,thr,ri,thr,ri,ri,thr)<<endl;
}
