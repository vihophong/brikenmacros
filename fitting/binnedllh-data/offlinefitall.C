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
void offlinefitall(Int_t rebin, Int_t niter, char* ri,Double_t neueff=0.62,Double_t neueff_err=0.053)
{
  gROOT->ProcessLine(".L fitDecayLL3hists_syserr.C");
  gROOT->ProcessLine(Form("fitDecayLL3hists_wrndcoin(\"%s\",\"../../../../results/fithistograms/%slowin.root\",\"parms/parms%s.txt\",%d,\"fitresults/fitresult%slowin.root\",%d,%f,%f)",ri,ri,ri,niter,ri,rebin,neueff,neueff_err));
  //cout<<Form("fitDecayLL3hists_wrndcoin(\"%s_%s\",\"../../../brikenmacros-rootfiles/fitting/fithistslowin/histlowin%s_%s.root\",\"parms/parms%slowin.txt\",0,\"fitresults/%s_%skeV.root\")",ri,thr,ri,thr,ri,ri,thr)<<endl;
}
