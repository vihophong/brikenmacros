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
void offlineimprate(char* input, char* ri, Int_t selectz)
{
    gROOT->ProcessLine(".L checkimprate.C");
      gROOT->ProcessLine(Form("checkimprate o(\"%s\",\"%s\");",input,ri));
      gROOT->ProcessLine(Form("o.Loop(%d);",selectz));

}
