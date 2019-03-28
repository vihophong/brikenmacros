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
void offlinecheckruntimeclover(char* input)
{
    gROOT->ProcessLine(".L checkruntimeclover.C");

      gROOT->ProcessLine(Form("checkruntimeclover o(\"%s\");",input));
      gROOT->ProcessLine("o.Loop();");

}
