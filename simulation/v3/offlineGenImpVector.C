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
void offlineGenImpVector(char* infile, char* riname, char* outfile,Int_t layer)
{
    gROOT->ProcessLine(".L genImpVector.C");
    gROOT->ProcessLine(Form("genImpVector o(\"%s\",\"%s\");",infile,riname));
    gROOT->ProcessLine(Form("o.Loop(%d,\"%s\");",layer,outfile));
}
