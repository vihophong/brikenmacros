#include "TSpectrum.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"

#include "TCut.h"
#include "TCutG.h"
#include "TLatex.h"
#include "TSpectrum.h"
#include "TF1.h"
#include <algorithm>

void makespec(Int_t chmin,Int_t chmax)
{
    TChain ch("group1");
    ch.Add("preparation/prerun00042.root");
    ch.Add("preparation/prerun00043.root");
    ch.Add("preparation/prerun00044.root");
    TH1F *h[100];

    for (Int_t ich=chmin;ich<=chmax;ich++){
        ch.Draw(Form("dgtz_e[%d]>>h%d(500,10,2010)",ich,ich),"","goff");
        h[ich]=(TH1F*) gDirectory->Get(Form("h%d",ich));
    }

}
