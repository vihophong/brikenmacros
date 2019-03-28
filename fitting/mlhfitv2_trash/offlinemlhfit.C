#include "stdio.h"
#include "string.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TTree.h"
void offlinemlhfit(char* fitname,char* infile,char* parmsfile,char* outfile,Int_t opt)
{
    

    gROOT->ProcessLine(".L mlhfitv2.C");
    
    gROOT->ProcessLine(Form("mlhfitv2(\"%s\",\"%s\",\"%s\",\"%s\",%d)",fitname,infile,parmsfile,outfile,opt));
        
}
