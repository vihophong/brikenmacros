#include "stdio.h"
#include "string.h"
void offLine(Int_t dssd,Int_t part) {
  //TString filenamein=Form("/home/ur16/rawdata/bigrips/run%d.ridf",infile);
  //TString filenameout1=Form("/home/ur16/treefiles/bigrips/bigrips_run%d.root",outfile);
  gROOT->ProcessLine(".L fitchannel.C");
  gROOT->ProcessLine("fitchannel ee;");
  if (part==1)
    gROOT->ProcessLine(Form("ee.Loop(%d,0,3,\"out1dssd%d\",0,5461)",dssd,dssd));
  else if (part==2)
    gROOT->ProcessLine(Form("ee.Loop(%d,0,3,\"out2dssd%d\",5461,10922)",dssd,dssd));
  else if (part==3)
    gROOT->ProcessLine(Form("ee.Loop(%d,0,3,\"out3dssd%d\",10922,16384)",dssd,dssd));
  else cout<<"wrong selection"<<endl;
}
