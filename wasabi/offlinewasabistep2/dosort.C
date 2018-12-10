#include "stdio.h"
#include "string.h"
void dosort(char* input,char* bripsin, char* output) {
  gROOT->ProcessLine(".L aidaclass.h");
  gROOT->ProcessLine(".L wasabisorter2.C");
  gROOT->ProcessLine(Form("sorter2(\"%s\",\"%s\",\"%s\",%d);",input,bripsin,output,1));
}
