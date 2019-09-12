#!/bin/bash                                                                   
root -b -q 'fitterUnbin.C("'$1'","'$2'","'$3'",'$4','$5')'

#void fitterBin(char* infilename,char* parmsfilename,char* outfilename,Double_t irejectrange=0.075, Double_t ineueff=0.62, Int_t ninterations=0,Int_t rebin=1)

