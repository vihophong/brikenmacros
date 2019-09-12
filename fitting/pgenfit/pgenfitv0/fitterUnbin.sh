#!/bin/bash                                                                   
root -b -q 'fitterUnbin.C("'$1'","'$2'","'$3'",'$4','$5','$6')'

#void fitterUnbin(char* infile, char* parmsfile,char* outfilename,Double_t deadtimeinput=0.08, Double_t neueffinput=0.62,Int_t rebin=1)


