#!/bin/bash                                                                   
root -b -q 'fitUnbin.C("'$1'","'$2'","'$3'",'$4','$5','$6')'

#void fitUnbin(char* inputfile,char* parmsfile,char*outfile,Double_t inputneueff =0.62,Double_t inputneueff_err=0.01,Double_t beginfit=0.08,Int_t fitopt=0)// fit option 0: no constrains(central value) fit option 1: constrains for sys error estimation; fit option 2: constrains for background only

#./fitUnbin.sh fittrees/Cd130.root Cd130full.txt Cd130fullnoisomer.root 0.6441523116 0.032151846784528 0.05 0
#./fitUnbin.sh fittrees/Cd130.root Cd130full_isomer.txt Cd130fullisomer.root 0.6441523116 0.032151846784528 0.05 0
#./fitUnbin.sh outhist.root In134full.txt test.root 0.62 0.01 0

