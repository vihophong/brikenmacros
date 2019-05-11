#!/bin/bash                                                                   
root -b -q 'offlinefitall.C('$3','$2',"'$1'",'$4','$5','$6','$7')'

#example ./offlinefitall.sh Ag130 0 10 0.05 0.62 0.053 38232
# ./offlinefitall.sh ri ninteration rebin start range efficiency efficiency_error nimplant
