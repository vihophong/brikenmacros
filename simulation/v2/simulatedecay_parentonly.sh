#!/bin/bash                                                                  

exedir=/home/phong/briken17/brikenmacrosupdate/brikenmacros/simulation/v2

echo "usage: ./simulatedecay.sh decayparms expparms outputfiles"           

if [ $# -eq 3 ]; 
then
    root -b -q ''$exedir'/simulatedecay_parentonly.C("'$1'","'$2'","'$3'.tree")'
    root -b -q ''$exedir'/checksimulation.C("'$3'.tree","'$3'")'
fi