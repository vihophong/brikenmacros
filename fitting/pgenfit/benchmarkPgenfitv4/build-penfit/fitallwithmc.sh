#!/bin/bash

#In isotoptes (low to high statistics)
./runsinglefit.sh parms/In136parms.txt lowin/In136.root outFit/In136out.root 0.05 1000
./runsinglefit.sh parms/In132parms.txt lowin/In132.root outFit/In132out.root 0.08 1000
./runsinglefit.sh parms/In131parms.txt lowin/In131.root outFit/In131out.root 0.08 500
./runsinglefit.sh parms/In135parms.txt lowin/In135.root outFit/In135out.root 0.08 500
./runsinglefit.sh parms/In134parms.txt lowin/In134.root outFit/In134out.root 0.08 300
./runsinglefit.sh parms/In133parms.txt lowin/In133.root outFit/In133out.root 0.08 100
