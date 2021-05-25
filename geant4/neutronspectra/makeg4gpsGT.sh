#!/bin/bash
#Ag isotopes
root -b -q 'makeg4gpsGT.C("neuspecs/minatospec/z47/129/nspec.out","neuspecs/minatospec/Ag129.mac",4)'
root -b -q 'makeg4gpsGT.C("neuspecs/minatospec/z47/130/nspec.out","neuspecs/minatospec/Ag130.mac",4)'

#Cd isotopes
root -b -q 'makeg4gpsGT.C("neuspecs/minatospec/z48/131/nspec.out","neuspecs/minatospec/Cd131.mac",4)'
root -b -q 'makeg4gpsGT.C("neuspecs/minatospec/z48/132/nspec.out","neuspecs/minatospec/Cd132.mac",4)'

#In isotopes
root -b -q 'makeg4gpsGT.C("neuspecs/minatospec/z49/133/nspec.out","neuspecs/minatospec/In133.mac",4)'
root -b -q 'makeg4gpsGT.C("neuspecs/minatospec/z49/134/nspec.out","neuspecs/minatospec/In134.mac",4)'
root -b -q 'makeg4gpsGT.C("neuspecs/minatospec/z49/135/nspec.out","neuspecs/minatospec/In135.mac",4)'
