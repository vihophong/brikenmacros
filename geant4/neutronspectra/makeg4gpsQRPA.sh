#!/bin/bash
#Ag isotopes
root -b -q 'makeg4gpsQRPA.C("neuspecs/qrpahfspec/original/bspec_47_129.dat","neuspecs/qrpahfspec/Ag129.mac",1)'
root -b -q 'makeg4gpsQRPA.C("neuspecs/qrpahfspec/original/bspec_47_130.dat","neuspecs/qrpahfspec/Ag130.mac",1)'

#Cd isotopes
root -b -q 'makeg4gpsQRPA.C("neuspecs/qrpahfspec/original/bspec_48_131.dat","neuspecs/qrpahfspec/Cd131.mac",1)'
root -b -q 'makeg4gpsQRPA.C("neuspecs/qrpahfspec/original/bspec_48_132.dat","neuspecs/qrpahfspec/Cd132.mac",1)'

#In isotopes
root -b -q 'makeg4gpsQRPA.C("neuspecs/qrpahfspec/original/bspec_49_131.dat","neuspecs/qrpahfspec/In131.mac",1)'
root -b -q 'makeg4gpsQRPA.C("neuspecs/qrpahfspec/original/bspec_49_132.dat","neuspecs/qrpahfspec/In132.mac",1)'
root -b -q 'makeg4gpsQRPA.C("neuspecs/qrpahfspec/original/bspec_49_133.dat","neuspecs/qrpahfspec/In133.mac",1)'
root -b -q 'makeg4gpsQRPA.C("neuspecs/qrpahfspec/original/bspec_49_134.dat","neuspecs/qrpahfspec/In134.mac",1)'
root -b -q 'makeg4gpsQRPA.C("neuspecs/qrpahfspec/original/bspec_49_135.dat","neuspecs/qrpahfspec/In135.mac",1)'
root -b -q 'makeg4gpsQRPA.C("neuspecs/qrpahfspec/original/bspec_49_136.dat","neuspecs/qrpahfspec/In136.mac",1)'

#Sn isotopes
root -b -q 'makeg4gpsQRPA.C("neuspecs/qrpahfspec/original/bspec_50_136.dat","neuspecs/qrpahfspec/Sn136.mac",1)'
root -b -q 'makeg4gpsQRPA.C("neuspecs/qrpahfspec/original/bspec_50_137.dat","neuspecs/qrpahfspec/Sn137.mac",1)'
