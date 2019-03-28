#!/bin/bash                                                                                                                                                             
root -b -q 'simulatedecay.C("'$1'","'$2'")'
root -b -q 'checksimulation.C("'$3'")'

