#!/bin/bash
FILES=decaynew_lowin/*
for f in $FILES
do
	echo "Processing $f file.."
	root -b -q 'checkdeadtime.C("'$f'")'
done
