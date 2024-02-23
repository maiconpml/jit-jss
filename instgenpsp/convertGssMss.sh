#!/usr/bin/env bash

exe=./main
#mod=MapPssToMss
mod=MapPssToGss
instDir=../insts/PSSP/PspNasiriStyleAllTaiJsp10copies/

FILE=$instDir/*.psi

for f in $FILE
do
	echo -n "`basename $f .ppsi` "
	#$exe $mod $f ./MssFromPss/"`basename $f _toMSS.ppsi` "
	$exe $mod $f ./GssFromPss/"`basename $f _toGSS.ppsi` "
done
