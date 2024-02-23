#!/usr/bin/env bash                                                                                   

######################################################################                                
# test runner script                                                                                  
######################################################################

#example of exe use: ./mainParallel ./parallelInst/AllJsp/k2/abz5Jsp_k2.ppsi 
#example output: 1221 0.25 0.332
#$1 is k expects in [2,5]


base=./
exe=${base}/mainParallel
instDir=$base/parallelInst/AllJsp

function runTest() {


	FILE=$instDir/k$1/*

	for f in $FILE
	do
		#"`basename ${instances[i]} .psi` "
		echo -n "`basename $f .ppsi` "
		$exe $f
	done
}

runTest $1
