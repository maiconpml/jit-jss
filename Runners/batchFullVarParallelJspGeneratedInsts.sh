#!/usr/bin/env bash                                                                                   

######################################################################                                
# test runner script                                                                                  
######################################################################

#example of exe use: ./mainParallel ./parallelInst/AllJsp/k2/abz5Jsp_k2.ppsi 
#example output: 1221 0.25 0.332
#$1 is k expects in [2,5]


base=./
exe=${base}/mainParallel
instDir=$base/GenParallelInsts/GeneratedInsts

function runTest() {

	echo Inst	Makes	LB	Dev	Time Evals

	FILE=$instDir/*

	for f in $FILE
	do
		#"`basename ${instances[i]} .psi` "
		echo -n "`basename $f .ppsi` " 
		$exe $f --tenure 8 --iters 2905 --iterDecrease 346 --bjSize 6 --learn true --trimmSize 0.2 --minLearnKW 0.6 --minLearnPeriod 30
		#--tenure 8 --iters 2905 --iterDecrease 346 --bjSize 6
		#--learn true --trimmSize 0.2 --minLearnKW 0.6 --minLearnPeriod 30
		#--startType jobTail
		#--allDelayNeigh false --critVanillaNeigh true
		#--dispatch true
	done
}

runTest $1
