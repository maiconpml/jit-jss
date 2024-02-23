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

	echo Inst Beta W Makes LB Dev Time Evals

	FILE=$instDir/*

	for f in $FILE
	do
		echo -n "`basename $f .ppsi` " 0.0 0.6 " "
		$exe $f --tenure 8 --iters 2905 --iterDecrease 346 --bjSize 6 --learn true --trimmSize 0.0 --minLearnKW 0.6 --minLearnPeriod 30

		echo -n "`basename $f .ppsi` " 0.1 0.6 " "
		$exe $f --tenure 8 --iters 2905 --iterDecrease 346 --bjSize 6 --learn true --trimmSize 0.1 --minLearnKW 0.6 --minLearnPeriod 30
		
		echo -n "`basename $f .ppsi` " 0.2 0.6 " "
		$exe $f --tenure 8 --iters 2905 --iterDecrease 346 --bjSize 6 --learn true --trimmSize 0.2 --minLearnKW 0.6 --minLearnPeriod 30
		
		echo -n "`basename $f .ppsi` " 0.3 0.6 " "
		$exe $f --tenure 8 --iters 2905 --iterDecrease 346 --bjSize 6 --learn true --trimmSize 0.3 --minLearnKW 0.6 --minLearnPeriod 30
		
		echo -n "`basename $f .ppsi` " 0.4 0.6 " "
		$exe $f --tenure 8 --iters 2905 --iterDecrease 346 --bjSize 6 --learn true --trimmSize 0.4 --minLearnKW 0.6 --minLearnPeriod 30
		
		echo -n "`basename $f .ppsi` " 0.5 0.6 " "
		$exe $f --tenure 8 --iters 2905 --iterDecrease 346 --bjSize 6 --learn true --trimmSize 0.5 --minLearnKW 0.6 --minLearnPeriod 30
		
		echo -n "`basename $f .ppsi` " 0.6 0.6 " "
		$exe $f --tenure 8 --iters 2905 --iterDecrease 346 --bjSize 6 --learn true --trimmSize 0.6 --minLearnKW 0.6 --minLearnPeriod 30
		
		echo -n "`basename $f .ppsi` " 0.7 0.6 " "
		$exe $f --tenure 8 --iters 2905 --iterDecrease 346 --bjSize 6 --learn true --trimmSize 0.7 --minLearnKW 0.6 --minLearnPeriod 30
		
		echo -n "`basename $f .ppsi` " 0.8 0.6 " "
		$exe $f --tenure 8 --iters 2905 --iterDecrease 346 --bjSize 6 --learn true --trimmSize 0.8 --minLearnKW 0.6 --minLearnPeriod 30
		
		echo -n "`basename $f .ppsi` " 0.9 0.6 " "
		$exe $f --tenure 8 --iters 2905 --iterDecrease 346 --bjSize 6 --learn true --trimmSize 0.9 --minLearnKW 0.6 --minLearnPeriod 30
		
		echo -n "`basename $f .ppsi` " 1.0 0.6 " "
		$exe $f --tenure 8 --iters 2905 --iterDecrease 346 --bjSize 6 --learn true --trimmSize 0.9 --minLearnKW 0.6 --minLearnPeriod 30
	done
}

runTest $1
