#!/usr/bin/env bash                                                                                   

######################################################################                                
# test runner script                                                                                  
######################################################################

#example: execBacthFullParallelJspAdjTime.sh 2 
#$1: k  expects in [2,5]
#example of exe use: ./main ./insts/all/ft10_3.psi ft10_3 --searchType tabuSearch
#example output: 1221 0.25 0.332


baseDir=./
exe=${baseDir}main
instListDir=${baseDir}Runners/JSP/
instBaseDir=${baseDir}/parallelInst/AllJsp/k1/

function runTest() {
    for ((index=4; index<11; index++)) {
		instList=Table${index}.dat
		mapfile -t allInsts < ${instListDir}${instList}
		for instInfo in "${allInsts[@]}"
		do
			:
			instName=`echo ${instInfo} | cut -d " " -f 1`
			maxSecs=`echo ${instInfo} | cut -d " " -f 2`
			#maxMilliSecs=`echo "${maxSecs} * 1000" | bc -l`
			smallNameS="`basename $instName .psi` "
			smallName="`basename $instName .psi`"

			echo -n $smallNameS

			${exe} ${instBaseDir}${smallName}_k1.ppsi ${smallName} --maxSecs ${maxSecs} --searchType tabuSearch --timeLog false
		done
	}
}

runTest

#base=./
#exe=${base}/mainParallel
#instDir=$base/parallelInst/AllJsp

#function runTest() {


#	FILE=$instDir/k$1/*

#	for f in $FILE
#	do
#		#"`basename ${instances[i]} .psi` "
#		echo -n "`basename $f .ppsi` "
#		$exe $f
#	done
#}

#runTest $1
