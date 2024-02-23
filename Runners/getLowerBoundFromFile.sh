#!/usr/bin/env bash                                                                                   

######################################################################                                
# test runner script                                                                                  
######################################################################

#usage: ./Runners/getLowerBoundFromFile.sh ./insts/MssFromPss/

function run() {


	FILE=$1/*

	for f in $FILE
	do
		#"`basename ${instances[i]} .psi` "
		##echo -n "`basename $f .psi` "
		#echo ./main $f  xxx --onlyLowerBound true
		./main $f  "`basename $f .psi` " --onlyLowerBound true 
	done
}

run $1
