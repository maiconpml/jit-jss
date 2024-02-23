#!/usr/bin/env bash                                                                                   

######################################################################                                
# test runner script                                                                                  
######################################################################

#example of c exec: ./main JspToParallel ../insts/JSP/abz5Jsp.psi ./test.txt --parallelK 2
#example of bash exec: ./execBacthJspToParallel.sh ../insts/DemirkolJsp/
#$1 is dir of source insts

exe=./main
targetDir=../parallelInst/AllJsp/

function runInstTransf() {

	#instances=(`ls -1 $1`)
	FILE=$1/*

	for f in $FILE
	do
		#for ((k=2; k<=5; k++)) {
		for ((k=1; k<=1; k++)) {
				echo "Processing $f"
				#echo JspToParallel $f $targetDir/"`basename ${f} .psi`"_k$k.ppsi --parallelK $k
				$exe JspToParallel $f $targetDir/k$k/"`basename ${f} .psi`"_k$k.ppsi --parallelK $k
		}
		
	done

}

runInstTransf $1
