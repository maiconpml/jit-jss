#!/usr/bin/env bash                                                                                   

######################################################################                                
# test runner script                                                                                  
######################################################################                                

#$EXE $INSTANCE_PATH $CANDIDATE --timeLog false --maxSecs 180.0 ${CAND_PARAMS}
#--tenure 9 --initialjumpLimit 1237 --bjSize 11 --acceptAlpha 0.4629 --perturbSize 22 --lowerBoundOrder true --perturbType mix --bubbleP 0.6012 --mixProb 0.2513

# basic configuration                                                                                 

#base=${HOME}/Home/psp/psp                                                                            
base=./
exe=${base}/main


function runTest() {
    instance_dir=$1
    lot=$2
    nlo=$3

    instances=(`ls -1 ${instance_dir}`)
    nin=${#instances[@]}

    ini=$(( (${lot}-1)*${nin}/${nlo}  ))
    fin=$(( (${lot})*${nin}/${nlo}  ))

	for ((i=$ini; i<$fin; i++)) {
		${exe} ${instance_dir}/${instances[i]} "`basename ${instances[i]} .psi` " --maxSecs 300.0  --seed 0 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 --timeLog false --lowerBoundOrder false
	}
}

runTest ${base}/insts/CompAnalysis $1 $2
