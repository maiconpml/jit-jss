#!/usr/bin/env bash                                                                                   

######################################################################                                
# test runner script                                                                                  
######################################################################

#usage: ./Runners/execPssToMssGss.sh ./insts/MssFromPss/ 1 5
# $1 inst dir  ----  batch $2 of a total $3



function run() {
	
	instance_dir=$1
    lot=$2
    nlo=$3

    instances=(`ls -1 ${instance_dir}`)
    nin=${#instances[@]}
	
    ini=$(( (${lot}-1)*${nin}/${nlo}  ))
    fin=$(( (${lot})*${nin}/${nlo}  ))

	for ((i=$ini; i<$fin; i++)) {
		for ((index=1; index<11; index++)) {
			./main ${instance_dir}/${instances[i]} "`basename ${instances[i]} .psi` " --maxSecs 300.0 --scheduleFile schedules$2_$3.txt --seed ${index} --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33
		}
	}

}

run $1 $2 $3





#function run() {


#	FILE=$1/*

#	for f in $FILE
#	do
		#./main $f  "`basename $f .psi` " --onlyLowerBound true 
		#./main $f  "`basename $f .psi` " --maxSecs 300.0 --scheduleFile schedules$2_$3.txt --seed 0 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 
	#done
#}

#run $1
