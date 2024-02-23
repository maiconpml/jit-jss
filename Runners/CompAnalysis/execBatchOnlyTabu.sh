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

#old: --tenure 5 --initialjumpLimit 8343 --jumpLimitDecrease 711 --bjSize 4 --maxD 168 --maxC 4 --acceptAlpha 0.2755 --sizePInsa 27 --sizePBubble 3 --perturbType mix --bubbleP 0.6217 --mixProb 0.2144
#new: --tenure 10 --initialjumpLimit 2763 --jumpLimitDecrease 401 --bjSize 4 --maxD 85 --maxC 2 --acceptAlpha 0.0624 --sizePInsa 27 --sizePBubble 3 --perturbType mix --bubbleP 0.6217 --mixProb 0.6399
#after light: ${exe} ${instance_dir}/${instances[i]} "`basename ${instances[i]} .psi` " --timeLog true --maxSecs 600.0 --scheduleFile schedules$2_$3.txt --tenure 10 --initialjumpLimit 2763 --jumpLimitDecrease 401 --bjSize 4 --maxD 15 --maxC 3 --acceptAlpha 0.0624 --sizePInsa 27 --sizePBubble 3 --lowerBoundOrder true --perturbType mix --bubbleP 0.6217 --mixProb 0.6399 --seed 0
#evolve: 1  --tenure 7 --initialjumpLimit 1618 --decreaseDivisor 6 --increaseIters 700 --bjSize 3 --maxD 90 --acceptAlpha 0.6829 --sizePInsa 18 --sizePBubble 6 --perturbType mix --bubbleP 0.5468 --mixProb 0.3359
#old:       --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 

function runTest() {
    instance_dir=$1
    lot=$2
    nlo=$3

    instances=(`ls -1 ${instance_dir}`)
    nin=${#instances[@]}

    ini=$(( (${lot}-1)*${nin}/${nlo}  ))
    fin=$(( (${lot})*${nin}/${nlo}  ))

	for ((i=$ini; i<$fin; i++)) {
		${exe} ${instance_dir}/${instances[i]} "`basename ${instances[i]} .psi` " --maxSecs 300.0 --scheduleFile schedules$2_$3.txt --seed 0 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 --searchType tabuSearch

		${exe} ${instance_dir}/${instances[i]} "`basename ${instances[i]} .psi` " --maxSecs 300.0 --scheduleFile schedules$2_$3.txt --seed 1111 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 --searchType tabuSearch

		${exe} ${instance_dir}/${instances[i]} "`basename ${instances[i]} .psi` " --maxSecs 300.0 --scheduleFile schedules$2_$3.txt --seed 2222 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 --searchType tabuSearch

		${exe} ${instance_dir}/${instances[i]} "`basename ${instances[i]} .psi` " --maxSecs 300.0 --scheduleFile schedules$2_$3.txt --seed 3333 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 --searchType tabuSearch

		${exe} ${instance_dir}/${instances[i]} "`basename ${instances[i]} .psi` " --maxSecs 300.0 --scheduleFile schedules$2_$3.txt --seed 4444 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 --searchType tabuSearch

		${exe} ${instance_dir}/${instances[i]} "`basename ${instances[i]} .psi` " --maxSecs 300.0 --scheduleFile schedules$2_$3.txt --seed 5555 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 --searchType tabuSearch

		${exe} ${instance_dir}/${instances[i]} "`basename ${instances[i]} .psi` " --maxSecs 300.0 --scheduleFile schedules$2_$3.txt --seed 6666 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 --searchType tabuSearch

		${exe} ${instance_dir}/${instances[i]} "`basename ${instances[i]} .psi` " --maxSecs 300.0 --scheduleFile schedules$2_$3.txt --seed 7777 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 --searchType tabuSearch

		${exe} ${instance_dir}/${instances[i]} "`basename ${instances[i]} .psi` " --maxSecs 300.0 --scheduleFile schedules$2_$3.txt --seed 8888 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 --searchType tabuSearch

		${exe} ${instance_dir}/${instances[i]} "`basename ${instances[i]} .psi` " --maxSecs 300.0 --scheduleFile schedules$2_$3.txt --seed 9999 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 --searchType tabuSearch
	}
}

runTest ${base}/insts/CompAnalysis $1 $2
