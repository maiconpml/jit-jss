#!/usr/bin/env bash                                                                                   
######################################################################                                
# test runner script                                                                                  
######################################################################               


#EXAMPLE: ./Runners/runTimedJspActualTime.sh Table4.dat 


exeDir=./
#exe=${exeDir}mainTimedIg
exe=${exeDir}main
instListDir=${exeDir}Runners/JSP/
instDir=${exeDir}insts/all/

#		${exe} ${instDir}${instName} ${instName} --timeLog true --maxSecs ${maxSecs} --scheduleFile schedules.txt --tenure 9 --initialjumpLimit 1237 --bjSize 11 --acceptAlpha 0.4629 --perturbSize 22 --lowerBoundOrder true --perturbType mix --bubbleP 0.6012 --mixProb 0.2513 --seed 0

function foo() {
	instList=$1
	mapfile -t allInsts < ${instListDir}${instList}
	for instInfo in "${allInsts[@]}"
	do
		:
		instName=`echo ${instInfo} | cut -d " " -f 1`
		maxSecs=`echo ${instInfo} | cut -d " " -f 2`

		maxSecs=`echo "${maxSecs} * 5" | bc -l`


		
		${exe} ${instDir}${instName} ${instName} --timeLog true --maxSecs ${maxSecs} --scheduleFile schedules$1.txt --seed 0 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 
		${exe} ${instDir}${instName} ${instName} --timeLog true --maxSecs ${maxSecs} --scheduleFile schedules$1.txt --seed 1111 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 
		${exe} ${instDir}${instName} ${instName} --timeLog true --maxSecs ${maxSecs} --scheduleFile schedules$1.txt --seed 2222 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 
		${exe} ${instDir}${instName} ${instName} --timeLog true --maxSecs ${maxSecs} --scheduleFile schedules$1.txt --seed 3333 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 
		${exe} ${instDir}${instName} ${instName} --timeLog true --maxSecs ${maxSecs} --scheduleFile schedules$1.txt --seed 4444 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 
		${exe} ${instDir}${instName} ${instName} --timeLog true --maxSecs ${maxSecs} --scheduleFile schedules$1.txt --seed 5555 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 
		${exe} ${instDir}${instName} ${instName} --timeLog true --maxSecs ${maxSecs} --scheduleFile schedules$1.txt --seed 6666 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 
		${exe} ${instDir}${instName} ${instName} --timeLog true --maxSecs ${maxSecs} --scheduleFile schedules$1.txt --seed 7777 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 
		${exe} ${instDir}${instName} ${instName} --timeLog true --maxSecs ${maxSecs} --scheduleFile schedules$1.txt --seed 8888 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 
		${exe} ${instDir}${instName} ${instName} --timeLog true --maxSecs ${maxSecs} --scheduleFile schedules$1.txt --seed 9999 --maxC 2 --tenure 6 --initialjumpLimit 2809 --decreaseDivisor 6 --bjSize 4 --maxD 95 --acceptAlpha 0.2169 --perturbType insa --sizePInsa 33 

	done
}

foo $1
