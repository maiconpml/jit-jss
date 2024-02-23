#!/usr/bin/env bash

#example (from root of fastPsp): executeTester.sh <path to inst dir> --- ./InstChars/executeTester.sh ./insts/PSSP/PspNasiriStyleAllTaiJsp10copies

exe=./mainTester

function run() {

	

	for inst in $1/*
	do
		:
		#echo $inst
		echo -n "`basename $inst .psi` "
		$exe $inst
	done

}

run $1
