/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#pragma once

//#define NDEBUG
using namespace std;

#include <assert.h>
#include<vector>
#include<string>

//#define COUNT_STEPS
#ifdef COUNT_STEPS
unsigned numSteps;
unsigned shownSeconds;
unsigned currentShownSeconds;
#endif

#define UINFTY UINT_MAX
#define ANALYZE_MODE

//#define RFILE													//print schedule in format wich is used to create gantt chart
//#define PRINT_ONLY_RESULT							//print only the penalties result
//#define PRINT_SCHEDULE								//print a detailed representaton of the schedule
//#define PRINT_ONLY_IS									//print only the initial solution
//#define PRINT_NEIGHBOURS_NB							//print the mean of the number of neighbours
//#define PRINT_DEFAULT										//print the result in the format: instance-path makes-lower-bound obtained-makes penalties earliness-penalties tardiness-penalties neighbours-number(if exists)
#define PRINT_VERIFY_OUTPUT
#define SHIFT_OPERS											

//only one of this 3
#define VANILLA_PSS

// ALIASES
#define UNDEFINED 11

#define FORTH 12
#define BACK 13

#define JOB 14
#define MACH 15

// CONSTS
#define MAX_MACHS 100
//#define MAX_JOBS 200
#define MAX_OPERS 2000


//FLOW OPTIONS
#define ITALIANO_WITH_COUNTER

//#define USE_DELTA_SWAP

#define INSA_START  1001
#define RAND_START 1002
#define GT 1003

#define LOCAL_SEARCH_INTENS 2001 //TODO
#define TABU_INTENS 2002

#define INSA_PERTURB 3001
#define BUBBLE_PERTURB 3002
#define MIX_PERTURB 3003
//#define CARLIER_PERTURB 3004 //TODO

#define METROPOLIS 4001
#define BEST_JUMP 4002

vector<string> resultList;

#ifdef PRINT_NEIGHBOURS_NB
	vector<double> neigh;
#endif // PRINT_NEIGHBOURS_NB