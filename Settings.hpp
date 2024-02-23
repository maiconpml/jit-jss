/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#pragma once

#define NDEBUG

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

//#define REOPT_GAPS


//#define PARALLEL_MACHS
//#define VAR_PARALLEL_MACHS
//#define EXTRA_DECODE

//only one of this 3
#define VANILLA_PSS
//#define EXTENDED_JSS
//#define EXTENDED_PSS

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

#define LOCAL_SEARCH_INTENS 2001 //TODO
#define TABU_INTENS 2002

#define INSA_PERTURB 3001
#define BUBBLE_PERTURB 3002
#define MIX_PERTURB 3003
//#define CARLIER_PERTURB 3004 //TODO

#define METROPOLIS 4001
#define BEST_JUMP 4002

//Parallel neighbouhood options
#define CRIT_VANILLA_PARALLEL_NEIGH 5001
#define ALL_DELAY_PARALLEL_NEIGH 5002

#define SMALLER_GAP 6001
#define LARGER_GAP 6002

#define PARALLEL_INSA_START 7001
#define PARALLEL_JOB_TAIL_START 7002


vector<string> resultList;


//#define HALF_MAKES

// #define PRINT_SCHEDULE "./schedules.txt"

double sumKW = 0.0;
unsigned iterKW = 0;
double at30KW = 0.0;
