/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#pragma once

//#define NDEBUG

using namespace std;

#include <assert.h>
#include<vector>
#include<string>

#define MAX_ELEMENTS 100

//#define PARALLEL_MACHS
//#define EXTRA_DECODE


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

vector<string> resultList;


//#define HALF_MAKES

// #define PRINT_SCHEDULE "./schedules.txt"
