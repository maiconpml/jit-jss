/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#pragma once

//#define NDEBUG
using namespace std;

#include <assert.h>
#include <climits>
#include <vector>
#include <string>

#define UINFTY UINT_MAX

// select only onr of these 5
//#define RFILE													//print schedule in the format that is used to create gantt chart
//#define PRINT_ONLY_RESULT							//print only the penalties result
//#define PRINT_SCHEDULE								//print a detailed representaton of the schedule
//#define PRINT_ONLY_IS									//print only the initial solution
#define PRINT_DEFAULT										//print the result in the format: instance-path makes-lower-bound obtained-makes time-spent(ms) penalties earliness-penalties tardiness-penalties neighbors-number(if defined)

//#define PRINT_NEIGHBOURS_NB						//print the mean of the number of neighbors concatenated with PRINT_DEFAULT; should be defined only with PRINT_DEFAULT

//#define PRINT_VERIFY_OUTPUT						//print solution in format that is used for a external verifier

#define SHIFT_OPERS											//use shift operations on solution evaluation 

#define VANILLA_PSS

// CONSTS
#define MAX_MACHS 20
//#define MAX_JOBS 200
//#define MAX_OPERS 300

//FLOW OPTIONS
#define ITALIANO_WITH_COUNTER

extern vector<string> resultList;

#ifdef PRINT_NEIGHBOURS_NB
extern vector<double> neigh;
#endif // PRINT_NEIGHBOURS_NB