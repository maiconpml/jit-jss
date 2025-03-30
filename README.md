# Just in Time Job Shop with Tabu Search

## Compile

### Debug build

```bashsession
user~$: mkdir Debug && cd Debug
user~$: cmake DCMAKE_BUILD=Debug ..
user~$: make
```

### Release build

```bashsession
user~$: mkdir Release && cd Release
user~$: cmake DCMAKE_BUILD=Release ..
user~$: make
```
## Run

```bashsession
./main [INSTANCE_PATH] [LOG_NAME] [Program Options]
```
### Program options

help - show help\
instPath - instance path\
name - name for log\
maxSecs - maximum time in seconds (default: 300)\
seed - random number generator seed (default: 13)\
tenure - tabu tenure (default: 8)\
initialjumpLimit - initial max iterations without improvement before backjump (default: 2500)\
decreaseDivisor - decrease in jumLimit each time jumps without improve (default: 7)\
bjSize - back jump list size (default: 20)\
maxD - cycle size to stop search and force backjump (default: 100)\
maxC -  number of repeats of cycle to force backjump (default: 2)\
timeLog -  Show time of each new best sollution found (default: true)\
scaleTime -  Change received time, for testing. (JSP) (default: 1.0)\
onlyMakesLowerBound -  To get lower bound values (default: false)\
schedulerType -  Type of scheduler used: 1 - early as possible (fastest), 2 - delaying when possible (better schedule but slower), 3 - 1 and 2 mixed (1 when tard penalties are dominating and 2 when earl penalties are dominating, 4 - cplex (slowest but optimal) (default: 1)\
useSwapAllNIter - Number of iterations without improvement to use neighborhood Swap All (default: 1)\
useCplexRelaxRatio - Ratio between earliness penalty and tardiness penalty to use neighborhood Cplex Relax (default: 0.25)\

## Tabu Search

The current tabu search uses 3 neighborhoods.\
Swap all - all possible swaps between two adjacent operations in the same machine\
Swap Critical oper - swaps between two firts and lasts operations of each block inside each late operation critical path. A critical path of a late operation O is defined by all operations that prevents O of being scheduled earlier.\
Cplex relax- cplex relax machine precedence of a block of earl operations to find the optimal sequence of the operations in that block.

Swap all is used when the search don't find a better solution in an x number of iterations\
Cplex relax is used when the earliness penalty dominates the tardiness operation by a factor x\
Swap critical oper is used by default, that is, when the other two are not used

## Evaluation of a sequence

The order of operations in each machine is called by sequence and a better sequence is found by the Tabu Search. To get the objective value of a sequence is necessary to set the start times of each operation. There are three ways to do it.

**Schedule as early as possible** - As the name sugests, given a sequence schedule each operation as early as possible. This scheduler is good for solutions where tardiness penalties dominate. If the flag SHIFT_OPERATIONS is defined, a algorithm to delay operations without incresing penalties are executed at the end of each schedule to minimize earliness penalties.

**Schedule delaying when advantageous** - Schedule operations starting of last to first. For each operation schedule at 0. If delaying this operation reduces penalties delay it until its deadline or the start of next operation in job or mach. If delaying the operation and the other operations next to it reduces penalties delay until some operation reach its deadline or the start of other operation. Repeat this step until delaying the operations don't reduces the penalties.

**Schedule with cplex** - Given a sequence cplex find the optimal schedule but is slower than the previous methods.

It's possible to choose wich method to use. The hybrid method uses the two first together. When the current state has more tardiness penalties use ScheduleAsEarly and when state has more earliness penalties use ScheduleDelaying.

## Current results

The current best results are obtained by the following configuration:\
Tabu tenure: 29\
Scheduler type: Hybrid\
Iterations to use swap all: 1\
Factor to use cplex relax: 1/4

```bashsession
  ./main instance name --tenure 29 --maxSecs 60 --schedulerType 3
```

