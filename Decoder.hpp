/*
  @author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/
#pragma once

#include <boost/multi_array.hpp>
#include <climits>
#include <algorithm> 
#include <chrono>

#include "Settings.hpp"
#include "Utilities.hpp"
#include "Instance.hpp"
#include "RandomNumberGenerator.hpp"

using namespace boost;

class Decode {
public:
	Decode() : makes(UINT_MAX) {
		assert(inst.O > 0);
		jobUsage.resize(inst.J);
		machUsage.resize(inst.M);
		schedule.resize(inst.O);
		indeg = inst.posetIndeg;
	}
	~Decode() {}

	void reset() {
		fill(jobUsage.begin(), jobUsage.end(),0);
		fill(machUsage.begin(), machUsage.end(),0);
		fill(schedule.begin(), schedule.end(),0);
	    indeg = inst.posetIndeg;
		makes = UINT_MAX;
		omega.clear();
	}

	
	void simpleDecode(const double delayWeight, const vector<unsigned> & priority, vector<unsigned> & indeg) {
#ifndef NDEBUG
		assert(priority.size() == inst.O);
		unsigned inserted = 0;
		assert(delayWeight>=0   && delayWeight<=1.0);
#endif //NDEBUG

		//fprintf(stderr, "\n\ta\n");

		reset();
		vector<unsigned> start(inst.O);
		unsigned bestStart;
		vector<unsigned> finish(inst.O);
		unsigned bestFinish;
		unsigned slack;
		unsigned chosenOp;
		unsigned chosenOpIndex;
		unsigned childOp;

		makes = 0;

		//fprintf(stderr, "\n\tb\n");

		vector<unsigned> omega = inst.posetRoots;
		indeg = inst.posetIndeg;

		start[0] = UINT_MAX;

		//fprintf(stderr, "\n\tc\n");

		while( ! omega.empty()) {
#ifndef NDEBUG
			inserted++;
			assert(inserted<inst.O);
#endif //NDEBUG

			bestStart = UINT_MAX;
			bestFinish = UINT_MAX;

			//fprintf(stderr, "\n\ttc - a\n");

			//finding start and finish -- best start and finish
			for(unsigned o : omega) {
				assert(o != 0);
				assert(o<inst.O);
				assert(indeg[o] == 0);
				start[o] = max(jobUsage[inst.operToJ[o]], machUsage[inst.operToM[o]]);
				bestStart = min(bestStart, start[o]);

				finish[o] = start[o] + inst.P[o];
				bestFinish = min(bestFinish, finish[o]);
			}
			assert(bestFinish>=bestStart);

			//fprintf(stderr, "\n\ttc - b\n");

			slack = delayWeight*(bestFinish - bestStart);

			//filter the set of omega such that the start is not delay by more then slack -- select chosenOp
			chosenOp = 0;
			assert(start[chosenOp] == UINT_MAX);
			for(unsigned index=0; index<omega.size(); index++) {
				//fprintf(stderr, "\n\ttc - b - a\n");
			//for(unsigned o : omega) {
				unsigned o = omega[index];
				assert(o!=0 && o<inst.O);
				//is in slack interval
				//fprintf(stderr, "\n\ttc - b - b\n");
				if(start[o] <= bestStart+slack) {
					//has precedence accoring to state
					if(chosenOp==0   ||   priority[chosenOp]>priority[o]) {
						chosenOp = o;
						chosenOpIndex = index;
					}
				}
				//fprintf(stderr, "\n\ttc - b - c\n");
			}
			assert(chosenOp != 0);

			
			//schedule chosenOp
			assert(jobUsage[inst.operToJ[chosenOp]] <= start[chosenOp]);
			assert(machUsage[inst.operToM[chosenOp]] <= start[chosenOp]);
			schedule[chosenOp] = start[chosenOp];
			makes = max(makes, start[chosenOp] + inst.P[chosenOp]);

			//fprintf(stderr, "\n\ttc - c\n");

			//update job and machine usage
			jobUsage[inst.operToJ[chosenOp]] = start[chosenOp] + inst.P[chosenOp];
			machUsage[inst.operToM[chosenOp]] = start[chosenOp] + inst.P[chosenOp];

			//remove inserted op from omega
			assert(omega[chosenOpIndex] == chosenOp);
			omega.erase(omega.begin()+chosenOpIndex);

			//fprintf(stderr, "\n\ttc - d\n");

			//insert freed operations to omega
			for(unsigned childMachine : inst.posets[inst.operToJ[chosenOp]].children[inst.operToM[chosenOp]]) {
				childOp = inst.jmToIndex[inst.operToJ[chosenOp]][childMachine];
				assert(indeg[childOp] > 0);
				indeg[childOp]--;
				if(indeg[childOp] == 0) {
					omega.push_back(childOp);
				}
			}
		}
		assert(inserted==inst.O-1);
		//fprintf(stderr, "\n\ttc - e\n");
	}

	string toString() const {
		string str = ""; 
		str.append("\njobUsage: ");
		for(unsigned o=0; o<inst.J; o++)
		    str.append(unsigStr(o) + "[" + unsigStr(jobUsage[o]) + "] - ");
		str.append("\nmachUsage: ");
		for(unsigned o=0; o<inst.M; o++)
		    str.append(unsigStr(o) + "[" + unsigStr(machUsage[o]) + "] - ");
		str.append("\nschedule: ");
		for(unsigned o=0; o<inst.O; o++)
		    str.append(unsigStr(o) + "[" + unsigStr(schedule[o]) + "] - ");
		str.append("\nomega: ");
		for(unsigned o=0; o<omega.size(); o++)
		    str.append(unsigStr(o) + "[" + unsigStr(omega[o]) + "] - ");
		str.append("\nindeg: ");
		for(unsigned o=0; o<inst.O; o++)
		    str.append(unsigStr(o) + "[" + unsigStr(indeg[o]) + "] - ");
		str.append("\nmakes: ");
		str.append(unsigStr(makes));

		return str;
	}

	
	vector<unsigned> jobUsage;
	vector<unsigned> machUsage;
	vector<unsigned> schedule;
	vector<unsigned> omega;
	vector<unsigned> indeg;
	unsigned makes;
};
