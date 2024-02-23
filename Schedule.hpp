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

#ifdef VAR_PARALLEL_MACHS







class Window {
public:
	Window() { 
		start = 0;
		finish = UINFTY;
	}
	Window(unsigned _start, unsigned _finish) : start(_start), finish(_finish) { assert(start < finish); }
	~Window() { }
	

	unsigned start;
	unsigned finish;
};











class Bucket {
public:
	Bucket() { vaccancy.push_back(Window()); }
	~Bucket() { }
	
	
	string toString() const {
		assert( ! vaccancy.empty());
		string str;

		str.append(unsigStr(vaccancy[0].start));
		for(unsigned u=1; u<vaccancy.size(); u++) {
			assert(vaccancy[u-1].start < vaccancy[u-1].finish);
			assert(vaccancy[u].start < vaccancy[u].finish);
			assert(vaccancy[u-1].finish < vaccancy[u].start);
			
			str.append("__");
			str.append(unsigStr(vaccancy[u-1].finish));
			str.append("##");
			str.append(unsigStr(vaccancy[u].start));
		}
		
		str.append("__");
		str.append(unsigStr(vaccancy.back().finish));
		
		return str;
	}
	

	
	
	//release is the job head, dispatching rule has to obey job order and release is start + duration of the predecessor
	//@return: Best possible spot to start the operation
	//gap is  return - release
	unsigned bestStart(unsigned release, unsigned duration) const {
		unsigned bestStart = UINFTY;
		
		for(const Window & w : vaccancy) {
			if(w.finish >= release+duration    &&    w.finish-w.start>=duration) {
				bestStart = w.start;
				break;
			}
		}
		
		
		assert(bestStart < UINFTY);
		return bestStart;
	}
	
	
	
	
	//expects the duration requested to be vancant
	//@return: new buket makespan
	unsigned occupy(const Window & block) {
		assert(block.start != UINFTY);
		assert(block.finish != UINFTY);
		
		
		unsigned buf;
		
		
		for(unsigned u=0; u<vaccancy.size(); u++) {
			
			if(vaccancy[u].finish >= block.finish ) {//will insert in this vacancy
				assert(vaccancy[u].start <= block.start);//else no place for block
				
				if(vaccancy[u].start == block.start     &&     vaccancy[u].finish == block.finish) {//we have a perfect fit, just delete vaccancy spot
					vaccancy.erase(vaccancy.begin()+u);
				} else if(vaccancy[u].start == block.start) { // no gap before
					vaccancy[u].start = block.finish;
				} else if(vaccancy[u].finish == block.finish) {// no gap after
					vaccancy[u].finish = block.start;
				} else {//gap on both sides
					buf = vaccancy[u].finish;
					vaccancy[u].finish = block.start;
					assert(vaccancy[u].start < vaccancy[u].finish);
					vaccancy.insert(vaccancy.begin()+u+1, Window(block.finish,buf));
				}
				
				
				break;//already insrted, can stop
			}
		}
		assert(vaccancy.back().finish == UINFTY);
		return vaccancy.back().start;
	}
	

	vector<Window> vaccancy;
};













class Machine {
public:
	Machine(unsigned reps) : makes(0) { 
		buckets.resize(reps);
	}
	~Machine() {}
	
	string toString() const {
		string str = "";
		for(unsigned k=0; k<buckets.size(); k++)
			str.append("\n\tk"+unsigStr(k)+" "+buckets[k].toString());
		return str;
	}
	
	//@return: start of the the operations
	//updates makes
	//Always schedule without delay if possible
	//@tiebreaker: smallestGap, largestGap
	unsigned insertOper(unsigned release, unsigned duration, unsigned tiebreaker) {
		unsigned bestBucket = UINFTY;
		unsigned bestBucketStart = UINFTY;
		unsigned aStart;
		bool foundNonDelay = false;
		unsigned buckMakes;
	
		for(unsigned u=0; u<buckets.size(); u++) {
			aStart = buckets[u].bestStart(release, duration);
			if( ! foundNonDelay) {//still have to delay
				if(aStart < bestBucketStart) {//starts as early as possible
					bestBucket = u;
					bestBucketStart = aStart;
				}
				if(aStart <= release)
					foundNonDelay = true;
			} else {//have already found nondelay
				assert(bestBucketStart < UINFTY);
				if(aStart <= release) { //found another nondelay
					if(tiebreaker == SMALLER_GAP) {
						if(aStart > bestBucketStart) {//new bucket start is thighter -> smaller gap
							bestBucket = u;
							bestBucketStart = aStart;
						} 
					} else if(tiebreaker == LARGER_GAP) {
						if(aStart < bestBucketStart) {//new bucket start is looser -> larger gap
							bestBucket = u;
							bestBucketStart = aStart;
						} 
					} else {
						throw errorText("Invalid option for tiebreaker","insertOper","Schedule.hpp");
					}
				}
			}
		}
		
		aStart = max(release, bestBucketStart);
		buckMakes = buckets[bestBucket].occupy(Window(aStart, aStart+duration));
		makes = max(makes, buckMakes);
		return aStart;
	}
	
	unsigned makes;
	vector<Bucket> buckets;
};










class Schedule {
public:
	Schedule() {  }
	~Schedule() {  }
	
	//clear and alloc
	void clear() {
		machines.clear();	
		starts.clear();
		alloc();
	}
	
	void alloc() {
		machines.clear();
		for(unsigned m=0; m<inst.M; m++)
			machines.push_back(Machine(inst.varK[m]));
		starts.resize(inst.O, 0);
		makes =0;
	}
	
	string toString() const {
		string str = "makes: " + unsigStr(makes);
		for(unsigned m=0; m<machines.size(); m++) {
			str.append("\nm"+unsigStr(m)+"\n");
			str.append(machines[m].toString());
		}
		
		assert(starts[0]==0);
		str.append("\n\nStarts o(t) :");
		for(unsigned u=1; u<starts.size(); u++) {
			str.append(unsigStr(u)+"("+unsigStr(starts[u])+")  ");
		}
		
		return str;
	}
	
	
	
	//schedules one by one - choses the job that the tp operation has LOWER priority
	void dipatchByPiority(const vector<unsigned> & priority, unsigned tiebreaker) {
		clear();
		
		vector<unsigned> freeOpers(inst.J, 0);
		vector<unsigned> release(inst.J, 0);
		unsigned curJob;
		unsigned curMach;
		unsigned curOp;
		unsigned highestPrioJob;
		unsigned highestPrio;
		unsigned curStart;
		
		for(unsigned aRoot : inst.roots) {
			curJob = inst.operToJ[aRoot];
			assert(freeOpers[curJob] == 0);
			
			freeOpers[curJob] = aRoot;
		}
		
		for(unsigned iter=0; iter<inst.O-1; iter++) {//-1 because oper 0 is not scheduled
			highestPrio = UINFTY;
			highestPrioJob = UINFTY;
			for(unsigned j=0; j<inst.J; j++) {
				if(freeOpers[j] != 0) {
					if(highestPrio > priority[freeOpers[j]]) {
						highestPrioJob = j;
						highestPrio = priority[freeOpers[j]];
					}
				}
			}
			assert(highestPrio < UINFTY);
			
			curOp = freeOpers[highestPrioJob];
			assert(curOp>0   &&   curOp<inst.O);
			curJob = inst.operToJ[curOp];
			curMach = inst.operToM[curOp];
			
			assert(curOp>0   &&   curOp<inst.O);
			if(inst.next[freeOpers[highestPrioJob]].size() == 1)
				freeOpers[highestPrioJob] = inst.next[freeOpers[highestPrioJob]][0];
			else
				freeOpers[highestPrioJob] = 0;
			
			curStart = machines[curMach].insertOper(release[curJob], inst.P[curOp], tiebreaker);
			release[curJob] = curStart + inst.P[curOp];
			starts[curOp] = curStart;
			makes = max(makes, machines[curMach].makes);		
		}	
	}
	
	
	

	void dispatchLongestJobTail(unsigned tiebreaker) {
		clear();
	
	
		vector<unsigned> jobTail(inst.J, 0);
		vector<unsigned> freeOpers(inst.J, 0);
		vector<unsigned> release(inst.J, 0);
		unsigned curJob;
		unsigned curMach;
		unsigned curOp;
		unsigned highestTailJob;
		unsigned highestTail;
		unsigned curStart;
		
		
		for(unsigned aRoot : inst.roots) {
			curJob = inst.operToJ[aRoot];
			
			assert(jobTail[curJob] == 0);
			assert(freeOpers[curJob] == 0);
			
			freeOpers[curJob] = aRoot;
			curOp = aRoot;
			for(unsigned u=0; u<inst.M; u++) {
				assert(curOp<inst.O);
				assert(u==inst.M-1    ||    inst.next[curOp].size()==1);
				assert(u!=inst.M-1    ||    inst.next[curOp].size()==0);
				assert(inst.operToJ[curOp] == curJob);
				
				jobTail[curJob] += inst.P[curOp];
				
				if(u!=inst.M-1)
					curOp = inst.next[curOp][0];
			}
		}
		
		for(unsigned iter=0; iter<inst.O-1; iter++) {//-1 because oper 0 is not scheduled
			highestTail = 0;
			highestTailJob = UINFTY;
			for(unsigned j=0; j<inst.J; j++) {
				//cout << "j: " << j << "    tail: " <<  jobTail[j] << "    ;    freeOper: " << freeOpers[j] << endl;
				if(highestTail < jobTail[j]) {
					highestTailJob = j;
					highestTail = jobTail[j];
				}
			}
			assert(highestTail > 0);
			
			curOp = freeOpers[highestTailJob];
			assert(curOp>0   &&   curOp<inst.O);
			curJob = inst.operToJ[curOp];
			curMach = inst.operToM[curOp];
			
			if(inst.next[freeOpers[highestTailJob]].size() == 1)
				freeOpers[highestTailJob] = inst.next[freeOpers[highestTailJob]][0];
			else
				freeOpers[highestTailJob] = 0;
			
			curStart = machines[curMach].insertOper(release[curJob], inst.P[curOp], tiebreaker);
			release[curJob] = curStart + inst.P[curOp];
			assert(jobTail[curJob] >= inst.P[curOp]);
			jobTail[curJob] -= inst.P[curOp];
			starts[curOp] = curStart;
			makes = max(makes, machines[curMach].makes);
			
		}
#ifndef NDEBUG
		for(unsigned u : jobTail) assert(u == 0);
#endif
	
	}
	
	//to set state with same machine orders
	void getMachOrders(vector<vector<unsigned>> & allMachOrder) const {
		assert(allMachOrder.size() == inst.M);
		vector<pair<unsigned,unsigned>> startsXoper;
		
		for(unsigned m=0; m<inst.M; m++) {
			startsXoper.clear();
			for(unsigned op : inst.machOpers[m])
				startsXoper.push_back(pair<unsigned, unsigned>(starts[op], op));
			sort(startsXoper.begin(), startsXoper.end());
			allMachOrder[m].clear();
			for(const pair<unsigned, unsigned> & p : startsXoper)
				allMachOrder[m].push_back(p.second);
		}
	
		
	}
	
	
	unsigned makes;
	vector<Machine> machines;	
	vector<unsigned> starts;
};


#endif //VAR_PARALLEL_MACHS
