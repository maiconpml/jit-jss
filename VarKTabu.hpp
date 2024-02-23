/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/
#pragma once

#include <algorithm> 
#include <chrono>

#include "Settings.hpp"
#include "Utilities.hpp"
#include "VarKDG.hpp"
//#include "Analyzer.hpp"

using namespace chrono;

#ifdef VAR_PARALLEL_MACHS



//tabu list - holds forbidden moves
//holds pairs of operations that can not be swapped. does not consider order
class TabuList {
public:
	TabuList() {  }
	
	TabuList(int _tenure) : curIter(0), tenure(_tenure) {
		if(tenure==0)
			return;
		tabu.resize(boost::extents[inst.O][inst.O]);
		fill(tabu.data(),tabu.data()+tabu.num_elements(),-1);
	}
	TabuList(TabuList && other) {
		this->swap(other);
	}
	TabuList(const TabuList& other) {
		curIter = other.curIter;
		tenure = other.tenure;
		tabu.resize(boost::extents[other.tabu.shape()[0]][other.tabu.shape()[1]]);
		tabu = other.tabu;
	}



	TabuList& operator=(TabuList other) {
		this->swap(other);
		return *this;
	}



	void swap(TabuList& other) {
		using std::swap;
		swap(curIter,other.curIter);
		swap(tenure,other.tenure);
		tabu.resize(boost::extents[other.tabu.shape()[0]][other.tabu.shape()[1]]);
		swap(tabu,other.tabu);
	}



	string toString() const {
		string str = "";
		
		str.append("curIter: "+unsigStr(curIter));
		str.append("  ;  tenure: "+unsigStr(tenure));
		str.append(" --- forbidden: ");
		str.append(listForbidden());

		return str;
	}



	void passTime(int time) {
		assert(time <= tenure);
		curIter += time;
	}



	string listForbidden() const {
		string str = "";
		for(unsigned o1=0; o1<inst.O; o1++) {
			for(unsigned o2=0; o2<inst.O; o2++) {
				if(o2 > o1) {
					assert(tabu[o1][o2] <= curIter+tenure);
					if( ! canMove(o1, o2)) {
						str.append(unsigStr(o1) + "-" + unsigStr(o2)+"["+unsigStr(age(o1,o2))+"]    ");
						assert( ! canMove(o2, o1));
					} else {
						assert(canMove(o2, o1));
					}
				}
			}
		}
		return str;
	}



	void insert(unsigned o1, unsigned o2) {
		assert(o1!=o2);
		assert(tenure != 0);
		assert(o1 < inst.O);
		assert(o2 < inst.O);
		
		curIter++;
		tabu[o2][o1] = curIter+tenure;
		tabu[o1][o2] = curIter+tenure;
	}



	bool canMove(unsigned o1, unsigned o2) const {
		assert(tenure > 0);
		assert(tabu[o1][o2] <= curIter + tenure);
		assert(tabu[o1][o2] == tabu[o2][o1]);
		return tabu[o1][o2] <= curIter;
	}



	int age(unsigned o1, unsigned o2) const {
		assert(tenure + curIter >= tabu[o1][o2]);
		assert(tabu[o1][o2] == tabu[o2][o1]);
		return tenure + curIter - tabu[o1][o2];
	}



	int timeToLeave(unsigned o1, unsigned o2) const {
		int delta = tabu[o1][o2] - curIter;
		assert(delta>0);
		assert(delta <= tenure);
		assert(tenure-age(o1, o2) == delta);
		return delta;
	}



	int curIter;
	int tenure;
	multi_array<int,2> tabu;
};








//Solution + neighbourhood unvisited + tabu list
//solution is only select (complete selection of machine orders
class TabuTrio {
public:
	TabuTrio() : dummy(true) {  }
	TabuTrio(const Select & _select, const vector<Cand> & _cands, const TabuList & _tabuList) : dummy(false), select(_select), cands(_cands), tabuList(_tabuList) {  }
	~TabuTrio() {  }



	string toString() const {
		string str;

		if(dummy)
			return "DUMMY TRIO";

		str.append("Cands size: " +unsigStr(cands.size()));

		return str;
	}



	bool dummy;
	Select select;
	vector<Cand> cands; //cand moves - removed best one which was used
	TabuList tabuList;
};








//elite list - conatins tabu trios - uses structure of nowicky smithnicky 2001
class JumpList {
public:
	JumpList(unsigned _maxL) : empty(true), basePos(0), topPos(0), maxL(_maxL) { ttList.resize(_maxL); }
	~JumpList() {  }

	string toString() const {
		string str = "";

		str.append(empty ? "empty" : "NOT empty");
		str.append(" -- basePos: "+unsigStr(basePos));
		str.append(" -- topPos: "+unsigStr(topPos));
		str.append(" -- maxL: "+unsigStr(maxL));
		//str.append(" -- TT ttList: "+unsigStr(ttList.size()));
		str.append(" -- TT ttList:");
		for(const TabuTrio & tt : ttList) {
			str.append("\n\t"+tt.toString());
		}

		return str;
	}



	void push(const TabuTrio & tt) {
		//cout << "\tPUSHED\n" << toString() << endl;

		assert(maxL>0);
		if(empty){
			empty = false;
			topPos = 0;
			basePos = 0;
		} else {
			topPos++;
			topPos  = topPos % maxL;
			if(topPos == basePos) {
				basePos++;
				basePos  = basePos % maxL;
			}
		}

		ttList[topPos] = tt;
	}


	
	void updateCands(const vector<Cand> & cands) {
		assert( ! empty);

#ifndef NDEBUG
		bool found;
		for(const Cand & p1 : cands) {
			found = false;
			for(const Cand & p2 : ttList[topPos].cands) {
				if(p1.bailiff == p2.bailiff    &&    p1.provost == p2.provost)
					found = true;
			}
			assert(found);
		}
#endif
		ttList[topPos].cands = cands;

		while(ttList[topPos].cands.empty()) {
			if(topPos == basePos) {
				empty = true;
				return;
			}

			if(topPos == 0)
				topPos = maxL-1;
			else
				topPos--;
		}
	}



	TabuTrio pop() {
		if(empty)
			return TabuTrio();

		while(ttList[topPos].cands.empty()) {
			if(topPos == basePos) {
				empty = true;
				return TabuTrio();
			}

			if(topPos == 0)
				topPos = maxL-1;
			else
				topPos--;
		}

		return ttList[topPos];
	}



	bool empty;
	unsigned basePos;
	unsigned topPos;
	const unsigned maxL;
	vector<TabuTrio> ttList;
};








class Params {
public:
	Params() { }
	~Params() { }



	string toString() const {
		string str = "";

		str.append("tenure: " + unsigStr(tenure)+ "  -  initialjumpLimit: "+ unsigStr(initialjumpLimit) + "  -  jumpLimitDecrease: " + unsigStr(jumpLimitDecrease) + "  -  bjSize: " + unsigStr(bjSize) + "  -  maxD: " + unsigStr(maxD) + "  -  maxC: " + unsigStr(maxC));

		str.append("\n" + (string)(learn ? "" : "NOT") + " learn   -  learnMinIters: " + unsigStr(learnMinIters) + "  -  learnMinKW: " + doubleStr(learnMinKW) + "  -  learnTrimm: " + doubleStr(learnTrimm)); 

		return str;
	}

	
	
	//pure tabu params
	unsigned tenure;
	unsigned initialjumpLimit;
	unsigned jumpLimitDecrease;
	unsigned bjSize;
	unsigned maxD;
	unsigned maxC;

	//learning params
	bool learn;
	unsigned learnMinIters;
	double learnMinKW;
	double learnTrimm;
};








//holds main functions of tabu search
namespace Tabu {

	//@return: makespan - UINFTY if no valid states
	//removes takes cand from cands, and adjusts select with the chosen move
	//--does not change cands if no valid move
	unsigned nsp(Select & select,  vector<Cand> & cands, TabuList & tabuList, Utility & utility, Stages & stages, vector<unsigned> & auxDists, unsigned aspiration) {
#ifndef NDEBUG
		assert( ! cands.empty());
		assert(utility.isAlloc());
		assert(stages.isAlloc());
		assert(auxDists.size() == inst.O);
	      const Select testSelect = select;
#endif

		unsigned bailiff;
		unsigned provost;
		unsigned direction;
		unsigned contraryDirection;
		unsigned pivot; //to undo the shift
		vector<bool> stillToTestCandPos(cands.size(), true);//controls which cands were tested - first goes trough unforbidden then forbidden

		//Unforbidden or aspiration here
		unsigned chosenMakes = UINFTY;
		unsigned chosenSwapPos;
		Cand chosenCand;

		//Forbidden Non Profitable
		int fnpStateAge = -1;
		Cand oldestCand;
		unsigned oldestFnpSwapPos;
		unsigned oldestMakes = UINFTY;;

		unsigned bufMakes;

		//unforbidden moves
		for(unsigned pPos=0; pPos<cands.size(); pPos++) {		
			bailiff = cands[pPos].bailiff;
			provost = cands[pPos].provost;
			direction = cands[pPos].direction;
			contraryDirection = (direction==FORTH ? contraryDirection==BACK : contraryDirection==FORTH);
			
			assert(bailiff != 0);
			assert(provost != 0);
			assert((direction==BACK  &&  contraryDirection == FORTH)    ||
				 (direction==FORTH  &&  contraryDirection == BACK));

			if(tabuList.canMove(bailiff, provost)) {
				stillToTestCandPos[pPos] = false;
				
				if(direction == FORTH) {
					pivot = select.mach[provost];
				} else  {// direction == BACK
					pivot = select._mach[provost];
				}
				assert(pivot != 0);
				
				select.shift(bailiff, provost, direction);
				
				bufMakes = select.fastMakes(utility, stages, auxDists);
				if(bufMakes < chosenMakes) {
					chosenMakes = bufMakes;
					chosenSwapPos = pPos;
					chosenCand = cands[pPos];
				}

				//undoing the shift
				select.shift(pivot, provost, contraryDirection);
				assert(select == testSelect);
			}//if tabuList.canMove(o1, o2)
		}//unforbinden for

		//forbidden moves
		for(unsigned pPos=0; pPos<cands.size(); pPos++) {
			if(stillToTestCandPos[pPos]) {
				bailiff = cands[pPos].bailiff;
				provost = cands[pPos].provost;
				direction = cands[pPos].direction;
				contraryDirection = (direction==FORTH ? contraryDirection==BACK : contraryDirection==FORTH);

				assert(bailiff != 0);
				assert(provost != 0);
				assert((direction==BACK  &&  contraryDirection == FORTH)    ||
					 (direction==FORTH  &&  contraryDirection == BACK));
				assert( ! tabuList.canMove(bailiff, provost));

				if(direction == FORTH) {
					pivot = select.mach[provost];
				} else  {// direction == BACK
					pivot = select._mach[provost];
				}
				assert(pivot != 0);

				select.shift(bailiff, provost, direction);
				
				bufMakes = select.fastMakes(utility, stages, auxDists);
				//aspiration
				if(bufMakes < aspiration   &&   bufMakes < chosenMakes) {
					chosenMakes = bufMakes;
					chosenSwapPos = pPos;
					chosenCand = cands[pPos];
				}

				//crap moves
				if(chosenMakes==UINT_MAX   &&   fnpStateAge < tabuList.age(bailiff, provost)   &&   bufMakes < UINFTY) {
					fnpStateAge = tabuList.age(bailiff, provost);
					oldestCand = cands[pPos];
					oldestFnpSwapPos = pPos;
					oldestMakes = bufMakes;
				}

				//undoing the shift
				select.shift(pivot, provost, contraryDirection);
				assert(select == testSelect);				
			}//still to test if
		}//forbiden for

		//perform move - remove used cand - update tabu list
		if(chosenMakes != UINT_MAX) {
			assert(chosenCand.isValid());
			assert(chosenSwapPos < cands.size());
			assert(chosenCand == cands[chosenSwapPos]);

			select.shift(chosenCand);
			assert(select.fastMakes(utility, stages, auxDists) == chosenMakes);
			
			tabuList.insert(chosenCand.bailiff, chosenCand.provost);
			cands.erase(cands.begin() + chosenSwapPos);

			return chosenMakes;
		} else {
			assert( ! chosenCand.isValid());
			//not a single valid solution in the neighbourhood
			if( ! oldestCand.isValid()) {
				assert(testSelect == select);
				return UINFTY;
			}

			assert(fnpStateAge >= 0);
			select.shift(oldestCand);

			tabuList.passTime(tabuList.timeToLeave(oldestCand.bailiff, oldestCand.provost));
			tabuList.insert(oldestCand.bailiff, oldestCand.provost);
			cands.erase(cands.begin() + oldestFnpSwapPos);

			return oldestMakes;
		}
	}

	

	void runTabu(DG & bestDg, vector<Cand> & cands, Utility & util, Stages & stages, const Params & params, double maxMillisecs, const high_resolution_clock::time_point & tpStart) {
		assert( ! bestDg.dirty);
		//(0) declare aux
		//cycle detector aux isCycle
		bool makesCycleDetected;
		vector<int> oldValues(params.maxD, -1);
		unsigned posCurMakes = 0;
		unsigned cycleLastPos = 0;
		unsigned cycleL = 0;
		const unsigned maxLen = params.maxD * params.maxC;

		//elite set jumping vars
		unsigned tabuIter =0;
		unsigned noImproveIters = 0;
		unsigned curJumpLimit = params.initialjumpLimit;
		bool jumped = false;
		bool newBestFound = true;

		DG curDg = bestDg;

		bool existViable = true; //did the last computed neighbourhood yielded viable neighbours?

		TabuTrio tt;
		TabuList tabuList(params.tenure);
		JumpList jumpList(params.bjSize);

		Select newBestDgSelect;
		TabuList newBestTL;

		vector<unsigned> auxDists(inst.O);

		//(1)main loop
		while(maxMillisecs > duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count()) {
			tabuIter++;
			noImproveIters++;

			makesCycleDetected = isCyle(oldValues, posCurMakes, cycleLastPos, cycleL, curDg.makes, newBestFound, params.maxD, maxLen);
			//jumping
			if((curJumpLimit<noImproveIters)   ||   makesCycleDetected   ||   ! existViable) {
				assert( ! newBestFound);
				existViable = true;

				jumped = true;
				if(curJumpLimit > params.jumpLimitDecrease) 
					curJumpLimit -= params.jumpLimitDecrease;
				noImproveIters = 0;

				tt = jumpList.pop();
				if( ! tt.dummy) {					
					curDg.select = tt.select;
#ifndef NDEBUG
					curDg.dirty = true;
#endif
					cands = tt.cands;
					tabuList = tt.tabuList;
					assert( ! cands.empty());
				} else { //Stop -> no more states in jump list
					return;
				}
			} else { //NOT JUMPING
				jumped = false;
				curDg.setMeta(cands, util, stages);
				if(cands.size() == 0) {
					return;
				}
			}
			assert( ! cands.empty());

			//stuff to store in elite list
			if(newBestFound) {
				assert( ! jumped);
				newBestDgSelect = curDg.select;
				newBestTL = tabuList;
			}

			curDg.makes  = nsp(curDg.select, cands, tabuList, util, stages, auxDists, bestDg.makes);
			existViable = curDg.makes < UINFTY;
#ifndef NDEBUG
			curDg.dirty = true;
#endif

			// last iter new best was found? store oldState with paths not taken
			if(newBestFound) {
				if( ! cands.empty())
					jumpList.push(TabuTrio(newBestDgSelect, cands, newBestTL));
				newBestFound = false;
			}

			//Removing the taken path from jump list
			if(jumped) {
				assert( ! newBestFound);
				jumpList.updateCands(cands);
			}

			//STEP new best ??
			if(curDg.makes < bestDg.makes) {
				noImproveIters = 0;
				bestDg = curDg;
				curJumpLimit = params.initialjumpLimit;
				newBestFound= true;

				bestDg.millisecs=duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
			}		
		}//while(maxMillisecs > duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())
	}
	


	//reads instance in instPath
	//generates initial solution with jobTailHeuristic()
	//calls runTabu with params
	//print results
	void fullTabu(const string & instPath, Params & params, double maxMillisecs) {
		//(0)initialize
		//read instance file
		inst.parse(instPath);

		high_resolution_clock::time_point tpStart = high_resolution_clock::now();

		//declare vars + allocate containers
		DG dg;
		vector<Cand> cands;
		Utility util;
		Stages stages;
		unsigned initialMakes;
		unsigned iniitalMillisecs;

		dg.alloc();
		util.alloc();
		stages.alloc();

		//(1)initial solution
		dg.jobTailHeuristic();
		dg.setMeta(cands, util, stages);
		initialMakes = dg.makes;
		assert(dg.makes < UINFTY);
		iniitalMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
		
		//(2)run tabu
		runTabu(dg, cands, util, stages, params, maxMillisecs, tpStart);
		dg.setMeta(cands, util, stages);

		//(3)print results
		cout << initialMakes << " " << iniitalMillisecs << " " << dg.makes << " " << dg.millisecs << " " << endl;

		//(4)verify solution is valid
		assert( ! dg.dirty);
	      inst.verifyPrallelSchedule(dg.heads, dg.makes);
	}
	
};

#endif //VAR_PARALLEL_MACHS
