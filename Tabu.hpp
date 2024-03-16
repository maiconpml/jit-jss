/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/
#pragma once

#include <algorithm> 
#include <chrono>

#include "Settings.hpp"
#include "Utilities.hpp"
#include "State.hpp"
#include "Analyzer.hpp"

using namespace chrono;


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




class TabuTrio {
public:
	TabuTrio() : dummy(true) {  }
	TabuTrio(const State & _state, const vector<pair<unsigned,unsigned>> & _cands, const TabuList & _tabuList) : dummy(false), state(_state), cands(_cands), tabuList(_tabuList) {  }
	~TabuTrio() {  }


	string toString() const {
		string str;

		if(dummy)
			return "DUMMY TRIO";

		str.append("Makes: " +unsigStr(state.makes)+"   ---   ");
		str.append("Cands size: " +unsigStr(cands.size()));

		return str;
	}


	bool dummy;
	State state;
	vector<pair<unsigned,unsigned>> cands; //cand moves - removed best one which was used
	TabuList tabuList;
};




class TabuTrioN6 {
public:
	TabuTrioN6() : dummy(true) {  }
	TabuTrioN6(const State & _state, const vector<tuple<unsigned, unsigned, unsigned,unsigned>> & _cands, const TabuList & _tabuList) : dummy(false), state(_state), cands(_cands), tabuList(_tabuList) {  }
	~TabuTrioN6() {  }


	string toString() const {
		string str;

		if(dummy)
			return "DUMMY TRIO";

		str.append("Makes: " +unsigStr(state.makes)+"   ---   ");
		str.append("Cands size: " +unsigStr(cands.size()));

		return str;
	}


	bool dummy;
	State state;
	vector<tuple<unsigned,unsigned, unsigned, unsigned>> cands; //cand moves - removed best one which was used
	TabuList tabuList;
};




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


	
	void updateCands(const vector<pair<unsigned, unsigned>> & cands) {
		assert( ! empty);

#ifndef NDEBUG
		bool found;
		for(const pair<unsigned, unsigned> & p1 : cands) {
			found = false;
			for(const pair<unsigned, unsigned> & p2 : ttList[topPos].cands) {
				if(p1.first == p2.first    &&    p1.second == p2.second)
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




class JumpListN6 {
public:
	JumpListN6(unsigned _maxL) : empty(true), basePos(0), topPos(0), maxL(_maxL) { ttList.resize(_maxL); }
	~JumpListN6() {  }

	string toString() const {
		string str = "";

		str.append(empty ? "empty" : "NOT empty");
		str.append(" -- basePos: "+unsigStr(basePos));
		str.append(" -- topPos: "+unsigStr(topPos));
		str.append(" -- maxL: "+unsigStr(maxL));
		//str.append(" -- TT ttList: "+unsigStr(ttList.size()));
		str.append(" -- TT ttList:");
		for(const TabuTrioN6 & tt : ttList) {
			str.append("\n\t"+tt.toString());
		}

		return str;
	}



	void push(const TabuTrioN6 & tt) {
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



	void updateCands(const vector<tuple<unsigned, unsigned, unsigned, unsigned>> & cands) {

		//cout << "\tUPDATED\n" << toString() << endl;

		assert( ! empty);
		assert(cands.size() == ttList[topPos].cands.size()-1); //removes exactly one for last move
#ifndef NDEBUG
		bool found;
		for(const tuple<unsigned, unsigned, unsigned, unsigned> & p1 : cands) {
			found = false;
			for(const tuple<unsigned, unsigned, unsigned, unsigned> & p2 : ttList[topPos].cands) {
				if(p1 == p2)
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



	TabuTrioN6 pop() {
		if(empty)
			return TabuTrioN6();

		while(ttList[topPos].cands.empty()) {
			if(topPos == basePos) {
				empty = true;
				return TabuTrioN6();
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
	vector<TabuTrioN6> ttList;
};




namespace Tabu {

	void nspN6(State & theState, unsigned & lastOp, TabuList & tabuList, vector<tuple<unsigned, unsigned, unsigned, unsigned>> & cands, unsigned aspiration, vector<unsigned> & dists, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, vector<unsigned> & heads, vector<unsigned> & tails) {
#ifndef NDEBUG
		assert( ! cands.empty());
		assert(dists.size() == inst.O);
		assert(prev.size() == inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);
		State testState = theState;
		bool cycle;
#endif
		/*
		cout << "calling nspN6 for: " << endl;
		cout << theState.toString() << endl;
		cout << "Cands are: " << endl;
		for(tuple<unsigned, unsigned, unsigned, unsigned> & t : cands) {
			assert(get<0>(t)==JOB || get<0>(t)==MACH);
			assert(get<1>(t)==FORTH || get<1>(t)==BACK);
			cout << (get<0>(t)==JOB ? "\tJOB " : "\tMACH ") << (get<1>(t)==FORTH ? "FORTH" : "BACK") << " bailiff: " << get<2>(t) << "  provost: " << get<3>(t) << endl;
		}
		*/
		
		unsigned type;
		unsigned direction;
		unsigned oBailiff;
		unsigned oProvost;
		vector<bool> stillToTestCandPos(cands.size(), true);

		//Unforbidden or aspiration here
		unsigned chosenMakes = UINT_MAX;
		unsigned chosenType;
		unsigned chosenDirection;
	    unsigned chosenBailiff = 0;
		unsigned chosenProvost = 0;
		unsigned chosenShiftPos;

		unsigned pivot; //to undo shift

		//Forbidden Non Profitable
		int fnpStateAge = -1;
		unsigned oldestType  = UINT_MAX;
		unsigned oldestBailiff = 0;
		unsigned oldestProvost = 0;
		unsigned oldestFnpShiftPos;
		unsigned oldestDirection  = UINT_MAX;

#ifndef NDEBUG
		cycle = 
#endif
			theState.compTailsComplete(tails, indeg, Q);
		assert( ! cycle);

#ifndef NDEBUG
		cycle = 
#endif
			theState.compHeadsComplete(heads, indeg, Q);
		assert( ! cycle);


		//unforbidden moves
		for(unsigned pPos=0; pPos<cands.size(); pPos++) {
			type = get<0>(cands[pPos]);
			direction = get<1>(cands[pPos]);
			oBailiff = get<2>(cands[pPos]);
			oProvost = get<3>(cands[pPos]);
			assert(oBailiff != 0   &&    oBailiff<inst.O);
			assert(oProvost != 0   &&    oProvost<inst.O);
			assert(type==JOB   ||   type==MACH);
			assert((type==JOB && inst.operToJ[oProvost] == inst.operToJ[oBailiff])
			   ||  (type==MACH && inst.operToM[oProvost] == inst.operToM[oBailiff]));

			if(tabuList.canMove(oBailiff, oProvost)) {
				stillToTestCandPos[pPos] = false;

				//get provost
				if(direction == FORTH) {
					if(type == JOB) {
						pivot = theState.job[oProvost];
					} else {
						pivot = theState.mach[oProvost];
					}
				} else {
					if(type == JOB) {
						pivot = theState._job[oProvost];
					} else {
						pivot = theState._mach[oProvost];
					}
				}

				//cout << "init unforbid shift" << endl;
				theState.shift(type, direction, oBailiff, oProvost);
				assert(theState.verify());

#ifndef NDEBUG
				cycle = 
#endif
					theState.setMeta(dists, lastOp, prev, indeg, Q);
				assert( ! cycle);



				if(theState.makes < chosenMakes) {
					chosenMakes = theState.makes;
					chosenType = type;
					chosenDirection = direction;
					chosenBailiff = oBailiff;
					chosenProvost = oProvost;
					chosenShiftPos = pPos;
				}
				
				//undoing the shift
				//cout << "back unforbid shift" << endl;
				theState.shift(type, (direction==FORTH ? BACK : FORTH), pivot, oProvost);
				assert(theState.verify());
				assert(theState == testState);
			}
		}

		//forbidden moves
		for(unsigned pPos=0; pPos<cands.size(); pPos++) {
			if(stillToTestCandPos[pPos]) {
				type = get<0>(cands[pPos]);
				direction = get<1>(cands[pPos]);
				oBailiff = get<2>(cands[pPos]);
				oProvost = get<3>(cands[pPos]);
				assert(type==JOB   ||   type==MACH);
				assert(direction==FORTH   ||   direction==BACK);
				assert(oBailiff != 0);
				assert(oProvost != 0);

				assert(! tabuList.canMove(oBailiff, oProvost));

				//get provost
				if(direction == FORTH) {
					if(type == JOB) {
						pivot = theState.job[oProvost];
					} else {
						pivot = theState.mach[oProvost];
					}
				} else {
					if(type == JOB) {
						pivot = theState._job[oProvost];
					} else {
						pivot = theState._mach[oProvost];
					}
				}
				//cout << "init forbid shift" << endl;
				theState.shift(type, direction, oBailiff, oProvost);
				assert(theState.verify());


#ifndef NDEBUG
				cycle = 
#endif
					theState.setMeta(dists, lastOp, prev, indeg, Q);
				assert( ! cycle);


				//aspiration
				if(theState.makes < aspiration   &&   theState.makes < chosenMakes) {
					chosenMakes = theState.makes;
					chosenType = type;
					chosenDirection = direction;
					chosenBailiff = oBailiff;
					chosenProvost = oProvost;
					chosenShiftPos = pPos;
				}

				//crap moves
				if(chosenMakes==UINT_MAX   &&   fnpStateAge < tabuList.age(oBailiff, oProvost)) {
					fnpStateAge = tabuList.age(oBailiff, oProvost);
					oldestType = type;
					oldestDirection = direction;
					oldestBailiff = oBailiff;
					oldestProvost = oProvost;
					oldestFnpShiftPos = pPos;
				}

			    //undoing the shift
				//cout << "back forbid shift" << endl;
				theState.shift(type, (direction==FORTH ? BACK : FORTH), pivot, oProvost);
				assert(theState.verify());
				assert(theState == testState);
			}
		}

		//perform move - remove used cand - update tabu list
		if(chosenMakes != UINT_MAX) {
			//good move found
			assert(chosenBailiff != 0   &&   chosenBailiff < inst.O);
			assert(chosenProvost != 0   &&   chosenProvost < inst.O);
			assert(chosenShiftPos < cands.size());
			//cout << "Best shift" << endl;
			theState.shift(chosenType, chosenDirection, chosenBailiff, chosenProvost);
			tabuList.insert(chosenBailiff, chosenProvost);
			cands.erase(cands.begin() + chosenShiftPos);
		} else {
			//only crap moves
			assert(chosenBailiff == 0);
			assert(chosenProvost == 0);
			assert(oldestBailiff != 0   &&   oldestBailiff < inst.O);
			assert(oldestProvost != 0   &&   oldestProvost < inst.O);
			assert(fnpStateAge >= 0);
			//cout << "Old shift" << endl;
			theState.shift(oldestType, oldestDirection, oldestBailiff, oldestProvost);
			tabuList.passTime(tabuList.timeToLeave(oldestBailiff, oldestProvost));
			tabuList.insert(oldestBailiff, oldestProvost);
			cands.erase(cands.begin() + oldestFnpShiftPos);
		}
#ifndef NDEBUG
		cycle = 
#endif
			theState.setMeta(dists, lastOp, prev, indeg, Q);
		assert( ! cycle);
	}




		//@return: is optimal?
	//printWhenBetter use 0 to not print. will print preString makes seconds (according to tpSTart) d 
	bool evolveTabuN6(State & theState, unsigned tenure, unsigned initialjumpLimit, unsigned decreaseDivisor, unsigned bjSize, const high_resolution_clock::time_point & tpStart, vector<unsigned> & dists, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, const unsigned maxD, const unsigned maxC, const unsigned lowerBound, unsigned maxMillisecs, unsigned printWhenBetter, const string & preString, vector<unsigned> & heads, vector<unsigned> & tails, bool timeLog)  {
#ifndef NDEBUG
		unsigned testMakes;
		bool cycle;
#endif
		unsigned lastOp;
#ifndef NDEBUG
		cycle = 
#endif
			theState.setMeta(dists, lastOp, prev, indeg, Q);
		assert( ! cycle);
		assert(theState.makes >= lowerBound);


		if(theState.makes < printWhenBetter) {
			if(timeLog) {
				resultList.push_back(preString + unsigStr(theState.makes) + " " +doubleStr((duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())/1000.0) + " d");
				//cout << preString << theState.makes << " " << (duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())/1000.0 << " d" << endl;
			}
		}

		if(theState.makes == lowerBound)
			return true;

		State curState = theState;
		State oldState;

#ifdef HALF_MAKES
		vector<unsigned> invertTable(inst.O);
#endif

		bool emptyNeighbourhood = false;

		vector<unsigned> jobBb;
		jobBb.reserve(inst.O);
		vector<unsigned> machBb;
		machBb.reserve(inst.O);
		vector<unsigned> critic;
		critic.reserve(inst.O);

		bool makesCycleDetected;
		vector<int> oldValues(maxD, -1);
		unsigned posCurMakes = 0;
		unsigned cycleLastPos = 0;
		unsigned cycleL = 0;
		const unsigned maxLen = maxD * maxC;

		vector<tuple<unsigned, unsigned, unsigned, unsigned>> cands;
		cands.reserve((inst.O)^2);
		TabuTrioN6 tt;

		unsigned tabuIter =0;
		unsigned noImproveIters = 0;
		unsigned curJumpLimit = initialjumpLimit;
		bool jumped = false;
		bool newBestFound = true;

		// 1 initiatizations
		TabuList oldTabuList;
		TabuList tabuList(tenure);
		//unsigned lowerBound = inst.lowerBoundTkz(indeg, Q);
		JumpListN6 jumpList(bjSize);

		// 2 main loop
		while(maxMillisecs > duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count()) {
			tabuIter++;

			noImproveIters++;

			//STEP Checking for cycles
			makesCycleDetected = State::detectRepeat(oldValues, posCurMakes, cycleLastPos, cycleL, curState.makes, newBestFound, maxD, maxLen);

			// STEP get trio from jump list or compute cands
			if((curJumpLimit < noImproveIters)   ||    makesCycleDetected   ||  emptyNeighbourhood) {//JUMPING
				assert( ! newBestFound);

				jumped = true;
				//if(curJumpLimit > jumpLimitDecrease) 
				//curJumpLimit -= jumpLimitDecrease;
				curJumpLimit -= (curJumpLimit/decreaseDivisor);
				//cout << curJumpLimit << endl;

				emptyNeighbourhood = false;
				
				noImproveIters = 0;

				tt = jumpList.pop();

				if( ! tt.dummy) {					
					curState = tt.state;
					cands = tt.cands;
					tabuList = tt.tabuList;
					assert( ! cands.empty());
				} else { //Stop -> no more states in jump list
					return false;
				}
			} else { //NOT JUMPING
				jumped = false;
				cands.clear();

				State::computeCritic(critic, lastOp, prev);
#ifndef NDEBUG
				assert( ! critic.empty());
				assert(critic.size() < inst.O);
				testMakes = 0;
				for(unsigned pos=0; pos<critic.size(); pos++)
					testMakes += inst.P[critic[pos]];
				assert(testMakes == curState.makes);
				for(unsigned pos=0; pos<critic.size()-1; pos++)
					assert(curState.job[critic[pos]] == critic[pos+1]    ||    curState.mach[critic[pos]] == critic[pos+1]);
#endif
				curState.compTailsComplete(tails, indeg, Q);
				curState.compHeadsComplete(heads, indeg, Q);
				curState.fillCandidatesN6(cands, jobBb, machBb, critic, heads,tails);
			}
			assert( ! cands.empty());

			if(newBestFound) {
				assert( ! jumped);
				oldState = curState;
				oldTabuList = tabuList;
			}

			//STEP Go to neighbour
			//nspN6(curState, lastOp, tabuList, cands, theState.makes, dists, prev, indeg, Q, heads, tails);
			if( ! cands.empty())
				nspN6(curState, lastOp, tabuList, cands, theState.makes, dists, prev, indeg, Q, heads, tails);
			else
				emptyNeighbourhood = true;

			//STEP last iter new best was found? store oldState with paths not taken
			if(newBestFound) {
				if( ! cands.empty())
					jumpList.push(TabuTrioN6(oldState, cands, oldTabuList));
				newBestFound = false;
			}

			if(jumped) {
				assert( ! newBestFound);
				jumpList.updateCands(cands);
			}

			//STEP new best ??
			if(curState.makes < theState.makes) {

				if(curState.makes < printWhenBetter)
					if(timeLog) {
						resultList.push_back(preString + unsigStr(curState.makes) + " " + doubleStr((duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())/1000.0) + " d");
						//cout << preString << curState.makes << " " << (duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())/1000.0 << " d" << endl;
					}

				noImproveIters = 0;
				theState = curState;
				curJumpLimit = initialjumpLimit;
				assert(curState.makes >= lowerBound);
				if(curState.makes == lowerBound) //Stop -> Optimal
					return true;

				newBestFound= true;

				theState.millisecsFound=duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
				curState.millisecsFound=duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
			}
		}//end while
		assert(maxMillisecs <= duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count());

		return false;
	}



	
	
	void nsp(State & theState, unsigned & lastOp, TabuList & tabuList, vector<pair<unsigned,unsigned>> & cands, unsigned aspiration, vector<unsigned> & dists, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, vector<unsigned> & heads, vector<unsigned> & tails, bool lbOrder) {
#ifndef NDEBUG
		assert( ! cands.empty());
		assert(dists.size() == inst.O);
		assert(prev.size() == inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);
		State testState = theState;
#endif
		bool cycle;
		unsigned o1;
		unsigned o2;
		vector<bool> stillToTestCandPos(cands.size(), true);
		vector<pair<unsigned, unsigned>> candVisitOrder(cands.size());

		unsigned pPos;

		//Unforbidden or aspiration here
		unsigned chosenMakes = UINT_MAX;
	    unsigned chosenO1 = 0;
		unsigned chosenO2 = 0;
		unsigned chosenSwapPos;

		//Forbidden Non Profitable
		int fnpStateAge = -1;
		unsigned oldestO1 = 0;
		unsigned oldestO2 = 0;
		unsigned oldestFnpSwapPos;

#ifndef NDEBUG
		cycle = 
#endif
			theState.compTailsComplete(tails, indeg, Q);
		//theState.compTails(tails, indeg, Q);
		assert( ! cycle);

#ifndef NDEBUG
		cycle = 
#endif
			theState.compHeadsComplete(heads, indeg, Q);
		assert( ! cycle);

		//preparing candVisitOrder
		for(unsigned curCandIndex=0; curCandIndex<cands.size(); curCandIndex++) {
			candVisitOrder[curCandIndex].first = theState.swapLowerBound(cands[curCandIndex].first, cands[curCandIndex].second, heads, tails);
			candVisitOrder[curCandIndex].second = curCandIndex;
		}
		if(lbOrder)
			sort(candVisitOrder.begin(), candVisitOrder.end());

		//unforbidden moves
		//for(unsigned pPos=0; pPos<cands.size(); pPos++) {
		for(unsigned curVisitOrderI=0; curVisitOrderI<candVisitOrder.size(); curVisitOrderI++) {

			if(lbOrder    &&    candVisitOrder[curVisitOrderI].first >= chosenMakes)
				break;

			pPos = candVisitOrder[curVisitOrderI].second;
			// changed 
			o1 = cands[pPos].first;
			o2 = cands[pPos].second;
			assert(o1 != 0);
			assert(o2 != 0);

			if(tabuList.canMove(o1, o2)) {
				stillToTestCandPos[pPos] = false;

				theState.swap(o1, o2);
				assert(theState.verify());

//#ifndef NDEBUG
				cycle = 
//#endif
					theState.setMeta(dists, lastOp, prev, indeg, Q);
				//assert( ! cycle);



				if(theState.penalties < chosenMakes && !cycle) { //PENAL //alt makes
					chosenMakes = theState.penalties;
					chosenO1 = o1;
					chosenO2 = o2;
					chosenSwapPos = pPos;
				}
				//undoing the swap
				theState.swap(o2, o1);
				assert(theState.verify());
				assert(theState == testState);
			}
		}

		//forbidden moves
		//for(unsigned pPos=0; pPos<cands.size(); pPos++) {
		for(unsigned curVisitOrderI=0; curVisitOrderI<candVisitOrder.size(); curVisitOrderI++) {

			if(lbOrder  &&  candVisitOrder[curVisitOrderI].first >= chosenMakes)
				break;

			pPos = candVisitOrder[curVisitOrderI].second;

			if(stillToTestCandPos[pPos]) {
				o1 = cands[pPos].first;
				o2 = cands[pPos].second;
				assert(o1 != 0);
				assert(o2 != 0);

				assert(! tabuList.canMove(o1, o2));

				theState.swap(o1, o2);
				assert(theState.verify());


//#ifndef NDEBUG
				cycle = 
//#endif
					theState.setMeta(dists, lastOp, prev, indeg, Q);
				//assert( ! cycle);


				//aspiration
				if(theState.penalties < aspiration   &&   theState.penalties < chosenMakes && !cycle) {
					chosenMakes = theState.penalties;
					chosenO1 = o1;
					chosenO2 = o2;
					chosenSwapPos = pPos;
				}

				//crap moves
				if(chosenMakes==UINT_MAX   &&   fnpStateAge < tabuList.age(o1, o2)) {
				    oldestO1 = o1;
					oldestO2 = o2;
					oldestFnpSwapPos = pPos;
					fnpStateAge = tabuList.age(o1, o2);
				}

				theState.swap(o2, o1);
				assert(theState.verify());
				assert(theState == testState);
			}
		}

		//perform move - remove used cand - update tabu list
		if(chosenMakes != UINT_MAX) {
			//good move found
			assert(chosenO1 != 0   &&   chosenO1 < inst.O);
			assert(chosenO2 != 0   &&   chosenO2 < inst.O);
			assert(chosenSwapPos < cands.size());
			theState.swap(chosenO1, chosenO2);
			tabuList.insert(chosenO1, chosenO2);
			cands.erase(cands.begin() + chosenSwapPos);
		} else {
			//only crap moves
			assert(chosenO1 == 0);
			assert(chosenO2 == 0);
			assert(oldestO1 != 0   &&   oldestO1 < inst.O);
			assert(oldestO2 != 0   &&   oldestO2 < inst.O);
			assert(fnpStateAge >= 0);
			theState.swap(oldestO1, oldestO2);
			tabuList.passTime(tabuList.timeToLeave(oldestO1, oldestO2));
			tabuList.insert(oldestO1, oldestO2);
			cands.erase(cands.begin() + oldestFnpSwapPos);
		}
#ifndef NDEBUG
		cycle = 
#endif
			theState.setMeta(dists, lastOp, prev, indeg, Q);
		assert( ! cycle);
	}




	//@return: is optimal?
	//printWhenBetter use 0 to not print. will print preString makes seconds (according to tpSTart) d 
	bool evolveTabu(State & theState, unsigned tenure, unsigned initialjumpLimit, unsigned decreaseDivisor, unsigned bjSize, const high_resolution_clock::time_point & tpStart, vector<unsigned> & dists, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, const unsigned maxD, const unsigned maxC, const unsigned lowerBound, unsigned maxMillisecs, unsigned printWhenBetter, const string & preString, vector<unsigned> & heads, vector<unsigned> & tails, bool timeLog, bool lbOrder)  {
#ifndef NDEBUG
		unsigned testMakes;
		bool cycle;
#endif
		unsigned lastOp;
#ifndef NDEBUG
		cycle = 
#endif
		theState.setMeta(dists, lastOp, prev, indeg, Q);
		assert( ! cycle);
		assert(theState.penalties >= lowerBound);
		//PENAL assert


		if(theState.penalties < printWhenBetter) { //PENAL alt makes
			if(timeLog) {
				resultList.push_back(preString + unsigStr(theState.makes) + " " +doubleStr((duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())/1000.0) + " d");
				//cout << preString << theState.makes << " " << (duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())/1000.0 << " d" << endl;
			}
		}

		if(theState.penalties == lowerBound) //PENAL alt makes
			return true;

		State curState = theState;
		State oldState;

#ifdef HALF_MAKES
		vector<unsigned> invertTable(inst.O);
#endif

		vector<unsigned> jobBb;
		jobBb.reserve(inst.O);
		vector<unsigned> machBb;
		machBb.reserve(inst.O);
		vector<unsigned> critic;
		critic.reserve(inst.O);

		bool makesCycleDetected;
		vector<int> oldValues(maxD, -1);
		unsigned posCurMakes = 0;
		unsigned cycleLastPos = 0;
		unsigned cycleL = 0;
		const unsigned maxLen = maxD * maxC;

		vector<pair<unsigned, unsigned>> cands;
		cands.reserve(inst.O);
		TabuTrio tt;

		unsigned tabuIter =0;
		unsigned noImproveIters = 0;
		unsigned curJumpLimit = initialjumpLimit;
		bool jumped = false;
		bool newBestFound = true;

		bool emptyNeighbourhood = false;

		// 1 initiatizations
		TabuList oldTabuList;
		TabuList tabuList(tenure);
		//unsigned lowerBound = inst.lowerBoundTkz(indeg, Q);
		JumpList jumpList(bjSize);

		// 2 main loop
		while(maxMillisecs > duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count()) {
			tabuIter++;

			noImproveIters++;

			//STEP Checking for cycles
			makesCycleDetected = State::detectRepeat(oldValues, posCurMakes, cycleLastPos, cycleL, curState.penalties, newBestFound, maxD, maxLen); //PENAL alt makes

			// STEP get trio from jump list or compute cands
			if((curJumpLimit < noImproveIters)   ||    makesCycleDetected   ||   emptyNeighbourhood) {//JUMPING
				assert( ! newBestFound);

				jumped = true;
				//if(curJumpLimit > jumpLimitDecrease) 
				//curJumpLimit -= jumpLimitDecrease;
				curJumpLimit -= (curJumpLimit/decreaseDivisor);
				//cout << curJumpLimit << endl;

				emptyNeighbourhood = false;
				
				noImproveIters = 0;

				tt = jumpList.pop();

				if( ! tt.dummy) {					
					curState = tt.state;
					cands = tt.cands;
					tabuList = tt.tabuList;
					assert( ! cands.empty());
				} else { //Stop -> no more states in jump list
					return false;
				}
			} else { //NOT JUMPING
				jumped = false;
				cands.clear();

				/*State::computeCritic(critic, lastOp, prev);
#ifndef NDEBUG
				assert( ! critic.empty());
				assert(critic.size() < inst.O);
				testMakes = 0;
				for(unsigned pos=0; pos<critic.size(); pos++)
					testMakes += inst.P[critic[pos]];
				assert(testMakes == curState.makes);
				for(unsigned pos=0; pos<critic.size()-1; pos++)
					assert(curState.job[critic[pos]] == critic[pos+1]    ||    curState.mach[critic[pos]] == critic[pos+1]);
#endif
				State::fillCandidatesN5(cands, jobBb, machBb, critic);*/

				State::fillCandidatesTest1(cands, curState.mach, lastOp);
			}
			assert( ! cands.empty());

			if(newBestFound) {
				assert( ! jumped);
				oldState = curState;
				oldTabuList = tabuList;
			}

#ifdef COUNT_STEPS
			numSteps++;
			currentShownSeconds = (duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())/1000;
			if(currentShownSeconds > shownSeconds) {
				shownSeconds = currentShownSeconds;
				cout << "TIME" << preString << " " << shownSeconds << " " << numSteps << endl;
			}
#endif //COUNT_STEPS

			//STEP Go to neighbour
			//nsp(curState, lastOp, tabuList, cands, theState.makes, dists, prev, indeg, Q, heads, tails, lbOrder);
			
			if( ! cands.empty())
				
 				nsp(curState, lastOp, tabuList, cands, theState.penalties, dists, prev, indeg, Q, heads, tails, lbOrder);
			else
				emptyNeighbourhood = true;

			//STEP last iter new best was found? store oldState with paths not taken
			if(newBestFound) {
				if( ! cands.empty())
					jumpList.push(TabuTrio(oldState, cands, oldTabuList));
				newBestFound = false;
			}

			if(jumped) {
				assert( ! newBestFound);
				jumpList.updateCands(cands);
			}

			//STEP new best ??
			if(curState.penalties < theState.penalties) {

				if(curState.penalties < printWhenBetter)
					if(timeLog) {
						resultList.push_back(preString + unsigStr(curState.penalties) + " " + doubleStr((duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())/1000.0) + " d");
						//cout << preString << curState.makes << " " << (duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())/1000.0 << " d" << endl;
					}

				noImproveIters = 0;
				theState = curState;
				curJumpLimit = initialjumpLimit;
				assert(curState.penalties >= lowerBound);
				if(curState.penalties == lowerBound) //Stop -> Optimal
					return true;

				newBestFound= true;

				theState.millisecsFound=duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
				curState.millisecsFound=duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
			}
		}//end while
		assert(maxMillisecs <= duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count());

		return false;
	}
	



	void nosmuTabu(const string & instPath, const string & name, unsigned tenure, unsigned initialjumpLimit, unsigned jumpLimitDecrease, unsigned bjSize, unsigned maxD, unsigned maxC, unsigned startType, unsigned maxMillisecs, bool timeLog) {
		high_resolution_clock::time_point tpStart = high_resolution_clock::now();

		inst.parse(instPath);

		vector<unsigned> dists(inst.O);
		vector<unsigned> prev(inst.O);
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		vector<unsigned> heads(inst.O); 
		vector<unsigned> tails(inst.O);
		vector<unsigned> invertedQ(inst.O);
		vector<bool> reach(inst.O);


		State theState;
		theState.alloc();
		bool bigJobFirst = false;

		unsigned lowerBound = inst.lowerBoundTkz(indeg, Q);

		if(startType==INSA_START)
			theState.insaPsp(bigJobFirst, indeg, Q, heads, tails, invertedQ, reach);
		else if(startType==RAND_START)
			theState.makeRand(indeg, Q);
		else
			throw errorText("Invalid start option.", "Tabu.hpp","nosmuTabu");

		unsigned lastOp;
		theState.setMeta(dists, lastOp, prev, indeg, Q);

		//unsigned startMakes = theState.makes;
		//unsigned startMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
		//theState.millisecsFound = startMillisecs;

		theState.millisecsFound = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();;

		// changed lower bound to 0
		evolveTabu(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, 0, maxMillisecs, UINT_MAX, "", heads, tails, timeLog, false);
		//evolveTabu(newState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, bestState.makes, preString, heads, tails, timeLog, lowerBoundOrder);
		

		//unsigned totalMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		//Name    Seed         InitMakes                      InitTime                                        Makes                                         TimeToTarget                                            TotalTime
		//cout << name << " " << initialjumpLimit << " " << startMakes << " " << startMillisecs/1000.0 << " " << theState.makes << " " <<  theState.millisecsFound/1000.0 <<  " " <<   totalMillisecs/1000.0 << endl;
		/*for(string s : resultList) 
			cout << name << " " << initialjumpLimit << " " << s << endl;*/
		
		if( ! theState.verifySchedule())
			throw errorText(theState.toString() + "\n\t\tBad schedule !!!!!","","");

		cout << instPath << " " << lowerBound << " ";
		theState.printPenaltys();
		
		// called here
	}






#ifdef PARALLEL_MACHS


	//******************************
	//PARALLEL VVVVVVVVVVVVVV
	//******************************
	bool nspParallel(State & theState, unsigned & lastOp, TabuList & tabuList, vector<pair<unsigned,unsigned>> & cands, unsigned aspiration, vector<unsigned> & dists, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & bucketLimits, multi_array<unsigned, 2> & bucketTops) {
#ifndef NDEBUG
		assert( ! cands.empty());
		assert(dists.size() == inst.O);
		assert(prev.size() == inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);
		State testState = theState;
#endif
		unsigned o1;
		unsigned o2;
		unsigned x;
		vector<bool> stillToTestCandPos(cands.size(), true);

		//Unforbidden or aspiration here
		unsigned chosenMakes = UINT_MAX;
	    unsigned chosenO1 = 0;
		unsigned chosenO2 = 0;
		unsigned chosenSwapPos;

		//Forbidden Non Profitable
		int fnpStateAge = -1;
		unsigned oldestO1 = 0;
		unsigned oldestO2 = 0;
		unsigned oldestFnpSwapPos;

		//XXX
		//cout << "called NSP parallel" <<endl; 
		//XXX

		//cout << theState.easyString() << endl;

		//unforbidden moves
		for(unsigned pPos=0; pPos<cands.size(); pPos++) {

			//cout << "cand: " <<  cands[pPos].first << " - " <<  cands[pPos].second << endl;

			o1 = cands[pPos].first;
			o2 = cands[pPos].second;
			assert(o1 != 0);
			assert(o2 != 0);

			if(tabuList.canMove(o1, o2)) {
				stillToTestCandPos[pPos] = false;
				//theState.swap(o1, o2);
				x = theState._mach[o2];
				assert(x != 0);
				theState.shiftBack(o1,o2);
				assert(theState.verify());

				theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);

				//XXX
				//cout << "\tunforb: " << theState.makes <<  endl; 
				//XXX

				if(theState.makes < chosenMakes) {

					//XXX
					//cout << "\t\tis best: " << theState.makes << endl; 
					//XXX

					chosenMakes = theState.makes;
					chosenO1 = o1;
					chosenO2 = o2;
					chosenSwapPos = pPos;
				}
				//undoing the swap
				//theState.swap(o2, o1);
				theState.shiftForth(o2, x);
				assert(theState.verify());
				assert(theState == testState);
			}
		}//unforbinden for 

		//forbidden moves
		for(unsigned pPos=0; pPos<cands.size(); pPos++) {
			if(stillToTestCandPos[pPos]) {
				o1 = cands[pPos].first;
				o2 = cands[pPos].second;
				assert(o1 != 0);
				assert(o2 != 0);

				assert(! tabuList.canMove(o1, o2));

				//theState.swap(o1, o2);
				x = theState._mach[o2];
				assert(x != 0);
				theState.shiftBack(o1,o2);
				assert(theState.verify());
				theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);

				//XXX
				//cout << "\tforb: " << theState.makes <<  endl; 
				//XXX

				//aspiration
				if(theState.makes < aspiration   &&   theState.makes < chosenMakes) {
					//XXX
					//cout << "\t\taspire: " << theState.makes <<  endl; 
					//XXX

					chosenMakes = theState.makes;
					chosenO1 = o1;
					chosenO2 = o2;
					chosenSwapPos = pPos;
				}

				//crap moves
				if(chosenMakes==UINT_MAX   &&   fnpStateAge < tabuList.age(o1, o2)) {
					//XXX
					//cout << "\t\tcrap: "  <<endl; 
					//XXX
				    oldestO1 = o1;
					oldestO2 = o2;
					oldestFnpSwapPos = pPos;
					fnpStateAge = tabuList.age(o1, o2);
				}

				//theState.swap(o2, o1);
				theState.shiftForth(o2, x);
				assert(theState.verify());
				assert(theState == testState);
			}//still to test if
		}//forbiden for		

		//perform move - remove used cand - update tabu list
		if(chosenMakes != UINT_MAX) {
			//XXX
			//cout << "chose best makes: " << chosenMakes <<  endl; 
			//XXX
			//good move found
			assert(chosenO1 != 0   &&   chosenO1 < inst.O);
			assert(chosenO2 != 0   &&   chosenO2 < inst.O);
			assert(chosenSwapPos < cands.size());
			//theState.swap(chosenO1, chosenO2);
			theState.shiftBack(chosenO1,chosenO2);
			tabuList.insert(chosenO1, chosenO2);
			cands.erase(cands.begin() + chosenSwapPos);
		} else {
			//only crap moves
			//XXX
			//cout << "crap move chosen: " << chosenMakes <<  endl; 
			//XXX
			assert(chosenO1 == 0);
			assert(chosenO2 == 0);

			if(oldestO1 == 0) {
				cands.clear();
				assert(oldestO2 ==0);
				theState.makes = UINT_MAX;
				return true;
			}

			assert(oldestO1 != 0   &&   oldestO1 < inst.O);
			assert(oldestO2 != 0   &&   oldestO2 < inst.O);
			assert(fnpStateAge >= 0);
			//theState.swap(oldestO1, oldestO2);
			theState.shiftBack(oldestO1,oldestO2);
			tabuList.passTime(tabuList.timeToLeave(oldestO1, oldestO2));
			tabuList.insert(oldestO1, oldestO2);
			cands.erase(cands.begin() + oldestFnpSwapPos);
		}

		theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);


		//cout << "final makes: " << theState.makes << endl;

		assert(theState.makes < UINT_MAX);

		return false;
	}







//printWhenBetter use 0 to not print. will print preString makes seconds (according to tpSTart) d 
	void evolveTabuParallel(State & theState, unsigned tenure, unsigned initialjumpLimit, unsigned jumpLimitDecrease, unsigned bjSize, const high_resolution_clock::time_point & tpStart, vector<unsigned> & dists, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & bucketLimits, multi_array<unsigned, 2> & bucketTops, const unsigned maxD, const unsigned maxC,  double maxMillisecs, unsigned printWhenBetter, const string & preString)  {

		unsigned lastOp;

		theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);

		if(theState.makes < printWhenBetter)
			cout << preString << theState.makes << " " << (duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())/1000.0 << " d" << endl;

		State curState = theState;
		State oldState;

		vector<unsigned> jobBb;
		jobBb.reserve(inst.O);
		vector<unsigned> machBb;
		machBb.reserve(inst.O);
		vector<unsigned> critic;
		critic.reserve(inst.O);

		bool makesCycleDetected;
		vector<int> oldValues(maxD, -1);
		unsigned posCurMakes = 0;
		unsigned cycleLastPos = 0;
		unsigned cycleL = 0;
		const unsigned maxLen = maxD * maxC;

		vector<pair<unsigned, unsigned>> cands;
		cands.reserve(inst.O);
		TabuTrio tt;

		unsigned tabuIter =0;
		unsigned noImproveIters = 0;
		unsigned curJumpLimit = initialjumpLimit;
		bool jumped = false;
		bool newBestFound = true;

		bool noViable = false;

		// 1 initiatizations
		TabuList oldTabuList;
		TabuList tabuList(tenure);
		JumpList jumpList(bjSize);

		// 2 main loop
		while(maxMillisecs > duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count()) {
			tabuIter++;

			noImproveIters++;

			//STEP Checking for cycles
			makesCycleDetected = State::detectRepeat(oldValues, posCurMakes, cycleLastPos, cycleL, curState.makes, newBestFound, maxD, maxLen);

			// STEP get trio from jump list or compute cands
			if((curJumpLimit < noImproveIters)   ||    makesCycleDetected || noViable) {//JUMPING
				assert( ! newBestFound);

				noViable = false;

				jumped = true;
				if(curJumpLimit > jumpLimitDecrease) 
					curJumpLimit -= jumpLimitDecrease;
				noImproveIters = 0;

				tt = jumpList.pop();

				if( ! tt.dummy) {					
					curState = tt.state;
					cands = tt.cands;
					tabuList = tt.tabuList;
					assert( ! cands.empty());
				} else { //Stop -> no more states in jump list
					return;
				}
			} else { //NOT JUMPING
				jumped = false;
				cands.clear();

				State::computeCritic(critic, lastOp, prev);
				//State::fillCandidatesN5(cands, jobBb, machBb, critic);
				State::fillCandAllInCritic(cands, critic);
				//curState.fillCandAllMachs(cands);
				if(cands.size() <=1)
					return;
			}
			assert( ! cands.empty());

			if(newBestFound) {
				assert( ! jumped);
				oldState = curState;
				oldTabuList = tabuList;
			}
			
			//STEP Go to neighbour
			noViable = nspParallel(curState, lastOp, tabuList, cands, theState.makes, dists, prev, indeg, Q, bucketLimits, bucketTops);
			//noViable = nspParallelGaps(curState, lastOp, tabuList, cands, theState.makes, dists, prev, indeg, Q, allGaps);

			//STEP last iter new best was found? store oldState with paths not taken
			if(newBestFound) {
				if( ! cands.empty())
					jumpList.push(TabuTrio(oldState, cands, oldTabuList));
				newBestFound = false;
			}

			if(jumped) {
				assert( ! newBestFound);
				jumpList.updateCands(cands);
			}

			//STEP new best ??
			if(curState.makes < theState.makes) {

				if(curState.makes < printWhenBetter)
					cout << preString << curState.makes << " " << (duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())/1000.0 << " d" << endl;

				noImproveIters = 0;
				theState = curState;
				curJumpLimit = initialjumpLimit;

				newBestFound= true;

				theState.millisecsFound=duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
				curState.millisecsFound=duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
			}
		} //end while
		assert(maxMillisecs <= duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count());
		return;
	}

	//******************************
	//PARALLEL^^^^^^^^^^^
	//******************************




    void nosmuTabuParallel(const string & instPath, unsigned tenure, unsigned initialjumpLimit, unsigned jumpLimitDecrease, unsigned bjSize, unsigned maxD, unsigned maxC, double maxMillisecs) {
		high_resolution_clock::time_point tpStart = high_resolution_clock::now();

		inst.parse(instPath);

		vector<unsigned> dists(inst.O);
		vector<unsigned> prev(inst.O);
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		multi_array<unsigned, 2> bucketLimits;
		bucketLimits.resize(boost::extents[inst.M][inst.K]);
		multi_array<unsigned, 2> bucketTops;
		bucketTops.resize(boost::extents[inst.M][inst.K]);
		vector<vector<vector<pair<unsigned, unsigned>>>> allGaps(inst.M, vector<vector<pair<unsigned, unsigned>>>(inst.K));


		State theState;
		theState.alloc();

		theState.insaParallel(indeg, Q, dists, bucketLimits);


		unsigned lastOp;
		theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);



		//unsigned startMakes = theState.makes;
		unsigned startMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		theState.millisecsFound = startMillisecs;

		//evolveTabuParallel(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, buckets, maxD, maxC, maxMillisecs, theState.makes, "");
		evolveTabuParallel(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, bucketLimits, bucketTops, maxD, maxC, maxMillisecs, 0, "");
		//evolveTabuParallelGaps(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, allGaps, maxD, maxC, maxMillisecs, 0, "");

		theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);

		unsigned totalMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		//Name           Makes                                         TimeToTarget                                            TotalTime
		cout << theState.makes << " " <<  theState.millisecsFound/1000.0 <<  " " <<   totalMillisecs/1000.0 << endl;
		//cout << theState.makes;
		
		inst.verifyPrallelSchedule(dists, theState.makes);
	}

#endif //PARALLEL_MACHS























#ifdef VAR_PARALLEL_MACHS


	//******************************
	//PARALLEL VVVVVVVVVVVVVV
	//******************************
	bool nspParallel(State & theState, unsigned & lastOp, TabuList & tabuList, vector<pair<unsigned,unsigned>> & cands, unsigned aspiration, vector<unsigned> & dists, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & bucketLimits, multi_array<unsigned, 2> & bucketTops) {
#ifndef NDEBUG
		assert( ! cands.empty());
		assert(dists.size() == inst.O);
		assert(prev.size() == inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);
		State testState = theState;
#endif
		unsigned o1;
		unsigned o2;
		unsigned x;
		vector<bool> stillToTestCandPos(cands.size(), true);

		//Unforbidden or aspiration here
		unsigned chosenMakes = UINT_MAX;
	    unsigned chosenO1 = 0;
		unsigned chosenO2 = 0;
		unsigned chosenSwapPos;

		//Forbidden Non Profitable
		int fnpStateAge = -1;
		unsigned oldestO1 = 0;
		unsigned oldestO2 = 0;
		unsigned oldestFnpSwapPos;

		//unforbidden moves
		for(unsigned pPos=0; pPos<cands.size(); pPos++) {

			o1 = cands[pPos].first;
			o2 = cands[pPos].second;
			assert(o1 != 0);
			assert(o2 != 0);

			if(tabuList.canMove(o1, o2)) {
				stillToTestCandPos[pPos] = false;
				//theState.swap(o1, o2);
				x = theState._mach[o2];
				assert(x != 0);
				theState.shiftBack(o1,o2);
				assert(theState.verify());

				theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);

				if(theState.makes < chosenMakes) {

					chosenMakes = theState.makes;
					chosenO1 = o1;
					chosenO2 = o2;
					chosenSwapPos = pPos;
				}
				//undoing the swap
				//theState.swap(o2, o1);
				theState.shiftForth(o2, x);
				assert(theState.verify());
				assert(theState == testState);
			}
		}//unforbinden for 

		//forbidden moves
		for(unsigned pPos=0; pPos<cands.size(); pPos++) {
			if(stillToTestCandPos[pPos]) {
				o1 = cands[pPos].first;
				o2 = cands[pPos].second;
				assert(o1 != 0);
				assert(o2 != 0);

				assert(! tabuList.canMove(o1, o2));

				//theState.swap(o1, o2);
				x = theState._mach[o2];
				assert(x != 0);
				theState.shiftBack(o1,o2);
				assert(theState.verify());
				theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);

				//aspiration
				if(theState.makes < aspiration   &&   theState.makes < chosenMakes) {

					chosenMakes = theState.makes;
					chosenO1 = o1;
					chosenO2 = o2;
					chosenSwapPos = pPos;
				}

				//crap moves
				if(chosenMakes==UINT_MAX   &&   fnpStateAge < tabuList.age(o1, o2)) {
				    oldestO1 = o1;
					oldestO2 = o2;
					oldestFnpSwapPos = pPos;
					fnpStateAge = tabuList.age(o1, o2);
				}

				//theState.swap(o2, o1);
				theState.shiftForth(o2, x);
				assert(theState.verify());
				assert(theState == testState);
			}//still to test if
		}//forbiden for		

		//perform move - remove used cand - update tabu list
		if(chosenMakes != UINT_MAX) {
			assert(chosenO1 != 0   &&   chosenO1 < inst.O);
			assert(chosenO2 != 0   &&   chosenO2 < inst.O);
			assert(chosenSwapPos < cands.size());
			
			//theState.swap(chosenO1, chosenO2);
			theState.shiftBack(chosenO1,chosenO2);
			tabuList.insert(chosenO1, chosenO2);
			cands.erase(cands.begin() + chosenSwapPos);
			
			
		} else {
			assert(chosenO1 == 0);
			assert(chosenO2 == 0);
			


			if(oldestO1 == 0) {
				cands.clear();
				assert(oldestO2 ==0);
				theState.makes = UINT_MAX;
				return true;
			}

			assert(oldestO1 != 0   &&   oldestO1 < inst.O);
			assert(oldestO2 != 0   &&   oldestO2 < inst.O);
			assert(fnpStateAge >= 0);
			//theState.swap(oldestO1, oldestO2);
			theState.shiftBack(oldestO1,oldestO2);
			tabuList.passTime(tabuList.timeToLeave(oldestO1, oldestO2));
			tabuList.insert(oldestO1, oldestO2);
			cands.erase(cands.begin() + oldestFnpSwapPos);
		}

		theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);


		//cout << "final makes: " << theState.makes << endl;

		assert(theState.makes < UINT_MAX);

		return false;
	}
	
	
	
	
#ifdef ANALYZE_MODE
		bool nspParallelAllDelay(State & theState, unsigned & lastOp, TabuList & tabuList, vector<pair<unsigned,unsigned>> & cands, unsigned aspiration, vector<unsigned> & dists, vector<unsigned> & prev, vector<vector<unsigned>> & allDelay, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & bucketLimits, multi_array<unsigned, 2> & bucketTops, vector<OperPairQuality> & qualityV, double trimSize) 
#else
		bool nspParallelAllDelay(State & theState, unsigned & lastOp, TabuList & tabuList, vector<pair<unsigned,unsigned>> & cands, unsigned aspiration, vector<unsigned> & dists, vector<unsigned> & prev, vector<vector<unsigned>> & allDelay, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & bucketLimits, multi_array<unsigned, 2> & bucketTops) 
#endif
		{

#ifndef NDEBUG
		assert( ! cands.empty());
		assert(dists.size() == inst.O);
		assert(prev.size() == inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);
		State testState = theState;
#endif

#ifdef ANALYZE_MODE
		qualityV.clear();
		const unsigned neighCut = (unsigned)(trimSize*(cands.size())+1.0);
		//XXX
		//cout << "\t\t\tneighCut: " << neighCut << endl;
		unsigned evaluated = 0; 
		bool stop = false;
#endif

		unsigned o1;
		unsigned o2;
		unsigned x;
		vector<bool> stillToTestCandPos(cands.size(), true);

		//Unforbidden or aspiration here
		unsigned chosenMakes = UINT_MAX;
	    unsigned chosenO1 = 0;
		unsigned chosenO2 = 0;
		unsigned chosenSwapPos;

		//Forbidden Non Profitable
		int fnpStateAge = -1;
		unsigned oldestO1 = 0;
		unsigned oldestO2 = 0;
		unsigned oldestFnpSwapPos;
		
		//unforbidden moves
		for(unsigned pPos=0; pPos<cands.size(); pPos++) {
#ifdef ANALYZE_MODE
			if(stop)
				break;
#endif //ANALYZE_MODE			
		
			o1 = cands[pPos].first;
			o2 = cands[pPos].second;
			assert(o1 != 0);
			assert(o2 != 0);
			
#ifdef ANALYZE_MODE
			observer.analyzer.allOperData[o1].inNeigh++;
			observer.analyzer.allOperData[o2].inNeigh++;
			observer.analyzer.allOperData[o1].inNeighAsBailiff++;
			observer.analyzer.allOperData[o2].inNeighAsProvost++;
#endif //ANALYZE_MODE

			if(tabuList.canMove(o1, o2)) {
				stillToTestCandPos[pPos] = false;
				//theState.swap(o1, o2);
				x = theState._mach[o2];
				assert(x != 0);
				theState.shiftBack(o1,o2);
				assert(theState.verify());

				theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);
#ifdef ANALYZE_MODE
				qualityV.push_back(OperPairQuality(o1, o2, theState.makes));
				if(theState.makes < UINFTY)
					evaluated++;
				if(evaluated >= neighCut)
					stop = true;
#endif //ANALYZE_MODE
				
				if(theState.makes < chosenMakes) {
					chosenMakes = theState.makes;
					chosenO1 = o1;
					chosenO2 = o2;
					chosenSwapPos = pPos;
				}
				//undoing the swap
				//theState.swap(o2, o1);
				theState.shiftForth(o2, x);
				assert(theState.verify());
				assert(theState == testState);
			}
		}//unforbinden for 

		//forbidden moves
		for(unsigned pPos=0; pPos<cands.size(); pPos++) {
#ifdef ANALYZE_MODE
			if(stop)
				break;
#endif //ANALYZE_MODE

			if(stillToTestCandPos[pPos]) {
				o1 = cands[pPos].first;
				o2 = cands[pPos].second;
				assert(o1 != 0);
				assert(o2 != 0);

				assert(! tabuList.canMove(o1, o2));

				//theState.swap(o1, o2);
				x = theState._mach[o2];
				assert(x != 0);
				theState.shiftBack(o1,o2);
				assert(theState.verify());
				theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);
#ifdef ANALYZE_MODE
				qualityV.push_back(OperPairQuality(o1, o2, theState.makes));
				if(theState.makes < UINFTY)
					evaluated++;
				if(evaluated >= neighCut)
					stop = true;
#endif //ANALYZE_MODE
				
				//aspiration
				if(theState.makes < aspiration   &&   theState.makes < chosenMakes) {
					chosenMakes = theState.makes;
					chosenO1 = o1;
					chosenO2 = o2;
					chosenSwapPos = pPos;
				}

				//crap moves
				if(chosenMakes==UINT_MAX   &&   fnpStateAge < tabuList.age(o1, o2)) {
				    oldestO1 = o1;
					oldestO2 = o2;
					oldestFnpSwapPos = pPos;
					fnpStateAge = tabuList.age(o1, o2);
				}

				//theState.swap(o2, o1);
				theState.shiftForth(o2, x);
				assert(theState.verify());
				assert(theState == testState);
			}//still to test if
		}//forbiden for		

		//perform move - remove used cand - update tabu list
		if(chosenMakes != UINT_MAX) {
			assert(chosenO1 != 0   &&   chosenO1 < inst.O);
			assert(chosenO2 != 0   &&   chosenO2 < inst.O);
			assert(chosenSwapPos < cands.size());
			
			//theState.swap(chosenO1, chosenO2);
			theState.shiftBack(chosenO1,chosenO2);
			tabuList.insert(chosenO1, chosenO2);
			cands.erase(cands.begin() + chosenSwapPos);
			
#ifdef ANALYZE_MODE
			observer.analyzer.allOperData[chosenO1].bestInNSP++;
			observer.analyzer.allOperData[chosenO2].bestInNSP++;
			observer.analyzer.allOperData[chosenO1].bestInNSPAsBailiff++;
			observer.analyzer.allOperData[chosenO2].bestInNSPAsProvost++;
			if(chosenMakes < aspiration) {
				observer.analyzer.allOperData[chosenO1].improvedCurrent++;
				observer.analyzer.allOperData[chosenO2].improvedCurrent++;
				observer.analyzer.allOperData[chosenO1].improvedCurrentAsBailiff++;
				observer.analyzer.allOperData[chosenO2].improvedCurrentAsProvost++;
			}
#endif //ANALYZE_MODE
		} else {
			assert(chosenO1 == 0);
			assert(chosenO2 == 0);

			if(oldestO1 == 0) {
				//cands.clear();
				assert(oldestO2 == 0);
				theState.makes = UINT_MAX;
				assert(testState == theState);
				return true;
			}

			assert(oldestO1 != 0   &&   oldestO1 < inst.O);
			assert(oldestO2 != 0   &&   oldestO2 < inst.O);
			assert(fnpStateAge >= 0);
			//theState.swap(oldestO1, oldestO2);
			theState.shiftBack(oldestO1,oldestO2);
			tabuList.passTime(tabuList.timeToLeave(oldestO1, oldestO2));
			tabuList.insert(oldestO1, oldestO2);
			cands.erase(cands.begin() + oldestFnpSwapPos);
		}

		theState.makes = theState.setMetaParallelAllDelay(dists, lastOp, prev, allDelay, indeg, Q, bucketLimits, bucketTops);

		assert(theState.makes < UINT_MAX);

		return false;
	}







	
//printWhenBetter use 0 to not print. will print preString makes seconds (according to tpSTart) d 	
#ifdef ANALYZE_MODE
	void evolveTabuParallel(State & theState, unsigned tenure, unsigned initialjumpLimit, unsigned jumpLimitDecrease, unsigned bjSize, const high_resolution_clock::time_point & tpStart, vector<unsigned> & dists, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & bucketLimits, multi_array<unsigned, 2> & bucketTops, const unsigned maxD, const unsigned maxC,  double maxMillisecs, unsigned printWhenBetter, const string & preString, bool critVannilaNeigh, bool allDelayNeigh, bool learn) {
#else  
    void evolveTabuParallel(State & theState, unsigned tenure, unsigned initialjumpLimit, unsigned jumpLimitDecrease, unsigned bjSize, const high_resolution_clock::time_point & tpStart, vector<unsigned> & dists, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & bucketLimits, multi_array<unsigned, 2> & bucketTops, const unsigned maxD, const unsigned maxC,  double maxMillisecs, unsigned printWhenBetter, const string & preString, bool critVannilaNeigh, bool allDelayNeigh)  {
#endif //ANALYZE_MODE 

	
#ifdef ANALYZE_MODE
#ifndef NDEBUG
		unsigned candSize;
#endif
		vector<OperPairQuality> qualityV;
		vector<vector<unsigned>> prioActual(inst.O, vector<unsigned>(inst.O, 0));
		vector<vector<unsigned>> prioPredict(inst.O, vector<unsigned>(inst.O, 0));
		vector<vector<unsigned>> R(inst.O, vector<unsigned>(inst.O));
		bool isTrimming = false;//learned enough and is now cutting newighbourhood
		double aKendalW;
#endif //ANALYZE_MODE

		vector<vector<unsigned>> allDelay(inst.O, vector<unsigned>());
		unsigned lastOp;

		theState.makes = theState.setMetaParallelAllDelay(dists, lastOp, prev, allDelay, indeg, Q, bucketLimits, bucketTops);

		if(theState.makes < printWhenBetter)
			cout << preString << theState.makes << " " << (duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())/1000.0 << " d" << endl;

		State curState = theState;
		State oldState;

		vector<unsigned> jobBb;
		jobBb.reserve(inst.O);
		vector<unsigned> machBb;
		machBb.reserve(inst.O);
		vector<unsigned> critic;
		critic.reserve(inst.O);

		bool makesCycleDetected;
		vector<int> oldValues(maxD, -1);
		unsigned posCurMakes = 0;
		unsigned cycleLastPos = 0;
		unsigned cycleL = 0;
		const unsigned maxLen = maxD * maxC;

		vector<pair<unsigned, unsigned>> cands;
		cands.reserve(inst.O);
		TabuTrio tt;

		unsigned tabuIter =0;
		unsigned noImproveIters = 0;
		unsigned curJumpLimit = initialjumpLimit;
		bool jumped = false;
		bool newBestFound = true;

		bool noViable = false;

		// 1 initiatizations
		TabuList oldTabuList;
		TabuList tabuList(tenure);
		JumpList jumpList(bjSize);

		// 2 main loop
		while(maxMillisecs > duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count()) {
			tabuIter++;

			noImproveIters++;

			//STEP Checking for cycles
			makesCycleDetected = State::detectRepeat(oldValues, posCurMakes, cycleLastPos, cycleL, curState.makes, newBestFound, maxD, maxLen);

			// STEP get trio from jump list or compute cands
			if((curJumpLimit < noImproveIters)   ||    makesCycleDetected || noViable) {//JUMPING
				assert( ! newBestFound);

#ifdef ANALYZE_MODE			
				//XXX
				observer.forget();
				isTrimming = false;
#endif //#ifdef ANALYZE_MODE
				

				noViable = false;

				jumped = true;
				if(curJumpLimit > jumpLimitDecrease) 
					curJumpLimit -= jumpLimitDecrease;
				noImproveIters = 0;

				tt = jumpList.pop();

				if( ! tt.dummy) {					
					curState = tt.state;
					cands = tt.cands;
					tabuList = tt.tabuList;
					assert( ! cands.empty());
				} else { //Stop -> no more states in jump list
					return;
				}
			} else { //NOT JUMPING
				jumped = false;
				cands.clear();

				State::computeCritic(critic, lastOp, prev);
				cands.clear();
				if(critVannilaNeigh) {
					State::fillCandAllInCritic(cands, critic);
				}
				if(allDelayNeigh) {
					State::fillCandCritGapFiller(cands, critic, allDelay);
					//State::fillCandFullDelaySweep(cands, allDelay);
					//State::fillCandParallelN5(cands, critic, allDelay);
				}
				if(cands.size() == 0) {
					return;
				}
			}
			assert( ! cands.empty());

			if(newBestFound) {
				assert( ! jumped);
				oldState = curState;
				oldTabuList = tabuList;
			}
			
			

			//STEP Go to neighbour
			if(! allDelayNeigh   &&   critVannilaNeigh) {
				noViable = nspParallel(curState, lastOp, tabuList, cands, theState.makes, dists, prev, indeg, Q, bucketLimits, bucketTops);
			} 
			if(allDelayNeigh) {
#ifdef ANALYZE_MODE
#ifndef NDEBUG
				candSize = cands.size();
#endif
				if(learn) {
					
					
					if(isTrimming) {
										
						observer.analyzer.sortCandsByPredict(cands);
						assert(cands.size() == candSize);
						noViable = nspParallelAllDelay(curState, lastOp, tabuList, cands, theState.makes, dists, prev, allDelay, indeg, Q, bucketLimits, bucketTops,  qualityV, observer.trimmSize);
						if(noViable) {//here was trimming but found nothing - revert to normal neighbourhood and forget observer
							assert(cands.size() == candSize);
							noViable = nspParallelAllDelay(curState, lastOp, tabuList, cands, theState.makes, dists, prev, allDelay, indeg, Q, bucketLimits, bucketTops,  qualityV, 1.0);
							isTrimming = false;
							if(! jumped   &&   noViable) {
								curState.setMetaParallelAllDelay(dists, lastOp, prev, allDelay, indeg, Q, bucketLimits, bucketTops);
								assert(cands.size() == candSize);
								curState.fillExtendedNeigh(cands, critic, allDelay);
								noViable = nspParallelAllDelay(curState, lastOp, tabuList, cands, theState.makes, dists, prev, allDelay, indeg, Q, bucketLimits, bucketTops,  qualityV, 1.0);
								if(noViable)
									cands.clear();
							} else if (noViable) {
								assert(jumped);
								cands.pop_back();
							}
						}
					} else {//not trimming
						observer.analyzer.predictInPosBestNSP(prioPredict, cands);
						//calling nsp normally
						noViable = nspParallelAllDelay(curState, lastOp, tabuList, cands, theState.makes, dists, prev, allDelay, indeg, Q, bucketLimits, bucketTops,  qualityV, 1.0);
						if(! jumped   &&   noViable) {
							curState.setMetaParallelAllDelay(dists, lastOp, prev, allDelay, indeg, Q, bucketLimits, bucketTops);
							curState.fillExtendedNeigh(cands, critic, allDelay);
							noViable = nspParallelAllDelay(curState, lastOp, tabuList, cands, theState.makes, dists, prev, allDelay, indeg, Q, bucketLimits, bucketTops,  qualityV, 1.0);
							if(noViable)
								cands.clear();
						}  else if(noViable) {
							assert(jumped);
							cands.pop_back();
						}
						
						if ( ! noViable ) {//not trimming and found viable - update observer
							//cout << "\tSomethign viable" << endl;
							//computing kendal w
							sort(qualityV.begin(), qualityV.end());
							for(unsigned u=0; u<qualityV.size(); u++)
								prioActual[qualityV[u].bailiff][qualityV[u].provost] = u+1;
							aKendalW = Analyzer::pairKendallW(prioActual, prioPredict, R, qualityV);	
							//update observer
							isTrimming = observer.addKendallW(aKendalW);

							//XXX
							sumKW += aKendalW;
							iterKW ++;
							if(iterKW == 30)
								at30KW = sumKW/30.0;
							//XXX
						}
					}//finished not trimming nsp call and updated staistics
				} else {//not learning
					noViable = nspParallelAllDelay(curState, lastOp, tabuList, cands, theState.makes, dists, prev, allDelay, indeg, Q, bucketLimits, bucketTops,  qualityV, 1.0);
					if(! jumped   &&   noViable) {
						//cout << "\n\t\t\t\tEXTENDED\n";
						curState.setMetaParallelAllDelay(dists, lastOp, prev, allDelay, indeg, Q, bucketLimits, bucketTops);
						curState.fillExtendedNeigh(cands, critic, allDelay);
						noViable = nspParallelAllDelay(curState, lastOp, tabuList, cands, theState.makes, dists, prev, allDelay, indeg, Q, bucketLimits, bucketTops,  qualityV, 1.0);
						if(noViable)
							cands.clear();
					}   else if(noViable) {
						assert(jumped);
						cands.pop_back();
					}
				}
				
#else
				noViable = nspParallelAllDelay(curState, lastOp, tabuList, cands, theState.makes, dists, prev, allDelay, indeg, Q, bucketLimits, bucketTops);
				if(! jumped   &&   noViable) {
					curState.setMetaParallelAllDelay(dists, lastOp, prev, allDelay, indeg, Q, bucketLimits, bucketTops);
					curState.fillExtendedNeigh(cands, critic, allDelay);
				}
#endif //ANALYZE_MODE

			}
			
			//noViable = nspParallelGaps(curState, lastOp, tabuList, cands, theState.makes, dists, prev, indeg, Q, allGaps);
		

			//STEP last iter new best was found? store oldState with paths not taken
			if(newBestFound) {
				if( ! cands.empty())
					jumpList.push(TabuTrio(oldState, cands, oldTabuList));
				newBestFound = false;
			}

			if(jumped) {
				assert( ! newBestFound);
				jumpList.updateCands(cands);
			}

			//STEP new best ??
			if(curState.makes < theState.makes) {

				if(curState.makes < printWhenBetter)
					cout << preString << curState.makes << " " << (duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count())/1000.0 << " d" << endl;

				noImproveIters = 0;
				theState = curState;
				curJumpLimit = initialjumpLimit;

				newBestFound= true;

				theState.millisecsFound=duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
				curState.millisecsFound=duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
			}
		} //end while
		assert(maxMillisecs <= duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count());
		
		
		return;
	}



#ifdef ANALYZE_MODE
    void nosmuTabuParallel(const string & instPath, unsigned tenure, unsigned initialjumpLimit, unsigned jumpLimitDecrease, unsigned bjSize, unsigned maxD, unsigned maxC, double maxMillisecs, bool critVannilaNeigh, bool allDelayNeigh, unsigned startType, bool learn, unsigned minLearnPeriod, double minLearnKW, double trimmSize) 
#else  
    void nosmuTabuParallel(const string & instPath, unsigned tenure, unsigned initialjumpLimit, unsigned jumpLimitDecrease, unsigned bjSize, unsigned maxD, unsigned maxC, double maxMillisecs, bool critVannilaNeigh, bool allDelayNeigh, unsigned startType) 
#endif //ANALYZE_MODE
    {
    
		high_resolution_clock::time_point tpStart = high_resolution_clock::now();

		inst.parse(instPath);
#ifdef ANALYZE_MODE
			observer.set(minLearnPeriod, minLearnKW, trimmSize);
#endif //ANALYZE_MODE

		vector<unsigned> dists(inst.O);
		vector<unsigned> prev(inst.O);
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		multi_array<unsigned, 2> bucketLimits;
		bucketLimits.resize(boost::extents[inst.M][inst.maxK]);
		multi_array<unsigned, 2> bucketTops;
		bucketTops.resize(boost::extents[inst.M][inst.maxK]);
		vector<vector<vector<pair<unsigned, unsigned>>>> allGaps(inst.M, vector<vector<pair<unsigned, unsigned>>>(inst.maxK));
		//unsigned startMakes;


		State theState;
		theState.alloc();

		if(startType == PARALLEL_INSA_START)
			theState.insaParallel(indeg, Q, dists, bucketLimits);
		else if(startType == PARALLEL_JOB_TAIL_START) {

			Schedule sch;
			sch.alloc();
			sch.dispatchLongestJobTail(SMALLER_GAP);
			vector<vector<unsigned>> allMachOrder(inst.M);
			sch.getMachOrders(allMachOrder);
			theState.importFromMachOrder(allMachOrder);
			
		} else {
			throw errorText("invalid start typer for parallel nosmu tabu","","");
		}

		unsigned lastOp;
		theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);
		inst.verifyPrallelSchedule(dists, theState.makes);
		//startMakes = theState.makes;

		//unsigned startMakes = theState.makes;
		unsigned startMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		theState.millisecsFound = startMillisecs;

		//evolveTabuParallel(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, buckets, maxD, maxC, maxMillisecs, theState.makes, "");
#ifdef ANALYZE_MODE
   		evolveTabuParallel(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, bucketLimits, bucketTops, maxD, maxC, maxMillisecs, 0, "", critVannilaNeigh, allDelayNeigh, learn);
#else  
		evolveTabuParallel(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, bucketLimits, bucketTops, maxD, maxC, maxMillisecs, 0, "", critVannilaNeigh, allDelayNeigh);
#endif //ANALYZE_MODE  

		
		//evolveTabuParallelGaps(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, allGaps, maxD, maxC, maxMillisecs, 0, "");

		theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);

		unsigned totalMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		//Name           Makes                                         TimeToTarget                                            TotalTime

		//XXX
		//cout << theState.makes << " " <<  /*startMakes << " " << theState.millisecsFound/1000.0 <<  " " <<*/   totalMillisecs/1000.0 << " " << observer.setMetaTimes << endl;
		//cout << maxKW << endl;
		double x = 0.0;
		if(iterKW > 0)
			x =  sumKW/(double)iterKW;
		cout << at30KW << " " << x << endl;

	
		inst.verifyPrallelSchedule(dists, theState.makes);
	}
	
	//******************************
	//PARALLEL^^^^^^^^^^^
	//******************************

#endif //VAR_PARALLEL_MACHS

}


