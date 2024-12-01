/*
  @author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/
#pragma once

#include "Settings.hpp"
#include "Utilities.hpp"
#include "Instance.hpp"

using namespace boost;

class State {
public:
	State() {  }
	~State() {  }


	bool operator==(const State& other) const;

	void alloc();

	void clear();

	string toString() const;

	//expects fully defined state
	string easyString() const;

	//x and y in same job 
	//max(head(Pm[y]), head(Pj[x]))+procT[X]+procT[Y]+max(tail(Sj]y]), tail(Sm[x]))
	//max(head(Pj[x], head(Pm[y])) +procT[y] + tail(Sm[x])
	//head(Pm[x]) + procT[x] +max(tail(Sj[y], tail(Sm[y])
	unsigned swapLowerBound(const unsigned x, const unsigned y, const vector<unsigned>& heads, const vector<unsigned>& tails) const;

	/*void swap(unsigned o1, unsigned o2) {

		assert(o1 != 0);
		assert(o2 != 0);
		assert(o1 != o2);

		unsigned prevO1 = _mach[o1];
		unsigned postO1 = mach[o1];
		unsigned prevO2 = _mach[o2];
		unsigned postO2 = mach[o2];

		mach[o1] = postO2;
		_mach[o2] = prevO1;

		if (prevO1 != 0) 
			mach[prevO1] = o2;

		if (postO2 != 0) 
			_mach[postO2] = o1;

		if (postO1 == o2) {
			_mach[o1] = o2;

			mach[o2] = o1;

			return;
		}

		_mach[o1] = prevO2;
		mach[o2] = postO1;

		if (postO1 != 0) 
			_mach[postO1] = o2;
		
		if (prevO2 != 0) 
			mach[prevO2] = o1;
	}*/

	void swap(unsigned o1, unsigned o2);

	//partial state has cycles?
	bool hasCycle(vector<unsigned>& indeg, vector<unsigned>& Q) const;

	void removeOper(unsigned targOp);

	void removeOperFromMach(unsigned targOp);

	//insert newOp right after prevOp and before postOp --- needs both in case prev is 0 needs to know cur root
	//care if using roots and leafs - must adjust outside
	void insertOper(unsigned newOp, unsigned prevJobOp, unsigned postJobOp, unsigned prevMachOp, unsigned postMachOp);

	void insertJobOper(unsigned newOp, unsigned prevJobOp, unsigned postJobOp);

	void insertMachOper(unsigned newOp, unsigned prevMachOp, unsigned postMachOp);

	bool compHeadsComplete(vector<unsigned>& heads, vector<unsigned>& indeg, vector<unsigned>& Q) const;

	bool compHeads(vector<unsigned>& heads, vector<unsigned>& indeg, vector<unsigned>& Q) const;

	bool compTailsComplete(vector<unsigned>& tails, vector<unsigned>& indeg, vector<unsigned>& Q) const;

	bool compTails(vector<unsigned>& tails, vector<unsigned>& indeg, vector<unsigned>& Q) const;

	//can reach fromIndex toIndex in Q?
	bool propagateReach(const vector<unsigned>& Q, vector<bool>& reach, unsigned fromIndexA, unsigned fromIndexB, unsigned toIndexA, unsigned toIndexB) const;

#ifdef VANILLA_PSS
	void insertInsaPos(unsigned curOp, vector<unsigned>& indeg, vector<unsigned>& Q, vector<bool>& inOpers, vector<unsigned>& jobRoot, vector<unsigned>& jobLeaf, vector<unsigned>& machRoot, vector<unsigned>& machLeaf, vector<unsigned>& heads, vector<unsigned>& tails, vector<unsigned>& invertedQ, vector<bool>& reach);
#endif //VANILLA_PSS

	void insaPsp(bool bigJobFirst, vector<unsigned>& indeg, vector<unsigned>& Q, vector<unsigned>& heads, vector<unsigned>& tails, vector<unsigned>& invertedQ, vector<bool>& reach);

	void gifflerThompson();

	static void computeCritic(vector<unsigned>& critic, unsigned lastOp, const vector<unsigned>& prev);

	static void blockBegins(vector<unsigned>& jobBb, vector<unsigned>& machBb, const vector<unsigned>& critic);

	static void fillCandidatesAllSwaps(vector<pair<unsigned, unsigned>>& cands, vector<unsigned>& mach);

	static void findCriticOper(vector<vector<unsigned>>& critic, vector<unsigned>& starts, vector<unsigned> _job, vector<unsigned> _mach);

	static void fillCandidatesCriticOper(vector<pair<unsigned, unsigned>>& cands, vector<vector<unsigned>>& criticOper);

	static void fillCandidatesCriticTotal(vector<pair<unsigned, unsigned>>& cands, vector<unsigned>& starts, vector<unsigned> _job, vector<unsigned> _mach);

	static void fillCandidatesCriticEarl(vector<pair<unsigned, unsigned>>& cands, vector<unsigned>& starts, vector<unsigned> _job, vector<unsigned> _mach);

	static void fillCandidatesN5(vector<pair<unsigned, unsigned>>& cand, vector<unsigned>& jobBb, vector<unsigned>& machBb, const vector<unsigned>& critic);

	//N5 neighbourhood first improvement - randomizes
//	void localSearch( vector<unsigned> & dists, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, vector<unsigned> & jobBb, vector<unsigned> & machBb, vector<pair<unsigned, unsigned>> & cands, vector<unsigned> & critic, unsigned lowerBound) {
//#ifndef NDEBUG
//		bool cycle;
//		assert(dists.size() == inst.O);
//		assert(prev.size() == inst.O);
//		assert(prev.size() == inst.O);
//		assert(Q.size() == inst.O);
//		assert(cands.capacity() == inst.O);
//		assert(critic.capacity() == inst.O);
//		assert(jobBb.capacity() == inst.O);
//		assert(machBb.capacity() == inst.O);
//#endif
//		unsigned lastOp;
//		unsigned thisMakes;
//
//		setMeta(dists, lastOp, prev, indeg, Q);
//
//		while(true) {
//			assert(makes >= lowerBound);
//			if(makes == lowerBound)
//				break;
//			thisMakes = makes;
//			computeCritic(critic, lastOp, prev);
//			assert( ! critic.empty());
//			assert(critic.size() < inst.O);
//			fillCandidatesN5(cands, jobBb, machBb, critic);
//			random_shuffle(cands.begin(), cands.end());
//
//			for(const pair<unsigned, unsigned> & p : cands) {
//				swap(p.first, p.second);
//#ifndef NDEBUG
//				cycle = 
//#endif
//					setMeta(dists, lastOp, prev, indeg, Q);
//				assert( ! cycle);
//
//				if(thisMakes > makes) {
//					break; //first improvement
//				} else {
//					swap(p.second, p.first); //undoing swap
//				}
//			}
//			if(thisMakes <= makes)
//				break;
//		}
//		setMeta(dists, lastOp, prev, indeg, Q);
//		assert(makes == thisMakes);
//	}

	//maxLen = maxD * maxC
	static bool detectRepeat(vector<int>& oldValues, unsigned& posCurPenalties, unsigned& cycleLastPos, unsigned& cycleL, double newPenalties, bool isNewBest, unsigned maxD, unsigned maxLen);

	//shift early operations making its completion time closer of its due date
	void shiftOperations(vector<unsigned>& starts, vector<unsigned>& Q);

	void getEarlBlock(vector<unsigned>& earlBlock) const;

	void cplexSolve(const unsigned maxSecs);

	/* Schedule operations relaxing machine precedence of operations in relaxBlock*/
	void schedulerCplexRelax(vector<unsigned>& relaxBlock);

	/* Schedule operations*/ 
	void schedulerCplex(vector<unsigned>& dists);

	//set makes with shift to minimize earliness penalties
	bool setMeta(vector<unsigned>& dists, unsigned& lastOp, vector<unsigned>& prev, vector<unsigned>& indeg, vector<unsigned>& Q, bool cplex);
	

//	//does not set critic but do set makes
//	//@return: cycle?
//	bool setMeta(vector<unsigned> & dists, unsigned & lastOp, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q)  {
//
//#ifdef SHIFT_OPERS
//		return setMetaWithShift(dists, lastOp, prev, indeg, Q);
//#endif //SHIFT_OPERS
//
//		assert(isAlloced());
//		assert(dists.size() == inst.O);
//		assert(prev.size() == inst.O);
//		assert(indeg.size() == inst.O);
//		assert(Q.size() == inst.O);
//
//		unsigned qInsert = 0;
//		unsigned qAccess = 0;
//
//		unsigned curOp;
//		unsigned newOp;
//
//		unsigned newMax;
//
//		makes = 0;
//
//		fill(dists.begin(), dists.end(), 0);
//		fill(indeg.begin(), indeg.end(), 0);
//
//		for(unsigned o=1; o<inst.O; o++) {
//			if(_job[o] != 0)
//				indeg[o]++;
//			if(_mach[o] != 0)
//				indeg[o]++;
//			if(indeg[o] == 0) {
//				prev[o] = 0;
//				Q[qInsert++] = o;
//			}
//		}
//
//		//assert(qInsert>0);
//
//		while(qAccess < qInsert) {
//			assert(qAccess<Q.size());
//			curOp = Q[qAccess++];
//
//			assert(indeg[curOp] == 0);
//
//			newMax = dists[curOp] + inst.P[curOp];
//			if(makes < newMax) {
//				makes = newMax;
//				lastOp = curOp;
//			}
//
//			//from JOB
//			newOp = job[curOp];
//			if(newOp != 0) {
//				assert(indeg[newOp] > 0);
//				indeg[newOp]--;
//				if(indeg[newOp] == 0) {
//					assert(qInsert<Q.size()-1);
//					Q[qInsert++] = newOp;
//				}
//				if(dists[newOp] < newMax) {
//					dists[newOp] = newMax;	
//					prev[newOp] = curOp;
//				}			
//			}
//
//			//from MACH
//			newOp = mach[curOp];
//			if(newOp != 0) {
//				assert(indeg[newOp] > 0);
//				indeg[newOp]--;
//				if(indeg[newOp] == 0) {
//					assert(qInsert<Q.size()-1);
//					Q[qInsert++] = newOp;
//				}
//				if(dists[newOp] < newMax) {
//					dists[newOp] = newMax; 
//					prev[newOp] = curOp;
//				}			
//			}
//		}
//		inst.calcPenalties(dists, ePenalty,lPenalty, operPenalties);
//		penalties = ePenalty + lPenalty;
//		startTime = dists;
//		
//		return qAccess<inst.O-1;
//	}

		//@return: starts
	vector<unsigned> genSchedule(bool cplex);

	bool verifySchedule(bool cplex);

	void printPenalties();

#ifndef NDEBUG
	bool isAlloced() const;

	bool isClear() const;

	//not all opers in state yet
	bool partialVerify(const vector<bool>& inOpers) const;


	bool verify() const;
#endif //NDEBUG


		vector<unsigned> job;							// job[o] is next operation of o in its job
		vector<unsigned> _job;						// _job[o] is previous operation of o in its job
		vector<unsigned> mach;						// mach[o] is next operation of o in its machine
		vector<unsigned> _mach;						// _mach[o] is previous operation of o in its machine
		vector<unsigned> startTime;				// startTime[o] is the startTime of operation o
		vector<double> operPenalties;			
		unsigned makes;
		unsigned millisecsFound;
		double lPenalty;
		double ePenalty;
		double penalties;
	
	};
