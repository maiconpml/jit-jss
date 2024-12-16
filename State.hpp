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

	static void fillCandidatesAllSwaps(vector<pair<unsigned, unsigned>>& cands, vector<unsigned>& starts, vector<unsigned>& _job, vector<unsigned>& _mach, vector<unsigned>& mach);

	static void fillCandidatesCriticTotal(vector<pair<unsigned, unsigned>>& cands, vector<unsigned>& starts, vector<unsigned>& _job, vector<unsigned>& _mach, vector<unsigned>& mach);

	//maxLen = maxD * maxC
	static bool detectRepeat(vector<int>& oldValues, unsigned& posCurPenalties, unsigned& cycleLastPos, unsigned& cycleL, double newPenalties, bool isNewBest, unsigned maxD, unsigned maxLen);

	//shift early operations making its completion time closer of its due date
	void shiftOperations(vector<unsigned>& starts, vector<unsigned>& Q);

	void getEarlBlock(vector<unsigned>& earlBlock) const;

	void cplexSolve(const unsigned maxSecs);

	/* Schedule operations relaxing machine precedence of operations in relaxBlock*/
	void schedulerCplexRelax(vector<unsigned>& relaxBlock);

	/* Schedule operations*/ 
	void scheduleCplex(vector<unsigned>& dists, bool cplex);

	//set makes with shift to minimize earliness penalties
	bool scheduleAsEarly(vector<unsigned>& starts, bool cplex);

		//@return: starts
	vector<unsigned> genSchedule(unsigned schedulerType, bool cplex);

	bool verifySchedule(unsigned schedulerType, bool cplex);

	void printPenalties(unsigned schedulerType);

	bool topoWalk( vector<unsigned>& indeg, vector<unsigned>& Q );

	void updateStrength(vector<unsigned>& lateCands, double & pS, double & hS );

	unsigned calcDelayTime(vector<unsigned> & starts,vector<unsigned>& lateCands,vector<unsigned>& limited);

	void delay(vector<unsigned> & starts,vector<unsigned>& lateCands, unsigned & t);

	void update(vector<unsigned>& lateCands, vector<unsigned> & limited, vector<unsigned> & heads, vector<unsigned> & starts);
		
	void forcedDelay(vector<unsigned> & starts,vector<unsigned>& lateCands, unsigned & op);

	bool scheduleDelaying(vector<unsigned>& starts, bool cplex);
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
