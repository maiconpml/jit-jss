/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/
#pragma once

#include <chrono>

#include "Settings.hpp"
#include "Utilities.hpp"
#include "State.hpp"
#include "Parameters.hpp"

using namespace chrono;

class TabuList {
public:
	TabuList();
	
	TabuList(int _tenure);

	TabuList(TabuList&& other);

	TabuList(const TabuList& other);

	TabuList& operator=(TabuList other);

	void swap(TabuList& other);

	string toString() const;

	void passTime(int time);

	string listForbidden() const;

	void insert(unsigned o1, unsigned o2);

	bool canMove(unsigned o1, unsigned o2) const;

	int age(unsigned o1, unsigned o2) const;

	int timeToLeave(unsigned o1, unsigned o2) const;

	int curIter;
	int tenure;
	multi_array<int,2> tabu;
};

class TabuTrio {
public:
	TabuTrio();

	TabuTrio(const State& _state, const vector<pair<unsigned, unsigned>>& _cands, const TabuList& _tabuList);

	~TabuTrio();

	string toString() const;

	bool dummy;
	State state;
	vector<pair<unsigned,unsigned>> cands; //cand moves - removed best one which was used
	TabuList tabuList;
};

class JumpList {
public:
	JumpList(unsigned _maxL);

	~JumpList();

	string toString() const;

	void push(const TabuTrio& tt);

	void updateCands(const vector<pair<unsigned, unsigned>>& cands);

	TabuTrio pop();

	bool empty;
	unsigned basePos;
	unsigned topPos;
	const unsigned maxL;
	vector<TabuTrio> ttList;
};

namespace Tabu {

	void nsp(State& theState, TabuList& tabuList, vector<pair<unsigned, unsigned>>& cands, double aspiration, vector<unsigned>& dists, vector<unsigned>& indeg, vector<unsigned>& Q, unsigned schedulerType, bool cplex);

	//@return: is optimal?
	//printWhenBetter use 0 to not print. will print preString makes seconds (according to tpSTart) d 
	bool evolveTabu(State& theState, const Parameters& param, const high_resolution_clock::time_point& tpStart, vector<unsigned>& dists, vector<unsigned>& indeg, vector<unsigned>& Q);
	
	void tabu(Parameters& param);
}


