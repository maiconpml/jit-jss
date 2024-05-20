/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/
#pragma once

#include <algorithm> 
#include <chrono>

#include "Settings.hpp"
#include "Utilities.hpp"
#include "Tabu.hpp"

using namespace chrono;


namespace IG {

	/*
	  Numerical Params:
	  @instPath:
	  @name: name of inst for log
	  @maxSecs:
	  @seed: 
	  @tenure: tabu tenure
	  @initialjumpLimit: iterations without improvements before backjump
	  @jumpLimitDecrease: decrease in iteration at each backjump without improvements
	  @bjSize: backjump list size
	  @maxD: size of cycle in tabu to be considered a cycle
	  @maxC: number of repeat cycles in tabu to backjump early
	  @acceptAlpha: for metropolis
	  @perturbSize: size of perturbation: number of operation or items removed and reinserted
	  @scheduleFile: file to print best found schedule. empty string means no logging
	*/
	void iteratedGreedy(const string & instPath, string & name, double maxSecs, unsigned seed, unsigned tenure, unsigned initialjumpLimit, unsigned decreaseDivisor, unsigned increaseIters, unsigned bjSize, unsigned maxD, unsigned maxC, double acceptAlpha, unsigned sizePBubble, unsigned sizePInsa, bool timeLog, bool lowerBoundOrder, unsigned perturbType, const string & scheduleFile, double bubbleP, double mixProb) {
		high_resolution_clock::time_point tpStart = high_resolution_clock::now();

	    RandomNumberGenerator<>::instance().seed(seed);
		srand(seed);

		inst.parse(instPath);
#ifdef COUNT_STEPS
		numSteps = 0;
		shownSeconds = 0;
#endif //COUNT_STEPS
		
		//trivial initializtions + allocs
		vector<unsigned> dists(inst.O);
		vector<unsigned> prev(inst.O);
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		vector<unsigned> jobBb;
		jobBb.reserve(inst.O);
		vector<unsigned> machBb;
		machBb.reserve(inst.O);
		vector<pair<unsigned, unsigned>> cands;
		cands.reserve(inst.O);
		vector<unsigned> critic;
		critic.reserve(inst.O);
		vector<bool> inOpers(inst.O, true);
		inOpers[0]=false;
		vector<unsigned> jobRoot(inst.J);
		vector<unsigned> jobLeaf(inst.J);
		vector<unsigned> machRoot(inst.M);
		vector<unsigned> machLeaf(inst.M);
		const unsigned maxMillisecs = maxSecs* 1000;
		unsigned curMillisecs;
		vector<unsigned> heads(inst.O);
		vector<unsigned> tails(inst.O);
		vector<unsigned> invertedQ(inst.O);
		vector<bool> reach(inst.O);
		const unsigned lowerBound = inst.lowerBoundTkz(indeg, Q);
		State theState;
		theState.alloc();
		State newState;
		State bestState;
		vector<vector<unsigned>> enforce(inst.O);
		for(vector<unsigned> & v : enforce)
			v.reserve(max(inst.J, inst.M));
		vector<unsigned> enfIndeg(inst.O);
		vector<unsigned> operPos(inst.O);
		double aRand;
		//unsigned bubblePertSize = (perturbSize/5);


		ElementProb ep;
		if(perturbType==BUBBLE_PERTURB   ||   perturbType==MIX_PERTURB)
			ep.setMap();

		string preString = name+" "+unsigStr(seed)+" ";

		//computing T for metropolis
		double meanProcTime = 0.0;
		for(unsigned p : inst.P) {
			meanProcTime += p;
		}
		meanProcTime /= (double)(inst.O-1);
		const double T = (double)acceptAlpha*(meanProcTime/10.0);

		//STEP generate first solution
		theState.insaPsp(false, indeg, Q, heads, tails, invertedQ, reach);

#ifdef COUNT_STEPS
	    tpStart = high_resolution_clock::now(); //removing insa time when counting steps per second
#endif //COUNT_STEPS

		//STEP 1st intensification
		Tabu::evolveTabu(theState, tenure, initialjumpLimit, decreaseDivisor, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, UINT_MAX, preString, heads, tails, timeLog, lowerBoundOrder);
		bestState = theState;

		curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		//IG body
		while(curMillisecs<maxMillisecs    &&    bestState.makes>lowerBound) {
			assert(theState.makes>=bestState.makes);

			newState = theState;

			//STEP perturb
			if(perturbType==BUBBLE_PERTURB)
				newState.bubblePerturb(bubbleP, sizePBubble, enforce, heads, indeg, Q, tails, reach, ep);//sematics not same of some params - reused some for efficiency
			else if(perturbType==INSA_PERTURB) 
				newState.insaPerturb(sizePInsa, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails, invertedQ, reach);
			else if(perturbType==MIX_PERTURB) {
				aRand = random_number<double>(0.0, 1.0);
				if(mixProb >= aRand) {
					newState.bubblePerturb(bubbleP, sizePBubble, enforce, heads, indeg, Q, tails, reach, ep);
				} else {
					newState.insaPerturb(sizePInsa, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails, invertedQ, reach);
				}
			} else
				throw errorText("Invalid perturb option","Iterated","iteratedGreedy");

			
			//STEP intensify
			Tabu::evolveTabu(newState, tenure, initialjumpLimit, decreaseDivisor, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, bestState.makes, preString, heads, tails, timeLog, lowerBoundOrder);
			assert(newState.verifySchedule());

			initialjumpLimit += increaseIters;

			//STEP new best ?
			if(bestState.makes > newState.makes) {
				bestState = newState;
			}

			//STEP accept ?
			if(metropolis(theState.makes, newState.makes, T))
				theState = newState;

			curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
		}
		//curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		if(timeLog) {
			//cout << name << " " << seed << " " << bestState.makes << " " << curMillisecs/1000.0 << " " << (string) (bestState.makes == lowerBound ? "o" : "t") << endl;
			for(string s : resultList) 
				cout << s << endl;
		} else {
			cout << bestState.makes << endl;
		}

		if( ! bestState.verifySchedule()   ||   bestState.makes < lowerBound)
			throw errorText(theState.toString() + "\n\t\tBad schedule !!!!!","","");
		
		if(scheduleFile.size() != 0) {
#ifndef NDEBUG
			throw errorText("Recording schedules while in DEBUG mode?","","");
#endif
			ofstream toFile(scheduleFile, fstream::app);
			if( ! toFile.is_open()) throw errorText("Failed to open file.","","");
			vector<unsigned> sched = bestState.genSchedule();
			toFile << instPath << " seed: "  << seed  << " makes: " << bestState.makes << endl << "sched:" ;
			for(unsigned s : sched) toFile << " " << s;
			toFile << endl << endl;
			toFile.close();
		}
	}
}
