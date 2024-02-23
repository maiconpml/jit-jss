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
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
#ifdef VAR_PARALLEL_MACHS
		void parllelIteratedGreedy(const string & instPath, string & name, double maxSecs, unsigned seed, unsigned tenure, unsigned initialjumpLimit, unsigned decreaseDivisor, unsigned increaseIters, unsigned bjSize, unsigned maxD, unsigned maxC, double acceptAlpha, unsigned sizePBubble, unsigned sizePInsa, bool timeLog, unsigned perturbType, const string & scheduleFile, double bubbleP, double mixProb) {
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
		multi_array<unsigned, 2> bucketLimits;
		bucketLimits.resize(boost::extents[inst.M][inst.maxK]);
		multi_array<unsigned, 2> bucketTops;
		bucketTops.resize(boost::extents[inst.M][inst.maxK]);
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
		//theState.insaPsp(false, indeg, Q, heads, tails, invertedQ, reach);
		theState.insaParallel(indeg, Q, heads, bucketLimits);
		
		unsigned lastOp;
		theState.makes = theState.setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);

#ifdef COUNT_STEPS
	    tpStart = high_resolution_clock::now(); //removing insa time when counting steps per second
#endif //COUNT_STEPS

		//STEP 1st intensification
		//Tabu::evolveTabu(theState, tenure, initialjumpLimit, decreaseDivisor, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, UINT_MAX, preString, heads, tails, timeLog, lowerBoundOrder);
		//XXX UNCOMMENT TO USE
		//unsigned jumpLimitDecrease = initialjumpLimit/decreaseDivisor;
		//XXX TODO ADD neighbs to use
		//Tabu::evolveTabuParallel(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, bucketLimits, bucketTops, maxD, maxC, maxMillisecs, 0, "");
		bestState = theState;

		curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		//IG body
		while(curMillisecs<maxMillisecs    &&    bestState.makes>lowerBound) {
			assert(theState.makes>=bestState.makes);

			newState = theState;

			//STEP perturb
			if(perturbType==BUBBLE_PERTURB)
				//TODO XXX XXX CHANGE
				newState.bubblePerturb(bubbleP, sizePBubble, enforce, heads, indeg, Q, tails, reach, ep);//sematics not same of some params - reused some for efficiency
			else if(perturbType==INSA_PERTURB)
				//TODO XXX XXX CHANGE
				newState.insaPerturb(sizePInsa, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails, invertedQ, reach);
			else if(perturbType==MIX_PERTURB) {
				aRand = random_number<double>(0.0, 1.0);
				if(mixProb >= aRand) {
				//TODO XXX XXX CHANGE
					newState.bubblePerturb(bubbleP, sizePBubble, enforce, heads, indeg, Q, tails, reach, ep);
				} else {
				//TODO XXX XXX CHANGE
					newState.insaPerturb(sizePInsa, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails, invertedQ, reach);
				}
			} else
				throw errorText("Invalid perturb option","Iterated","iteratedGreedy");

			
			//STEP intensify
			//Tabu::evolveTabu(newState, tenure, initialjumpLimit, decreaseDivisor, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, bestState.makes, preString, heads, tails, timeLog, lowerBoundOrder);
			//XXX TODO ADD neighbs to use
			//Tabu::evolveTabuParallel(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, bucketLimits, bucketTops, maxD, maxC, maxMillisecs, 0, "");
			assert(inst.verifyPrallelSchedule(dists, theState.makes));

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

		//if( ! bestState.verifySchedule()   ||   bestState.makes < lowerBound)
			//throw errorText(theState.toString() + "\n\t\tBad schedule !!!!!","","");
		
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
#endif //VAR_PARALLEL_MACHS
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	/*
	  Numerical Params:
	  @instPath:
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
	  @name: name of inst for log
	*/
	//print: name - seed - makes - millisecs - situation: d/t/o (during exec, timed out, optimal by lower bound)
	/*
	void timedIgDeltaSwap(const string & instPath, string & name, unsigned maxSecs, unsigned seed, unsigned tenure, unsigned initialjumpLimit, unsigned jumpLimitDecrease, unsigned bjSize, unsigned maxD, unsigned maxC, double acceptAlpha, unsigned perturbSize) {
		high_resolution_clock::time_point tpStart = high_resolution_clock::now();

	    RandomNumberGenerator<>::instance().seed(seed);
		srand(seed);

		inst.parse(instPath);

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

		vector<unsigned> dirtyDists(inst.O);
		vector<bool> dirty(inst.O);

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


		//cout << "FINISHED FIRST INTENSIFICATION" << endl;

		//STEP 1st intensification
		//Tabu::evolveTabu(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, UINT_MAX, preString);
		Tabu::evolveTabuDeltaSwap(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, UINT_MAX, preString, heads, tails, dirtyDists, dirty);

		bestState = theState;

		curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
		//if(bestState.makes > lowerBound) {//esle gonna print on end
			//name - seed - makes - millisecs - situation: d/t/o (during exec, timed out, optimal by lower bound)
			//cout << name << " " << seed << " " << bestState.makes << " " << curMillisecs/1000.0 << " " << "d" << endl;
		//}

		//IG body
		while(curMillisecs<maxMillisecs    &&    bestState.makes>lowerBound) {
			assert(theState.makes>=bestState.makes);

			//cout << "AAAAA: " << theState.makes << endl;

			newState = theState;

			//STEP perturb
			newState.insaPerturb(perturbSize, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails, invertedQ, reach);

			//cout << "BBBBB" << newState.makes<< endl;

			//STEP intensify
			Tabu::evolveTabuLbOrder(newState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, bestState.makes, preString, heads, tails);
			assert(newState.verifySchedule());

			//cout << "CCC"  << newState.makes<< endl;

			//STEP new best ?
			if(bestState.makes > newState.makes) {
				bestState = newState;
				//if(bestState.makes > lowerBound) {
				//curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
					//name - seed - makes - millisecs - situation: d/t/o (during exec, timed out, optimal by lower bound)
					//cout << name << " " << seed << " " << bestState.makes << " " << curMillisecs/1000.0 << " " << "d" << endl;
				//}
			}

			//STEP accept ?
			if(metropolis(theState.makes, newState.makes, T))
				theState = newState;

			curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

			//cout << "DDDD" << endl;
		}
		//curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		cout << name << " " << seed << " " << bestState.makes << " " << curMillisecs/1000.0 << " " << (string) (bestState.makes == lowerBound ? "o" : "t") << endl;

		if( ! bestState.verifySchedule()   ||   bestState.makes < lowerBound)
			throw errorText(theState.toString() + "\n\t\tBad schedule !!!!!","","");
		
#ifdef PRINT_SCHEDULE
#ifndef NDEBUG
		throw errorText("Recording schedules while in DEBUG mode?","","");
#endif
		ofstream toFile(PRINT_SCHEDULE, fstream::app);
		if( ! toFile.is_open()) throw errorText("Failed to open file.","","");
		vector<unsigned> sched = bestState.genSchedule();
		toFile << instPath << " makes: " << bestState.makes << endl << "sched:" ;
		for(unsigned s : sched) toFile << " " << s;
		toFile << endl << endl;
		toFile.close();
#endif
	}
	*/




	/*
	  Numerical Params:
	  @instPath:
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
	  @name: name of inst for log
	*/
	/*
	//print: name - seed - makes - millisecs - situation: d/t/o (during exec, timed out, optimal by lower bound)
	void timedIgLbOrdered(const string & instPath, string & name, double maxSecs, unsigned seed, unsigned tenure, unsigned initialjumpLimit, unsigned jumpLimitDecrease, unsigned bjSize, unsigned maxD, unsigned maxC, double acceptAlpha, unsigned perturbSize) {
		high_resolution_clock::time_point tpStart = high_resolution_clock::now();

	    RandomNumberGenerator<>::instance().seed(seed);
		srand(seed);

		inst.parse(instPath);

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

		//STEP 1st intensification
		//Tabu::evolveTabu(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, UINT_MAX, preString);
		Tabu::evolveTabuLbOrder(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, UINT_MAX, preString, heads, tails);

		bestState = theState;

		curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
		//if(bestState.makes > lowerBound) {//esle gonna print on end
			//name - seed - makes - millisecs - situation: d/t/o (during exec, timed out, optimal by lower bound)
			//cout << name << " " << seed << " " << bestState.makes << " " << curMillisecs/1000.0 << " " << "d" << endl;
		//}

		//IG body
		while(curMillisecs<maxMillisecs    &&    bestState.makes>lowerBound) {
			assert(theState.makes>=bestState.makes);

			newState = theState;

			//STEP perturb
			//XXX PERTURB BUBBLE HERE
			
		//double p = 0.5;
		//unsigned perturbSize = 4;
		//vector<unsigned> indeg(inst.O);
		//vector<unsigned> Q(inst.O);
		//vector<unsigned> operPos(inst.O);
		//vector<bool> reach(inst.O);
		//vector<unsigned> dists(inst.O);
		//bool bigJobFirst = false;
		//vector<vector<unsigned>> enforce(inst.O);
		//for(vector<unsigned> & v : enforce)
	//		v.reserve(max(inst.J, inst.M));
//		vector<unsigned> enfIndeg(inst.O);
//		unsigned lastOp;
//		vector<unsigned> prev(inst.O);
//		vector<unsigned> heads(inst.O); 
//		vector<unsigned> tails(inst.O);
//		vector<unsigned> invertedQ(inst.O);
			 
			//newState..bubblePerturb(p, perturbSize, enforce, enfIndeg, indeg, Q, operPos, reach);
			newState.insaPerturb(perturbSize, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails, invertedQ, reach);

			//STEP intensify
			Tabu::evolveTabuLbOrder(newState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, bestState.makes, preString, heads, tails);
			assert(newState.verifySchedule());

			//STEP new best ?
			if(bestState.makes > newState.makes) {
				bestState = newState;
				//if(bestState.makes > lowerBound) {
				//curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
					//name - seed - makes - millisecs - situation: d/t/o (during exec, timed out, optimal by lower bound)
					//cout << name << " " << seed << " " << bestState.makes << " " << curMillisecs/1000.0 << " " << "d" << endl;
				//}
			}

			//STEP accept ?
			if(metropolis(theState.makes, newState.makes, T))
				theState = newState;

			curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
		}
		//curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		cout << name << " " << seed << " " << bestState.makes << " " << curMillisecs/1000.0 << " " << (string) (bestState.makes == lowerBound ? "o" : "t") << endl;

		if( ! bestState.verifySchedule()   ||   bestState.makes < lowerBound)
			throw errorText(theState.toString() + "\n\t\tBad schedule !!!!!","","");
		
#ifdef PRINT_SCHEDULE
#ifndef NDEBUG
		throw errorText("Recording schedules while in DEBUG mode?","","");
#endif
		ofstream toFile(PRINT_SCHEDULE, fstream::app);
		if( ! toFile.is_open()) throw errorText("Failed to open file.","","");
		vector<unsigned> sched = bestState.genSchedule();
		toFile << instPath << " makes: " << bestState.makes << endl << "sched:" ;
		for(unsigned s : sched) toFile << " " << s;
		toFile << endl << endl;
		toFile.close();
#endif
	}
	*/





	/*
	  Numerical Params:
	  @instPath:
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
	  @name: name of inst for log
	*/
	/*
	//print: name - seed - makes - millisecs - situation: d/t/o (during exec, timed out, optimal by lower bound)
	void timedIg(const string & instPath, string & name, double maxSecs, unsigned seed, unsigned tenure, unsigned initialjumpLimit, unsigned jumpLimitDecrease, unsigned bjSize, unsigned maxD, unsigned maxC, double acceptAlpha, unsigned perturbSize) {
		high_resolution_clock::time_point tpStart = high_resolution_clock::now();

	    RandomNumberGenerator<>::instance().seed(seed);
		srand(seed);

		inst.parse(instPath);

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

		//STEP 1st intensification
		Tabu::evolveTabu(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, UINT_MAX, preString);

		bestState = theState;

		curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
		//if(bestState.makes > lowerBound) {//esle gonna print on end
			//name - seed - makes - millisecs - situation: d/t/o (during exec, timed out, optimal by lower bound)
			//cout << name << " " << seed << " " << bestState.makes << " " << curMillisecs/1000.0 << " " << "d" << endl;
		//}

		//IG body
		while(curMillisecs<maxMillisecs    &&    bestState.makes>lowerBound) {
			assert(theState.makes>=bestState.makes);

			newState = theState;

			//STEP perturb
			newState.insaPerturb(perturbSize, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails, invertedQ, reach);

			//STEP intensify
			Tabu::evolveTabu(newState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, bestState.makes, preString);
			assert(newState.verifySchedule());

			//STEP new best ?
			if(bestState.makes > newState.makes) {
				bestState = newState;
				//if(bestState.makes > lowerBound) {
				//curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
					//name - seed - makes - millisecs - situation: d/t/o (during exec, timed out, optimal by lower bound)
					//cout << name << " " << seed << " " << bestState.makes << " " << curMillisecs/1000.0 << " " << "d" << endl;
				//}
			}

			//STEP accept ?
			if(metropolis(theState.makes, newState.makes, T))
				theState = newState;

			curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
		}
		//curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		cout << name << " " << seed << " " << bestState.makes << " " << curMillisecs/1000.0 << " " << (string) (bestState.makes == lowerBound ? "o" : "t") << endl;

		if( ! bestState.verifySchedule()   ||   bestState.makes < lowerBound)
			throw errorText(theState.toString() + "\n\t\tBad schedule !!!!!","","");
		
#ifdef PRINT_SCHEDULE
#ifndef NDEBUG
		throw errorText("Recording schedules while in DEBUG mode?","","");
#endif
		ofstream toFile(PRINT_SCHEDULE, fstream::app);
		if( ! toFile.is_open()) throw errorText("Failed to open file.","","");
		vector<unsigned> sched = bestState.genSchedule();
		toFile << instPath << " makes: " << bestState.makes << endl << "sched:" ;
		for(unsigned s : sched) toFile << " " << s;
		toFile << endl << endl;
		toFile.close();
#endif
	}
	*/







	/*
	  Numerical Params:
	  @instPath:
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
	  @name: name of inst for log
	*/
	/*
	//print: name - seed - makes - millisecs - situation: d/t/o (during exec, timed out, optimal by lower bound)
	void igRace(const string & instPath, unsigned maxSecs, unsigned seed, unsigned tenure, unsigned initialjumpLimit, unsigned jumpLimitDecrease, unsigned bjSize, unsigned maxD, unsigned maxC, double acceptAlpha, unsigned perturbSize) {
		high_resolution_clock::time_point tpStart = high_resolution_clock::now();

	    RandomNumberGenerator<>::instance().seed(seed);
		srand(seed);

		inst.parse(instPath);

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

		//computing T for metropolis
		double meanProcTime = 0.0;
		for(unsigned p : inst.P) {
			meanProcTime += p;
		}
		meanProcTime /= (double)(inst.O-1);
		const double T = (double)acceptAlpha*(meanProcTime/10.0);

		//STEP generate first solution
		theState.insaPsp(false, indeg, Q, heads, tails, invertedQ, reach);

		//STEP 1st intensification
		Tabu::evolveTabu(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, 0, "");

		bestState = theState;

		curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		//IG body
		while(curMillisecs<maxMillisecs    &&    bestState.makes>lowerBound) {
			assert(theState.makes>=bestState.makes);

			newState = theState;

			//STEP perturb
			newState.insaPerturb(perturbSize, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails, invertedQ, reach);

			//STEP intensify
			Tabu::evolveTabu(newState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, 0, "");
			assert(newState.verifySchedule());

			//STEP new best ?
			if(bestState.makes > newState.makes) {
				bestState = newState;
				if(bestState.makes > lowerBound) {
					curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
				}
			}

			//STEP accept ?
			if(metropolis(theState.makes, newState.makes, T))
				theState = newState;

			curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
		}
		curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		cout << bestState.makes;

		if( ! bestState.verifySchedule()   ||   bestState.makes < lowerBound)
			throw errorText(theState.toString() + "\n\t\tBad schedule !!!!!","","");
		
#ifdef PRINT_SCHEDULE
#ifndef NDEBUG
		throw errorText("Recording schedules while in DEBUG mode?","","");
#endif
		ofstream toFile(PRINT_SCHEDULE, fstream::app);
		if( ! toFile.is_open()) throw errorText("Failed to open file.","","");
		vector<unsigned> sched = bestState.genSchedule();
		toFile << instPath << " makes: " << bestState.makes << endl << "sched:" ;
		for(unsigned s : sched) toFile << " " << s;
		toFile << endl << endl;
		toFile.close();
#endif
	}
	*/






	/*
	  Numerical Params:
	  @instPath:
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
	  @bubbleP:
	  Categorical params:
	  @startType: first solution generator: INSA_START -  RAND_START
	  @intensificType: intenification step of ig: LOCAL_SEARCH_INTENS  -  TABU_INTENS
	  @perturbType: perturbation type: INSA_PERTURB  --  BUBBLE_PERTURB  --  CARLIER_PERTURB
	*/
	/*
	void igGeneric(const string & instPath, unsigned maxSecs, unsigned seed, unsigned tenure, unsigned initialjumpLimit, unsigned jumpLimitDecrease, unsigned bjSize, unsigned maxD, unsigned maxC, double acceptAlpha, unsigned perturbSize, double bubbleP, unsigned startType, unsigned intensificType, unsigned perturbType, unsigned acceptType) {
		high_resolution_clock::time_point tpStart = high_resolution_clock::now();

	    RandomNumberGenerator<>::instance().seed(seed);
		srand(seed);

		inst.parse(instPath);

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
		unsigned lastOp;
		const unsigned maxMillisecs = maxSecs* 1000;
		unsigned curMillisecs;
		unsigned startMillisecs;
		unsigned totalMillisecs;
		unsigned startMakes;
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

		//computing T for metropolis
		double meanProcTime = 0.0;
		for(unsigned p : inst.P) {
			meanProcTime += p;
		}
		meanProcTime /= (double)(inst.O-1);
		const double T = (double)acceptAlpha*(meanProcTime/10.0);

		//STEP generate first solution
		if(startType==INSA_START) {
			theState.insaPsp(false, indeg, Q, heads, tails, invertedQ, reach);
		} else {// if(startType==RAND_START)
			theState.makeRand(indeg, Q);
		}
		theState.setMeta(dists, lastOp, prev, indeg, Q);

		theState.millisecsFound = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		//STEP 1st intensification
		if(intensificType == LOCAL_SEARCH_INTENS) {
			theState.localSearch(dists, prev, indeg, Q, jobBb, machBb, cands, critic, lowerBound);
			theState.millisecsFound = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
		} else {// (intensificType == TABU_INTENS) 
			//millisecsFound updated in Tabu
			Tabu::evolveTabu(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, 0, "");
	    }

		//cout << "first instesification done:" <<  theState.makes << " at " << theState.millisecsFound/1000.0 << endl;

		startMakes = theState.makes;
		startMillisecs = theState.millisecsFound;
		curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		bestState = theState;

		//IG body
		while(curMillisecs<maxMillisecs    &&    bestState.makes>lowerBound) {
			assert(theState.makes>=bestState.makes);

			newState = theState;

			//STEP perturb
			if(perturbType == INSA_PERTURB) {
				newState.insaPerturb(perturbSize, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails, invertedQ, reach);
			} else { //(perturbType == BUBBLE_PERTURB)
				newState.bubblePerturb(bubbleP, perturbSize, enforce, enfIndeg, indeg, Q, operPos, reach);
			}

			//STEP intensify
			if(intensificType == LOCAL_SEARCH_INTENS) {
				newState.localSearch(dists,prev, indeg, Q, jobBb, machBb, cands, critic, lowerBound);
				newState.millisecsFound = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
			} else {// (intensificType == TABU_INTENS) 
				Tabu::evolveTabu(newState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, lowerBound, maxMillisecs, 0, "");
			}
			assert(newState.verifySchedule());

			//STEP new best ?
			if(bestState.makes > newState.makes) 
				bestState = newState;

			//STEP accept ?
			if(acceptType == METROPOLIS) {
				if(metropolis(theState.makes, newState.makes, T)) {
					theState = newState;
				}
			} else if(acceptType == BEST_JUMP) {
				curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
				if(((curMillisecs/(double)maxMillisecs)*(curMillisecs/(double)maxMillisecs)) > random_number<double>(0.0, 1.0)) {
					theState = bestState;
				} else {
					theState = newState;
				}
			}

			curMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
		}
		totalMillisecs = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

		//Name              Seed                      InitMakes                      InitT                                        Makes                                         TTT                                            TotalTime
		cout << " " << seed << " "  << startMakes << " " << startMillisecs/1000.0 << " " << bestState.makes << " " <<  bestState.millisecsFound/1000.0 <<  " " <<   totalMillisecs/1000.0 << endl;

		if( ! bestState.verifySchedule())
			throw errorText(theState.toString() + "\n\t\tBad schedule !!!!!","","");
		
#ifdef PRINT_SCHEDULE
#ifndef NDEBUG
		throw errorText("Recording schedules while in DEBUG mode?","","");
#endif
		ofstream toFile(PRINT_SCHEDULE, fstream::app);
		if( ! toFile.is_open()) throw errorText("Failed to open file.","","");
		vector<unsigned> sched = bestState.genSchedule();
		toFile << instPath << " makes: " << bestState.makes << endl << "sched:" ;
		for(unsigned s : sched) toFile << " " << s;
		toFile << endl << endl;
		toFile.close();
#endif
	}
	*/
}//end IG namespace
