#include <string>
#include <iostream>

#include "Tabu.hpp"

void Tabu::nsp(State& theState, TabuList& tabuList, vector<pair<unsigned, unsigned>>& cands, double aspiration, vector<unsigned>& starts, vector<unsigned>& indeg, vector<unsigned>& Q, unsigned schedulerType) {
#ifndef NDEBUG
	assert(!cands.empty());
	assert(starts.size() == inst.O);
	assert(indeg.size() == inst.O);
	assert(Q.size() == inst.O);
	State testState = theState;
#endif
	high_resolution_clock::time_point tpStart;
	bool cycle;
	unsigned o1;
	unsigned o2;
	vector<bool> stillToTestCandPos(cands.size(), true);
	vector<pair<unsigned, unsigned>> candVisitOrder(cands.size());
	vector<unsigned> heads(inst.O);
	vector<unsigned> tails(inst.O);

	unsigned pPos;

	//Unforbidden or aspiration here
	double chosenPenalties = UINT_MAX;
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
	assert(!cycle);

#ifndef NDEBUG
	cycle =
#endif
		theState.compHeadsComplete(heads, indeg, Q);
	assert(!cycle);

	//preparing candVisitOrder
	for (unsigned curCandIndex = 0; curCandIndex < cands.size(); curCandIndex++) {
		candVisitOrder[curCandIndex].first = theState.swapLowerBound(cands[curCandIndex].first, cands[curCandIndex].second, heads, tails);
		candVisitOrder[curCandIndex].second = curCandIndex;
	}

	//unforbidden moves
	//for(unsigned pPos=0; pPos<cands.size(); pPos++) {
	for (unsigned curVisitOrderI = 0; curVisitOrderI < candVisitOrder.size(); curVisitOrderI++) {

		pPos = candVisitOrder[curVisitOrderI].second;
		// changed 
		o1 = cands[pPos].first;
		o2 = cands[pPos].second;
		assert(o1 != 0);
		assert(o2 != 0);

		if (tabuList.canMove(o1, o2)) {
			stillToTestCandPos[pPos] = false;

			theState.swap(o1, o2);
			assert(theState.verify());

			if (schedulerType == 4) {
				cycle = theState.scheduleCplex(starts);
			}
			else if (schedulerType == 2) {
				cycle = theState.scheduleDelaying(starts);
			}
			else if (schedulerType == 1) {
				cycle = theState.scheduleAsEarly(starts);
			}

			if (theState.penalties < chosenPenalties && !cycle) {
				chosenPenalties = theState.penalties;
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
	for (unsigned curVisitOrderI = 0; curVisitOrderI < candVisitOrder.size(); curVisitOrderI++) {

		pPos = candVisitOrder[curVisitOrderI].second;

		if (stillToTestCandPos[pPos]) {
			o1 = cands[pPos].first;
			o2 = cands[pPos].second;
			assert(o1 != 0);
			assert(o2 != 0);

			assert(!tabuList.canMove(o1, o2));

			theState.swap(o1, o2);
			assert(theState.verify());

			if (schedulerType == 4) {
				cycle = theState.scheduleCplex(starts);
			}
			else if (schedulerType == 2) {
				cycle = theState.scheduleDelaying(starts);
			}
			else if (schedulerType == 1) {
				cycle = theState.scheduleAsEarly(starts);
			}

			//aspiration
			if (theState.penalties < aspiration && theState.penalties < chosenPenalties && !cycle) {
				chosenPenalties = theState.penalties;
				chosenO1 = o1;
				chosenO2 = o2;
				chosenSwapPos = pPos;
			}

			//crap moves
			if (chosenPenalties == UINT_MAX && fnpStateAge < tabuList.age(o1, o2) && !cycle) {
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
	if (chosenPenalties != UINT_MAX) {
		//good move found
		assert(chosenO1 != 0 && chosenO1 < inst.O);
		assert(chosenO2 != 0 && chosenO2 < inst.O);
		assert(chosenSwapPos < cands.size());
		theState.swap(chosenO1, chosenO2);
		tabuList.insert(chosenO1, chosenO2);
		cands.erase(cands.begin() + chosenSwapPos);
	}
	else {
		//only crap moves
		assert(chosenO1 == 0);
		assert(chosenO2 == 0);
		assert(oldestO1 != 0 && oldestO1 < inst.O);
		assert(oldestO2 != 0 && oldestO2 < inst.O);
		assert(fnpStateAge >= 0);
		theState.swap(oldestO1, oldestO2);
		tabuList.passTime(tabuList.timeToLeave(oldestO1, oldestO2));
		tabuList.insert(oldestO1, oldestO2);
		cands.erase(cands.begin() + oldestFnpSwapPos);
	}
	if (schedulerType == 4) {
#ifndef NDEBUG
		cycle =
#endif
			theState.scheduleCplex(starts);
	}
	else if (schedulerType == 2) {
#ifndef NDEBUG
		cycle =
#endif
			theState.scheduleDelaying(starts);
	}
	else if (schedulerType == 1) {
#ifndef NDEBUG
		cycle =
#endif
			theState.scheduleAsEarly(starts);
	}
	assert(!cycle);
}

//@return: is optimal?
//printWhenBetter use 0 to not print. will print preString makes seconds (according to tpSTart) d 
bool Tabu::evolveTabu(State& theState, const Parameters& param, const high_resolution_clock::time_point& tpStart, vector<unsigned>& starts, vector<unsigned>& indeg, vector<unsigned>& Q) {
#ifndef NDEBUG
	bool cycle;
#endif
#ifndef NDEBUG
	cycle =
#endif
	theState.scheduleCplex(starts);
	assert(!cycle);
	assert(theState.penalties >= 0);

	if (theState.penalties == 0)
		return true;

	State curState = theState;
	State oldState;

	vector<vector<unsigned>> criticOper(inst.O);

	bool penaltiesCycleDetected;
	vector<int> oldValues(param.maxD, -1);
	unsigned posCurPenalties = 0;
	unsigned cycleLastPos = 0;
	unsigned cycleL = 0;
	const unsigned maxLen = param.maxD * param.maxC;

	vector<pair<unsigned, unsigned>> cands;
	cands.reserve(inst.O);
	TabuTrio tt;

	unsigned tabuIter = 0;
	unsigned noImproveIters = 0;
	unsigned trySwapNeigh = 0;
	unsigned curJumpLimit = param.initialjumpLimit;
	bool jumped = false;
	bool newBestFound = true;
	bool cplex = false;
	bool emptyneighborhood = false;

	// 1 initiatizations
	TabuList oldTabuList;
	TabuList tabuList(param.tenure);
	JumpList jumpList(param.bjSize);

	// 2 main loop
	while (param.maxMillisecs > duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count()) {
		tabuIter++;
		noImproveIters++;

		//STEP Checking for cycles
		penaltiesCycleDetected = State::detectRepeat(oldValues, posCurPenalties, cycleLastPos, cycleL, curState.penalties, newBestFound, param.maxD, maxLen); //PENAL alt makes

		// STEP get trio from jump list or compute cands
		if ((curJumpLimit < noImproveIters) || penaltiesCycleDetected || emptyneighborhood) {//JUMPING
			assert(!newBestFound);

			jumped = true;
			curJumpLimit -= (curJumpLimit / param.decreaseDivisor);

			emptyneighborhood = false;

			noImproveIters = 0;

			tt = jumpList.pop();

			if (!tt.dummy) {
				curState = tt.state;
				cands = tt.cands;
				tabuList = tt.tabuList;
				assert(!cands.empty());
			}
			else { //Stop -> no more states in jump list
				return false;
			}
		}
		else { //NOT JUMPING
			jumped = false;
			cands.clear();

			if (trySwapNeigh == 1) {
				State::fillCandidatesAllSwaps(cands, curState.startTime, curState._job, curState._mach, curState.mach);
			}
			else if (curState.lPenalty <= theState.penalties / 4) {
				vector<unsigned> earlBlock;
				curState.getEarlBlock(earlBlock);
				curState.schedulerCplexRelax(earlBlock);
				cplex = true;
				State::fillCandidatesCriticTotal(cands, curState.startTime, curState._job, curState._mach, curState.mach);
			}
			else {
				State::fillCandidatesCriticTotal(cands, curState.startTime, curState._job, curState._mach, curState.mach);
			}
#ifdef PRINT_NEIGHBOURS_NB
			neigh.push_back(cands.size());
#endif // PRINT_NEIGHBOURS_NB
		}
		assert(!cands.empty());

		if (newBestFound) {
			assert(!jumped);
			oldState = curState;
			oldTabuList = tabuList;
		}

		//STEP Go to neighbour			

		if(!cplex){
			if (!cands.empty()) {
				if(param.schedulerType == 3){
					if (curState.lPenalty >= curState.ePenalty) {
						nsp(curState, tabuList, cands, theState.penalties, starts, indeg, Q, 1);
					}
					else {
						nsp(curState, tabuList, cands, theState.penalties, starts, indeg, Q, 2);
					}
				}
				else{
					nsp(curState, tabuList, cands, theState.penalties, starts, indeg, Q, param.schedulerType);
				}
			}
			else {
				emptyneighborhood = true;
			}
		}

		//STEP last iter new best was found? store oldState with paths not taken
		if (newBestFound) {
			if (!cands.empty())
				jumpList.push(TabuTrio(oldState, cands, oldTabuList));
			newBestFound = false;
		}

		if (jumped) {
			assert(!newBestFound);
			jumpList.updateCands(cands);
		}

		//STEP new best ??
		if (curState.penalties < theState.penalties) {

			trySwapNeigh = 0;
			noImproveIters = 0;
			theState = curState;
			curJumpLimit = param.initialjumpLimit;
			assert(curState.penalties >= 0);
			if (curState.penalties == 0) //Stop -> Optimal
				return true;

			newBestFound = true;

			theState.millisecsFound = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
			curState.millisecsFound = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
		}
		else {
			if (trySwapNeigh == 1) trySwapNeigh = 0;
			else trySwapNeigh += 1;
		}

		if (cplex) cplex = false;
	}//end while

	cout << tabuIter << endl;

	assert(param.maxMillisecs <= duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count());

	return false;
}

void Tabu::tabu(Parameters& param) {

	high_resolution_clock::time_point tpStart = high_resolution_clock::now();

	inst.parse(param.instPath);

	vector<unsigned> starts(inst.O);
	vector<unsigned> indeg(inst.O);
	vector<unsigned> Q(inst.O);

	State theState;
	theState.alloc();

	unsigned lowerBound = inst.lowerBoundTkz(indeg, Q);

	theState.gifflerThompson();

	if (param.schedulerType == 4) {
		theState.scheduleCplex(starts);
	}
	else if (param.schedulerType == 2 || param.schedulerType == 3) {
		theState.scheduleDelaying(starts);
	}
	else {
		theState.scheduleAsEarly(starts);
	}
	
	
#ifdef PRINT_ONLY_IS
	cout << param.instPath << " " << lowerBound << " ";
	theState.printPenalties();
	return;
#endif //PRINT_ONLY_IS

	theState.millisecsFound = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

	evolveTabu(theState, param, tpStart, starts, indeg, Q);

	theState.scheduleCplex(starts);

	theState.millisecsFound = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

	if (!theState.verifySchedule(4))
		throw errorText(theState.toString() + "\n\t\tBad schedule !!!!!", "", "");

#ifdef PRINT_DEFAULT
	cout << param.instPath << " " << lowerBound << " ";
#endif //PRINT_ONLY_RESULT
	theState.printPenalties(4);
}