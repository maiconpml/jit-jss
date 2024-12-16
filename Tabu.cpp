#include <string>
#include <iostream>

#include "Tabu.hpp"

void Tabu::nsp(State& theState, unsigned& lastOp, TabuList& tabuList, vector<pair<unsigned, unsigned>>& cands, double aspiration, vector<unsigned>& dists, vector<unsigned>& prev, vector<unsigned>& indeg, vector<unsigned>& Q, vector<unsigned>& heads, vector<unsigned>& tails, bool lbOrder, bool cplex) {
#ifndef NDEBUG
	assert(!cands.empty());
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
	if (lbOrder)
		sort(candVisitOrder.begin(), candVisitOrder.end());

	//unforbidden moves
	//for(unsigned pPos=0; pPos<cands.size(); pPos++) {
	for (unsigned curVisitOrderI = 0; curVisitOrderI < candVisitOrder.size(); curVisitOrderI++) {

		if (lbOrder && candVisitOrder[curVisitOrderI].first >= chosenPenalties)
			break;

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

			cycle = theState.setMeta(dists, lastOp, prev, indeg, Q, cplex);
			//assert( ! cycle);



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

		if (lbOrder && candVisitOrder[curVisitOrderI].first >= chosenPenalties)
			break;

		pPos = candVisitOrder[curVisitOrderI].second;

		if (stillToTestCandPos[pPos]) {
			o1 = cands[pPos].first;
			o2 = cands[pPos].second;
			assert(o1 != 0);
			assert(o2 != 0);

			assert(!tabuList.canMove(o1, o2));

			theState.swap(o1, o2);
			assert(theState.verify());

			cycle = theState.setMeta(dists, lastOp, prev, indeg, Q, cplex);
			//assert( ! cycle);

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
#ifndef NDEBUG
	cycle =
#endif
		theState.setMeta(dists, lastOp, prev, indeg, Q, cplex);
	assert(!cycle);
}

//@return: is optimal?
//printWhenBetter use 0 to not print. will print preString makes seconds (according to tpSTart) d 
bool Tabu::evolveTabu(State& theState, unsigned tenure, unsigned initialjumpLimit, unsigned decreaseDivisor, unsigned bjSize, const high_resolution_clock::time_point& tpStart, vector<unsigned>& dists, vector<unsigned>& prev, vector<unsigned>& indeg, vector<unsigned>& Q, const unsigned maxD, const unsigned maxC, const double lowerBound, unsigned maxMillisecs, unsigned printWhenBetter, const string& preString, vector<unsigned>& heads, vector<unsigned>& tails, bool timeLog, bool lbOrder) {
#ifndef NDEBUG
	bool cycle;
#endif
	unsigned lastOp;
#ifndef NDEBUG
	cycle =
#endif
		theState.setMeta(dists, lastOp, prev, indeg, Q, true);
	assert(!cycle);
	assert(theState.penalties >= lowerBound);

	if (theState.penalties < printWhenBetter) {
		if (timeLog) {
			resultList.push_back(preString + unsigStr(theState.penalties) + " " + doubleStr((duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count()) / 1000.0) + " d");
		}
	}

	if (theState.penalties == lowerBound)
		return true;

	State curState = theState;
	State oldState;

	vector<unsigned> jobBb;
	jobBb.reserve(inst.O);
	vector<unsigned> machBb;
	machBb.reserve(inst.O);
	vector<unsigned> critic;
	critic.reserve(inst.O);
	vector<vector<unsigned>> criticOper(inst.O);

	bool penaltiesCycleDetected;
	vector<int> oldValues(maxD, -1);
	unsigned posCurPenalties = 0;
	unsigned cycleLastPos = 0;
	unsigned cycleL = 0;
	const unsigned maxLen = maxD * maxC;

	vector<pair<unsigned, unsigned>> cands;
	cands.reserve(inst.O);
	TabuTrio tt;

	unsigned tabuIter = 0;
	unsigned noImproveIters = 0;
	unsigned trySwapNeigh = 0;
	unsigned curJumpLimit = initialjumpLimit;
	bool jumped = false;
	bool newBestFound = true;

	bool emptyNeighbourhood = false;

	// 1 initiatizations
	TabuList oldTabuList;
	TabuList tabuList(tenure);
	//unsigned lowerBound = inst.lowerBoundTkz(indeg, Q);
	JumpList jumpList(bjSize);
	bool cplex = true;

	// 2 main loop
	while (maxMillisecs > duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count()) {
		tabuIter++;
		noImproveIters++;

		//STEP Checking for cycles
		penaltiesCycleDetected = State::detectRepeat(oldValues, posCurPenalties, cycleLastPos, cycleL, curState.penalties, newBestFound, maxD, maxLen); //PENAL alt makes

		// STEP get trio from jump list or compute cands
		if ((curJumpLimit < noImproveIters) || penaltiesCycleDetected || emptyNeighbourhood) {//JUMPING
			assert(!newBestFound);

			jumped = true;
			//if(curJumpLimit > jumpLimitDecrease) 
			//curJumpLimit -= jumpLimitDecrease;
			curJumpLimit -= (curJumpLimit / decreaseDivisor);
			//cout << curJumpLimit << endl;

			emptyNeighbourhood = false;

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



			//				State::computeCritic(critic, lastOp, prev);
			////#ifndef NDEBUG
			////				assert( ! critic.empty());
			////				assert(critic.size() < inst.O);
			////				testMakes = 0;
			////				for(unsigned pos=0; pos<critic.size(); pos++)
			////					testMakes += inst.P[critic[pos]];
			////				assert(testMakes == curState.makes);
			////				for(unsigned pos=0; pos<critic.size()-1; pos++)
			////					assert(curState.job[critic[pos]] == critic[pos+1]    ||    curState.mach[critic[pos]] == critic[pos+1]);
			////#endif
			//				State::fillCandidatesN5(cands, jobBb, machBb, critic);

							//State::fillCandidatesAllSwaps(cands, curState.mach);

			if (trySwapNeigh == 1) {
				State::fillCandidatesAllSwaps(cands, curState.mach);
			}
			else if (curState.lPenalty <= theState.penalties / 4) {
				vector<unsigned> earlBlock;
				curState.getEarlBlock(earlBlock);
				curState.schedulerCplexRelax(earlBlock);
				cplex = true;
				State::fillCandidatesCriticTotal(cands, curState.startTime, curState._job, curState._mach);
			}
			else {
				//State::findCriticOper(criticOper, curState.startTime, curState._job, curState._mach);
				//State::fillCandidatesCriticOper(cands, criticOper);
				State::fillCandidatesCriticTotal(cands, curState.startTime, curState._job, curState._mach);
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

#ifdef COUNT_STEPS
		numSteps++;
		currentShownSeconds = (duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count()) / 1000;
		if (currentShownSeconds > shownSeconds) {
			shownSeconds = currentShownSeconds;
			cout << "TIME" << preString << " " << shownSeconds << " " << numSteps << endl;
		}
#endif //COUNT_STEPS

		//STEP Go to neighbour			
		if (!cands.empty())
			nsp(curState, lastOp, tabuList, cands, theState.penalties, dists, prev, indeg, Q, heads, tails, lbOrder, cplex);
		else
			emptyNeighbourhood = true;

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

			if (curState.penalties < printWhenBetter)
				if (timeLog) {
					resultList.push_back(preString + unsigStr(curState.penalties) + " " + doubleStr((duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count()) / 1000.0) + " d");
				}

			trySwapNeigh = 0;
			noImproveIters = 0;
			theState = curState;
			curJumpLimit = initialjumpLimit;
			assert(curState.penalties >= lowerBound);
			if (curState.penalties == lowerBound) //Stop -> Optimal
				return true;

			newBestFound = true;

			theState.millisecsFound = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
			curState.millisecsFound = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();
		}
		else {
			if (trySwapNeigh == 1) trySwapNeigh = 0;
			else trySwapNeigh += 1;
		}

		//if (cplex) cplex = false;
	}//end while

	assert(maxMillisecs <= duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count());

	return false;
}

void Tabu::nosmuTabu(const string& instPath, const string& name, unsigned tenure, unsigned initialjumpLimit, unsigned jumpLimitDecrease, unsigned bjSize, unsigned maxD, unsigned maxC, unsigned startType, unsigned maxMillisecs, bool timeLog) {

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

	if (startType == GT)
		theState.gifflerThompson();
	else if (startType == INSA_START)
		theState.insaPsp(bigJobFirst, indeg, Q, heads, tails, invertedQ, reach);
	else
		throw errorText("Invalid start option.", "Tabu.hpp", "nosmuTabu");

	unsigned lastOp;
	theState.setMeta(dists, lastOp, prev, indeg, Q, true);

#ifdef PRINT_ONLY_IS

	cout << instPath << " " << lowerBound << " ";
	theState.printPenalties();
	return;
#endif //PRINT_ONLY_IS

	theState.millisecsFound = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

	evolveTabu(theState, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, tpStart, dists, prev, indeg, Q, maxD, maxC, 0, maxMillisecs, UINT_MAX, "", heads, tails, timeLog, false);

	theState.setMeta(dists, lastOp, prev, indeg, Q, true);

	theState.millisecsFound = duration_cast<milliseconds>(high_resolution_clock::now() - tpStart).count();

	if (!theState.verifySchedule(true))
		throw errorText(theState.toString() + "\n\t\tBad schedule !!!!!", "", "");

#ifdef PRINT_DEFAULT
	cout << instPath << " " << lowerBound << " ";
#endif //PRINT_ONLY_RESULT
	theState.printPenalties();
}