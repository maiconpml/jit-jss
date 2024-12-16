#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include <boost/multi_array.hpp>

#include "Instance.hpp"

using namespace std;

string Inst::toString() const {
	string str = "";

	str.append("J: " + unsigStr(J));
	str.append(" ; M: " + unsigStr(M));
	str.append(" ; O: " + unsigStr(O));

	str.append("\nP:");
	for (unsigned u = 1; u < O; u++)
		str.append(unsigStr(u) + "(" + unsigStr(P[u]) + ")  ");

	str.append("\njSize:");
	for (unsigned n : jSize)
		str.append(" " + unsigStr(n));
	str.append("\nmSize:");
	for (unsigned n : mSize)
		str.append(" " + unsigStr(n));

	str.append("\njobOpers: ");
	for (const vector<unsigned>& aJob : jobOpers) {
		str.append("\n\t");
		for (unsigned o : aJob)
			str.append(" " + unsigStr(o));
	}
	str.append("\nmachOpers: ");
	for (const vector<unsigned>& aMach : machOpers) {
		str.append("\n\t");
		for (unsigned o : aMach)
			str.append(" " + unsigStr(o));
	}

	str.append("\noperToJ: ");
	for (unsigned pos : operToJ)
		str.append(" " + unsigStr(pos));
	str.append("\noperToM: ");
	for (unsigned pos : operToM)
		str.append(" " + unsigStr(pos));

	str.append("\njmToIndex: ");
	assert(jmToIndex.size() == J);
	for (unsigned j = 0; j < J; j++) {
		assert(jmToIndex[j].size() == M);
		str.append("\n\t");
		for (unsigned m = 0; m < M; m++) {
			str.append(" " + unsigStr(jmToIndex[j][m]));
		}
	}

	str.append("\nroots: ");
	for (unsigned o : roots)
		str.append(" " + unsigStr(o));
	str.append("\nleafs: ");
	for (unsigned o : leafs)
		str.append(" " + unsigStr(o));

	str.append("\nnext: ");
	assert(next.size() == O);
	for (unsigned o = 0; o < O; o++) {
		str.append("    " + unsigStr(o) + ":");
		for (unsigned nextO : next[o])
			str.append("," + unsigStr(nextO));
	}

	str.append("\nprev: ");
	assert(prev.size() == O);
	for (unsigned o = 0; o < O; o++) {
		str.append("    " + unsigStr(o) + ":");
		for (unsigned prevO : prev[o])
			str.append("," + unsigStr(prevO));
	}

	assert(posets.size() == J);
	for (const Poset& aPoset : posets)
		str.append("\n" + aPoset.toString());

	return str;
}

void Inst::setHeads(vector<unsigned>& heads, vector<unsigned>& indeg, vector<unsigned>& Q) const {
#ifndef NDEBUG
	assert(heads.size() == O);
	assert(indeg.size() == O);
	assert(Q.size() == O);
	bool found;
	fill(Q.begin(), Q.end(), 0);
#endif
	fill(heads.begin(), heads.end(), 0);

	unsigned insertQ = 0;
	unsigned removeQ = 0;
	unsigned curOp;

	for (unsigned o = 1; o < O; o++)
		indeg[o] = prev[o].size();

	for (unsigned r : roots) {
		if (indeg[r] == 0) {
			Q[insertQ++] = r;
		}
	}

#ifndef NDEBUG
	for (unsigned o = 1; o < O; o++) {
		if (indeg[o] == 0) {
			found = false;
			for (unsigned r : roots) {
				if (o == r) {
					found = true;
					break;
				}
			}
			assert(found);
		}
	}
#endif

	while (removeQ < insertQ) {
		curOp = Q[removeQ++];

		assert(curOp != 0);
		assert(indeg[curOp] == 0);

		for (unsigned p : next[curOp]) {
			assert(indeg[p] > 0);
			indeg[p]--;
			if (indeg[p] == 0)
				Q[insertQ++] = p;

			heads[p] = max(heads[p], heads[curOp] + P[curOp]);
		}
	}
#ifndef NDEBUG
	assert(insertQ == O - 1);//no cycle
	for (unsigned o = 1; o < O; o++) {
		if (prev[o].size() == 0)
			assert(heads[o] == 0);//roots
		else
			assert(heads[o] != 0);//not roots
	}
#endif
}

void Inst::setTails(vector<unsigned>& tails, vector<unsigned>& indeg, vector<unsigned>& Q) const {
#ifndef NDEBUG
	assert(tails.size() == O);
	assert(indeg.size() == O);
	assert(Q.size() == O);
	bool found;
	fill(Q.begin(), Q.end(), 0);
#endif
	fill(tails.begin(), tails.end(), 0);

	unsigned insertQ = 0;
	unsigned removeQ = 0;
	unsigned curOp;

	for (unsigned o = 1; o < O; o++)
		indeg[o] = next[o].size();

	for (unsigned l : leafs) {
		if (indeg[l] == 0) {
			Q[insertQ++] = l;
		}
	}

#ifndef NDEBUG
	for (unsigned o = 1; o < O; o++) {
		if (indeg[o] == 0) {
			found = false;
			for (unsigned l : leafs) {
				if (o == l) {
					found = true;
					break;
				}
			}
			assert(found);
		}
	}
#endif

	while (removeQ < insertQ) {
		curOp = Q[removeQ++];

		assert(curOp != 0);
		assert(indeg[curOp] == 0);

		for (unsigned p : prev[curOp]) {
			assert(indeg[p] > 0);
			indeg[p]--;
			if (indeg[p] == 0)
				Q[insertQ++] = p;

			tails[p] = max(tails[p], tails[curOp] + P[curOp]);
		}
	}

#ifndef NDEBUG
	assert(insertQ == O - 1);//no cycle
	for (unsigned o = 1; o < O; o++) {
		if (next[o].size() == 0)
			assert(tails[o] == 0);//leafs
		else
			assert(tails[o] != 0);//not leafs
	}
#endif
}

//max for each machine: of min head + min tail + procTs of that machine
unsigned Inst::lowerBoundNasiri(vector<unsigned>& indeg, vector<unsigned>& Q) const {
	assert(O > 0);
	assert(indeg.size() == O);
	assert(Q.size() == O);

	vector<unsigned> heads(O);
	vector<unsigned> tails(O);
	//vector<unsigned> indeg(O);
	//vector<unsigned> Q(O);

	unsigned minHead;
	unsigned minTail;
	unsigned machT;
	unsigned lb = 0;

	setHeads(heads, indeg, Q);
	setTails(tails, indeg, Q);

	for (const vector<unsigned>& aMach : machOpers) {
		minHead = UINT_MAX;
		minTail = UINT_MAX;
		machT = 0;

		for (unsigned o : aMach) {
			assert(o != 0);
			assert(o < O);
			minHead = min(minHead, heads[o]);
			minTail = min(minTail, tails[o]);
			machT += P[o];
		}
		lb = max(lb, minHead + machT + minTail);
	}
	return lb;
}

//max of lowerBoundnasiri and job lower bound
unsigned Inst::lowerBoundTkz(vector<unsigned>& indeg, vector<unsigned>& Q) const {
	unsigned lb = lowerBoundNasiri(indeg, Q);
	unsigned jobSum;

	for (unsigned j = 0; j < J; j++) {
		jobSum = 0;
		for (unsigned o : jobOpers[j])
			jobSum += P[o];
		lb = max(lb, jobSum);
	}

	return lb;
}

//max of lowerBoundnasiri and job lower bound
unsigned Inst::lowerBoundTkz() const {
	vector<unsigned> indeg(O);
	vector<unsigned> Q(O);
	return lowerBoundTkz(indeg, Q);
}

unsigned Inst::lowerBoundMss() const {
	unsigned sum;
	unsigned lb = 0;

	for (unsigned j = 0; j < J; j++) {
		sum = 0;
		for (unsigned o : jobOpers[j])
			sum += P[o];
		lb = max(lb, sum);
	}
	for (unsigned m = 0; m < M; m++) {
		sum = 0;
		for (unsigned o : machOpers[m])
			sum += P[o];
		lb = max(lb, sum);
	}

	return lb;
}

// Verify solution given by starts
bool Inst::verifySchedule(const vector<unsigned>& starts, unsigned expecMakes) const {
	unsigned foundMakes = 0;

	for (unsigned o1 = 1; o1 < O; o1++) {
		for (unsigned o2 = 1; o2 < O; o2++) {
			if (o1 != o2) {
				//previous were ok?
				if (hasPrec(o1, o2)) {
					if (starts[o1] + P[o1] > starts[o2]) {
						cout << "Prev not OK: " << o1 << "->" << o2 << endl;
						return false;
					}
				}

				//Job lock
				if (operToJ[o1] == operToJ[o2]) {
					if (starts[o1] < starts[o2]) {
						if (starts[o1] + P[o1] > starts[o2]) {
							cout << "job lock not OK: " << o1 << "-" << o2 << endl;
							return false;
						}
					}
					else {
						if (starts[o2] + P[o2] > starts[o1]) {
							cout << "job lock not OK: " << o1 << "-" << o2 << endl;
							return false;
						}
					}
				}

				//mach lock
				if (operToM[o1] == operToM[o2]) {
					if (starts[o1] < starts[o2]) {
						if (starts[o1] + P[o1] > starts[o2]) {
							cout << "mach lock not OK: " << o1 << "-" << o2 << endl;
							return false;
						}
					}
					else {
						if (starts[o2] + P[o2] > starts[o1]) {
							cout << "mach lock not OK: " << o1 << "-" << o2 << endl;
							return false;
						}
					}
				}
			}
		}
	}

	for (unsigned o = 1; o < O; o++)
		foundMakes = max(foundMakes, starts[o] + P[o]);

	if (foundMakes != expecMakes) {
		cout << "makes not OK: " << endl;
		return false;
	}

	return true;
}

// Generate format file to be consumed by r library to generate gantt chart
void Inst::generateGantt(const vector<unsigned>& starts) {

	cout << "machines tasks timeStarts timeEnds penaltie\n";

	for (unsigned i = 1; i < O; ++i) {
		cout << operToM[i] << " " << operToJ[i] << " " << starts[i] << " " << starts[i] + P[i] << " ";
		if (starts[i] + P[i] < deadlines[i]) cout << "1";
		else if (starts[i] + P[i] > deadlines[i]) cout << "2";
		else cout << "0";
		cout << endl;
	}
}

// Calculate earliness and tardiness of schedule given by starts
void Inst::calcPenalties(const vector<unsigned>& starts, double& ePenalty, double& lPenalty, vector<double>& operPenalties) {
	ePenalty = 0;
	lPenalty = 0;
	int curDueDate;
	int curStart;
	double curTardiness;
	double curEarliness;

	operPenalties.resize(O, 0);

	for (unsigned o = 1; o < O; ++o) {

		curStart = starts[o];
		curDueDate = deadlines[o];

		curEarliness = max(curDueDate - (curStart + (int)P[o]), 0) * earlPenalties[o];
		curTardiness = max((curStart + (int)P[o]) - curDueDate, 0) * tardPenalties[o];

		assert(curEarliness >= 0);
		assert(curTardiness >= 0);
		ePenalty += curEarliness;
		lPenalty += curTardiness;
		operPenalties[o] = curEarliness > 0 ? curEarliness : curTardiness;
	}
}

// Print in format for external verifier
void Inst::printOutForTest(const vector<unsigned>& starts) {
	cout << M << " " << J << endl;
	for (unsigned a = 0; a < J; a++) {
		for (unsigned b = 0; b < M; b++) {
			int o = -1;
			unsigned auxI = jmToIndex[a][b];
			for (unsigned i = 0; i < M; i++) {
				if (jobOpers[a][i] == auxI) {
					o = i;
					break;
				}
			}
			cout << b << " " << a << " " << P[auxI] << " " << starts[auxI] << " " << o << endl;
		}
	}
}

// Print solution
void Inst::printPenalties(const vector<unsigned>& starts, const unsigned& makes) {

	double sumTardPenalties = 0;
	double sumEarlPenalties = 0;
	double sumPenalties = 0;
	int curDueDate;
	int curStart;
	double curTardiness;
	double curEarliness;

#ifdef RFILE
	generateGantt(starts);
	return;
#endif //RFILE

	for (unsigned o = 1; o < O; ++o) {

		curStart = starts[o];
		curDueDate = deadlines[o];

		curEarliness = max(curDueDate - (curStart + (int)P[o]), 0);
		curTardiness = max((curStart + (int)P[o]) - curDueDate, 0);

		assert(curEarliness >= 0);
		assert(curTardiness >= 0);

		sumEarlPenalties += curEarliness * earlPenalties[o];
		sumTardPenalties += curTardiness * tardPenalties[o];
	}

	sumPenalties = sumEarlPenalties + sumTardPenalties;

#ifdef PRINT_SCHEDULE

	for (unsigned j = 0; j < J; ++j) {

		printf("%02d - ", j);

		for (unsigned m = 0; m < M; ++m) {
			int i = jmToIndex[j][m];
			curStart = starts[i];
			curDueDate = deadlines[i];
			curEarliness = max(curDueDate - (curStart + (int)P[i]), 0);
			curTardiness = max((curStart + (int)P[i]) - curDueDate, 0);
			printf("[%03u|%03u|%03u|%2.0f|%2.0f] ", starts[i], P[i], deadlines[i], curEarliness, curTardiness);
		}
		cout << endl;
	}

	return;
#endif //PRINT_SCHEDULE

#ifdef PRINT_DEFAULT
	cout << sumPenalties << " " << sumEarlPenalties << " " << sumTardPenalties << " ";
#endif //PRINT_ONLY_RESULT
#ifdef PRINT_ONLY_RESULT
	cout << sumPenalties;
#endif // PRINT_ONLY_RESULT


#ifdef PRINT_NEIGHBOURS_NB
	double mean = 0;
	for (double n : neigh) {
		mean += n;
	}
	if (neigh.size()) {
		mean = mean / neigh.size();
	}
	cout << mean << " ";
#endif // PRINT_NEIGHBOURS_NB


#ifdef PRINT_VERIFY_OUTPUT 
	printOutForTest(starts);

#endif
	cout << endl;
}

//o1 preceedes o2?
bool Inst::hasPrec(unsigned o1, unsigned o2) const {
	assert(o1 < O && o1 != 0);
	assert(o2 < O && o2 != 0);
	assert(o1 != 0 && o2 != 0);

	if (operToJ[o1] != operToJ[o2])
		return false;
	return posets[operToJ[o1]].cg.path(operToM[o1], operToM[o2]);
}

// reads file instance
void Inst::parse(const string& filePath) {
	P.clear();
	operToJ.clear();
	operToM.clear();
	jSize.clear();
	mSize.clear();
	jobOpers.clear();
	machOpers.clear();
	deadlines.clear();
	earlPenalties.clear();
	tardPenalties.clear();
	string bufferStr;
	double bufferT;
	vector<pair<unsigned, unsigned>> order;
	unsigned fstM, sndM;
	unsigned fstO, sndO;
	vector<unsigned> fromDummy;
	vector<unsigned> toDummy;

	ifstream streamFromFile;
	streamFromFile.open(filePath);
	if (!streamFromFile.is_open())
		throw errorText("Could not open instance file: " + filePath, "Instance.hpp", "Inst::parse()");

	streamFromFile >> J;

	streamFromFile >> M;

	P.push_back(0);
	deadlines.push_back(0);
	earlPenalties.push_back(0);
	tardPenalties.push_back(0);
	operToJ.push_back(UINT_MAX);
	operToM.push_back(UINT_MAX);

	//initializing
	O = 1;
	jSize.resize(J, 0);
	mSize.resize(M, 0);

	jmToIndex.resize(extents[J][M]);
	jobOpers.resize(J);
	machOpers.resize(M);

	unsigned auxM;

	for (unsigned j = 0; j < J; j++) {
		for (unsigned m = 0; m < M; m++) {
			streamFromFile >> auxM;
			streamFromFile >> bufferT;

			assert(bufferT != 0);

			P.push_back(bufferT);
			assert(P[O] == bufferT);
			jmToIndex[j][auxM] = O;
			jobOpers[j].push_back(O);
			machOpers[auxM].push_back(O);

			streamFromFile >> bufferT;
			deadlines.push_back(bufferT);
			assert(deadlines[O] == bufferT);
			streamFromFile >> bufferT;
			earlPenalties.push_back(bufferT);
			assert(earlPenalties[O] == bufferT);
			streamFromFile >> bufferT;
			tardPenalties.push_back(bufferT);
			assert(tardPenalties[O] == bufferT);

			O++;
			operToJ.push_back(j);
			operToM.push_back(auxM);
			jSize[j]++;
			mSize[auxM]++;
		}
	}

	assert(O <= J * M + 1);
	assert(operToJ.size() == O);
	assert(operToM.size() == O);

	streamFromFile.close();

	assert(!streamFromFile.is_open());

	streamFromFile.open(filePath);

	assert(streamFromFile.is_open());

	//Meta that needs O and P
	next.resize(O);
	prev.resize(O);

	streamFromFile >> bufferT;
	streamFromFile >> bufferT;

	//Reading precedences
	for (unsigned j = 0; j < J; j++) {

		order.clear();

		fromDummy.clear();
		toDummy.clear();

		streamFromFile >> fstM;
		streamFromFile >> bufferT;
		streamFromFile >> bufferT;
		streamFromFile >> bufferT;
		streamFromFile >> bufferT;
		for (unsigned u = 1; u < M; u++) {
			assert(fstM < M);
			fstO = jmToIndex[j][fstM];
			assert(fstO < O);

			streamFromFile >> sndM;
			assert(sndM < M);
			sndO = jmToIndex[j][sndM];
			assert(sndO < O);

			if (fstO != 0 && sndO != 0) {
				next[fstO].push_back(sndO);
				prev[sndO].push_back(fstO);
			}

			order.push_back(pair<unsigned, unsigned>(fstM, sndM));
			fstM = sndM;
			streamFromFile >> bufferT;
			streamFromFile >> bufferT;
			streamFromFile >> bufferT;
			streamFromFile >> bufferT;
		}

		posets.push_back(Poset(order, M));
	}

	for (unsigned j = 0; j < J; j++) {
		assert(!posets[j].roots.empty());
		for (unsigned m : posets[j].roots) {
			if (jmToIndex[j][m] != 0)
				roots.push_back(jmToIndex[j][m]);
		}
	}
	for (unsigned j = 0; j < J; j++) {
		assert(!posets[j].leafs.empty());
		for (unsigned m : posets[j].leafs) {
			if (jmToIndex[j][m] != 0)
				leafs.push_back(jmToIndex[j][m]);
		}
	}

	streamFromFile.close();
}