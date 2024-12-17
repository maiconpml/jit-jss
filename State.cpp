#include <iostream>
#include <queue>
#include <set>

#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

#include "State.hpp"

bool State::operator==(const State& other) const {
	assert(verify());
	assert(other.verify());

	for (unsigned pos = 0; pos < inst.O; pos++) {
		if (job[pos] != other.job[pos]) return false;
		if (_job[pos] != other._job[pos]) return false;
		if (mach[pos] != other.mach[pos]) return false;
		if (_mach[pos] != other._mach[pos]) return false;
	}

	return true;
}

void State::alloc() {
	assert(inst.O != 0);
	assert(!isAlloced());
	job.resize(inst.O);
	_job.resize(inst.O);
	mach.resize(inst.O);
	_mach.resize(inst.O);
}

void State::clear() {
	assert(inst.O != 0);
	assert(mach.size() == inst.O);
	assert(job.size() == inst.O);

	//roots.clear();
	//leafs.clear();
	fill(job.begin(), job.end(), 0);
	fill(_job.begin(), _job.end(), 0);
	fill(mach.begin(), mach.end(), 0);
	fill(_mach.begin(), _mach.end(), 0);
}

string State::toString() const {
	assert(job[0] == 0);
	assert(_job[0] == 0);
	assert(mach[0] == 0);
	assert(_mach[0] == 0);

	string str = "makes: " + unsigStr(makes);

	str.append("\njob: ");
	for (unsigned o = 1; o < job.size(); o++) {
		str.append("  " + unsigStr(o) + "->" + unsigStr(job[o]));
		assert(job[o] == 0 || _job[job[o]] == o);
	}

	str.append("\nmach:");
	for (unsigned o = 1; o < mach.size(); o++) {
		str.append("  " + unsigStr(o) + "->" + unsigStr(mach[o]));
		assert(mach[o] == 0 || _mach[mach[o]] == o);
	}

	return str;
}

//expects fully defined state
string State::easyString() const {
	assert(job[0] == 0);
	assert(_job[0] == 0);
	assert(mach[0] == 0);
	assert(_mach[0] == 0);

	unsigned root;
	unsigned curOp;

	string str = "makes: " + unsigStr(makes);

	str.append("\njob:");
	for (unsigned j = 0; j < inst.J; j++) {
		root = 0;
		for (unsigned o : inst.jobOpers[j]) {
			if (_job[o] == 0) {
				assert(root == 0);
				root = o;
			}
		}
		curOp = root;
		str.append("\n\t");
		for (unsigned u = 0; u < inst.jSize[j]; u++) {
			assert(curOp != 0);
			str.append(unsigStr(curOp) + " ");
			curOp = job[curOp];
		}
	}

	str.append("\nmach:");
	for (unsigned m = 0; m < inst.M; m++) {
		root = 0;
		for (unsigned o : inst.machOpers[m]) {
			if (_mach[o] == 0) {
				assert(root == 0);
				root = o;
			}
		}
		curOp = root;
		str.append("\n\t");
		for (unsigned u = 0; u < inst.mSize[m]; u++) {
			assert(curOp != 0);
			str.append(unsigStr(curOp) + " ");
			curOp = mach[curOp];
		}
	}

	return str;
}

//x and y in same job 
//max(head(Pm[y]), head(Pj[x]))+procT[X]+procT[Y]+max(tail(Sj]y]), tail(Sm[x]))
//max(head(Pj[x], head(Pm[y])) +procT[y] + tail(Sm[x])
//head(Pm[x]) + procT[x] +max(tail(Sj[y], tail(Sm[y])
unsigned State::swapLowerBound(const unsigned x, const unsigned y, const vector<unsigned>& heads, const vector<unsigned>& tails) const {
	assert(heads.size() == inst.O);
	assert(tails.size() == inst.O);
	assert(x != 0);
	assert(y != 0);
	assert(inst.operToJ[x] == inst.operToJ[y] || inst.operToM[x] == inst.operToM[y]);
	assert((job[x] == y && _job[y] == x) || (mach[x] == y && _mach[y] == x));

	unsigned hA1;
	unsigned hA2;
	unsigned tB1;
	unsigned tB2;
	unsigned hAlpha;
	unsigned tBetha;
	unsigned hGamma;
	unsigned tDelta;
	unsigned lb;

	if (inst.operToJ[x] == inst.operToJ[y]) {
		hA1 = (_job[x] != 0 ? heads[_job[x]] + inst.P[_job[x]] : 0);
		hA2 = (_mach[y] != 0 ? heads[_mach[y]] + inst.P[_mach[y]] : 0);
		tB1 = (job[y] != 0 ? tails[job[y]] + inst.P[job[y]] : 0);
		tB2 = (mach[x] != 0 ? tails[mach[x]] + inst.P[mach[x]] : 0);
		hGamma = (_mach[x] != 0 ? heads[_mach[x]] + inst.P[_mach[x]] : 0);
		tDelta = (mach[y] != 0 ? tails[mach[y]] + inst.P[mach[y]] : 0);
	}
	else {//same mach
		hA1 = (_mach[x] != 0 ? heads[_mach[x]] + inst.P[_mach[x]] : 0);
		hA2 = (_job[y] != 0 ? heads[_job[y]] + inst.P[_job[y]] : 0);
		tB1 = (mach[y] != 0 ? tails[mach[y]] + inst.P[mach[y]] : 0);
		tB2 = (job[x] != 0 ? tails[job[x]] + inst.P[job[x]] : 0);
		hGamma = (_job[x] != 0 ? heads[_job[x]] + inst.P[_job[x]] : 0);
		tDelta = (job[y] != 0 ? tails[job[y]] + inst.P[job[y]] : 0);
	}
	hAlpha = max(hA1, hA2);
	tBetha = max(tB1, tB2);

	lb = hAlpha + inst.P[y] + inst.P[x] + tBetha;
	lb = max(lb, hAlpha + inst.P[y] + tDelta);
	lb = max(lb, hGamma + inst.P[x] + tBetha);

	return lb;
}

void State::swap(unsigned o1, unsigned o2) {
	assert(o1 != 0);
	assert(o2 != 0);
	assert(o1 != o2);

	unsigned prev;
	unsigned post;

#ifdef VANILLA_PSS
	if (inst.operToJ[o1] == inst.operToJ[o2]) {
		//JOB swap
		assert(job[o1] == o2);
		assert(_job[o2] == o1);
		prev = _job[o1];
		post = job[o2];
		//forward
		if (prev != 0)
			job[prev] = o2;
		job[o1] = post;
		job[o2] = o1;
		//backward
		_job[o1] = o2;
		_job[o2] = prev;
		if (post != 0)
			_job[post] = o1;
	}
	else {
		assert(inst.operToM[o1] == inst.operToM[o2]);
		//MACH swap
		assert(mach[o1] == o2);
		assert(_mach[o2] == o1);
		prev = _mach[o1];
		post = mach[o2];
		//forward
		if (prev != 0)
			mach[prev] = o2;
		mach[o1] = post;
		mach[o2] = o1;
		//backward
		_mach[o1] = o2;
		_mach[o2] = prev;
		if (post != 0)
			_mach[post] = o1;
	}
#endif //VANILLA_PSS
}

//partial state has cycles?
bool State::hasCycle(vector<unsigned>& indeg, vector<unsigned>& Q) const {

	assert(isAlloced());
	assert(indeg.size() == inst.O);
	assert(Q.size() == inst.O);

	unsigned qInsert = 0;
	unsigned qAccess = 0;

	unsigned curOp;
	unsigned newOp;

	fill(indeg.begin(), indeg.end(), 0);

	for (unsigned o = 1; o < inst.O; o++) {
		indeg[o] = inst.prev[o].size();
		if (_job[o] != 0)
			indeg[o]++;
		if (_mach[o] != 0)
			indeg[o]++;
		if (indeg[o] == 0)
			Q[qInsert++] = o;
	}


	while (qAccess < qInsert) {
		assert(qAccess < Q.size());
		curOp = Q[qAccess++];
		assert(indeg[curOp] == 0);

		//job order
		newOp = job[curOp];
		if (newOp != 0) {
			assert(indeg[newOp] > 0);
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = newOp;
			}
		}
		//mach order
		newOp = mach[curOp];
		if (newOp != 0) {
			assert(indeg[newOp] > 0);
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = newOp;
			}
		}
		//job precs
		for (unsigned n : inst.next[curOp]) {
			assert(indeg[n] > 0);
			indeg[n]--;
			if (indeg[n] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = n;
			}
		}
	}
	assert(qAccess <= inst.O - 1);
	return qAccess < inst.O - 1;
}

void State::removeOper(unsigned targOp) {
	assert(targOp != 0);

	const unsigned prevJobOp = _job[targOp];
	assert(prevJobOp == 0 || job[prevJobOp] == targOp);
	const unsigned postJobOp = job[targOp];
	assert(postJobOp == 0 || _job[postJobOp] == targOp);
	const unsigned prevMachOp = _mach[targOp];
	assert(prevMachOp == 0 || mach[prevMachOp] == targOp);
	const unsigned postMachOp = mach[targOp];
	assert(postMachOp == 0 || _mach[postMachOp] == targOp);

	job[targOp] = 0;
	_job[targOp] = 0;
	mach[targOp] = 0;
	_mach[targOp] = 0;
	if (prevJobOp != 0) job[prevJobOp] = postJobOp;
	if (postJobOp != 0) _job[postJobOp] = prevJobOp;
	if (prevMachOp != 0) mach[prevMachOp] = postMachOp;
	if (postMachOp != 0) _mach[postMachOp] = prevMachOp;
}

void State::removeOperFromMach(unsigned targOp) {
	assert(targOp != 0);

	const unsigned prevMachOp = _mach[targOp];
	assert(prevMachOp == 0 || mach[prevMachOp] == targOp);
	const unsigned postMachOp = mach[targOp];
	assert(postMachOp == 0 || _mach[postMachOp] == targOp);

	mach[targOp] = 0;
	_mach[targOp] = 0;
	if (prevMachOp != 0) mach[prevMachOp] = postMachOp;
	if (postMachOp != 0) _mach[postMachOp] = prevMachOp;
}

//insert newOp right after prevOp and before postOp --- needs both in case prev is 0 needs to know cur root
//care if using roots and leafs - must adjust outside
void State::insertOper(unsigned newOp, unsigned prevJobOp, unsigned postJobOp, unsigned prevMachOp, unsigned postMachOp) {
	insertJobOper(newOp, prevJobOp, postJobOp);
	insertMachOper(newOp, prevMachOp, postMachOp);
}

void State::insertJobOper(unsigned newOp, unsigned prevJobOp, unsigned postJobOp) {
	assert(newOp != 0);
	//assert(prevJobOp==0   ||   job[prevJobOp]==postJobOp);
	//assert(postJobOp==0   ||   _job[postJobOp]==prevJobOp);

	job[newOp] = postJobOp;
	if (postJobOp != 0) _job[postJobOp] = newOp;
	_job[newOp] = prevJobOp;
	if (prevJobOp != 0) job[prevJobOp] = newOp;
}

void State::insertMachOper(unsigned newOp, unsigned prevMachOp, unsigned postMachOp) {
	assert(newOp != 0);
	//assert(prevMachOp==0   ||   mach[prevMachOp]==postMachOp);
	//assert(postMachOp==0   ||   _mach[postMachOp]==prevMachOp);

	mach[newOp] = postMachOp;
	if (postMachOp != 0) _mach[postMachOp] = newOp;
	_mach[newOp] = prevMachOp;
	if (prevMachOp != 0) mach[prevMachOp] = newOp;
}

bool State::compHeadsComplete(vector<unsigned>& heads, vector<unsigned>& indeg, vector<unsigned>& Q) const {
	assert(isAlloced());
	assert(heads.size() == inst.O);
	assert(indeg.size() == inst.O);
	assert(Q.size() == inst.O);

	unsigned qInsert = 0;
	unsigned qAccess = 0;

	unsigned curOp;
	unsigned newOp;

	unsigned newMax;

	fill(heads.begin(), heads.end(), 0);
	fill(indeg.begin(), indeg.end(), 0);

	for (unsigned o = 1; o < inst.O; o++) {
		if (_job[o] != 0)
			indeg[o]++;
		if (_mach[o] != 0)
			indeg[o]++;
		if (indeg[o] == 0)
			Q[qInsert++] = o;
	}

	assert(qInsert > 0);

	while (qAccess < qInsert) {
		assert(qAccess < Q.size());
		curOp = Q[qAccess++];
		assert(indeg[curOp] == 0);
		newMax = heads[curOp] + inst.P[curOp];

		//job order
		newOp = job[curOp];
		if (newOp != 0) {
			assert(indeg[newOp] > 0);
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = newOp;
			}
			heads[newOp] = max(heads[newOp], newMax);
		}
		//mach order
		newOp = mach[curOp];
		if (newOp != 0) {
			assert(indeg[newOp] > 0);
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = newOp;
			}
			heads[newOp] = max(heads[newOp], newMax);
		}
	}

	assert(qAccess <= inst.O - 1);
	return qAccess < inst.O - 1;
}

bool State::compHeads(vector<unsigned>& heads, vector<unsigned>& indeg, vector<unsigned>& Q) const {
	assert(isAlloced());
	assert(heads.size() == inst.O);
	assert(indeg.size() == inst.O);
	assert(Q.size() == inst.O);

	unsigned qInsert = 0;
	unsigned qAccess = 0;

	unsigned curOp;
	unsigned newOp;

	unsigned newMax;

	fill(heads.begin(), heads.end(), 0);

	for (unsigned o = 1; o < inst.O; o++) {
		indeg[o] = inst.prev[o].size();
		if (_job[o] != 0)
			indeg[o]++;
		if (_mach[o] != 0)
			indeg[o]++;
		if (indeg[o] == 0)
			Q[qInsert++] = o;
	}

	//assert(qInsert>0);//cycle may reate no first op in topo

	while (qAccess < qInsert) {
		assert(qAccess < Q.size());
		curOp = Q[qAccess++];
		assert(indeg[curOp] == 0);
		newMax = heads[curOp] + inst.P[curOp];

		//job order
		newOp = job[curOp];
		if (newOp != 0) {
			assert(indeg[newOp] > 0);
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = newOp;
			}
			heads[newOp] = max(heads[newOp], newMax);
		}
		//mach order
		newOp = mach[curOp];
		if (newOp != 0) {
			assert(indeg[newOp] > 0);
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = newOp;
			}
			heads[newOp] = max(heads[newOp], newMax);
		}
		//job precs
		for (unsigned n : inst.next[curOp]) {
			assert(indeg[n] > 0);
			indeg[n]--;
			if (indeg[n] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = n;
			}
			heads[n] = max(heads[n], newMax);
		}
	}

	assert(qAccess <= inst.O - 1);
	return qAccess < inst.O - 1;
}

bool State::compTailsComplete(vector<unsigned>& tails, vector<unsigned>& indeg, vector<unsigned>& Q) const {
	assert(isAlloced());
	assert(tails.size() == inst.O);
	assert(indeg.size() == inst.O);
	assert(Q.size() == inst.O);
	unsigned qInsert = 0;
	unsigned qAccess = 0;

	unsigned curOp;
	unsigned newOp;

	unsigned newMax;

	fill(tails.begin(), tails.end(), 0);
	fill(indeg.begin(), indeg.end(), 0);

	for (unsigned o = 1; o < inst.O; o++) {
		if (job[o] != 0)
			indeg[o]++;
		if (mach[o] != 0)
			indeg[o]++;
		if (indeg[o] == 0)
			Q[qInsert++] = o;
	}

	assert(qInsert > 0);

	while (qAccess < qInsert) {
		assert(qAccess < Q.size());
		curOp = Q[qAccess++];
		assert(indeg[curOp] == 0);
		newMax = tails[curOp] + inst.P[curOp];

		//job order
		newOp = _job[curOp];
		if (newOp != 0) {
			assert(indeg[newOp] > 0);
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = newOp;
			}
			tails[newOp] = max(tails[newOp], newMax);
		}
		//mach order
		newOp = _mach[curOp];
		if (newOp != 0) {
			assert(indeg[newOp] > 0);
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = newOp;
			}
			tails[newOp] = max(tails[newOp], newMax);
		}
	}

	assert(qAccess <= inst.O - 1);
	return qAccess < inst.O - 1;
}

bool State::compTails(vector<unsigned>& tails, vector<unsigned>& indeg, vector<unsigned>& Q) const {
	assert(isAlloced());
	assert(tails.size() == inst.O);
	assert(indeg.size() == inst.O);
	assert(Q.size() == inst.O);
	unsigned qInsert = 0;
	unsigned qAccess = 0;

	unsigned curOp;
	unsigned newOp;

	unsigned newMax;

	fill(tails.begin(), tails.end(), 0);

	for (unsigned o = 1; o < inst.O; o++) {
		indeg[o] = inst.next[o].size();
		if (job[o] != 0)
			indeg[o]++;
		if (mach[o] != 0)
			indeg[o]++;
		if (indeg[o] == 0)
			Q[qInsert++] = o;
	}

	//assert(qInsert>0);

	while (qAccess < qInsert) {
		assert(qAccess < Q.size());
		curOp = Q[qAccess++];
		assert(indeg[curOp] == 0);
		newMax = tails[curOp] + inst.P[curOp];

		//job order
		newOp = _job[curOp];
		if (newOp != 0) {
			assert(indeg[newOp] > 0);
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = newOp;
			}
			tails[newOp] = max(tails[newOp], newMax);
		}
		//mach order
		newOp = _mach[curOp];
		if (newOp != 0) {
			assert(indeg[newOp] > 0);
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = newOp;
			}
			tails[newOp] = max(tails[newOp], newMax);
		}
		//job precs
		for (unsigned n : inst.prev[curOp]) {
			assert(indeg[n] > 0);
			indeg[n]--;
			if (indeg[n] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = n;
			}
			tails[n] = max(tails[n], newMax);
		}
	}

	assert(qAccess <= inst.O - 1);
	return qAccess < inst.O - 1;
}

//can reach fromIndex toIndex in Q?
bool State::propagateReach(const vector<unsigned>& Q, vector<bool>& reach, unsigned fromIndexA, unsigned fromIndexB, unsigned toIndexA, unsigned toIndexB) const {
	unsigned curOp;
	unsigned newOp;
	unsigned minIndex = min(fromIndexA, fromIndexB);
	unsigned maxIndex = max(toIndexA, toIndexB);

	unsigned targA = Q[toIndexA];
	unsigned targB = Q[toIndexB];

	//cout << "\t\t\t\ttargA: " << targA << "  ;  targB: " << targB << endl;

	fill(reach.begin(), reach.end(), false);
	if (Q[fromIndexA] != 0)
		reach[Q[fromIndexA]] = true;
	if (Q[fromIndexB] != 0)
		reach[Q[fromIndexB]] = true;

	for (unsigned topoIndex = minIndex; topoIndex < maxIndex; topoIndex++) {
		curOp = Q[topoIndex];
		//assert(curOp != targA   && curOp != targB);
		if (reach[curOp]) {
			//job order
			newOp = job[curOp];
			if (newOp != 0) {
				if (newOp == targA || newOp == targB) {
					//tcout << "\t\t\t\t\tX - newOp" << newOp << " ; curOp: " << curOp  << endl;
					return true;
				}
				reach[newOp] = true;
			}
			//mach order
			newOp = mach[curOp];
			if (newOp != 0) {
				if (newOp == targA || newOp == targB) {
					//cout << "\t\t\t\t\tY - newOp"  << newOp << " ; curOp: " << curOp  << endl;
					return true;
				}
				reach[newOp] = true;
			}
			//job precs
			for (unsigned n : inst.next[curOp]) {
				assert(n != 0);
				if (n == targA || n == targB) {
					//cout << "\t\t\t\t\tZ - newOp"  << newOp << " ; curOp: " << curOp  << endl;
					return true;
				}
				reach[n] = true;
			}
		}
	}

	return false;
}

#ifdef VANILLA_PSS

void State::insertInsaPos(unsigned curOp, vector<unsigned>& indeg, vector<unsigned>& Q, vector<bool>& inOpers, vector<unsigned>& jobRoot, vector<unsigned>& jobLeaf, vector<unsigned>& machRoot, vector<unsigned>& machLeaf, vector<unsigned>& heads, vector<unsigned>& tails, vector<unsigned>& invertedQ, vector<bool>& reach) {
#ifndef NDEBUG
	bool cycle;
	vector<unsigned> auxQ(inst.O);
	vector<unsigned> auxIndeg(inst.O);
	assert(curOp != 0 && curOp < inst.O);
	assert(indeg.size() == inst.O);
	assert(Q.size() == inst.O);
	assert(inOpers.size() == inst.O);
	assert(!inOpers[curOp]);
	assert(jobRoot.size() == inst.J);
	assert(jobLeaf.size() == inst.J);
	assert(machRoot.size() == inst.M);
	assert(machLeaf.size() == inst.M);
	for (unsigned jr : jobRoot) {
		assert(jr == 0 || inOpers[jr]);
		assert(_job[jr] == 0);
	}
	for (unsigned jl : jobLeaf) {
		assert(jl == 0 || inOpers[jl]);
		assert(job[jl] == 0);
	}
	for (unsigned mr : machRoot) {
		assert(mr == 0 || inOpers[mr]);
		assert(_mach[mr] == 0);
	}
	for (unsigned ml : machLeaf) {
		assert(ml == 0 || inOpers[ml]);
		assert(mach[ml] == 0);
	}
	bool found;
	for (unsigned j = 0; j < inst.J; j++) {
		assert((jobRoot[j] && jobLeaf[j]) || (!jobRoot[j] && !jobLeaf[j]));
		found = false;
		for (unsigned o : inst.jobOpers[j]) found |= inOpers[o];
		if (found) {
			assert(jobRoot[j] != 0);
			assert(jobLeaf[j] != 0);
		}
		else {
			assert(jobRoot[j] == 0);
			assert(jobLeaf[j] == 0);
		}
	}
	for (unsigned m = 0; m < inst.M; m++) {
		assert((machRoot[m] && machLeaf[m]) || (!machRoot[m] && !machLeaf[m]));
		found = false;
		for (unsigned o : inst.machOpers[m]) found |= inOpers[o];
		if (found) {
			assert(machRoot[m] != 0);
			assert(machLeaf[m] != 0);
		}
		else {
			assert(machRoot[m] == 0);
			assert(machLeaf[m] == 0);
		}
	}
#endif

	unsigned curCost;
	const unsigned curJob = inst.operToJ[curOp];
	const unsigned curMach = inst.operToM[curOp];

	unsigned bestCost = UINT_MAX;
	unsigned bufOp;
	unsigned jobHead;
	unsigned jobTail;

	//opers that delimit where we can insert curOp in job
	unsigned fromHere = 0;
	unsigned toHere = 0;
	bool leftPrecFound = false;
	unsigned nextOper;

	unsigned bestJobPrev = 0;
	unsigned bestJobNext = 0;
	unsigned bestMachPrev = 0;
	unsigned bestMachNext = 0;

	assert(!inOpers[curOp]);
	inOpers[curOp] = true;

	//checking limits in job for curOp
	bufOp = jobRoot[curJob];
	while (bufOp != 0) {
		if (inst.hasPrec(bufOp, curOp)) {
			fromHere = bufOp;
			leftPrecFound = true;
		}
		if (inst.hasPrec(curOp, bufOp)) {
			toHere = bufOp;
			break;
		}
		bufOp = job[bufOp];
	}
	if (!leftPrecFound) {//no precs -> fromHere is null and next is root (may be also null if empty)
		assert(fromHere == 0);
		nextOper = jobRoot[curJob];
	}
	else {// precs found
		assert(jobRoot[curJob] != 0);//not empty
		assert(fromHere != 0);//something is actually prec to curOp
		nextOper = job[fromHere];
	}

#ifndef NDEBUG
	cycle =
#endif
		compTails(tails, indeg, Q);
	assert(!cycle);
#ifndef NDEBUG
	cycle =
#endif
		compHeads(heads, indeg, Q);//important that heads is after tails - Q contains a topo sort
	assert(!cycle);

	//finding inverted table of Q
	for (unsigned u = 0; u < inst.O; u++)
		invertedQ[Q[u]] = u;
	//evaluating position pairs
	insertJobOper(curOp, fromHere, nextOper);//insert in 1st valid job position
	while (true) {//job runthrough
		//finding job heads and tails --- either from order or from precs
		jobHead = heads[_job[curOp]] + inst.P[_job[curOp]];
		for (unsigned p : inst.prev[curOp])
			jobHead = max(jobHead, heads[p] + inst.P[p]);
		jobTail = tails[job[curOp]] + inst.P[job[curOp]];
		for (unsigned n : inst.next[curOp])
			jobTail = max(jobTail, tails[n] + inst.P[n]);

		insertMachOper(curOp, 0, machRoot[curMach]);//insert in 1st mach position
		while (true) {//mach runthrough
			curCost = jobHead + jobTail;
			curCost = max(curCost, jobHead + tails[mach[curOp]] + inst.P[mach[curOp]]);
			curCost = max(curCost, heads[_mach[curOp]] + inst.P[_mach[curOp]] + jobTail);
			curCost = max(curCost, heads[_mach[curOp]] + inst.P[_mach[curOp]] + tails[mach[curOp]] + inst.P[mach[curOp]]);

			if (bestCost > curCost) {
				assert(Q[invertedQ[job[curOp]]] == job[curOp]);
				assert(Q[invertedQ[_job[curOp]]] == _job[curOp]);
				assert(Q[invertedQ[mach[curOp]]] == mach[curOp]);
				assert(Q[invertedQ[_mach[curOp]]] == _mach[curOp]);
				if (!propagateReach(Q, reach, invertedQ[mach[curOp]], invertedQ[job[curOp]], invertedQ[_job[curOp]], invertedQ[_mach[curOp]])) {
					bestCost = curCost;
					bestJobPrev = _job[curOp];
					bestJobNext = job[curOp];
					bestMachPrev = _mach[curOp];
					bestMachNext = mach[curOp];
					assert(!hasCycle(auxIndeg, auxQ));
				}
				else {
					assert(hasCycle(auxIndeg, auxQ));
				}
			}

			if (mach[curOp] == 0)
				break;
			swap(curOp, mach[curOp]);
		}//mach runthrough
		removeOperFromMach(curOp);

		if (job[curOp] == toHere)
			break;
		assert(job[curOp] != 0);
		swap(curOp, job[curOp]);
	}//job runthrough
	assert(bestCost < UINT_MAX);

	//remove from last postion and reinsert in best place
	removeOper(curOp);
	insertOper(curOp, bestJobPrev, bestJobNext, bestMachPrev, bestMachNext);

	assert(partialVerify(inOpers));

	//adjust roots and leafs
	if (bestJobNext == jobRoot[curJob])
		jobRoot[curJob] = curOp;
	if (bestJobPrev == jobLeaf[curJob])
		jobLeaf[curJob] = curOp;
	if (bestMachNext == machRoot[curMach])
		machRoot[curMach] = curOp;
	if (bestMachPrev == machLeaf[curMach])
		machLeaf[curMach] = curOp;
}

#endif //VANILLA_PSS

void State::insaPsp(bool bigJobFirst, vector<unsigned>& indeg, vector<unsigned>& Q, vector<unsigned>& heads, vector<unsigned>& tails, vector<unsigned>& invertedQ, vector<bool>& reach) {
	assert(isAlloced());
	assert(isClear());
	assert(indeg.size() == inst.O);
	assert(Q.size() == inst.O);

	unsigned curOp;

	unsigned curCost;
	unsigned bigJob = UINT_MAX;
	unsigned bigJobCost = 0;

	vector<bool> inOpers(inst.O, false);
	vector<unsigned> jobRoot(inst.J, 0);
	vector<unsigned> jobLeaf(inst.J, 0);
	vector<unsigned> machRoot(inst.M, 0);
	vector<unsigned> machLeaf(inst.M, 0);

	vector<pair<unsigned, unsigned>> costVoper;
	costVoper.reserve(inst.O);


	//inserting big job first
	if (bigJobFirst) {
		for (unsigned j = 0; j < inst.J; j++) {
			curCost = 0;
			for (unsigned o : inst.jobOpers[j])
				curCost += inst.P[o];
			if (curCost > bigJobCost) {
				bigJobCost = curCost;
				bigJob = j;
			}
		}
		assert(bigJob < inst.J);
		for (unsigned o : inst.jobOpers[bigJob]) {
			insertInsaPos(o, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails, invertedQ, reach);
		}
	}

	//ordering the rest of the opers
	for (unsigned o = 1; o < inst.O; o++) {
		if (!inOpers[o])
			costVoper.push_back(pair<unsigned, unsigned>(inst.P[o], o));
	}
	sort(costVoper.rbegin(), costVoper.rend());

	assert((bigJobFirst && costVoper.size() == (inst.O - 1 - inst.jSize[bigJob]))
		|| (!bigJobFirst && costVoper.size() == (inst.O - 1)));

	//inserting missing ops
	for (unsigned pos = 0; pos < costVoper.size(); pos++) {
		curOp = costVoper[pos].second;
		insertInsaPos(curOp, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails, invertedQ, reach);
	}
}

void State::computeCritic(vector<unsigned>& critic, unsigned lastOp, const vector<unsigned>& prev) {
	assert(lastOp != 0);
	assert(lastOp < inst.O);
	assert(critic.capacity() == inst.O);

	critic.clear();

	while (lastOp != 0) {

		assert(critic.size() < inst.O);
		critic.push_back(lastOp);
		lastOp = prev[lastOp];
	}
	reverse(critic.begin(), critic.end());
}

void State::blockBegins(vector<unsigned>& jobBb, vector<unsigned>& machBb, const vector<unsigned>& critic) {
	assert(!critic.empty());
	assert(critic.size() < inst.O);
	assert(jobBb.capacity() == inst.O);
	assert(machBb.capacity() == inst.O);

	jobBb.clear();
	machBb.clear();
	jobBb.push_back(0);
	machBb.push_back(0);

	for (unsigned pos = 0; pos < critic.size() - 1; pos++) {
		assert(critic[pos] != 0);
		assert(inst.operToJ[critic[pos]] == inst.operToJ[critic[pos + 1]]
			|| inst.operToM[critic[pos]] == inst.operToM[critic[pos + 1]]);
#ifdef VANILLA_PSS
		//job BBs
		if (inst.operToJ[critic[pos]] != inst.operToJ[critic[pos + 1]])
			jobBb.push_back(pos + 1);
#endif //VANILLA_PSS 
		//mach BBs
		if (inst.operToM[critic[pos]] != inst.operToM[critic[pos + 1]])
			machBb.push_back(pos + 1);
	}
}

bool State::topoWalk(vector<unsigned>& indeg, vector<unsigned>& Q ){
	unsigned qInsert = 0;
	unsigned qAccess = 0;
	unsigned curOp;
	unsigned newOp;
	//unsigned newMax;
	makes = 0;

	fill(indeg.begin(), indeg.end(), 0);

	for (unsigned o = 1; o < inst.O; o++) {
		if (_job[o] != 0)
			indeg[o]++;
		if (_mach[o] != 0)
			indeg[o]++;
		if (indeg[o] == 0) {
			Q[qInsert++] = o;
		}
	}
		
	while(qAccess < qInsert){
		curOp = Q[qAccess++];
		newOp = job[curOp];
		
		if (newOp != 0) {
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				Q[qInsert++] = newOp;
			}
		}

		newOp = mach[curOp];
		if (newOp != 0) {
			assert(indeg[newOp] > 0);
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = newOp;
			}
		}
	}

	return qAccess < inst.O - 1;
}

void State::updateStrength(vector<unsigned>& lateCands, double & pS, double & hS ){
	pS = 0;
	hS = 0;
	for(unsigned i = 1; i < lateCands.size();i++){
		if(lateCands[i] == 1) pS += inst.earlPenalties[i];
		else if(lateCands[i] == 2) hS += inst.tardPenalties[i];
	}
}

unsigned State::calcDelayTime(vector<unsigned> & starts,vector<unsigned>& lateCands,vector<unsigned>& limited){
	unsigned t = UINT_MAX;
	unsigned m;
	vector<unsigned> limitedAux;
	for(unsigned i = 1; i< lateCands.size(); i++){
		if(!lateCands[i]) continue;
		m = UINT_MAX;
		unsigned aux = starts[i] + inst.P[i];
		if(inst.deadlines[i] > aux && inst.deadlines[i] - aux < m) m = inst.deadlines[i] - aux;
		if(mach[i] && (starts[mach[i]] > aux || !lateCands[mach[i]]) && starts[mach[i]] - aux < m) m = starts[mach[i]] - aux;
		if (job[i] && (starts[job[i]] > aux || !lateCands[job[i]]) && starts[job[i]] - aux < m) m = starts[job[i]] - aux;
		if(m < t){
			t = m;
			limitedAux.clear();
			limitedAux.push_back(i);
		} else if(m == t) limitedAux.push_back(i);;
	}
	fill(limited.begin(), limited.end(), 0);
	for (unsigned i = 0; i < limitedAux.size(); ++i) {
		limited[limitedAux[i]] = 1;
	}
	return t;
}

void State::delay(vector<unsigned> & starts,vector<unsigned>& lateCands, unsigned & t){
	for(unsigned i = 1; i < lateCands.size(); i++){
		if(lateCands[i]){
			starts[i] = starts[i] + t;
		}
	}
}

void State::update(vector<unsigned>& lateCands, vector<unsigned> & limited, vector<unsigned> & heads, vector<unsigned> & starts){
	queue<unsigned> remove;
	queue<unsigned> auxL;
	vector<unsigned> headsCpy(heads);
	vector<unsigned> headsRmv(inst.O, 0);
	unsigned iterations = heads.size();
	heads.clear();
	for (unsigned o = 0; o < iterations; ++o) {
		if(inst.P[headsCpy[o]]+ starts[headsCpy[o]] == inst.deadlines[headsCpy[o]]){
			headsRmv[headsCpy[o]] = 1;
			remove.push(o);

			while(!remove.empty()){
				unsigned c = remove.front();
				remove.pop();

				if(lateCands[mach[c]] && inst.P[mach[c]]+ starts[mach[c]] >= inst.deadlines[mach[c]])  remove.push(mach[c]);
				else if (lateCands[mach[c]]) headsCpy.push_back(mach[c]);

				if(lateCands[job[c]] && inst.P[job[c]]+ starts[job[c]] >= inst.deadlines[job[c]])  remove.push(job[c]);
				else if (lateCands[job[c]]) headsCpy.push_back(job[c]);

				if(limited[c]) limited[c] = 0;
				lateCands[c] = 0;

			}
		}
	}
	for (unsigned i = 0; i < headsCpy.size(); ++i) {
		if (!headsRmv[headsCpy[i]]) heads.push_back(headsCpy[i]);
	}

	for(unsigned i = 1; i< inst.O; i++){
		if(limited[i]) auxL.push(i);
	}

	while(!auxL.empty()){
		unsigned o = auxL.front();
		auxL.pop();
		if(inst.P[o]+ starts[o] == inst.deadlines[o] && lateCands[o]==1) lateCands[o]  = 2;
		if(job[o] && inst.P[o]+ starts[o] == starts[job[o]]){
			auxL.push(job[o]);
			if(inst.P[job[o]]+ starts[job[o]] >= inst.deadlines[job[o]]) lateCands[job[o]] = 2;
			else lateCands[job[o]] = 1;
		}

		if(mach[o] && inst.P[o]+ starts[o] == starts[mach[o]]){
			auxL.push(mach[o]);
			if(inst.P[mach[o]]+ starts[mach[o]] >= inst.deadlines[mach[o]]) lateCands[mach[o]] = 2;
			else lateCands[mach[o]] = 1;
		}
	}
}

void State::forcedDelay(vector<unsigned> & starts,vector<unsigned>& lateCands, unsigned & op){
	vector<unsigned> auxP = inst.P;
	unsigned co = UINT_MAX;
	unsigned objDelay;
	vector<unsigned> heads;
	vector<unsigned> limited(inst.O, 0);
	if(job[op]) co =  starts[job[op]];
	if(mach[op] && starts[mach[op]] < co) co =  starts[mach[op]];

	objDelay = inst.P[op] -  co ;
	starts[op] = 0;
	inst.P[op] = co;
	heads.push_back(op);
	limited[op] = 1;
	update(lateCands, limited, heads, starts);

	while(objDelay > 0){
		
		unsigned t = calcDelayTime(starts,lateCands,limited);
		if(objDelay < t) t = objDelay;

		delay(starts, lateCands, t);
		update(lateCands, limited, heads, starts);
		objDelay -= t;
	}

	starts[op] = 0;
	inst.P = auxP;
}

	//maxLen = maxD * maxC
bool State::detectRepeat(vector<int>& oldValues, unsigned& posCurPenalties, unsigned& cycleLastPos, unsigned& cycleL, double newPenalties, bool isNewBest, unsigned maxD, unsigned maxLen) {
	assert(oldValues.size() == maxD);
	assert(maxD != 0);


	// new best solution: reset cycle detector
	if (isNewBest) {
		posCurPenalties = 0;
		cycleLastPos = 0;
		cycleL = 0;
		for (unsigned u = 0; u < maxD; u++)
			oldValues[u] = -1;
	}
	// store current makespan in ring buffer

	++posCurPenalties;
	posCurPenalties = posCurPenalties % maxD;
	oldValues[posCurPenalties] = newPenalties; //posCurPenalties is latest makes pos
	//fprintf(stderr, "cur makes: %d.\n", newPenalties);
	// if current makespan occurred earlier
	if (cycleLastPos) {
		// if current make is also the same
		if (!(oldValues[(posCurPenalties + cycleLastPos) % maxD] - newPenalties)) {//maxD is lenght of ring buffer
			// if this has been the case for maxc periods: we have cycle
			if (++cycleL >= maxLen) {//cycleL is number of times it was already found to be equal
				//printf(    "oldValues d [%3d] (%8ld)   %6.2f'\n",cycleLastPos,iter,showtime(timep));
				//getchar();
				return true;
			}
			// otherwise: keep watching
			return false;// it is false but is finidng a cycle just not big enough yet
		}
	}
	// not detecting a cycle, or current cycle is not a cycle: reset, find last occurrence of current makespan
	cycleL = cycleLastPos = 0;
	for (unsigned u = 1; u < maxD; u++) {
		if (!(oldValues[(posCurPenalties + u) % maxD] - newPenalties))
			cycleLastPos = u; //cycleLastPos is position of the of last makespan - if cycleLastPos is 0 then it wasnt in list
	}
	return false;
}

//shift early operations making its completion time closer of its due date
void State::shiftOperations(vector<unsigned>& starts, vector<unsigned>& Q) {

	vector<unsigned> newStarts(inst.O, 0);
	vector<unsigned> startJobSuccessor(inst.J, UINT_MAX);
	vector<unsigned> startMachSuccessor(inst.M, UINT_MAX);

	unsigned qInsert = Q.size();
	unsigned curOp;
	unsigned newMakes = 0;

	while (qInsert > 0) {

		curOp = Q[--qInsert];

		if (!curOp) continue;

		unsigned sMS = startMachSuccessor[inst.operToM[curOp]];
		unsigned sJS = startJobSuccessor[inst.operToJ[curOp]];

		unsigned shiftLimit = min(sJS, sMS);

		unsigned shiftTarget = min(inst.deadlines[curOp], shiftLimit);
		shiftTarget = shiftTarget - inst.P[curOp];

		newStarts[curOp] = max(starts[curOp], shiftTarget);

		newMakes = max(newStarts[curOp] + inst.P[curOp], newMakes);

		startJobSuccessor[inst.operToJ[curOp]] = newStarts[curOp];
		startMachSuccessor[inst.operToM[curOp]] = newStarts[curOp];
	}

	assert(newMakes >= makes);
	makes = newMakes;

	starts = newStarts;
}

void State::getEarlBlock(vector<unsigned>& earlBlock) const {
	for (unsigned i = 1; i < inst.O; ++i) {
		if (startTime[i] + inst.P[i] >= inst.deadlines[i]) continue;

		unsigned auxOp = i;
		vector<unsigned> earlBlockAux;
		earlBlockAux.push_back(i);
		while (mach[auxOp] && startTime[mach[auxOp]] < startTime[job[auxOp]]) {
			earlBlockAux.push_back(mach[auxOp]);
			auxOp = mach[auxOp];
		}
		if (earlBlockAux.size() > earlBlock.size()) {
			earlBlock = earlBlockAux;
		}
	}
}

void heapify(vector<unsigned>& arr, vector<unsigned>& dists, int n, int i) {

	// Initialize largest as root
	int largest = i;

	// left index = 2*i + 1
	int l = 2 * i + 1;

	// right index = 2*i + 2
	int r = 2 * i + 2;

	// If left child is larger than root
	if (l < n && dists[arr[l]] > dists[arr[largest]])
		largest = l;

	// If right child is larger than largest so far
	if (r < n && dists[arr[r]] > dists[arr[largest]])
		largest = r;

	// If largest is not root
	if (largest != i) {
		swap(arr[i], arr[largest]);

		// Recursively heapify the affected sub-tree
		heapify(arr, dists, n, largest);
	}
}

// Main function to do heap sort
void heapSort(vector<unsigned>& arr, vector<unsigned>& dists) {
	int n = arr.size();

	// Build heap (rearrange vector)
	for (int i = n / 2 - 1; i >= 0; i--)
		heapify(arr, dists, n, i);

	// One by one extract an element from heap
	for (int i = n - 1; i > 0; i--) {

		// Move current root to end
		swap(arr[0], arr[i]);

		// Call max heapify on the reduced heap
		heapify(arr, dists, i, 0);
	}
}

void State::cplexSolve(const unsigned maxSecs) {

	//vector<unsigned> auxJobs = inst.roots;
	//vector<vector<unsigned>> machOrder(inst.M, vector<unsigned>(0));
	vector<unsigned> dists;

	dists.clear();
	dists.push_back(0);

	IloEnv jitEnv;
	IloModel jitModel(jitEnv);

	IloNumVarArray completionTimes(jitEnv, inst.O - 1, 0, IloInfinity, ILOINT);
	try {
		IloExpr objExpr(jitEnv);
		for (unsigned i = 1; i < inst.O; ++i) {
			/* Objective*/
			objExpr += (IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) * inst.earlPenalties[i]) + (IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) * inst.tardPenalties[i]);

			/* Job precedence constraints*/
			if (inst.next[i].size()) {
				jitModel.add(completionTimes[i - 1] <= completionTimes[inst.next[i][0] - 1] - inst.P[inst.next[i][0]]);
			}

			/* Penalties constraints*/
			/*jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= inst.deadlines[i] - completionTimes[i - 1]);
			jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= completionTimes[i - 1] - inst.deadlines[i]);
			jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= 0);
			jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= 0);*/
		}

		for (unsigned i = 0; i < inst.J; ++i) {
			jitModel.add(completionTimes[inst.roots[i] - 1] - inst.P[inst.roots[i]] >= 0);
		}

		for (unsigned i = 0; i < inst.M; ++i) {
			for (unsigned j = 0; j < inst.J; ++j) {
				for (unsigned k = j + 1; k < inst.J; ++k) {
					jitModel.add(completionTimes[inst.machOpers[i][j] - 1] <= completionTimes[inst.machOpers[i][k] - 1] - inst.P[inst.machOpers[i][k]] || completionTimes[inst.machOpers[i][k] - 1] <= completionTimes[inst.machOpers[i][j] - 1] - inst.P[inst.machOpers[i][j]]);
				}
			}
		}

		jitModel.add(IloMinimize(jitEnv, objExpr));

		IloCplex jitCplex(jitEnv);
		jitCplex.setParam(IloCplex::Param::TimeLimit, maxSecs);
		jitCplex.setParam(IloCplex::Param::WorkMem, 8192);
		jitCplex.setParam(IloCplex::Param::MIP::Strategy::File, 3);
		jitCplex.setOut(jitEnv.getNullStream());

		jitCplex.extract(jitModel);

		jitCplex.solve();


		for (unsigned i = 1; i < inst.O; ++i) {
			dists.push_back(round(jitCplex.getValue(completionTimes[i - 1])) - inst.P[i]);
		}

		penalties = jitCplex.getObjValue();
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}

	double penal = 0;

	cout << endl;

	for (unsigned i = 1; i < dists.size(); ++i) {
		penal += (int)inst.deadlines[i] - (int)(dists[i] + inst.P[i]) >= 0 ? ((int)inst.deadlines[i] - (int)(dists[i] + inst.P[i])) * inst.earlPenalties[i] : ((int)(dists[i] + inst.P[i]) - (int)inst.deadlines[i]) * inst.tardPenalties[i];
	}

	cout << penal << endl <<  endl;

	jitEnv.end();

	cout << penalties << endl << endl;

	for (unsigned i = 0; i < inst.O; ++i) {
		cout << dists[i] << " ";
	}
	cout << endl << endl;
		
	vector<unsigned> opersAux(inst.O - 1, 0);

	for (unsigned i = 1; i <= inst.O; ++i) {
		opersAux[i - 1] = i;
	}

	heapSort(opersAux, dists);

	vector<vector<unsigned>> machOrderAux(inst.M, vector<unsigned>());

	for (unsigned i = 0; i < inst.O-1; ++i) {
		machOrderAux[inst.operToM[opersAux[i]]].push_back(opersAux[i]);
	}

	for (unsigned i = 0; i < inst.M; ++i) {
		for (unsigned j = 0; j < inst.J - 1; ++j) {
			mach[machOrderAux[i][j]] = machOrderAux[i][j + 1];
		}
	}

	for (unsigned i = 0; i < inst.O; ++i) {
		cout << mach[i] << " ";
	}

	cout << endl;
}

		//@return: starts
vector<unsigned> State::genSchedule(unsigned schedulerType, bool cplex) {

	vector<unsigned> starts(inst.O, 0);
#ifndef NDEBUG
	bool cycle;
#endif

	if (cplex || schedulerType == 3) {
#ifndef NDEBUG
		cycle =
#endif
			scheduleCplex(starts);
	}
	else if (schedulerType == 2) {
#ifndef NDEBUG
		cycle =
#endif
			scheduleDelaying(starts);
	}
	else {
#ifndef NDEBUG
		cycle =
#endif
			scheduleAsEarly(starts);
	}

	assert(!cycle);

	return starts;
}

bool State::verifySchedule(unsigned schedulerType, bool cplex) {

	unsigned testMakes = makes;//genSchedule may change makes - store here to be sure
	vector<unsigned> schedule = genSchedule(schedulerType, cplex);

	return inst.verifySchedule(schedule, testMakes);
}

void State::printPenalties(unsigned schedulerType) {
	vector<unsigned> schedule = genSchedule(schedulerType, true);
#ifdef PRINT_DEFAULT
	cout << makes << " " << millisecsFound << " ";
#endif //PRINT_ONLY_RESULT

	inst.printPenalties(schedule, makes);
}

#ifndef NDEBUG
bool State::isAlloced() const {
	assert(inst.O != 0);
	if (mach.size() != inst.O) return false;
	if (job.size() != inst.O) return false;

	return true;
}

bool State::isClear() const {
	assert(inst.O != 0);
	assert(isAlloced());
	//if( ! roots.empty()) return false;
	for (unsigned o : job)
		if (o != 0) return false;

	for (unsigned o : _job)
		if (o != 0) return false;

	for (unsigned o : mach)
		if (o != 0) return false;

	for (unsigned o : _mach)
		if (o != 0) return false;

	return true;
}

//not all opers in state yet
bool State::partialVerify(const vector<bool>& inOpers) const {
	assert(isAlloced());
	assert(!inOpers[0]);
	vector<unsigned> jobLeafs(inst.J, 0);
	vector<unsigned> machLeafs(inst.M, 0);
	vector<unsigned> jobRoots(inst.J, 0);
	vector<unsigned> machRoots(inst.M, 0);
	unsigned j;
	unsigned m;

	if (job[0] != 0 || _job[0] != 0) {
		cout << "\n\t\tERROR !\n" << toString() << endl << "job[0] != 0    ||   _job[0] !=0" << endl;
		return false;
	}
	if (mach[0] != 0 || _mach[0] != 0) {
		cout << "\n\t\tERROR !\n" << toString() << endl << "mach[0] != 0    ||   _mach[0] !=0" << endl;
		return false;
	}

	for (unsigned o = 1; o < inst.O; o++) {
		j = inst.operToJ[o];
		m = inst.operToM[o];
#ifdef VANILLA_PSS
		//not in must have no next or prev
		if (!inOpers[o]) {
			if (job[o] != 0 || _job[o] != 0) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "job[o] != 0   ||    _job[o] != 0" << " ; o:" << o << endl;
				return false;
			}

			if (mach[o] != 0 || _mach[o] != 0) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "mach[o] != 0   ||    _mach[o] != 0" << " ; o:" << o << endl;
				return false;
			}
		}
#endif //#ifdef VANILLA_PSS
		//one leaf per job
		if (inOpers[o] && job[o] == 0) {
			jobLeafs[j]++;
			if (jobLeafs[j] > 1) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "jobLeafs[inst.operToJ[o]] > 1" << " ; o:" << o << " ; j:" << j << endl;
				return false;
			}
		}
		//one root per job
		if (inOpers[o] && _job[o] == 0) {
			jobRoots[j]++;
			if (jobRoots[j] > 1) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "jobRoots[inst.operToJ[o]] > 1" << " ; o:" << o << " ; j:" << j << endl;
				return false;
			}
		}
		//one leaf per mach
		if (inOpers[o] && mach[o] == 0) {
			machLeafs[m]++;
			if (machLeafs[m] > 1) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "machLeafs[m] > 1" << " ; o:" << o << " ; m:" << m << endl;
				return false;
			}
		}
		//one root per mach
		if (inOpers[o] && _mach[o] == 0) {
			machRoots[m]++;
			if (machRoots[m] > 1) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "machRoots[m] > 1" << " ; o:" << o << " ; m:" << m << endl;
				return false;
			}
		}

		//opers in range
		if (job[o] >= inst.O) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "job[o] >= inst.O" << " ; o:" << o << endl;
			return false;
		}
		if (_job[o] >= inst.O) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "_job[o]  >= inst.O" << " ; o:" << o << endl;
			return false;
		}
		if (mach[o] >= inst.O) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "mach[o] >= inst.O" << " ; o:" << o << endl;
			return false;
		}
		if (_mach[o] >= inst.O) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "_mach[o] >= inst.O" << " ; o:" << o << endl;
			return false;
		}

		//same job
		if (job[o] != 0 && inst.operToJ[o] != inst.operToJ[job[o]]) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "job[o]!=0   &&   inst.operToJ[o] != inst.operToJ[job[o]]" << " ; o:" << o << endl;
			return false;
		}
		if (_job[o] != 0 && inst.operToJ[o] != inst.operToJ[_job[o]]) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "_job[o]!=0   &&   inst.operToJ[o] != inst.operToJ[_job[o]]" << " ; o:" << o << endl;
			return false;
		}

		//same mach
		if (mach[o] != 0 && inst.operToM[o] != inst.operToM[mach[o]]) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "mach[o]!=0   &&   inst.operToM[o] != inst.operToM[mach[o]]" << " ; o:" << o << endl;
			return false;
		}
		if (_mach[o] != 0 && inst.operToM[o] != inst.operToM[_mach[o]]) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "_mach[o]!=0   &&   inst.operToM[o] != inst.operToM[_mach[o]]" << " ; o:" << o << endl;
			return false;
		}

		//check precs
		if (job[o] != 0 && inst.hasPrec(job[o], o)) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "job[o]!=0   &&    inst.hasPrec(job[o], o)" << " ; o:" << o << endl;
			return false;
		}
		if (_job[o] != 0 && inst.hasPrec(o, _job[o])) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "_job[o]!=0   &&    inst.hasPrec(o, _job[o])" << " ; o:" << o << endl;
			return false;
		}

		//invert compatible
		if (_job[o] != 0 && job[_job[o]] != o) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "_job[o]!=0   &&   ! job[_job[o]] == o" << " ; o:" << o << endl;
			return false;
		}
		if (job[o] != 0 && _job[job[o]] != o) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "job[o]!=0   &&   ! _job[job[o]] == o" << " ; o:" << o << endl;
			return false;
		}
		if (_mach[o] != 0 && mach[_mach[o]] != o) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "_mach[o]!=0   &&   ! mach[_mach[o]] == o" << " ; o:" << o << endl;
			return false;
		}
		if (mach[o] != 0 && _mach[mach[o]] != o) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "mach[o]!=0   &&   ! _mach[mach[o]] == o" << " ; o:" << o << endl;
			return false;
		}

		for (unsigned provost = o + 1; provost < inst.O; provost++) {
			//no repeat of opers
			if (job[o] != 0 && job[o] == job[provost]) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "job[o]!=0   &&   job[o] == job[provost]" << " ; o:" << o << " ; provost:" << provost << endl;
				return false;
			}
			if (_job[o] != 0 && _job[o] == _job[provost]) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "_job[o]!=0   &&   _job[o] == _job[provost]" << " ; o:" << o << " ; provost:" << provost << endl;
				return false;
			}
			if (mach[o] != 0 && mach[o] == mach[provost]) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "mach[o]!=0   &&   mach[o] == mach[provost]" << " ; o:" << o << " ; provost:" << provost << endl;
				return false;
			}
			if (_mach[o] != 0 && _mach[o] == _mach[provost]) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "_mach[o]!=0   &&   _mach[o] == _mach[provost]" << " ; o:" << o << " ; provost:" << provost << endl;
				return false;
			}
		}
	}

	return true;
}

bool State::verify() const {
	vector<bool> inOpers(inst.O, true);
	inOpers[0] = false;
	return partialVerify(inOpers);
}
#endif //NDEBUG