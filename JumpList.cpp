#include "Tabu.hpp"

JumpList::JumpList(unsigned _maxL) : empty(true), basePos(0), topPos(0), maxL(_maxL) { ttList.resize(_maxL); }

JumpList::~JumpList() {  }

string JumpList::toString() const {
	string str = "";

	str.append(empty ? "empty" : "NOT empty");
	str.append(" -- basePos: " + unsigStr(basePos));
	str.append(" -- topPos: " + unsigStr(topPos));
	str.append(" -- maxL: " + unsigStr(maxL));
	//str.append(" -- TT ttList: "+unsigStr(ttList.size()));
	str.append(" -- TT ttList:");
	for (const TabuTrio& tt : ttList) {
		str.append("\n\t" + tt.toString());
	}

	return str;
}

void JumpList::push(const TabuTrio& tt) {
	//cout << "\tPUSHED\n" << toString() << endl;

	assert(maxL > 0);
	if (empty) {
		empty = false;
		topPos = 0;
		basePos = 0;
	}
	else {
		topPos++;
		topPos = topPos % maxL;
		if (topPos == basePos) {
			basePos++;
			basePos = basePos % maxL;
		}
	}

	ttList[topPos] = tt;
}

void JumpList::updateCands(const vector<pair<unsigned, unsigned>>& cands) {
	assert(!empty);

#ifndef NDEBUG
	bool found;
	for (const pair<unsigned, unsigned>& p1 : cands) {
		found = false;
		for (const pair<unsigned, unsigned>& p2 : ttList[topPos].cands) {
			if (p1.first == p2.first && p1.second == p2.second)
				found = true;
		}
		assert(found);
	}
#endif
	ttList[topPos].cands = cands;

	while (ttList[topPos].cands.empty()) {
		if (topPos == basePos) {
			empty = true;
			return;
		}

		if (topPos == 0)
			topPos = maxL - 1;
		else
			topPos--;
	}
}

TabuTrio JumpList::pop() {
	if (empty)
		return TabuTrio();

	while (ttList[topPos].cands.empty()) {
		if (topPos == basePos) {
			empty = true;
			return TabuTrio();
		}

		if (topPos == 0)
			topPos = maxL - 1;
		else
			topPos--;
	}

	return ttList[topPos];
}