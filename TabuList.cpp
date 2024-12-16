#include "Tabu.hpp"

TabuList::TabuList() {  }

TabuList::TabuList(int _tenure) : curIter(0), tenure(_tenure) {
	if (tenure == 0)
		return;
	tabu.resize(boost::extents[inst.O][inst.O]);
	fill(tabu.data(), tabu.data() + tabu.num_elements(), -1);
}

TabuList::TabuList(TabuList&& other) {
	this->swap(other);
}

TabuList::TabuList(const TabuList& other) {
	curIter = other.curIter;
	tenure = other.tenure;
	tabu.resize(boost::extents[other.tabu.shape()[0]][other.tabu.shape()[1]]);
	tabu = other.tabu;
}

TabuList& TabuList::operator=(TabuList other) {
	this->swap(other);
	return *this;
}

void TabuList::swap(TabuList& other) {
	using std::swap;
	swap(curIter, other.curIter);
	swap(tenure, other.tenure);
	tabu.resize(boost::extents[other.tabu.shape()[0]][other.tabu.shape()[1]]);
	swap(tabu, other.tabu);
}

string TabuList::toString() const {
	string str = "";

	str.append("curIter: " + unsigStr(curIter));
	str.append("  ;  tenure: " + unsigStr(tenure));
	str.append(" --- forbidden: ");
	str.append(listForbidden());

	return str;
}

void TabuList::passTime(int time) {
	assert(time <= tenure);
	curIter += time;
}

string TabuList::listForbidden() const {
	string str = "";
	for (unsigned o1 = 0; o1 < inst.O; o1++) {
		for (unsigned o2 = 0; o2 < inst.O; o2++) {
			if (o2 > o1) {
				assert(tabu[o1][o2] <= curIter + tenure);
				if (!canMove(o1, o2)) {
					str.append(unsigStr(o1) + "-" + unsigStr(o2) + "[" + unsigStr(age(o1, o2)) + "]    ");
					assert(!canMove(o2, o1));
				}
				else {
					assert(canMove(o2, o1));
				}
			}
		}
	}
	return str;
}

void TabuList::insert(unsigned o1, unsigned o2) {
	assert(o1 != o2);
	assert(tenure != 0);
	assert(o1 < inst.O);
	assert(o2 < inst.O);

	curIter++;
	tabu[o2][o1] = curIter + tenure;
	tabu[o1][o2] = curIter + tenure;
}

bool TabuList::canMove(unsigned o1, unsigned o2) const {
	assert(tenure > 0);
	assert(tabu[o1][o2] <= curIter + tenure);
	assert(tabu[o1][o2] == tabu[o2][o1]);
	return tabu[o1][o2] <= curIter;
}

int TabuList::age(unsigned o1, unsigned o2) const {
	assert(tenure + curIter >= tabu[o1][o2]);
	assert(tabu[o1][o2] == tabu[o2][o1]);
	return tenure + curIter - tabu[o1][o2];
}

int TabuList::timeToLeave(unsigned o1, unsigned o2) const {
	int delta = tabu[o1][o2] - curIter;
	assert(delta > 0);
	assert(delta <= tenure);
	assert(tenure - age(o1, o2) == delta);
	return delta;
}