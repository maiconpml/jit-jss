#include "Instance.hpp"

Poset::Poset(const vector<pair<unsigned, unsigned>>& precs, unsigned _nEl) : nEl(_nEl), children(_nEl), parents(_nEl), cg(_nEl) {
	for (unsigned u = 0; u < nEl; u++) {
		leafs.push_back(u);
		roots.push_back(u);
	}
	for (const pair<unsigned, unsigned>& p : precs) {
		addPrec(p.first, p.second);
	}
}

Poset::Poset(unsigned _nEl) : nEl(_nEl), children(_nEl), parents(_nEl), cg(_nEl) {
	for (unsigned u = 0; u < nEl; u++) {
		leafs.push_back(u);
		roots.push_back(u);
	}
}

Poset::~Poset() {  }

bool Poset::operator!=(const Poset& p) { return !((*this) == p); }

bool Poset::operator==(const Poset& p) {

	if (nEl != p.nEl)
		return false;

	assert(children.size() == nEl);
	assert(children.size() == p.children.size());
	for (unsigned a = 0; a < nEl; a++) {
		if (children[a].size() != p.children[a].size())
			return false;
		for (unsigned b = 0; b < children[a].size(); b++) {
			if (children[a][b] != p.children[a][b])
				return false;
		}
	}

	return true;
}

//@return: created cycle doing so?
bool Poset::addAll(const vector<pair<unsigned, unsigned>>& precs) {
	bool anyCycle = false;
	for (const pair<unsigned, unsigned>& p : precs) {
		anyCycle |= addPrec(p.first, p.second);
	}
	return anyCycle;
}

//@return: created cycle doing so?
bool Poset::addPrec(unsigned from, unsigned to) {
	assert(!cg.hasCycle());

	unsigned u = 0;
	if (!cg.path(from, to)) {

		cg.add(from, to);

		children[from].push_back(to);
		parents[to].push_back(from);

		if (children[from].size() == 1) {
			for (u = 0; u < leafs.size(); u++) {
				if (leafs[u] == from)
					break;
			}
			assert(leafs[u] == from);//or else u==leafs.size() because from is not in leafs
			leafs.erase(leafs.begin() + u);
			assert(!leafs.empty());
		}
		if (parents[to].size() == 1) {
			for (u = 0; u < roots.size(); u++) {
				if (roots[u] == to)
					break;
			}
			assert(roots[u] == to);//or else u==leafs.size() because from is not in leafs
			roots.erase(roots.begin() + u);
			assert(!roots.empty());
		}
	}
	return cg.hasCycle();
}

string Poset::toString() const {
	string str = (cg.hasCycle() ? "\t\tCycle" : "\t\tNo Cycle");
	str.append("   ;   roots{ ");
	for (unsigned r : roots)
		str.append(unsigStr(r) + " ");
	str.append("}");

	str.append(" leafs{ ");
	for (unsigned l : leafs)
		str.append(unsigStr(l) + " ");
	str.append("}");

	str.append("\n\tForward: ");
	for (unsigned o = 0; o < children.size(); o++) {
		str.append(unsigStr(o) + "->{ ");
		for (unsigned c : children[o]) {
			str.append(unsigStr(c) + " ");
		}
		str.append("} ");
	}

	str.append(" \n\tbackward: ");
	for (unsigned o = 0; o < parents.size(); o++) {
		str.append("{ ");
		for (unsigned p : parents[o]) {
			str.append(unsigStr(p) + " ");
		}
		str.append("}->" + unsigStr(o) + " ");
	}

	return str;
}

double Poset::orderStr() const {
	double os = 0.0;

	for (unsigned from = 0; from < nEl; from++) {
		for (unsigned to = 0; to < nEl; to++) {
			if (from != to) {
				if (cg.path(from, to)) {
					os++;
				}
			}
		}
	}

	unsigned comb = (nEl * (nEl - 1)) / 2;

	os /= (double)comb;

	return os;
}

