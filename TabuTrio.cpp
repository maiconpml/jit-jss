#include "Tabu.hpp"

TabuTrio::TabuTrio() : dummy(true) {  }

TabuTrio::TabuTrio(const State& _state, const vector<pair<unsigned, unsigned>>& _cands, const TabuList& _tabuList) : dummy(false), state(_state), cands(_cands), tabuList(_tabuList) {  }

TabuTrio::~TabuTrio() {  }

string TabuTrio::toString() const {
	string str;

	if (dummy)
		return "DUMMY TRIO";

	str.append("Makes: " + unsigStr(state.makes) + "   ---   ");
	str.append("Cands size: " + unsigStr(cands.size()));

	return str;
}