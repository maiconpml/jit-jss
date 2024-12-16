#include "State.hpp"

#include <set>

void State::gifflerThompson() {

	set<unsigned> ready;
	vector<unsigned> ready0;
	vector<unsigned> ready1;
	vector<unsigned> jobStartTimes(inst.J, 0);
	vector<unsigned> machStartTimes(inst.M, 0);
	vector<unsigned> jobLeafs(inst.J, 0);
	vector<unsigned> jobLeafs2(inst.J, 0);
	vector<unsigned> machLeafs(inst.M, 0);
	vector<unsigned> machLeafs2(inst.M, 0);

	unsigned completionTime;
	unsigned earlCompletion;
	unsigned mach = 0;
	unsigned auxStartTime;

	for (unsigned op : inst.roots) {
		ready.insert(op);
	}

	while (!ready.empty()) {

		earlCompletion = UINT_MAX;

		ready0.clear();
		ready1.clear();

		for (unsigned op : ready) {
			assert(op < inst.O);
			completionTime = max(jobStartTimes[inst.operToJ[op]], machStartTimes[inst.operToM[op]]) + inst.P[op];
			if (completionTime < earlCompletion) {
				earlCompletion = completionTime;
				mach = inst.operToM[op];
			}
		}

		for (unsigned op : ready) {
			if (inst.operToM[op] == mach) {
				ready0.push_back(op);
			}
		}

		for (unsigned op : ready0) {
			auxStartTime = max(jobStartTimes[inst.operToJ[op]], machStartTimes[inst.operToM[op]]);
			if (auxStartTime < earlCompletion) {
				ready1.push_back(op);
			}
		}

		assert(ready1.size() > 0);
		unsigned op = ready1[0];

		for (unsigned i = 1; i < ready1.size(); ++i) {
			if (inst.deadlines[op] > inst.deadlines[ready1[i]]) {
				op = ready1[i];
			}
		}

		if (jobLeafs[inst.operToJ[op]]) {
			insertJobOper(jobLeafs[inst.operToJ[op]], jobLeafs2[inst.operToJ[op]], op);
			jobLeafs2[inst.operToJ[op]] = jobLeafs[inst.operToJ[op]];
		}

		jobLeafs[inst.operToJ[op]] = op;

		if (machLeafs[inst.operToM[op]]) {
			insertMachOper(machLeafs[inst.operToM[op]], machLeafs2[inst.operToM[op]], op);
			machLeafs2[inst.operToM[op]] = machLeafs[inst.operToM[op]];
		}

		machLeafs[inst.operToM[op]] = op;

		ready.erase(op);

		if (inst.next[op].size()) {
			ready.insert(inst.next[op][0]);
		}

		auxStartTime = max(jobStartTimes[inst.operToJ[op]], machStartTimes[inst.operToM[op]]);
		jobStartTimes[inst.operToJ[op]] = auxStartTime + inst.P[op];
		machStartTimes[inst.operToM[op]] = auxStartTime + inst.P[op];
	}
}