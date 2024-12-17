#include "State.hpp"

#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

#include <queue>

// set candidates neighbors with late operations critical path
void State::fillCandidatesCriticTotal(vector<pair<unsigned, unsigned>>& cands, vector<unsigned>& starts, vector<unsigned>& _job, vector<unsigned>& _mach, vector<unsigned>& mach) {

	for (unsigned op = 1; op < inst.O; ++op) {

		if (starts[op] + inst.P[op] > inst.deadlines[op]) {

			unsigned auxOp = op;
			vector<unsigned> opCritic;

			opCritic.push_back(op);
			while (_job[auxOp] != 0 || _mach[auxOp] != 0) {

				while (_mach[auxOp] && starts[_mach[auxOp]] + inst.P[_mach[auxOp]] > starts[_job[auxOp]] + inst.P[_job[auxOp]]) {
					opCritic.push_back(_mach[auxOp]);
					auxOp = _mach[auxOp];
				}
				if (opCritic.size() > 1) {
					cands.push_back(pair<unsigned, unsigned>(opCritic[1], opCritic[0]));
					if (opCritic.size() > 2)cands.push_back(pair<unsigned, unsigned>(opCritic[opCritic.size() - 1], opCritic[opCritic.size() - 2]));
				}
				auxOp = _job[auxOp];
				opCritic.clear();
			}
		}
	}

	unsigned op1Aux, op2Aux;
	for (unsigned i = 0; i < cands.size(); ++i) {
		op1Aux = cands[i].first;
		op2Aux = cands[i].second;

		for (unsigned j = i + 1; j < cands.size(); ++j) {
			if (cands[j].first == op1Aux && cands[j].second == op2Aux) {
				cands.erase(cands.begin() + j);
			}
		}
	}
}

// set neighbors candidates with all swaps between two adjacent operations
void State::fillCandidatesAllSwaps(vector<pair<unsigned, unsigned>>& cands, vector<unsigned>& starts, vector<unsigned>& _job, vector<unsigned>& _mach, vector<unsigned>& mach) {

	assert(mach.size() == inst.O);

	for (unsigned currentOp = 1; currentOp < inst.O; ++currentOp) {

		if (mach[currentOp]) {
			cands.push_back(pair<unsigned, unsigned>(currentOp, mach[currentOp]));
		}
	}
}

/* Schedule operations relaxing machine precedence of a block of early operations. Already selects the best neighbor*/
void State::schedulerCplexRelax(vector<pair<unsigned, unsigned>>& cands, vector<unsigned>& starts, vector<unsigned>& _job, vector<unsigned>& _mach, vector<unsigned>& mach) {

	vector<unsigned> auxJobs = inst.roots;
	vector<vector<unsigned>> machOrder(inst.M, vector<unsigned>(0));
	vector<unsigned> dists;

	dists.clear();
	dists.push_back(0);

	vector<unsigned> relaxBlock;
	getEarlBlock(relaxBlock);

	IloEnv jitEnv;
	IloModel jitModel(jitEnv);

	IloNumVarArray completionTimes(jitEnv, inst.O - 1, 0, IloInfinity, ILOINT);
	try {
		IloExpr objExpr(jitEnv);
		for (unsigned i = 1; i < inst.O; ++i) {
			/* Objective*/
			objExpr += (IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) * inst.earlPenalties[i]) + (IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) * inst.tardPenalties[i]);

			/* Start time of a job's first operation must be greater than zero*/
			if (!_job[i]) jitModel.add(completionTimes[i - 1] - inst.P[i] >= 0);

			/* Job precedence constraints*/
			if (job[i]) {
				jitModel.add(completionTimes[i - 1] <= completionTimes[job[i] - 1] - inst.P[job[i]]);
			}

			/* Penalties constraints*/
			/*jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= inst.deadlines[i] - completionTimes[i - 1]);
			jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= completionTimes[i - 1] - inst.deadlines[i]);
			jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= 0);
			jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= 0);*/
		}

		vector<unsigned> machHeads(inst.M, 0);

		/* Find machine heads*/
		for (unsigned i = 1; i < inst.O; ++i) {
			if (!_mach[i]) {
				machHeads[inst.operToM[i]] = i;
			}
		}

		/* Add or constraints between relaxed operations*/
		for (unsigned i = 0; i < relaxBlock.size(); ++i) {
			for (unsigned j = i + 1; j < relaxBlock.size(); ++j) {
				jitModel.add(completionTimes[relaxBlock[i] - 1] <= completionTimes[relaxBlock[j] - 1] - inst.P[relaxBlock[j]] || completionTimes[relaxBlock[j] - 1] <= completionTimes[relaxBlock[i] - 1] - inst.P[relaxBlock[i]]);
			}
		}

		/* For each machine, set machine precendence constraints*/
		for (unsigned i = 0; i < inst.M; ++i) {
			unsigned auxOp = machHeads[i];
			while (mach[auxOp]) {
				/* If current machine is the relaxed machine, set machine precedence constraints between non-relaxed and relaxed operations*/
				if (relaxBlock.size() && i == inst.operToM[relaxBlock[0]]) {
					/* If the next operations is a relaxed operation, set all precendence constraints between current and relaxed operations*/
					if (mach[auxOp] == relaxBlock[0]) {
						for (unsigned j = 0; j < relaxBlock.size(); ++j) {
							jitModel.add(completionTimes[auxOp - 1] <= completionTimes[relaxBlock[j] - 1] - inst.P[relaxBlock[j]]);
						}
						auxOp = mach[relaxBlock[relaxBlock.size() - 1]];
						continue;
					}
					/* If current operation is directly after relaxed operations, set all precedence constraints between relaxed and current operations*/
					else if (_mach[auxOp] == relaxBlock[relaxBlock.size() - 1]) {
						for (unsigned j = 0; j < relaxBlock.size(); ++j) {
							jitModel.add(completionTimes[relaxBlock[j] - 1] <= completionTimes[auxOp - 1] - inst.P[auxOp]);
						}
					}
				}
				jitModel.add(completionTimes[auxOp - 1] <= completionTimes[mach[auxOp] - 1] - inst.P[mach[auxOp]]);
				auxOp = mach[auxOp];
			}
		}


		jitModel.add(IloMinimize(jitEnv, objExpr));

		IloCplex jitCplex(jitEnv);
		//jitCplex.setParam(IloCplex::Param::TimeLimit, 1);
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

	jitEnv.end();

	makes = 0;
	for (unsigned i = 1; i < inst.O; ++i) {
		if (makes < dists[i] + inst.P[i]) {
			makes = dists[i] + inst.P[i];
		}
	}

	/* For each relaxed operation, verify if its order has changed*/
	for (unsigned i = 0; i < relaxBlock.size() - 1; ++i) {
		/* If the order of any relaxed operation has changed adjust order of relaxed machine*/
		if (dists[relaxBlock[i]] > dists[relaxBlock[i + 1]]) {

			priority_queue<pair<unsigned, unsigned>, vector<pair<unsigned, unsigned>>, greater<pair<unsigned, unsigned>>> q;

			for (unsigned j = 0; j < relaxBlock.size(); ++j) {
				q.push(pair<unsigned, unsigned>(dists[relaxBlock[j]], relaxBlock[j]));
			}

			unsigned prevOper = _mach[relaxBlock[0]];
			unsigned nextOper = mach[relaxBlock[relaxBlock.size() - 1]];
			unsigned curOper;
			while (!q.empty()) {
				curOper = q.top().second;
				q.pop();
				if (prevOper) mach[prevOper] = curOper;
				_mach[curOper] = prevOper;
				prevOper = curOper;
			}
			mach[curOper] = nextOper;
			if (nextOper)	_mach[nextOper] = curOper;
			break;
		}
	}

	inst.calcPenalties(dists, ePenalty, lPenalty, operPenalties);
	penalties = ePenalty + lPenalty;
	startTime = dists;
}