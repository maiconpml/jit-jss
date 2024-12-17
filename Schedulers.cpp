#include "State.hpp"

#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

// schedule operations as early as possible based on current sequence, return true if current sequence has a cycle
bool State::scheduleAsEarly(vector<unsigned>& starts) {
	assert(isAlloced());
	assert(starts.size() == inst.O);

	vector<unsigned> indeg(inst.O);
	vector<unsigned> Q(inst.O);

	unsigned qInsert = 0;
	unsigned qAccess = 0;

	unsigned curOp;
	unsigned newOp;

	unsigned newMax;

	makes = 0;

	fill(starts.begin(), starts.end(), 0);
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

	//assert(qInsert>0);

	while (qAccess < qInsert) {
		assert(qAccess < Q.size());
		curOp = Q[qAccess++];

		assert(indeg[curOp] == 0);

		newMax = starts[curOp] + inst.P[curOp];
		if (makes < newMax) {
			makes = newMax;
		}

		//from JOB
		newOp = job[curOp];
		if (newOp != 0) {
			assert(indeg[newOp] > 0);
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = newOp;
			}
			if (starts[newOp] < newMax) {
				starts[newOp] = newMax;
			}
		}

		//from MACH
		newOp = mach[curOp];
		if (newOp != 0) {
			assert(indeg[newOp] > 0);
			indeg[newOp]--;
			if (indeg[newOp] == 0) {
				assert(qInsert < Q.size() - 1);
				Q[qInsert++] = newOp;
			}
			if (starts[newOp] < newMax) {
				starts[newOp] = newMax;
			}
		}
	}

	if (qAccess == inst.O - 1) {

#ifdef SHIFT_OPERS
#ifndef NDEBUG
		inst.calcPenalties(starts, ePenalty, lPenalty, operPenalties);

		double testPenalties = ePenalty + lPenalty;
#endif //NDEBUG

		shiftOperations(starts, Q);
#endif //SHIFT_OPERS

		inst.calcPenalties(starts, ePenalty, lPenalty, operPenalties);
		penalties = ePenalty + lPenalty;
		startTime = starts;
#ifdef SHIFT_OPERS
		assert(penalties <= testPenalties);
#endif //SHIFT_OPERS
	}

	return qAccess < inst.O - 1;
}

// schedule operations delaying when is possible to decrease penalties based on current sequence, return true if current sequence has a cycle
bool State::scheduleDelaying(vector<unsigned>& starts) {
	unsigned delayTime = 0;
	vector<unsigned> limited(inst.O, 0);
	vector<unsigned> heads;
	vector<unsigned> ops;
	vector<unsigned> lateCands(inst.O, 0);
	vector<unsigned> indeg(inst.O);
	vector<unsigned> Q(inst.O);
	double pushStrength, holdStrength;
	fill(starts.begin(), starts.end(), 0);
	if (topoWalk(indeg, Q)) return true;
	ops = Q;
	reverse(ops.begin(), ops.end());
	for (unsigned i = 1; i < ops.size(); i++) {
		lateCands[ops[i]] = 1; // 1- late 2- early/on schedule
		heads.push_back(ops[i]); // first earlies operations 
		if ((job[ops[i]] && (starts[job[ops[i]]] < inst.P[ops[i]] || starts[job[ops[i]]] == 0)) || (mach[ops[i]] && (starts[mach[ops[i]]] < inst.P[ops[i]] || starts[mach[ops[i]]] == 0))) {
			forcedDelay(starts, lateCands, ops[i]);
		}
		updateStrength(lateCands, pushStrength, holdStrength);
		while (pushStrength > holdStrength) {
			delayTime = calcDelayTime(starts, lateCands, limited);
			delay(starts, lateCands, delayTime);
			update(lateCands, limited, heads, starts);
			updateStrength(lateCands, pushStrength, holdStrength);
		}

		fill(lateCands.begin(), lateCands.end(), 0);
		fill(limited.begin(), limited.end(), 0);
		heads.clear();
	}
	//verifySchedule(false);
	for (unsigned times = 0; times < 0; ++times) {
		for (unsigned i = 1; i < ops.size(); i++) {
			if (starts[ops[i]] + inst.P[ops[i]] >= inst.deadlines[ops[i]]) {
				lateCands[ops[i]] = 2;
			}
			else {
				lateCands[ops[i]] = 1; // 1- late 2- early/on schedule
			}
			heads.push_back(ops[i]); // first earlies operations 
			if ((job[ops[i]] && starts[job[ops[i]]] < inst.P[ops[i]]) || (mach[ops[i]] && starts[mach[ops[i]]] < inst.P[ops[i]])) {
				forcedDelay(starts, lateCands, ops[i]);
			}
			limited[ops[i]] = 1;
			update(lateCands, limited, heads, starts);
			updateStrength(lateCands, pushStrength, holdStrength);
			while (pushStrength > holdStrength) {
				delayTime = calcDelayTime(starts, lateCands, limited);
				delay(starts, lateCands, delayTime);
				update(lateCands, limited, heads, starts);
				updateStrength(lateCands, pushStrength, holdStrength);
			}

			fill(lateCands.begin(), lateCands.end(), 0);
			fill(limited.begin(), limited.end(), 0);
			heads.clear();
		}
	}
	makes = 0;
	for (unsigned i = 1; i < inst.O; ++i) {
		if (starts[i] + inst.P[i] > makes) makes = starts[i] + inst.P[i];
	}

	inst.calcPenalties(starts, ePenalty, lPenalty, operPenalties);
	penalties = ePenalty + lPenalty;
	startTime = starts;

	return false;
}

// schedule with cplex solver based on current sequence, return true if current sequence has a cycle
bool State::scheduleCplex(vector<unsigned>& starts) {

	vector<unsigned> auxJobs = inst.roots;
	vector<vector<unsigned>> machOrder(inst.M, vector<unsigned>(0));

	starts.clear();
	starts.push_back(0);

	vector<unsigned> indeg(inst.O);
	vector<unsigned> Q(inst.O);

	if (topoWalk(indeg, Q)) return true;

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

			/* Machine precedence constraints*/
			if (mach[i]) {
				jitModel.add(completionTimes[i - 1] <= completionTimes[mach[i] - 1] - inst.P[mach[i]]);
			}

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

		jitModel.add(IloMinimize(jitEnv, objExpr));

		IloCplex jitCplex(jitEnv);
		//jitCplex.setParam(IloCplex::Param::TimeLimit, 0.5);
		jitCplex.setOut(jitEnv.getNullStream());

		jitCplex.extract(jitModel);

		jitCplex.solve();



		for (unsigned i = 1; i < inst.O; ++i) {
			starts.push_back(round(jitCplex.getValue(completionTimes[i - 1])) - inst.P[i]);
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

	/* Adjust new makes*/
	makes = 0;
	for (unsigned i = 1; i < inst.O; ++i) {
		if (makes < starts[i] + inst.P[i]) {
			makes = starts[i] + inst.P[i];
		}
	}

	/* Calculate penalties*/
	inst.calcPenalties(starts, ePenalty, lPenalty, operPenalties);
	penalties = ePenalty + lPenalty;
	startTime = starts;

	return false;
}