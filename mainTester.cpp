/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

//#define NDEBUG

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace boost;

#include "Settings.hpp"
#include "Instance.hpp"
#include "State.hpp"
#include "Tabu.hpp"
//#include "Order.hpp"
#include "Iterated.hpp"
#include "Schedule.hpp"
#include "Analyzer.hpp"

#include <boost/multi_array.hpp>
#include <set>

using namespace boost;

int main(int argc, char *argv[]) {

			
	try {
		inst.parse("./dummyInstJobDom2nd.ppsi");
		//inst.parse("./GholamiInst/la20ParallelGholami.ppsi");
		//inst.parse("./dummyInstJobDom.ppsi");
		cout << inst.toString() << endl;
		
		vector<vector<unsigned>> prioActual(inst.O, vector<unsigned>(inst.O, UINFTY));
		vector<vector<unsigned>> prioPredict(inst.O, vector<unsigned>(inst.O, UINFTY));
		vector<vector<unsigned>> R(inst.O, vector<unsigned>(inst.O));
		vector<OperPairQuality> qualityV;
		qualityV.push_back(OperPairQuality(4,3,0));
		qualityV.push_back(OperPairQuality(7,12,0));
		qualityV.push_back(OperPairQuality(8,1,0));
		
		prioActual[4][3] = 1;
		prioActual[7][12] = 2;
		prioActual[8][1] = 3;

		
		prioPredict[4][3] = 3;
		prioPredict[7][12] = 1;
		prioPredict[8][1] = 2;

		
		cout << Analyzer::pairKendallW(prioActual, prioPredict, R, qualityV) << endl;
		
		/*Analyzer an;
		an.alloc();
		
		cout << an.toString() << endl;*/
		
		//Schedule sch;
		//sch.alloc();
		//cout << sch.toString() << endl;
		
		//unsigned tiebreaker = SMALLER_GAP;
		//unsigned curStart;
		
		
		/*
		set<pair<unsigned,unsigned>> s;
		s.insert(pair<unsigned,unsigned>(1,2));
		s.insert(pair<unsigned,unsigned>(13,2));
		s.insert(pair<unsigned,unsigned>(1,2));
		for(const pair<unsigned, unsigned> & p : s)cout << p.first << " " << p.second << endl;
		*/
		
		
		
		//curStart = s.machines[6].insertOper(2, 3, tiebreaker);
		//cout << "curStart: " << curStart << endl;
		//cout << s.toString() << endl;
		
		//curStart = s.machines[6].insertOper(1, 3, tiebreaker);
		//cout << "curStart: " << curStart << endl;
		//cout << s.toString() << endl;
		

		//sch.dispatchLongestJobTail(tiebreaker);
		//cout << sch.toString() << endl;

		
		//State st;
		//st.alloc();
		
		//s.dispatchLongestJobTail(SMALLER_GAP);
		
		//vector<vector<unsigned>> allMachOrder(inst.M);
		//sch.getMachOrders(allMachOrder);
		/*
		for(unsigned m=0; m<inst.M; m++) {
			cout << "m:  "; 
			for(unsigned op : allMachOrder[m]) {
				cout << " " << op;
			}
			cout << endl;
		}*/
		
		/*
		st.importFromMachOrder(allMachOrder);
		
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		unsigned lastOp;
		vector<unsigned> prev(inst.O);
		vector<unsigned> heads(inst.O); 
		//vector<unsigned> prev(inst.O);
		multi_array<unsigned, 2> bucketLimits;
		bucketLimits.resize(boost::extents[inst.M][inst.maxK]);
		multi_array<unsigned, 2> bucketTops;
		bucketTops.resize(boost::extents[inst.M][inst.maxK]);
		
		//(vector<unsigned> & dists, unsigned & lastOp, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & bucketLimits, multi_array<unsigned, 2> & bucketTops)
		st.makes = st.setMetaParallel(heads, lastOp, prev, indeg, Q,  bucketLimits,  bucketTops);
		cout << st.toString() << endl;
		
		cout << "schedule: " << sch.makes << endl;
		cout << "state: " << st.makes << endl;
		*/
		/*
		s.machines[0].buckets[0].occupy(Window(10,13));
		
		cout << s.machines[0].buckets[0].toString() << endl;
		
		s.machines[0].buckets[0].occupy(Window(0,3));
		
		cout << s.machines[0].buckets[0].toString() << endl;
		
		s.machines[0].buckets[0].occupy(Window(5,6));
		
		cout << s.machines[0].buckets[0].toString() << endl;
		
		cout << s.machines[0].buckets[0].bestStart(0,2) << endl;
		//s.machines[0].buckets[0].occupy(Window(6,10));
		
		//cout << s.machines[0].buckets[0].toString() << endl;
		*/
		
		/*
		inst.parse("./dummyParallelGholami.ppsi");
		cout << inst.toString() << endl;
		
		State s;
		s.alloc();
		
		vector<unsigned> order = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
		
		s.dispatchInOrder(order);
		*/
		
		//cout << s.toString() << endl;
/*
		inst.parse("./insts/GSP/Nasiri/ActualStageTai/actStage06.psi");
		cout << inst.toString() << endl;
		
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
*/
		/*
		cout << "lower bound nasiri: " << inst.lowerBoundNasiri(indeg, Q) << endl;
		
		unsigned sum;

		cout << "machine sums:" << endl;
		for(const vector<unsigned> & aMach : inst.machOpers) {
			sum = 0;
			for(unsigned o : aMach)
				sum+=inst.P[o];
			cout << "\t" << sum << endl;
		}

		cout << "job sums:" << endl;
		for(const vector<unsigned> & aJob : inst.jobOpers) {
			sum = 0;
			for(unsigned o : aJob)
				sum+=inst.P[o];
			cout << "\t" << sum << endl;
		}
		*/
//		cout << inst.lowerBoundAllBefore() << endl;
			

		/*
		// (1) process commandline options
		po::options_description desc("Options");
		desc.add_options()
			("help", "show help")
			("instPath", "instance")
			;
		po::positional_options_description pod;
		pod.add("instPath", 1); //instance is positional as well
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
		po::notify(vm);

		// (2) check necessaary options
		if (vm.count("help") || !vm.count("instPath")) {
			cout << desc << endl;
			return 0;
		}

		// (3) initialize params
		const string instPath = vm["instPath"].as<string>();

		inst.parse(instPath);
		cout << inst.toString() << endl;

		bool bigJobFirst = false;
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		vector<unsigned> heads(inst.O); 
		vector<unsigned> tails(inst.O);
		vector<unsigned> invertedQ(inst.O);
		vector<bool> reach(inst.O);
		//unsigned lastOp;
		vector<unsigned> prev(inst.O);
		
		State s;
		s.alloc();
		s.insaPsp(bigJobFirst, indeg, Q, heads, tails, invertedQ, reach);
		unsigned tenure = 8;
		unsigned initialjumpLimit = 1000;
		unsigned decreaseDivisor = 8;
		unsigned bjSize = 5;
		const high_resolution_clock::time_point tpStart = high_resolution_clock::now();
		vector<unsigned> dists(inst.O);
		//vector<unsigned> prev(inst.O);
		//vector<unsigned> indeg;
		//vector<unsigned> Q
		const unsigned maxD = 50;
		const unsigned maxC = 2;
		const unsigned lowerBound = 0;
		unsigned maxMillisecs = 10000;
		unsigned printWhenBetter = UINT_MAX;
		const string preString = "xxx";
		//vector<unsigned> & heads, vector<unsigned> & tails,
		bool timeLog = true;
		Tabu::evolveTabuN6(s, tenure, initialjumpLimit, decreaseDivisor, bjSize, tpStart,  dists, prev, indeg, Q,  maxD, maxC, lowerBound, maxMillisecs, printWhenBetter, preString, heads, tails, timeLog);
		//Tabu::evolveTabu(s, tenure, initialjumpLimit, decreaseDivisor, bjSize, tpStart,  dists, prev, indeg, Q,  maxD, maxC, lowerBound, maxMillisecs, printWhenBetter, preString, heads, tails, timeLog, true);
		cout << s.toString() << endl;
		*/
		
		/*
		bool cycle = s.setMeta(heads, lastOp, prev, indeg, Q);
		if(cycle)
			throw errorText("CYCLE!!","","");
		//cout << s.toString() << endl;
		cout << s.easyString() << endl;

		vector<unsigned> critic;
		critic.reserve(inst.O);
		State::computeCritic(critic, lastOp, prev);

		cout << "Critic:";
		for(unsigned o : critic) cout << "   " <<  "<o:" << o << ",j:" << inst.operToJ[o] << ",m:"<<inst.operToM[o] << ">";
		cout << endl;
		
		vector<tuple<unsigned, unsigned, unsigned, unsigned>> cands;
		cands.reserve((inst.O)^2);
		vector<unsigned> jobBb;
		jobBb.reserve(inst.O);
		vector<unsigned> machBb;
		machBb.reserve(inst.O);

		s.compTailsComplete(tails, indeg, Q);
		s.compHeadsComplete(heads, indeg, Q);
		
		//s.fillCandidatesN6(cands, jobBb, machBb, critic, heads, tails);


		cout << "cands:";
		for(tuple<unsigned, unsigned, unsigned, unsigned> & t : cands) cout << (get<1>(t)==FORTH ? " - F," : " - B,") << get<2>(t) << "," << get<3>(t);
		cout << endl;
		*/
		/*
		unsigned tenure = 8;
		unsigned initialjumpLimit = 2500;
		unsigned jumpLimitDecrease = 400;
		unsigned bjSize = 5;
		unsigned maxD = 100;
		unsigned maxC = 2;
		unsigned startType = INSA_START;
		unsigned maxMillisecs = 30000;
		bool timeLog = false;

		Tabu::nosmuTabu(instPath, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, maxD, maxC, startType, maxMillisecs, timeLog);
		*/
		// (5) execute
		//cout << instPath << endl;

		/*
		bool bigJobFirst = false;
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		vector<unsigned> heads(inst.O); 
		vector<unsigned> tails(inst.O);
		vector<unsigned> invertedQ(inst.O);
		vector<bool> reach(inst.O);
		
		State s;
		s.alloc();
		s.insaPsp(bigJobFirst, indeg, Q, heads, tails, invertedQ, reach);
		cout << s.toString() << endl;
		cout << s.easyString() << endl;

		unsigned o1 = 3;
		unsigned o2 = 7;
		unsigned direction = BACK;
		cout << "Moving " << o1 << (direction==FORTH ? " FORTH " : " BACK ") << "to " << o2 << endl;
		s.shift(o1,o2,direction);
		
		cout << s.toString() << endl;
		cout << s.easyString() << endl;
		*/

		/*
		inst.parse("./insts/all/admuStage05.psi");
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		cout << inst.lowerBoundTkz(indeg, Q) << endl;
		*/
		
		/*
		//const string instPath = "insts/Dummy/newDummy.psi";
		const string instPath = "insts/all/ft10_1.psi";
		//const string instPath = "insts/all/stageNas80_1.psi";
		//const string instPath = "insts/all/pspTa01_0.psi";
		inst.parse(instPath);

		double p = 0.5;
		unsigned perturbSize = 4;
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		vector<unsigned> operPos(inst.O);
		vector<bool> reach(inst.O);
		vector<unsigned> dists(inst.O);
		bool bigJobFirst = false;
		vector<vector<unsigned>> enforce(inst.O);
		for(vector<unsigned> & v : enforce)
			v.reserve(max(inst.J, inst.M));
		vector<unsigned> enfIndeg(inst.O);
		unsigned lastOp;
		vector<unsigned> prev(inst.O);
		vector<unsigned> heads(inst.O); 
		vector<unsigned> tails(inst.O);
		vector<unsigned> invertedQ(inst.O);

		State s;
		s.alloc();
		s.insaPsp(bigJobFirst, indeg, Q, heads, tails, invertedQ, reach);
		bool cycle = s.setMeta(dists, lastOp, prev, indeg, Q);
		if(cycle)
			throw errorText("CYCLE!!","","");
		//cout << s.toString() << endl;
		cout << s.makes << endl;
		//cout << s.easyString() << endl;

		if( ! s.verifySchedule())
			cout << "\n\t!! invalid schedule\n"; 
		cout << "BUBBLE" << endl;
		for(unsigned seed=0; seed<100; seed++) {
			fprintf(stderr, "\n\tseed: %d\n", seed);
			//unsigned seed = 0;
			RandomNumberGenerator<>::instance().seed(seed);
			srand(seed);
			s.bubblePerturb(p, perturbSize, enforce, enfIndeg, indeg, Q, operPos, reach);
		    cycle = s.setMeta(dists, lastOp, prev, indeg, Q);
			if(cycle)
				throw errorText("CYCLE!!","","");
			if( ! s.verifySchedule())
				cout << "\n\t!! invalid schedule\n"; 
			cout << "\tmakes: " << s.makes << endl;
		}
		*/
		

		/*
		const string instPath = "./insts/JSP/abz5Jsp.psi";
		//const string instPath = "./insts/PSSP/PspNasiriStyleAllTaiJsp10copies/pspTa01_0.psi";
		//const string instPath = "./insts/OSP/Tai/ta01Osp.psi";
		//const string instPath = "./insts/MSP/MixedLiu/mixed_base16_no5_rep0.psi";
		inst.parse(instPath);
		//cout << inst.toString() << endl;
		ElementProb ep;
		ep.setMap();

		double sumProb = 0.0;

		cout << "J: " << inst.J << "   M: " << inst.M << endl;

		for(unsigned j=0; j<inst.J; j++) {
			cout << "\t" << j << "-" << (1.0 - (inst.posets[j].orderStr())) << endl;
			sumProb += (1.0 - (inst.posets[j].orderStr()));
		}

		vector<unsigned> jobChoices(inst.J, 0);
		vector<unsigned> machChoices(inst.M, 0);

		int type;
		unsigned index;
		const unsigned multi = 1000000;
		const unsigned maxTests = (multi*sumProb)  +   inst.M*multi;
		cout << "maxTests: " << maxTests << endl;
		for(unsigned t=0; t<maxTests; t++) {
			tie(type, index) = ep.randElement();
			if(type==JOB)
				jobChoices[index]++;
			else
				machChoices[index]++;
		}
		
		cout << "Selected Jobs:" << endl;
		for(unsigned j=0; j<inst.J; j++)
			cout << "\t" << j << ":" << jobChoices[j] << endl;

		cout << "Selected Machs:" << endl;
		for(unsigned m=0; m<inst.M; m++)
			cout << "\t" << m << ":" << machChoices[m] << endl;
		*/
		/*
#ifdef PARALLEL_MACHS
		
		const string instPath = "./parallelInst/ParallelJspLawrence01_15/la01_k2.ppsi";
		inst.parse(instPath);

		State s;
		s.alloc();

		vector<unsigned> dists(inst.O); 
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		multi_array<unsigned, 2> buckets;
		buckets.resize(boost::extents[inst.M][inst.K]);
		vector<unsigned> machRoot(inst.O);
		bool cycle;

	    s.insaParallel(indeg, Q, dists, buckets);

		//cout << s.toString() << endl;

		unsigned lastOp;
		vector<unsigned> prev(inst.O);
		s.setMetaParallel(dists, lastOp, prev, indeg, Q, buckets);

		cout << "Start LS: " << s.makes << endl;

		vector<pair<unsigned, unsigned>> cands;
		cands.reserve(inst.O);
		vector<unsigned> critic;
		critic.reserve(inst.O);
		s.localSearchParallel(dists, prev, indeg, Q, cands, critic, buckets);

		cout << s.easyString() << endl;
		s.shiftForth(86, 96);
		cout << s.easyString() << endl;
		*/

		/*
		cout << "Finish LS: " << s.makes << endl;

		double p = 0.95;
		unsigned perturbSize = 1;
		vector<vector<unsigned>> enforce(inst.O);
		vector<unsigned> enfIndeg(inst.O);
		vector<unsigned> operPos(inst.O);
		vector<bool> reach(inst.O);
		for(unsigned rep=0; rep<10; rep++) {
			s.bubbleParallel(p, perturbSize, enforce, enfIndeg, indeg, Q, operPos, reach);
			cycle = s.setMetaParallel(dists, lastOp, prev, indeg, Q, buckets);
			if(cycle) s.makes = UINT_MAX;
			cout << "\tBubble: " << s.makes << endl;
			s.localSearchParallel(dists, prev, indeg, Q, cands, critic, buckets);
			cycle = s.setMetaParallel(dists, lastOp, prev, indeg, Q, buckets);
			if(cycle) s.makes = UINT_MAX;
			cout << "\tLS: " << s.makes << endl;
		}

		//inst.verifyPrallelSchedule(dists);

#endif //PARALLEL_MACHS
		*/


		/*
		const string instPath = "./parallelInst/dummy.ppsi";
		inst.parse(instPath);

		//cout << inst.toString() << endl;

		vector<pair<unsigned, unsigned>> v = {pair<unsigned, unsigned>(0, UINT_MAX-1)};
		vector<vector<pair<unsigned, unsigned>>> machGaps(2, v);
		unsigned jobHead = 0;
		unsigned oper = 4;

		cout << "inserting cost: " << inst.P[oper] << "  ;  head: " << jobHead << endl; 
		State::insertInGap(machGaps, jobHead, oper);
		cout << State::printBucket(machGaps); 

		jobHead = 10;
		oper = 9;
		cout << "inserting cost: " << inst.P[oper] << "  ;  head: " << jobHead << endl; 
		State::insertInGap(machGaps, jobHead, oper);
		cout << State::printBucket(machGaps); 

		jobHead = 5;
		oper = 4;
		cout << "inserting cost: " << inst.P[oper] << "  ;  head: " << jobHead << endl; 
		State::insertInGap(machGaps, jobHead, oper);
		cout << State::printBucket(machGaps); 
		*/

		/*
		vector<unsigned> v = {1, 2, 3, 4, 5};
		vector<unsigned>::iterator i = v.begin();
		while((*i) != 5) i++;
		v.insert(i, 777);
		for(unsigned u : v) cout << u << " ";
		cout << endl;

		i = v.begin();
		while((*i) != 777) i++;
		v.erase(i);
		for(unsigned u : v) cout << u << " ";
		cout << endl;
		*/


		/*
		vector<pair<unsigned, unsigned>> gapV = { pair<unsigned, unsigned>(1,8),
											   pair<unsigned, unsigned>(3,3),
											   pair<unsigned, unsigned>(4,5),
											   pair<unsigned, unsigned>(6,4),
											   pair<unsigned, unsigned>(7,UINT_MAX-1)
		};

		vector<pair<unsigned, unsigned>>operV;
		for(unsigned u=0; u<=8; u++)
			operV.push_back(pair<unsigned, unsigned>(u, UINT_MAX));




		cout << "gaps:";
		for(const pair<unsigned, unsigned> & gap : gapV)
			cout << "  ;  S:" << gap.first << "d:" << gap.second;
		cout << endl;

		vector<pair<unsigned, unsigned>>::iterator low;
		cout << "Positions:\n";
		for(const pair<unsigned, unsigned> & oper : operV) {
			low=upper_bound (gapV.begin(), gapV.end(), oper);
			if(low != gapV.begin()) low--;
			if(low == gapV.end())
				cout << "\toper: " << oper.first << "\tEND !" << endl;
			else
				cout << "\toper: " << oper.first << " in gap:   S:" << low->first << "d:" << low->second << endl;
		}
		*/

		/*
		//const string instPath = "insts/Dummy/newDummy.psi";
		const string instPath = "insts/all/ft10_1.psi";
		//const string instPath = "insts/all/stageNas25_1.psi";

		inst.parse(instPath);
		cout << inst.toString() << endl;

		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		vector<unsigned> operPos(inst.O);
		vector<bool> reach(inst.O);
		vector<unsigned> dists(inst.O);
		bool bigJobFirst = false;
		State s;
		s.alloc();
		s.insaPsp(bigJobFirst, dists, indeg, Q);
		cout << s.easyString() << endl;

		vector<vector<unsigned>> enforce(inst.O);
		for(vector<unsigned> & v : enforce)
			v.reserve(max(inst.J, inst.M));
		vector<unsigned> enfIndeg(inst.O);
		unsigned type = MACH;
		unsigned index = 7;
		s.mustEnforce(enforce, enfIndeg, type, index,  indeg, Q, operPos, reach);
		cout << "\nenforce\n";
		for(unsigned fromOp : (type==JOB  ?  inst.jobOpers[index] : inst.machOpers[index])) {
			cout << "\t" << fromOp << ":";
			for(unsigned toOp : enforce[fromOp]) {
				cout << " " << toOp;
			}
			cout << endl;
		}
		cout << "\nenfIndeg:";
		for(unsigned o : (type==JOB  ?  inst.jobOpers[index] : inst.machOpers[index])) {
			cout << o << "(" << enfIndeg[o] << ")   ";
		}
		cout << endl;
		*/
		/*
		vector<unsigned> indeg;
		vector<unsigned> Q;
		vector<Inst> iV(4);
		for(unsigned u=1; u<=3; u++) {
			iV[u].parse("insts/all/actStage0" + unsigStr(u)+".psi");
			indeg.resize(iV[u].O);
			Q.resize(iV[u].O);
			cout << "actStage0" << u << " " << iV[u].lowerBoundTkz(indeg, Q) << endl;
		}
		*/
		/*
		unsigned i1 = argv[1];
		unsigned i2 = argv[2];
		const string instName = "pspTa" + (string)(i1<10 ? "0" : "") + unsigStr(i1)+"_"+unsigStr(i2)+".psi";
		const string instPath = "insts/all/"+instName;
	    inst.parse(instPath);
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		cout << inst.toString() << endl;
		cout << inst.lowerBoundNasiri(indeg, Q);
		*/
		/*
		vector<unsigned> indeg;
		vector<unsigned> Q;
		vector<Inst> iV(800);
		unsigned ii = 0;
		for(unsigned i1=1; i1<=80;i1++) {
			for(unsigned i2=0; i2<10; i2++) {
				iV[ii].parse("insts/all/pspTa" + (string)(i1<10 ? "0" : "") + unsigStr(i1)+"_"+unsigStr(i2)+".psi");
				indeg.resize(iV[ii].O);
				Q.resize(iV[ii].O);
				cout << "pspTa" + (string)(i1<10 ? "0" : "") + unsigStr(i1)+"_"+unsigStr(i2) << " NA NA " << iV[ii].lowerBoundNasiri(indeg, Q) << endl;
				ii++;
			}
			}*/

		/*
		//const string instPath = "insts/Dummy/newDummy.psi";
		const string instPath = "insts/all/ft10_1.psi";
		//const string instPath = "insts/all/stageNas25_1.psi";

		inst.parse(instPath);
		unsigned seed = 13;
		RandomNumberGenerator<>::instance().seed(seed);
		srand(seed);

		vector<unsigned> dists(inst.O);
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		vector<unsigned> prev(inst.O);
		vector<unsigned> jobBb;
		jobBb.reserve(inst.O);
		vector<unsigned> machBb;
		machBb.reserve(inst.O);
		vector<pair<unsigned, unsigned>> cands;
		cands.reserve(inst.O);
		vector<unsigned> critic;
		critic.reserve(inst.O);
		vector<bool> inOpers(inst.O, true);
		inOpers[0]=false;
		vector<unsigned> jobRoot(inst.J);
		vector<unsigned> jobLeaf(inst.J);
		vector<unsigned> machRoot(inst.M);
		vector<unsigned> machLeaf(inst.M);


		State theState;
		theState.alloc();

		theState.makeRand(indeg, Q);
		unsigned lowerBound = inst.lowerBoundTkz(indeg, Q);
		cout << "lowerBound: " << lowerBound << endl;
		theState.localSearch(dists, prev, indeg, Q, jobBb, machBb, cands, critic, lowerBound);
		//const double bubbleP = 0.98;
		const unsigned perturbSize = 15;
		//theState.bubblePerturb(bubbleP, perturbSize);
		for(unsigned u=0; u<10; u++) {
			theState.insaPerturb(perturbSize, dists, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf);
			theState.localSearch(dists, prev, indeg, Q, jobBb, machBb, cands, critic, lowerBound);
			cout << theState.makes << endl;
		}
		*/
		
		
		/*
		const string instPath = "insts/Dummy/newDummy.psi";
		//const string instPath = "insts/all/ft10_1.psi";
		//const string instPath = "insts/all/stageNas80_1.psi";
		unsigned maxSecs = 10;
		unsigned seed = 13;
		unsigned tenure = 8;
		unsigned initialjumpLimit = 2500;
		unsigned jumpLimitDecrease = 400;
		unsigned bjSize = 5;
		unsigned maxD = 100;
		unsigned maxC = 2;
		double acceptAlpha = 0.95;
		unsigned perturbSize = 5;
		double bubbleP = 0.95;
		unsigned startType = INSA_START;
		//unsigned startType = RAND_START;
		unsigned perturbType = INSA_PERTURB;
		//unsigned perturbType = BUBBLE_PERTURB;
		unsigned intensificType = LOCAL_SEARCH_INTENS;
		//unsigned intensificType = TABU_INTENS;
		*/

		//RandomNumberGenerator<>::instance().seed(1);
		/*double r[10]; // 10 doubles randomicos entre 0 e 1
		for (int i = 0; i < 10; ++i) {
			r[i] = random_number<double>(0.0, 1.0);
			// se for int, poderia ser por exemplo
			// r[i] = random_number<int>(0, 25000)
			cout << r[i] << " ";
			}*/
		/*
		//string instPath = "insts/Dummy/newDummy.psi";
		string instPath = "insts/all/ft10_1.psi";
		//string instPath = "insts/all/stageNas80_1.psi";
		inst.parse(instPath);
		cout << inst.toString() << endl;
		*/
		/*
		srand(13);
		//string instPath = "insts/Dummy/newDummy.psi";
		//string instPath = "insts/all/ft10_1.psi";
		string instPath = "insts/all/stageNas80_1.psi";
		inst.parse(instPath);
		//cout << inst.toString() << endl;
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		vector<unsigned> dists(inst.O);
		unsigned lastOp;
		vector<unsigned> prev(inst.O);
		State s;
		s.alloc();
		s.makeRand(indeg, Q);
		bool cycle = s.setMeta(dists, lastOp, prev, indeg, Q);
		if(cycle)
			throw errorText("CYCLE!!","","");
		//cout << s.toString() << endl;
		cout << s.makes << endl;
		*/
		/*
		//string instPath = "insts/Dummy/newDummy.psi";
		string instPath = "insts/all/ft10_1.psi";
		//string instPath = "insts/all/stageNas80_1.psi";
		//string instPath = "insts/all/pspTa08_9.psi";


		unsigned tenure = 8;
		//unsigned initialjumpLimit = 2500;
		unsigned initialjumpLimit = 40000;
		unsigned jumpLimitDecrease = 400;
		unsigned bjSize = 5;
		const unsigned maxD = 100;
		const unsigned maxC = 2;
		Tabu::nosmuTabu(instPath, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, maxD, maxC, INSA_START);
		*/
		/*
		inst.parse("insts/all/stageNas80_1.psi");
		TabuList tl(8);
		for(unsigned o1=0; o1<3; o1++) {
			for(unsigned o2=10; o2<13; o2++) {
				tl.insert(o1, o2);
				cout << tl.toString() << endl << endl;
			}
		}
		*/
		/*
		//inst.parse("insts/Dummy/newDummy.psi");
		//inst.parse("insts/all/ft10_1.psi");
		//inst.parse("insts/all/abz7_1.psi");
		inst.parse("insts/all/stageNas80_1.psi");
		//cout << inst.toString() << endl;

		unsigned lastOp;
		
		State s;
		s.alloc();

		vector<unsigned> dists(inst.O);
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);

		bool bigJobFirst = false;

		s.insaPsp(bigJobFirst, dists, indeg, Q);

		vector<unsigned> prev(inst.O);
		unsigned cycle = s.setMeta(dists, lastOp, prev, indeg, Q);
		if(cycle)
			cout << "\n\t!! cycle\n";
		if( ! s.verifySchedule())
			cout << "\n\t!! invalid schedule\n"; 
		//cout << s.toString() << endl;
		cout << s.makes << endl;
		//vector<unsigned> schedule = s.genSchedule();
		//cout << "schedule: ";
		//for(unsigned o=1; o<inst.O; o++) cout << " o:" << o << "[" << schedule[o] << "]";
		//cout << endl;
		*/
		/*
		srand(7);
		inst.parse("insts/all/ft10_1.psi");
		//inst.parse("insts/all/stageNas80_1.psi");
		//inst.parse("insts/Dummy/basicDummy.psi");
		//cout << inst.toString() << endl;

		Order order;
		order.alloc();
		Adj adj;
		adj.alloc();

		vector<unsigned> dists(inst.O);
		unsigned makes;
		vector<unsigned> prev(inst.O);
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);

		bool cycle;

		for(unsigned u=0; u<10000; u++) {
			//order.genRandom();
			assert(order.verify());
			//cout << order.toString() << endl;
			adj.convert(order);

			cycle = adj.setMeta(dists, makes, prev, indeg, Q);
			if( ! cycle)
				cout << adj.toString() << endl;

			adj.clear();
		}
		*/
		//cout << inst.lowerBoundNasiri() << endl;

		/*

#ifndef NDEBUG
		cout << "DEBUG MODE" << endl;
		//assert(false);
#endif
		const string instPath = "./insts/all/ft10_1.psi";
		unsigned tenure = 8;
		unsigned initialjumpLimit = 40000;
		//unsigned initialjumpLimit = 5000;
		//unsigned initialjumpLimit = 2500;
		unsigned jumpLimitDecrease = 400;
		unsigned bjSize = 5;
		const unsigned maxD = 100;
		const unsigned maxC = 2;
		steady_clock::time_point tpStart = steady_clock::now();
		State s = Tabu::nosmuTabu(instPath, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, maxD, maxC);
		cout << "totalSecs: " << (duration_cast<chrono::nanoseconds>(steady_clock::now() - tpStart).count())/1000000000.0 << endl; 

		cout << s.makes << endl;

		cout << "XXX_GARBAGE_makesCall: " << XXX_GARBAGE_makesCall << endl;
		 */

		/*
		theInst.parse("./insts/all/ft10_1.psi");
		State s;
		multi_array<int, 2> maxDists;
		maxDists.resize(extents[theInst.nJobs][theInst.nMachs]);
		multi_array<Oper, 2> prev;
		prev.resize(extents[theInst.nJobs][theInst.nMachs]);
		multi_array<unsigned,2> indeg;
		indeg.resize(extents[theInst.nJobs][theInst.nMachs]);
		vector<Oper> q((theInst.nJobs)*(theInst.nMachs));
		s.makeInsa();
		bool cycle = s.setMeta(maxDists, prev, indeg, q);
		if(cycle)
			cout << "\n\n\t\tcycle ??\n";
		//cout << s.toString() << endl;
		const Oper o1 = Oper(3,4);
		const Oper o2 = Oper(7, 4);
		cout << "swaping: " << o1.toString() << o2.toString() << endl;
		s.swapOpers(o1, o2, maxDists, prev);

		cout << s.toString() << endl;
		*/


		//Inst i("./insts/all/ft10_1.psi");
		//theInst.parse("./insts/Dummy/basicDummy.psi");
		//theInst.parse("./insts/all/ft10_1.psi");
		//cout << theInst.toString() << endl;

		/*
		State s;
		multi_array<int, 2> maxDists;
		maxDists.resize(extents[theInst.nJobs][theInst.nMachs]);
		multi_array<Oper, 2> prev;
		prev.resize(extents[theInst.nJobs][theInst.nMachs]);
		multi_array<unsigned,2> indeg;
		indeg.resize(extents[theInst.nJobs][theInst.nMachs]);
		vector<Oper> q((theInst.nJobs)*(theInst.nMachs));
		s.makeInsa();
		bool cycle = s.setMeta(maxDists, prev, indeg, q);
		if(cycle)
			cout << "\n\n\t\tcycle ??\n";
		//cout << s.toString() << endl;

		cout << "insa state: " << endl;
		cout << s.toString() << endl;

		vector<pair<Oper,Oper>> cand;
		s.fillCandidatesN5(cand);
		for(const pair<Oper,Oper> & p : cand)
			cout << p.first.toString() << "-" << p.second.toString() << endl;

		unsigned tabuTenure = 8;
		TabuList tabuList(tabuTenure, theInst.nJobs, theInst.nMachs);

		unsigned aspiration = 930;
		Tabu::nsp(s, tabuList, cand, aspiration, maxDists, prev, indeg, q);
		cout << s.toString() << endl;
		*/
		/*
		multi_array<int, 2> heads;
		heads.resize(extents[theInst.nJobs][theInst.nMachs]);
		theInst.setHeads(heads);
		cout << "\nheads";
		for(unsigned j=0; j<theInst.nJobs; j++) {
			cout << "\n\t";
			for(unsigned m=0; m<theInst.nMachs; m++) {
				cout << heads[j][m] << " ";
			}
		}
		
		multi_array<int, 2> tails;
		tails.resize(extents[theInst.nJobs][theInst.nMachs]);
		theInst.setTails(tails);
		cout << "\ntails";
		for(unsigned j=0; j<theInst.nJobs; j++) {
			cout << "\n\t";
			for(unsigned m=0; m<theInst.nMachs; m++) {
				cout << tails[j][m] << " ";
			}
		}

		cout << "\nlb: " << theInst.lowerBoundNasiri() << endl;
		*/
		//PartialState ps = PartialState::makeInsaState(false);
		//cout << ps.toString() << endl;


		/*
		cout << "computed dists:";
		for(unsigned j=0; j<theInst.nJobs; j++) {
			cout << "\n\t";
			for(unsigned m=0; m<theInst.nMachs; m++) {
				cout << maxDists[j][m] << " ";
			}
			}*/
		/*
		s.setDists(maxDists, indeg, q);
		cout << "\ndirect dists:";
		for(unsigned j=0; j<theInst.nJobs; j++) {
			cout << "\n\t";
			for(unsigned m=0; m<theInst.nMachs; m++) {
				cout << maxDists[j][m] << " ";
			}
		}

		s.setInvertDists(maxDists, indeg, q);
		cout << "\ninvert dists:";
		for(unsigned j=0; j<theInst.nJobs; j++) {
			cout << "\n\t";
			for(unsigned m=0; m<theInst.nMachs; m++) {
				cout << maxDists[j][m] << " ";
			}
			}*/


		/*
		bool scheduleOk = s.verifySchedule();
		if( ! scheduleOk) {
			cout << s.toString() << endl;
			cout << theInst.toString() << endl;
			throw errorText("Bad schedule !!!!! .","","");
			}*/

		/*
		for(unsigned i1=1; i1<=80;i1++) {
				cout << i1 <<  " " << endl;
			}
		*/


	}catch(const string & e) {
		cout << e << endl;
	}
	return 0;
}


/*
typedef boost::multi_array<double, 3> array_type;
 array_type myarray[2][3][4];
  typedef boost::multi_array_types::index_range range;
  array_type::index_gen indices;
  array_type::array_view<2>::type myview =    myarray[ indices[range(0,2)][1][range(0,4,2)] ];

  for (array_type::index i = 0; i != 2; ++i)
    for (array_type::index j = 0; j != 2; ++j)
	assert(myview[i][j] == myarray[i][1][j*2]);
*/

	/*
	multi_array<int, 2> m1;
	m1.resize(boost::extents[2][3]);
	//fill(m1.data(),m1.data()+m1.num_elements(),-1);
	//cout << m1.shape()[0];

	int filler = 0;
	for(unsigned i=0; i<2; i++) {
		for(unsigned j=0; j<3; j++) {
			m1[i][j] = filler++;
		}
	}

	for(unsigned i=0; i<m1.shape()[0]; i++) {
		for(unsigned j=0; j<m1.shape()[1]; j++) {
			cout << m1[i][j] << " ";
		}
		cout << endl;
	}
	*/
	/*
	typedef boost::multi_array_types::index_range range;
    multi_array<unsigned, 2>::array_view<1>::type myView();
	for(unsigned i=0; i<2; i++) {
		myView = m1[ indices[range(0,2)][i]];
	}
	*/

	//getchar();
