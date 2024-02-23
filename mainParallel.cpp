/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace boost;

#include <chrono>

#include "Settings.hpp"
#include "State.hpp"
#include "Tabu.hpp"
#include "Schedule.hpp"


int main(int argc, char *argv[]) {

//#ifndef PARALLEL_MACHS
//	cout << "PARELE MAIN WITHOUT PARALLEL OPTION TURNED ON - FAIL" << endl;
//#endif
//#ifdef PARALLEL_MACHS
	try{

// (1) process commandline options
		po::options_description desc("Options");
		desc.add_options()
			("help", "show help")
			("ins", "instance")
			("tenure", po::value<unsigned>()->default_value(8), "tabu tenure")
			("iters", po::value<unsigned>()->default_value(2500), "initial max iterations withour improvement before backjump")
			("iterDecrease", po::value<unsigned>()->default_value(400), "decrease in jumLimit each time jumps without improve")
			("bjSize", po::value<unsigned>()->default_value(5), "back jump list size")
			("maxD", po::value<unsigned>()->default_value(100), "cycle size to stop search and force backjump")
			("maxC", po::value<unsigned>()->default_value(2), "number of repeats of cycle to force backjump")
			("maxMillisecs", po::value<double>()->default_value(DBL_MAX), "max millisecs to finish")
			("critVanillaNeigh", po::value<bool>()->default_value(false), "type of neighbouhood: 1-critVanillaNeigh")
			("allDelayNeigh", po::value<bool>()->default_value(true), "type of neighbouhood: 2-allDelayNeigh")
			("dispatch", po::value<bool>()->default_value(false), "use dispatch rule: longest tails first")
			("startType", po::value<string>()->default_value("jobTail"), "insa or jobTail")
			("learn", po::value<bool>()->default_value(false), "learn which opers improve and cut neighbouhood after reaching maxLearn kendall w.")
			("minLearnPeriod", po::value<unsigned>()->default_value(20), "min iterations before trimming")
			("minLearnKW", po::value<double>()->default_value(0.6), "Kendall w to start trimming - expected from 0.5 to 0.75")
			("trimmSize", po::value<double>()->default_value(0.5), "percent of the neighbourhood we evaluate 0.0 to 1.0 - iff too low it will just revert to normal neighbouhood")
			;

		po::positional_options_description pod;
		pod.add("ins", 1); //instance is positional as well
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
		po::notify(vm);

		// (2) check options
		if (vm.count("help") || !vm.count("ins")) {
			cout << desc << endl;
			return 0;
		}

		// (3) initialize params
		const string instPath = vm["ins"].as<string>();
		const unsigned tenure = vm["tenure"].as<unsigned>();
		const unsigned iters = vm["iters"].as<unsigned>();
		const unsigned iterDecrease = vm["iterDecrease"].as<unsigned>();
		const unsigned bjSize = vm["bjSize"].as<unsigned>();
		const unsigned maxD = vm["maxD"].as<unsigned>();
		const unsigned maxC = vm["maxC"].as<unsigned>();
		const double maxMillisecs = vm["maxMillisecs"].as<double>();
		const bool critVanillaNeigh = vm["critVanillaNeigh"].as<bool>();
		const bool allDelayNeigh = vm["allDelayNeigh"].as<bool>();
		const bool dispatch = vm["dispatch"].as<bool>();
		const string startTypeStr = vm["startType"].as<string>();
		unsigned startType;
		if(startTypeStr == "insa")
			startType = PARALLEL_INSA_START;
		else if(startTypeStr == "jobTail")
			startType = PARALLEL_JOB_TAIL_START;
		else
			throw errorText("no start type selected","","");
		
		if(!critVanillaNeigh   &&   !allDelayNeigh) {
			throw errorText("no neighbouhoood selected","","");
		}
		bool learn = vm["learn"].as<bool>();
		unsigned minLearnPeriod = vm["minLearnPeriod"].as<unsigned>();
		double minLearnKW = vm["minLearnKW"].as<double>();
		double trimmSize = vm["trimmSize"].as<double>();
		
		//cout << startType << endl;

#ifndef NDEBUG
		cout << "DEBUG MODE\n";
		cout << "instPath:" << instPath << endl;
		cout << "tenure:" << tenure << " iters:" << iters << " iterDecrease:" << iterDecrease << " bjSize:" << bjSize << " maxD:" << maxD << " maxC:" << maxC << endl;
#ifdef ANALYZE_MODE
		cout << "learn: " << learn << "  minLearnPeriod: " << minLearnPeriod <<   "  minLearnKW: " << minLearnKW << "  trimmSize: " << trimmSize << endl;
#endif //ANALYZE_MODE
#endif //NDEBUG

		// (4) execute code
	    //Tabu::nosmuTabu(instPath, tenure, iters, iterDecrease, bjSize, maxD, maxC, startType, maxMillisecs);
	    if(dispatch) {
	    	inst.parse(instPath);
			Schedule sch;
			//XXX
			//vector<unsigned> prio(inst.O);
			//for(unsigned u=0; u<inst.O; u++)
			//	prio[u] = u;
			//XXX
			sch.dispatchLongestJobTail(SMALLER_GAP);
			//s.dipatchByPiority(prio, SMALLER_GAP);
			//cout << sch.toString() << endl;
			//cout << sch.makes << endl;
		} else {
#ifdef ANALYZE_MODE
			Tabu::nosmuTabuParallel(instPath, tenure, iters, iterDecrease, bjSize, maxD, maxC, maxMillisecs, critVanillaNeigh, allDelayNeigh, startType, learn, minLearnPeriod, minLearnKW, trimmSize);
#else
			Tabu::nosmuTabuParallel(instPath, tenure, iters, iterDecrease, bjSize, maxD, maxC, maxMillisecs, critVanillaNeigh, allDelayNeigh, startType);
#endif //ANALYZE_MODE
		}
		//cout << inst.toString() << endl;
		/*
		const string instPath = "./parallelInst/ParallelJspLawrence01_15/la01_k2.ppsi";
		//const string instPath = "./parallelInst/la01_k1.ppsi";
		unsigned tenure = 8;
		unsigned initialjumpLimit = 2500;
		unsigned jumpLimitDecrease = 400;
		unsigned bjSize = 5;
		unsigned maxD = 100;
		unsigned maxC = 2;
		unsigned maxMillisecs = 5000;

		Tabu::nosmuTabuParallel(instPath, tenure, initialjumpLimit, jumpLimitDecrease, bjSize, maxD, maxC, maxMillisecs);
		*/
		
		/*
		const string instPath = "./parallelInst/ParallelJspLawrence01_15/la01_k2.ppsi";
		//const string instPath = "./parallelInst/la01_k1.ppsi";
		inst.parse(instPath);
		//cout << inst.toString() << endl;
		
		State s;
		s.alloc();

		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		vector<unsigned> dists(inst.O);
		multi_array<unsigned, 2> buckets;
		buckets.resize(boost::extents[inst.M][inst.kReps]);
		s.insaParallel(indeg, Q, dists, buckets);
		
		unsigned lastOp;
		vector<unsigned> prev(inst.O);
		s.setMetaParallel(dists, lastOp, prev, indeg, Q, buckets);

		cout << s.toString() << endl;
		*/

		/*		
		State s;
		s.alloc();
		unsigned targOp = 96;
		vector<unsigned> heads(inst.O);
		vector<unsigned> indeg(inst.O);
		vector<unsigned> Q(inst.O);
		multi_array<unsigned, 2> buckets;
		buckets.resize(boost::extents[inst.M][inst.kReps]);
		unsigned aDist = s.maxDistToTargParallel(targOp, heads, indeg, Q, buckets);
		cout << "forthDist: " << aDist << endl;
		aDist = s.maxDistFromTargParallel(targOp, heads, indeg, Q, buckets);
		cout << "backDist: " << aDist << endl;
		cout << "cost: " << inst.P[targOp] << endl;
		*/

		/*
		// (1) process commandline options
		po::options_description desc("Options");
		desc.add_options()
			("help", "show help")
			("instPath", "instance")
			("maxSecs", po::value<unsigned>()->default_value(180), "maximum time in seconds")
			("seed", po::value<unsigned>()->default_value(13), "random number generator seed")
			("tenure", po::value<unsigned>()->default_value(8), "tabu tenure")
			("initialjumpLimit", po::value<unsigned>()->default_value(2500), "initial max iterations withour improvement before backjump")
			("jumpLimitDecrease", po::value<unsigned>()->default_value(400), "decrease in jumLimit each time jumps without improve")
			("bjSize", po::value<unsigned>()->default_value(5), "back jump list size")
			("maxD", po::value<unsigned>()->default_value(100), "cycle size to stop search and force backjump")
			("maxC", po::value<unsigned>()->default_value(2), "number of repeats of cycle to force backjump")
			("acceptAlpha", po::value<double>()->default_value(0.1466), "For Metropolis acceptance criterion")
			("perturbSize", po::value<unsigned>()->default_value(10), "Perturbation size")
			("bubbleP", po::value<double>()->default_value(0.98), "Probablity of insertion in bubble perturb")
			("startType", po::value<string>()->default_value("INSA_START"), "type of start: INSA_START or RAND_START")
			("intensificType", po::value<string>()->default_value("TABU_INTENS"), "type of start: TABU_INTENS or LOCAL_SEARCH_INTENS")
			("perturbType", po::value<string>()->default_value("INSA_PERTURB"), "Perturbation type: CARLIER_PERTURB, BUBBLE_PERTURB, INSA_PERTURB")
			("acceptType", po::value<string>()->default_value("METROPOLIS"), "Jum type: METROPOLIS, BEST_JUMP")
			;
		po::positional_options_description pod;
		pod.add("instPath", 1); //instance is positional as well
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
		po::notify(vm);
		// (2) check options
		if (vm.count("help") || !vm.count("instPath")) {
			cout << desc << endl;
			return 0;
		}
		// (3) initialize params
		string instPath = vm["instPath"].as<string>();
		unsigned maxSecs = vm["maxSecs"].as<unsigned>();
		unsigned seed = vm["seed"].as<unsigned>();
		unsigned tenure = vm["tenure"].as<unsigned>();
		unsigned initialjumpLimit = vm["initialjumpLimit"].as<unsigned>();
		unsigned jumpLimitDecrease = vm["jumpLimitDecrease"].as<unsigned>();
		unsigned bjSize = vm["bjSize"].as<unsigned>();
		unsigned maxD = vm["maxD"].as<unsigned>();
		unsigned maxC = vm["maxC"].as<unsigned>();
		double acceptAlpha = vm["acceptAlpha"].as<double>();
		unsigned perturbSize = vm["perturbSize"].as<unsigned>();
		double bubbleP = vm["bubbleP"].as<double>();
		string startTypeStr = vm["startType"].as<string>();
		unsigned startType;
		if(startTypeStr == "RAND_START") {
			startType = RAND_START;
		}  else if(startTypeStr == "INSA_START") {
			startType = INSA_START;
		} else {
			throw errorText("startTypeStr: "+startTypeStr ,"mainIG.cpp","main()");
		}
		string intensificTypeStr = vm["intensificType"].as<string>();
		unsigned intensificType;
		if(intensificTypeStr == "LOCAL_SEARCH_INTENS") {
			intensificType = LOCAL_SEARCH_INTENS;
		}  else if(intensificTypeStr == "TABU_INTENS") {
			intensificType = TABU_INTENS;
		} else {
			throw errorText("intensificTypeStr: "+intensificTypeStr ,"mainIG.cpp","main()");
		}
		string perturbTypeStr = vm["perturbType"].as<string>();
		unsigned perturbType;
		if(perturbTypeStr == "INSA_PERTURB") {
			perturbType = INSA_PERTURB;
		} else {
			throw errorText("perturbTypeStr: "+perturbTypeStr ,"mainIG.cpp","main()");
		}

		string acceptTypeStr = vm["acceptType"].as<string>();
		unsigned acceptType;
		if(acceptTypeStr == "METROPOLIS") {
			acceptType = METROPOLIS;
		} else if(acceptTypeStr == "BEST_JUMP") {
			acceptType = BEST_JUMP;
		} else {
			throw errorText("acceptTypeStr: "+acceptTypeStr ,"mainIG.cpp","main()");
		}


#ifndef NDEBUG
		cout << "DEBUG MODE" << endl;
		cout << "instPath:" << instPath << endl;
		cout << "maxSecs: " << maxSecs << "   seed: " << seed << "   tenure: " << tenure << "   initialjumpLimit: " << initialjumpLimit  << "   jumpLimitDecrease: " << jumpLimitDecrease << "   bjSize: " << bjSize << "   maxD: " << maxD << "   maxC: " << maxC << "   acceptAlpha: " << acceptAlpha << "   perturbSize: " << perturbSize << endl;
		cout << "startTypeStr: " << startTypeStr << "   intensificTypeStr: " << intensificTypeStr << "   perturbTypeStr: " << perturbTypeStr << "   acceptTypeStr: " << acceptTypeStr << endl;
#endif

		IG::igGeneric(instPath, maxSecs, seed, tenure, initialjumpLimit, jumpLimitDecrease,bjSize, maxD, maxC, acceptAlpha, perturbSize, bubbleP,startType, intensificType, perturbType, acceptType);
		*/

	}catch(const string & e) {
		cout << e << endl;
	}		
//#endif //PARALLEL_MACHS
	return 0;
}


