/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace boost;

#include <iostream>

#include "Settings.hpp"
#include "State.hpp"
#include "Tabu.hpp"
#include "Parameters.hpp"

#ifdef PRINT_NEIGHBOURS_NB
vector<double> neigh;
#endif // PRINT_NEIGHBOURS_NB

vector<string> resultList;
Inst inst;

int main(int argc, char *argv[]) {
	
	try{
		// (1) process commandline options
		po::options_description desc("Options");
		desc.add_options()
			("help", "show help")
			("instPath", "instance")
			("name","name for log")
			("maxSecs", po::value<double>()->default_value(300.0), "maximum time in seconds")
			//("maxMilli", po::value<unsigned>()->default_value(300000), "maximum time in milliseconds")
			("seed", po::value<unsigned>()->default_value(13), "random number generator seed")
			("tenure", po::value<unsigned>()->default_value(8), "tabu tenure")
			("initialjumpLimit", po::value<unsigned>()->default_value(2500), "initial max iterations withour improvement before backjump")
			//("jumpLimitDecrease", po::value<unsigned>()->default_value(400), "decrease in jumLimit each time jumps without improve")
			("decreaseDivisor", po::value<unsigned>()->default_value(7), "decrease in jumLimit each time jumps without improve")
			("bjSize", po::value<unsigned>()->default_value(20), "back jump list size")
			("maxD", po::value<unsigned>()->default_value(100), "cycle size to stop search and force backjump")
			("maxC", po::value<unsigned>()->default_value(2), "number of repeats of cycle to force backjump")
			//("acceptAlpha", po::value<double>()->default_value(0.1466), "For Metropolis acceptance criterion")
			("timeLog", po::value<bool>()->default_value(true), "Show time of each new best sollution found")
			("scaleTime", po::value<double>()->default_value(1.0), "Change received time, for testing. (JSP)")
			("onlyMakesLowerBound", po::value<bool>()->default_value(false), "To get lower bound values")
			("schedulerType", po::value<unsigned>()->default_value(1), "Type of scheduler used: 1 - early as possible (fastest), 2 - delaying when possible (better schedule but slower), 3 - 1 and 2 mixed (1 when tard penalties are dominating and 2 when earl penalties are dominating, 4 - cplex (slowest but optimal)")
			("useSwapAllNIter", po::value<unsigned>()->default_value(1), "Number of iterations without improvement to use neighborhood Swap All")
			;
		po::positional_options_description pod;
		pod.add("instPath", 1); //instance is positional as well
		pod.add("name", 1); //name is positional as well
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
		po::notify(vm);

		// (2) check necessaary options
		if (vm.count("help") || !vm.count("instPath") || !vm.count("name")) {
			cout << "necessary fields not present   --   ";
			cout << desc << endl;
			return 0;
		}

		Parameters param;

		// (3) initialize params
		param.instPath = vm["instPath"].as<string>();
		param.name = vm["name"].as<string>();
		unsigned maxSecs = vm["maxSecs"].as<double>();
		param.seed = vm["seed"].as<unsigned>();
		param.tenure = vm["tenure"].as<unsigned>();
		param.initialjumpLimit = vm["initialjumpLimit"].as<unsigned>();
		param.decreaseDivisor = vm["decreaseDivisor"].as<unsigned>();
		param.bjSize = vm["bjSize"].as<unsigned>();
		param.maxD = vm["maxD"].as<unsigned>();
		param.maxC = vm["maxC"].as<unsigned>();
		//const double acceptAlpha = vm["acceptAlpha"].as<double>();
		param.timeLog = vm["timeLog"].as<bool>();
		param.scaleTime = vm["scaleTime"].as<double>();
		param.onlyMakesLowerBound = vm["onlyMakesLowerBound"].as<bool>();
		param.useSwapAllNIter = vm["useSwapAllNIter"].as<unsigned>();
		param.schedulerType = vm["schedulerType"].as<unsigned>();
		if (param.schedulerType < 1 || param.schedulerType > 4) {
			cout << "Invalid type for schedulerType" << endl;
			return 0;
		}

		if(param.onlyMakesLowerBound) {
			inst.parse(param.instPath);

			cout << param.name << " " << inst.lowerBoundTkz() << endl;
			return 0;
		}

		// (3.1) adjusting time
		maxSecs *= param.scaleTime;
		param.maxMillisecs = maxSecs * 1000;
		
		//(4) vesrbose if in DEBUG
#ifndef NDEBUG
		cout << "DEBUG MODE" << endl;
		cout << "instPath:" << param.instPath << endl;
		cout << "maxSecs: " << maxSecs << "   seed: " << param.seed << "   tenure: " << param.tenure << "   initialjumpLimit: " << param.initialjumpLimit  << "   decreaseDivisor: " << param.decreaseDivisor << "   bjSize: " << param.bjSize << "   maxD: " << param.maxD << "   maxC: " << param.maxC << "   timeLog: " << (param.timeLog ? "Y" : "n")  << "   scaleTime: " << param.scaleTime << endl;
#endif

		//(5) execute
		Tabu::tabu(param);
		
	}catch(const string & e) {
		cout << e << endl;
	}

	return 0;
}
