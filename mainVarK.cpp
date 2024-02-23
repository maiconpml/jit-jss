/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace boost;

#include <chrono>

#include "Settings.hpp"
#include "Utilities.hpp"
#include "VarKTabu.hpp"


int main(int argc, char *argv[]) {
	try{
	      //fullTabu(const string & instPath, Params & params, double maxMillisecs)

		// (1) process commandline options
		po::options_description desc("Options");
		desc.add_options()
			("help", "show help")
			
			("ins", "instance")
			("maxMillisecs", po::value<unsigned>()->default_value(300000), "Max millisecs to finish.")

			//tabu params
			("tenure", po::value<unsigned>()->default_value(8), "Tabu tenure")
			("initialjumpLimit", po::value<unsigned>()->default_value(2500), "Initial max iterations withour improvement before backjump.")
			("jumpLimitDecrease", po::value<unsigned>()->default_value(400), "Decrease in jumLimit each time jumps without improve.")
			("bjSize", po::value<unsigned>()->default_value(5), "Back jump list size. List of elite solutions.")
			("maxD", po::value<unsigned>()->default_value(100), "cycle size to stop search and force backjump.")
			("maxC", po::value<unsigned>()->default_value(2), "number of repeats of cycle to force backjump.")

			//learn params
			("learn", po::value<bool>()->default_value(false), "Use learning.")
			("learnMinIters", po::value<unsigned>()->default_value(30), "Minimum number of iterations learinng before triggering the trimming.")
			("learnMinKW", po::value<double>()->default_value(0.7), "Minimum Kendall-W before triggering the trimming.")
			("learnTrimm", po::value<double>()->default_value(0.2), "Size of the investigated neighbouhood.")
			;

		po::positional_options_description pod;
		pod.add("ins", 1); //instance is positional as well
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
		po::notify(vm);

		if (vm.count("help") || !vm.count("ins")) {
			cout << desc << endl;
			return 0;
		}

		// (2) initialize params
		Params params;
		const string instPath = vm["ins"].as<string>();
		const unsigned maxMillisecs = vm["maxMillisecs"].as<unsigned>();
		//tabu params
		params.tenure = vm["tenure"].as<unsigned>();
		params.initialjumpLimit = vm["initialjumpLimit"].as<unsigned>();
		params.jumpLimitDecrease = vm["jumpLimitDecrease"].as<unsigned>();
		params.bjSize = vm["bjSize"].as<unsigned>();
		params.maxD = vm["maxD"].as<unsigned>();
		params.maxC = vm["maxC"].as<unsigned>();
		//learn params
	      params.learn = vm["learn"].as<bool>();
		params.learnMinIters = vm["learnMinIters"].as<unsigned>();
		params.learnMinKW = vm["learnMinKW"].as<double>();
		params.learnTrimm = vm["learnTrimm"].as<double>();
		
#ifndef NDEBUG
		cout << "DEBUG MODE\n";
		cout << "instPath:" << instPath << "  -  maxMillisecs: " << maxMillisecs << endl;
		cout << params.toString() << endl;
#endif //NDEBUG

		// (3) execute
		Tabu::fullTabu(instPath, params, maxMillisecs);
		
	}catch(const string & e) {
		cout << e << endl;
	}
	
	return 0;
}
