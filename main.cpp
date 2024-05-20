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
#include "Iterated.hpp"

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
			("increaseIters", po::value<unsigned>()->default_value(0), "Increase the max number of iterations after each perturb")
			("bjSize", po::value<unsigned>()->default_value(5), "back jump list size")
			("maxD", po::value<unsigned>()->default_value(100), "cycle size to stop search and force backjump")
			("maxC", po::value<unsigned>()->default_value(2), "number of repeats of cycle to force backjump")
			("acceptAlpha", po::value<double>()->default_value(0.1466), "For Metropolis acceptance criterion")
			//("perturbSize", po::value<unsigned>()->default_value(20), "Perturbation size")
			("sizePBubble", po::value<unsigned>()->default_value(4), "Perturbation size for bubble")
			("sizePInsa", po::value<unsigned>()->default_value(20), "Perturbation size for insa")
			("timeLog", po::value<bool>()->default_value(true), "Show time of each new best sollution found")
			("lowerBoundOrder",  po::value<bool>()->default_value(true), "Order neighbourhood by lower bound")
			("perturbType", po::value<string>()->default_value("insa"), "Either bubble, insa, mix")
			("scaleTime", po::value<double>()->default_value(1.0), "Change received time, for testing. (JSP)")
			("scheduleFile",  po::value<string>()->default_value(""), "File to write best schedule. Empty string means no wirting.")
			("bubbleP", po::value<double>()->default_value(0.5), "Bubble perturb probability of insertion.")
			("mixProb", po::value<double>()->default_value(0.05), "Chance to bubble instead of insa.")
			("searchType", po::value<string>()->default_value("iteratedGreedy"), "Type of search: iteratedGreedy - tabuSearch")
			("onlyLowerBound", po::value<bool>()->default_value(false), "To get lower bound values")
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

		// (3) initialize params
		const string instPath = vm["instPath"].as<string>();
		string name = vm["name"].as<string>();
		double maxSecs = vm["maxSecs"].as<double>();
		const unsigned seed = vm["seed"].as<unsigned>();
		const unsigned tenure = vm["tenure"].as<unsigned>();
		const unsigned initialjumpLimit = vm["initialjumpLimit"].as<unsigned>();
		const unsigned decreaseDivisor = vm["decreaseDivisor"].as<unsigned>();
		const unsigned increaseIters = vm["increaseIters"].as<unsigned>();
		const unsigned bjSize = vm["bjSize"].as<unsigned>();
		const unsigned maxD = vm["maxD"].as<unsigned>();
		const unsigned maxC = vm["maxC"].as<unsigned>();
		const double acceptAlpha = vm["acceptAlpha"].as<double>();
		//const unsigned perturbSize = vm["perturbSize"].as<unsigned>();
		const unsigned sizePBubble = vm["sizePBubble"].as<unsigned>();
		const unsigned sizePInsa = vm["sizePInsa"].as<unsigned>();
		const bool timeLog = vm["timeLog"].as<bool>();
		const bool lowerBoundOrder = vm["lowerBoundOrder"].as<bool>();
		const string perturbTypeStr = vm["perturbType"].as<string>();
		/*unsigned perturbType;
		if (perturbTypeStr.compare("insa") == 0)
			perturbType = INSA_PERTURB;
		else if (perturbTypeStr.compare("bubble") == 0)
			perturbType = BUBBLE_PERTURB;
		else if (perturbTypeStr.compare("mix") == 0)
			perturbType = MIX_PERTURB;
		else if (perturbTypeStr.compare("gt") == 0)
			perturbType = GT;
		else
			throw errorText("Invalid opetion for perturbation","main","main");*/
		const double scaleTime = vm["scaleTime"].as<double>();
		const string scheduleFile = vm["scheduleFile"].as<string>();
		const double bubbleP = vm["bubbleP"].as<double>();
		const double mixProb = vm["mixProb"].as<double>();
		const string searchTypeStr = vm["searchType"].as<string>();
		const bool onlyLowerBound = vm["onlyLowerBound"].as<bool>();


		if(onlyLowerBound) {
			inst.parse(instPath);

			//cout << name << " " << initialjumpLimit << "  " << inst.lowerBoundMss() << " " << inst.lowerBoundTkz() << endl;
			cout << name << " " << inst.lowerBoundTkz() << endl;
			return 0;
		}
		

		// (3.1) adjusting time
		maxSecs *= scaleTime;
	
		//(4) vesrbose if in DEBUG
#ifndef NDEBUG
		cout << "DEBUG MODE" << endl;
		cout << "instPath:" << instPath << endl;
		cout << "maxSecs: " << maxSecs << "   seed: " << seed << "   tenure: " << tenure << "   initialjumpLimit: " << initialjumpLimit  << "   decreaseDivisor: " << decreaseDivisor << "   increaseIters: " << increaseIters <<"   bjSize: " << bjSize << "   maxD: " << maxD << "   maxC: " << maxC << "   acceptAlpha: " << acceptAlpha << "   sizePBubble: " << sizePBubble <<  "   sizePInsa: " << sizePInsa << "   timeLog: " << (timeLog ? "Y" : "n")  << "   lowerBoundOrder: " << (lowerBoundOrder ? "Y" : "n")  << "   perturbTypeStr: " << perturbTypeStr << "   scaleTime: " << scaleTime << "   scheduleFile: " << scheduleFile << "    bubbleP: " << bubbleP << "   mixProb: " << mixProb << endl;
#endif

		//(5) execute
		Tabu::nosmuTabu(instPath, name, tenure, initialjumpLimit, decreaseDivisor, bjSize, maxD, maxC, GT, 1000*maxSecs, timeLog);
		
	}catch(const string & e) {
		cout << e << endl;
	}

	return 0;
}
