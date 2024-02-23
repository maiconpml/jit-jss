/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace boost;


#include "../Settings.hpp"
#include "InstGen.hpp"


int main(int argc, char *argv[]) {
	try{
		
		// (1) process commandline options
		po::options_description desc("Options");
		desc.add_options()
			("help", "show help")
			("mode", "mode: JspToParallel, MapPssToMss, MapPssToGss") //options: JspToParallel, MapPssToMss
			("fromPath", "instance source")
			("targetPath", "instance target")
			("parallelK", po::value<unsigned>()->default_value(0), "number of parallel machines")
			;
		po::positional_options_description pod;
		pod.add("mode", 1);
		pod.add("fromPath", 1); 
		pod.add("targetPath", 1);
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
		po::notify(vm);

		// (2) check necessaary options
		if (vm.count("help") || !vm.count("mode")) {
			cout << desc << endl;
			return 0;
		}

		// (3) initialize params
		const string mode = vm["mode"].as<string>();

		PsiInst pi;
		
		if(mode.compare("JspToParallel")==0) {
			if (!vm.count("fromPath") || !vm.count("targetPath")  ||  !vm.count("parallelK")) {
				cout << desc << endl;
				return 0;
			}
			const string fromInstPath = vm["fromPath"].as<string>();
			const string toInstPath = vm["targetPath"].as<string>();
			const unsigned parallelK = vm["parallelK"].as<unsigned>();
			pi.readPsiInst(fromInstPath);
			pi.printParallel(toInstPath, parallelK);
		} else if(mode.compare("MapPssToMss")==0) {
			if (!vm.count("fromPath") || !vm.count("targetPath")) {
				cout << desc << endl;
				return 0;
			}
			const string fromInstPath = vm["fromPath"].as<string>();
			const string toInstPath = vm["targetPath"].as<string>();
			pi.mapPssToMSS(fromInstPath, toInstPath);
		} else if(mode.compare("MapPssToGss")==0) {
			if (!vm.count("fromPath") || !vm.count("targetPath")) {
				cout << desc << endl;
				return 0;
			}
			const string fromInstPath = vm["fromPath"].as<string>();
			const string toInstPath = vm["targetPath"].as<string>();
			pi.mapPssToGSS(fromInstPath, toInstPath);
		} else {
			cout << desc << endl;
			return 0;
		}

	}catch(const string & e) {
		cout << e << endl;
	}		

	return 0;
}
