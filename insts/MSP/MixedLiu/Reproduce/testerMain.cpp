/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#include "Settings.hpp"
//#include "Parser.hpp"
#include "InstGen.hpp"


using namespace std;


int main(int argc, char *argv[]) {

	//srand(17);
	try{

		vector<unsigned> noSizeVec = {5, 10, 15, 20};
		vector<unsigned> instIndVec = {16, 24, 29, 36};
		InstGen ig;
		
		srand(13);

		string inPath;
		string outPath;

		for(unsigned instInd : instIndVec) {
			inPath = "./SourceInsts/la"+unsigStr(instInd)+"_raw.txt";
			cout << "Generating from: " << inPath << endl;
			for(unsigned noSize : noSizeVec) {
				cout << "\tno " << noSize << endl;
				for(unsigned rep=0; rep<20;rep++) {
					if(rep>=10)
						outPath = "../mixed_base"+unsigStr(instInd)+"_no"+unsigStr(noSize)+"_rep"+unsigStr(rep)+".psi";
					else
						outPath = "../mixed_base"+unsigStr(instInd)+"_no"+unsigStr(noSize)+"_rep0"+unsigStr(rep)+".psi";

					InstGen::nonTaiJspToMixed(inPath, outPath, noSize);
					//ig =  InstGen::fromRawPsspToMixed("./JSSPinstNonTai/la"+unsigStr(instInd)+"_raw.txt", noSize);
					//ig.printPsp("./mixed_base"+unsigStr(instInd)+"_no"+unsigStr(noSize)+"_rep"+unsigStr(rep),true);
				}
			}
		}
		
		//InstGen::whizzkidToPsi("../Insts/SOURCES/onlyGSS_nativeFormat/whizzkids97.gss", "./APAGAR.txt");

		//InstGen::bloomGroupToPsiDir("../Insts/SOURCES/onlyGS_nativeFormat", "gss", "./TESTE");

		//InstGen::bloomGroupToPsi("../Insts/SOURCES/onlyGSS_nativeFormat/abz7_7.gss", "./APAGAR.txt");

		//InstGen::admuToStage("./AdmuJsp/jcmax.txt", "./NasStagePartitions", "./AdmuJsp");

		//InstGen::spToActualStage("./TaiJspTaiFormatOnly50", "./NasStagePartitions", "./ActualStage");

		//makeNasActualStage(const string & jspTaiFormatPath, const string & nasStagesDir)
		//InstGen ig = InstGen::makeNasActualStage("./TaiJspTaiFormat/ta01.txt", "./NasStagePartitions");
		//ig.printPsp("", false);

		/*
		vector<unsigned> noSizeVec = {5, 10, 15, 20};
		vector<unsigned> instIndVec = {16, 24, 29, 36};
		InstGen ig;
		
		for(unsigned instInd : instIndVec) {
			for(unsigned noSize : noSizeVec) {
				for(unsigned rep=0; rep<10;rep++) {
					ig =  InstGen::fromRawPsspToMixed("./JSSPinstNonTai/la"+unsigStr(instInd)+"_raw.txt", noSize);
					ig.printPsp("./mixed_base"+unsigStr(instInd)+"_no"+unsigStr(noSize)+"_rep"+unsigStr(rep),true);
				}
			}
			}*/

		//InstGen ig = InstGen::fromRawPsspToMixed("./JSSPinstNonTai/la16_raw.txt", 5);
		//ig.printPsp("", false);

		//unsigned reps = 10;
		//InstGen::batchJspToStage("./TaiJspTaiFormat","./StageShopNasInsts", reps);

		//InstGen ig = InstGen::jspToStage("./TaiJspTaiFormat/ta01.txt");
		//ig.printPsp("", false);

		/*	
		InstGen ig;
		ig = InstGen::jspTaiFormatToPsi("./TaiJspTaiFormat/ta01.txt");

		//printPsp(const string & pspPath, bool toFile)		
		ig.printPsp("./TaiJspTaiFormat/ta01.txt",false);

		//vector<unsigned> part = InstGen::genStagePartition(10);
		//for(unsigned u : part) cout << u << " ";
		//cout << endl;

		//ig.setNasiriStage();
		*/
		
	} catch(const string e) {
		cout << e << endl;
	}
	return 0;
}

	

