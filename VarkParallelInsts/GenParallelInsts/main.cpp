/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace boost;


#include "../Settings.hpp"
#include "InstGen.hpp"

int main(int argc, char *argv[]) {

	string fileName;

	srand(10);
	
	//parallel jss insts gholami new tab 3 
	unsigned maxK = 5;
	unsigned m;
	unsigned j;
	unsigned index = 0;
	for(unsigned x=3; x<=10; x++) {
		m = x;
		j = x;
		for(unsigned maxP = 10; maxP<=1000; maxP*=10) {
			index++;
			//for(unsigned rep=0; rep<10; rep++) { 
				fileName = "./GholamiExtraInstsNoRep/GhSo_"+unsigStr(index)+"noRep_jm_"+unsigStr(x)+"_maxP_"+unsigStr(maxP-1)+".ppsi";
				InstGen::genRandInst(fileName, j ,m, maxK, maxP-1);
			//}
		}
	}
	
	
	
	
	
	
	for(unsigned x=10; x<=160; x+=10) {
		m = x;
		j = x;
		index++;
		const unsigned maxP = 99;
		//for(unsigned rep=0; rep<10; rep++) {
			fileName = "./GholamiExtraInstsNoRep/GhSo_"+unsigStr(index)+"noRep_jm_"+unsigStr(x)+"_maxP_"+unsigStr(maxP)+".ppsi";
			InstGen::genRandInst(fileName, j ,m, maxK, maxP);
		//}
	}
	
	
	
	/*
	//parallel jss tkz new insts
	unsigned maxP = 99;
	for(unsigned j=5; j<=50; j+=5) {
		for(unsigned m=5; m<=30; m+=5) {
			for(unsigned maxK=2; maxK<=6; maxK++) {
				fileName = "./GeneratedInsts/genParallel_J_"+unsigStr(j)+"_M_"+unsigStr(m)+"_maxK_"+unsigStr(maxK)+".ppsi";
				InstGen::genRandInst(fileName, j ,m, maxK, maxP);
			}
		}
	}
	*/
	
	//fileName = "";
	//fileName = "./GeneratedInsts/test.pjss";
	//InstGen::genRandInst(fileName, 5 ,5, 3, 99);

	return 0;
}
