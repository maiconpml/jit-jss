/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#pragma once


#include "Settings.hpp"
#include "Utilities.hpp"
#include "italiano.hpp"

#include <vector>
#include <set>
#include <deque>

#include <tuple>

using namespace std;

//#########################################################################################################################################################################################
class Arrow {
public:
	Arrow() {  }
	Arrow(unsigned from, unsigned to) {
		this->from = from;
		this->to = to;
	}
	~Arrow() {  }
	
	string toString() const {
		return unsigStr(from)+"->"+unsigStr(to);
	}
	
	bool operator==(const Arrow & arrow) const {
		return this->from == arrow.from   &&   this->to == arrow.to;
	}
	
	bool operator<(const Arrow & arrow) const {
		if(this->from < arrow.from) {
			return true;
		} else if(arrow.from < this->from) {
			return false;
		} else {
			if(this->to < arrow.to)
				return true;
			else
				return false;
		}
	}
	
	unsigned from;
	unsigned to;
};
//#########################################################################################################################################################################################
class InstGen {
public:
	InstGen(unsigned _nJobs, unsigned _nMachs) : nJobs(_nJobs), nMachs(_nMachs) {
		assert(nJobs <= MAX_ELEMENTS);
		jobCGs.resize(nJobs, ClosedGraph<MAX_ELEMENTS>(nMachs));
		jobSets.resize(nJobs);
	}
	InstGen() {  }
	~InstGen() {  }

	//Liu instnace for MSP
	static void nonTaiJspToMixed(const string & inPath, const string & outPath, unsigned n0Size) {

		cout << "\t\t\t\t  AAAA" << endl;
		
	    InstGen ig;
		ifstream fromFile;

		fromFile.open(inPath);
		if( ! fromFile.is_open())
			throw errorText("Could not open partitionedraw nonTaillard's instance file: "+inPath, "Parser.hpp", "fromRawPsspToMixed");

		fromFile >> ig.nJobs;
		ig.nJobs+=n0Size;
		fromFile >> ig.nMachs;

		vector<unsigned> order;
		unsigned aOrder;
		vector<unsigned> costs;
		unsigned aCost;

		ig.procTs.resize(ig.nJobs, vector<unsigned>(ig.nMachs));
		ig.jobCGs.resize(ig.nJobs, ClosedGraph<MAX_ELEMENTS>(ig.nMachs));
		ig.jobSets.resize(ig.nJobs);
		
		//Setting NJ times
		for(unsigned j=0; j<ig.nJobs-n0Size; j++) {
			order.clear();
			costs.clear();
			for(unsigned m=0; m<ig.nMachs; m++) {
				fromFile >> aOrder;
				fromFile >> aCost;
				order.push_back(aOrder);
				costs.push_back(aCost);
			}

			for(unsigned index=0; index<ig.nMachs; index++) {
				ig.procTs[j][order[index]] = costs[index];
			}

			//Add arrows to NJ
	    	 for(unsigned index=0; index<ig.nMachs-1; index++) {
				 ig.addArrow(j, Arrow(order[index], order[index+1]));
			}
		}

		//Setting NO times
		for(unsigned j=ig.nJobs-n0Size; j<ig.nJobs; j++) {
			for(unsigned index=0; index<ig.nMachs; index++) {
				ig.procTs[j][index] = rand()%98+1;//CONFIRM LIMITS
			}
		}

		ig.printPsp(outPath, true);
	}

	static void bloomGroupToPsiDir(const string &  inDirPath, const string inFIleTerm, const string &  outDirPath) {
		vector<string> fileNames =  filesInDir(inDirPath, inFIleTerm);
		string newEndName;

		for(const string & fn : fileNames) {
			newEndName = fn;
			assert(newEndName.size()>3);
			newEndName[newEndName.size()-1] = 'i';
			newEndName[newEndName.size()-2] = 's';
			newEndName[newEndName.size()-3] = 'p';

			bloomGroupToPsi(inDirPath+"/"+fn, outDirPath+"/"+newEndName);
		}
	}


	static void whizzkidToPsi(const string &  inFilePath, const string &  outFilePath) {
		ifstream fromFile;
		string bufString;

		InstGen groupInst;

		unsigned nJobs;
		unsigned nMachs;

		unsigned curSize;
		unsigned baseI;
		unsigned fromOp;
		unsigned toOp;
		vector<vector<bool>> usedOps;
		vector<unsigned> usedOpsPerJob;
		unsigned nUsedOps = 0;
		unsigned bufU;

		vector<vector<unsigned>> readTimes;
		vector<vector<unsigned>> readOrder;
		vector<unsigned> readGroups;
		vector<unsigned> adjustedGroups;
		vector<vector<unsigned>> groupSizes;

		fromFile.open(inFilePath);
		if( ! fromFile.is_open())
			throw errorText("Could not open Bloom groups  instance file: "+inFilePath, "InstGen.hpp", "InstGen::bloomGroupToPsi()");
		
	    fromFile >> nJobs;
		assert(nJobs>0   &&   nJobs<1000);
		fromFile >> nMachs;
		assert(nMachs>0   &&   nMachs<1000);

		readTimes.resize(nJobs, vector<unsigned>(nMachs));
		readOrder.resize(nJobs, vector<unsigned>(nMachs));
		readGroups.resize(nJobs*nMachs);

		usedOps.resize(nJobs, vector<bool>(nMachs, false));
		usedOpsPerJob.resize(nJobs);

		groupInst.procTs.resize(nJobs, vector<unsigned>(nMachs, 0));		
		groupInst.jobCGs.resize(nJobs, ClosedGraph<MAX_ELEMENTS>(nMachs));
		groupInst.jobSets.resize(nJobs);

		adjustedGroups.resize(nJobs*nMachs);

	    for(unsigned j=0; j<nJobs; j++) {
			fromFile >> usedOpsPerJob[j];
			//assert(bufU == nMachs);
			for(unsigned m=0; m<usedOpsPerJob[j]; m++) {
				fromFile >> readOrder[j][m];
				usedOps[j][readOrder[j][m]] = true;
				nUsedOps++;
				fromFile >> readTimes[j][m];
			}
		}

		for(unsigned u=0; u<nUsedOps; u++) {
			fromFile >> readGroups[u];
		}

		fromFile.close();

		for(unsigned j=0; j<nJobs; j++) {
			bufU = usedOpsPerJob[j];
			for(unsigned m=0; m<nMachs; m++) {
				if( ! usedOps[j][m]) {
					assert(bufU<nMachs);
					readOrder[j][bufU] = m;
					readTimes[j][bufU] = 0;
					bufU++;
				}
			}
		}

		unsigned bailiff = 0;
		unsigned provost = 0;
		unsigned delta = 0;
		for(unsigned j=0; j<nJobs; j++) {
			for(unsigned u=0; u<usedOpsPerJob[j]; u++) {
				adjustedGroups[provost] = readGroups[bailiff]+delta;
				bailiff++;
				provost++;
			}
			for(unsigned u=usedOpsPerJob[j]; u<nMachs; u++) {
				assert(bailiff<readGroups.size());
				delta++;
				adjustedGroups[provost] = readGroups[bailiff-1]+delta;
				provost++;
			}
		}
		
		for(unsigned j=0; j<nJobs; j++) {
			for(unsigned u=0; u<nMachs; u++) {
				groupInst.procTs[j][readOrder[j][u]] = readTimes[j][u];
			}
		}
		/*
		cout << "\nreadTimes:\n";
		for(const vector<unsigned> & v : readTimes) {
			cout << "\t";
			for(unsigned u : v) {
				cout << u << " ";
			}
			cout << endl;
		}

		cout << "\nreadOrder:\n";
		for(const vector<unsigned> & v : readOrder) {
			cout << "\t";
			for(unsigned u : v) {
				cout << u << " ";
			}
			cout << endl;
		}

		cout << "\nadjustedGroups: ";
		for(unsigned u : adjustedGroups) {
			cout << u << " ";
		}
		getchar();*/

		
		curSize = 1;
		baseI = 0;
		groupSizes.resize(nJobs, vector<unsigned>());
		for(unsigned j=0; j<nJobs; j++) {
			for(unsigned m=0; m<nMachs; m++) {
				baseI++;
				if(adjustedGroups[baseI] == adjustedGroups[baseI-1]) {
					curSize++;
				} else if(baseI < nJobs*nMachs) {
					assert(adjustedGroups[baseI] == adjustedGroups[baseI-1]+1);
					groupSizes[j].push_back(curSize);
					curSize = 1;
				}
			}
		}
		groupSizes[nJobs-1].push_back(curSize);

	    for(unsigned j=0; j<nJobs; j++) {
			baseI = 0;
			for(unsigned bailiff=0; bailiff<groupSizes[j].size()-1; bailiff++) {
				for(unsigned prevStageI=0; prevStageI<groupSizes[j][bailiff]; prevStageI++) {
					fromOp = readOrder[j][baseI + prevStageI];
					for(unsigned nextStageI=0; nextStageI<groupSizes[j][bailiff+1]; nextStageI++) {
						toOp = readOrder[j][baseI + groupSizes[j][bailiff]  + nextStageI];
						//cout << "from: " << fromOp << " ; toOp: " << toOp << endl;
						//getchar();
						groupInst.addArrow(j, Arrow(fromOp, toOp));
					}
				}
				baseI += groupSizes[j][bailiff];
			}
			//cout << endl;
		}

		//groupInst.printPsp("", false);
		//getchar();
		//stageInst.printPsp(outDir+"/admuStage"+(string)(admuI<10 ? "0" : "")+unsigStr(admuI)+".psi", true);
		groupInst.printPsp(outFilePath, true);
	}


	static void bloomGroupToPsi(const string &  inFilePath, const string &  outFilePath) {

		ifstream fromFile;
		string bufString;

		InstGen groupInst;

		unsigned nJobs;
		unsigned nMachs;

		unsigned bufU;
		unsigned curSize;
		unsigned baseI;
		unsigned fromOp;
		unsigned toOp;

		vector<vector<unsigned>> readTimes;
		vector<vector<unsigned>> readOrder;
		vector<unsigned> readGroups;
		vector<vector<unsigned>> groupSizes;

		fromFile.open(inFilePath);
		if( ! fromFile.is_open())
			throw errorText("Could not open Bloom groups  instance file: "+inFilePath, "InstGen.hpp", "InstGen::bloomGroupToPsi()");
		
	    fromFile >> nJobs;
		assert(nJobs>0   &&   nJobs<1000);
		fromFile >> nMachs;
		assert(nMachs>0   &&   nMachs<1000);

		readTimes.resize(nJobs, vector<unsigned>(nMachs));
		readOrder.resize(nJobs, vector<unsigned>(nMachs));
		readGroups.resize(nJobs*nMachs);

		groupInst.procTs.resize(nJobs, vector<unsigned>(nMachs, 0));		
		groupInst.jobCGs.resize(nJobs, ClosedGraph<MAX_ELEMENTS>(nMachs));
		groupInst.jobSets.resize(nJobs);

	    for(unsigned j=0; j<nJobs; j++) {
			fromFile >> bufU;
			assert(bufU == nMachs);
			for(unsigned m=0; m<nMachs; m++) {
				fromFile >> readOrder[j][m];
				fromFile >> readTimes[j][m];
			}
		}

		for(unsigned u=0; u<nJobs*nMachs; u++) {
			fromFile >> readGroups[u];
		}

		fromFile.close();

		for(unsigned j=0; j<nJobs; j++) {
			for(unsigned u=0; u<nMachs; u++) {
				groupInst.procTs[j][readOrder[j][u]] = readTimes[j][u];
			}
		}

		curSize = 1;
		baseI = 0;
		groupSizes.resize(nJobs, vector<unsigned>());
		for(unsigned j=0; j<nJobs; j++) {
			for(unsigned m=0; m<nMachs; m++) {
				baseI++;
				if(readGroups[baseI] == readGroups[baseI-1]) {
					curSize++;
				} else if(baseI < nJobs*nMachs) {
					assert(readGroups[baseI] == readGroups[baseI-1]+1);
					groupSizes[j].push_back(curSize);
					curSize = 1;
				}
			}
		}
		groupSizes[nJobs-1].push_back(curSize);

		/*for(unsigned j=0; j<nJobs; j++) {
			for(unsigned gs : groupSizes[j]) {
				cout << gs << " ";
			}
			cout << endl;
			}*/

	    for(unsigned j=0; j<nJobs; j++) {
			baseI = 0;
			for(unsigned bailiff=0; bailiff<groupSizes[j].size()-1; bailiff++) {
				for(unsigned prevStageI=0; prevStageI<groupSizes[j][bailiff]; prevStageI++) {
					fromOp = readOrder[j][baseI + prevStageI];
					for(unsigned nextStageI=0; nextStageI<groupSizes[j][bailiff+1]; nextStageI++) {
						toOp = readOrder[j][baseI + groupSizes[j][bailiff]  + nextStageI];
						//cout << "from: " << fromOp << " ; toOp: " << toOp << endl;
						//getchar();
						groupInst.addArrow(j, Arrow(fromOp, toOp));
					}
				}
				baseI += groupSizes[j][bailiff];
			}
			//cout << endl;
		}

		//groupInst.printPsp("", false);
		//getchar();
		//stageInst.printPsp(outDir+"/admuStage"+(string)(admuI<10 ? "0" : "")+unsigStr(admuI)+".psi", true);
		groupInst.printPsp(outFilePath, true);
	}


	static void admuToStage(const string &  inFilePath, const string & nasStagesDir, const string & outDir) {
		
		ifstream fromFile;
		string bufString;

		InstGen stageInst;

		unsigned nJobs;
		unsigned nMachs;

		unsigned fromOp;
		unsigned toOp;
		unsigned baseIndex;

		unsigned admuI = 1;

	    vector<vector<unsigned>> readTimes;
		vector<vector<unsigned>> readOrder;

		vector<vector<unsigned>> allStagePart;
		
		fromFile.open(inFilePath);
		if( ! fromFile.is_open())
			throw errorText("Could not open Admu all instance file: "+inFilePath, "Parser.hpp", "Parser::sepAdmu()");

		while(true) {
			for(unsigned i=0; i<12; i++) {
				fromFile >> bufString;
				//cout << bufString << endl;
				//getchar();
			}
			if(fromFile.eof())
				break;

			fromFile >> nJobs;
			fromFile >> nMachs;

			readTimes.clear();
			readOrder.clear();
			readTimes.resize(nJobs, vector<unsigned>(nMachs));
			readOrder.resize(nJobs, vector<unsigned>(nMachs));

		   for(unsigned j=0; j<nJobs; j++) {
				for(unsigned m=0; m<nMachs; m++) {
					fromFile >> readOrder[j][m];
					fromFile >> readTimes[j][m];
				}
			}

			allStagePart = readNasPart(nasStagesDir, nJobs, nMachs);

			stageInst = InstGen();
			stageInst.procTs.resize(nJobs, vector<unsigned>(nMachs, 0));
			for(unsigned j=0; j<nJobs; j++) {
				for(unsigned u=0; u<nMachs; u++) {
					stageInst.procTs[j][readOrder[j][u]-1] = readTimes[j][u];
				}
			}
			stageInst.jobCGs.resize(nJobs, ClosedGraph<MAX_ELEMENTS>(nMachs));
			stageInst.jobSets.resize(nJobs);

			/*for(unsigned j=0; j<nJobs; j++) {
				for(unsigned m=0; m<nMachs; m++) {
				    cout << readOrder[j][m] << "\t";
				}
				cout << endl;
			}
			cout << endl;
			for(unsigned j=0; j<nJobs; j++) {
				for(unsigned m=0; m<nMachs; m++) {
				    cout << readTimes[j][m] << "\t";
				}
				cout << endl;
			}
			getchar();*/

			for(unsigned j=0; j<nJobs; j++) {
				baseIndex = 0;
				for(unsigned bailiff=0; bailiff<allStagePart[j].size()-1; bailiff++) {
					for(unsigned prevStageI=0; prevStageI<allStagePart[j][bailiff]; prevStageI++) {
						fromOp = readOrder[j][baseIndex + prevStageI]-1;
						for(unsigned nextStageI=0; nextStageI<allStagePart[j][bailiff+1]; nextStageI++) {
							toOp = readOrder[j][baseIndex + allStagePart[j][bailiff]  + nextStageI]-1;
							//cout << "from: " << fromOp << " ; toOp: " << toOp << endl;
							//getchar();
							stageInst.addArrow(j, Arrow(fromOp, toOp));
						}
					}
					baseIndex += allStagePart[j][bailiff];
				}
				//cout << endl;
			}

			//stageInst.printPsp("", false);
			//getchar();
			stageInst.printPsp(outDir+"/admuStage"+(string)(admuI<10 ? "0" : "")+unsigStr(admuI)+".psi", true);
			admuI++;
		}
	}



	//spToActualStage(const string & jspTaiFormatPath, const string & nasStagesDir)
	//expects jsp to be in tai format
	static void spToActualStage(const string & inJspDir, const string & nasStagesDir, const string & outDir) {
		vector<string> fileNames =  filesInDir(inJspDir, ".txt");
		sort(fileNames.begin(), fileNames.end(), lexiOrder);

		InstGen ig;
		unsigned taIndex=0;

		for(const string & fileName : fileNames) {
			taIndex++;
			
			ig = jspToActualStage(inJspDir+"/"+fileName, nasStagesDir);
			ig.printPsp(outDir+"/actStage"+(string)(taIndex<10 ? "0" : "")+unsigStr(taIndex)+".psi", true);
		}
	}


	

	//expects files in Tai format
	static void batchJspToStage(const string & inDirPath, const string & outDirPath, unsigned reps) {

		vector<string> fileNames =  filesInDir(inDirPath, ".txt");
		sort(fileNames.begin(), fileNames.end(), lexiOrder);

	    InstGen ig;
		unsigned taIndex=0;
		
		for(const string & fileName : fileNames) {
			taIndex++;

			for(unsigned r=0; r<reps; r++) {
				ig = jspToStage(inDirPath+"/"+fileName);
				ig.printPsp(outDirPath+"/stageNas"+(string)(taIndex<10 ? "0" : "")+unsigStr(taIndex)+"_"+unsigStr(r)+".psi", true);
			}
		}
	} 



	//will add if awwor is not in set, does not check transitive closure
	void addArrow(unsigned jobIndex, const Arrow & a) {
		assert(jobIndex < jobCGs.size());
		assert( ! jobCGs[jobIndex].path(a.to, a.from));//forming cycle?
		
		unsigned oldSize = jobSets[jobIndex].size();
		jobSets[jobIndex].insert(a);
		if(oldSize < jobSets[jobIndex].size()) {
			if( ! jobCGs[jobIndex].path(a.from, a.to))
				jobCGs[jobIndex].add(a.from, a.to);
		} else {
			assert(jobCGs[jobIndex].path(a.from, a.to));
		}
	}
	
	
	
	void clear() {
		nJobs = 0;
		nMachs = 0;
		procTs.clear();
		jobCGs.clear();
		jobSets.clear();
	}



	void printPsp(const string & pspPath, bool toFile) const {
		assert(jobSets.size() > 0);
		assert(procTs.size() == jobSets.size());
		
		string str ="#START#\n\nTotalJobs: ";
		
		str.append(unsigStr(procTs.size()));
		str.append(" TotalMachines: ");
		str.append(unsigStr(procTs[0].size()));
		
		str.append("\nCosts:");
		for(const vector<unsigned> & v : procTs) {
			str.append("\n\t");
			for(unsigned u : v) {
				str.append(unsigStr(u)+" ");
			}
		}
		str.append("\n");
		
		for(const set<Arrow> & sa : jobSets) {
			str.append("\nJob: "+unsigStr(sa.size())+"\n");
			for(const Arrow & a : sa) {
				str.append("\n\t"+unsigStr(a.from)+" -> "+unsigStr(a.to));
			}
			str.append("\n");
		}
		
		if(toFile) {
			ofstream toFile;
			toFile.open(pspPath);
			if( ! toFile.is_open())
				throw errorText("Could not open file to write psp instance: "+pspPath, "Parser.hpp", "Parser::printPsp()");
			toFile << str; 
			toFile.close();
		} else {
			cout << str << endl;
		}
	}


	//@return nJobs, nMachs, times, orders
	static tuple<unsigned, unsigned, vector<vector<unsigned>>,vector<vector<unsigned>>> readTaiJspInfo(const string & filePath) {

		unsigned nJobs;
		unsigned nMachs;
		vector<vector<unsigned>> readTimes;
		vector<vector<unsigned>> readOrder;

		ifstream fromFile;
		string bufStr;
		
		fromFile.open(filePath);
		if( ! fromFile.is_open())
			throw errorText("Could not open partitioned Taillard's instance file: "+filePath, "Parser.hpp", "Parser::jspTaiFormatToPsi()");
			
		fromFile >> bufStr;
		assert(bufStr.compare("Jobs:")==0);
		fromFile >> nJobs;
		
		fromFile >> bufStr;
		assert(bufStr.compare("Machs:")==0);
		fromFile >> nMachs;
		
		assert(nJobs>0);
		assert(nJobs<10000);
		assert(nMachs>0);
		assert(nMachs<10000);
		
		fromFile >> bufStr;
		assert(bufStr.compare("Times:")==0);
		readTimes.resize(nJobs, vector<unsigned>(nMachs));
		for(unsigned j=0; j<nJobs; j++) {
			for(unsigned m=0; m<nMachs; m++) {
				fromFile >> readTimes[j][m];
			}
		}
		
		fromFile >> bufStr;
		assert(bufStr.compare("Order:")==0);
		readOrder.resize(nJobs, vector<unsigned>(nMachs));
		for(unsigned j=0; j<nJobs; j++) {
			for(unsigned m=0; m<nMachs; m++) {
				fromFile >> readOrder[j][m];
			}
		}
		fromFile.close();

		return  tuple<unsigned, unsigned, vector<vector<unsigned>>,vector<vector<unsigned>>> (nJobs, nMachs, readTimes, readOrder);
	}


	
	static InstGen jspTaiFormatToPsi(const string & filePath) {
		unsigned nJobs;
		unsigned nMachs;
		vector<vector<unsigned>> readTimes;
		vector<vector<unsigned>> readOrder;
		
		tie(nJobs, nMachs, readTimes, readOrder) = readTaiJspInfo(filePath);
		
		InstGen gen(nJobs, nMachs);
		
		gen.procTs.resize(nJobs, vector<unsigned>(nMachs, 0));
		for(unsigned j=0; j<nJobs; j++) {
			for(unsigned u=0; u<nMachs; u++) {
				gen.procTs[j][readOrder[j][u]-1] = readTimes[j][u];
			}
		}
		//cout <<"\n!!!\n";
		//getchar();
		
		gen.jobCGs.clear(); 
		gen.jobSets.clear();
		
		gen.jobCGs.resize(nJobs, ClosedGraph<MAX_ELEMENTS>(nMachs));
		gen.jobSets.resize(nJobs);
		
		for(unsigned j=0; j<nJobs; j++) {
			for(unsigned u=0; u<nMachs-1; u++) {
				gen.addArrow(j, Arrow(readOrder[j][u]-1, readOrder[j][u+1]-1));
			}
		}

		
		return gen;
	}



	static vector<unsigned> genStagePartition(unsigned size) {
	
		unsigned nParts;
		vector<unsigned> bloat;

		nParts = 1+(rand() % size);
		bloat = vector<unsigned>(nParts,1);

		assert(nParts <= size);

		for(unsigned insertI=0; insertI<size-nParts; insertI++) {
			bloat[rand() % nParts]++;
		}

	    return bloat;
	}



	static InstGen jspToActualStage(const string & jspTaiFormatPath, const string & nasStagesDir) {
		
		InstGen stageInst;
		vector<vector<unsigned>> allStagePart;

		unsigned nJobs;
		unsigned nMachs;
		vector<vector<unsigned>> readTimes;
		vector<vector<unsigned>> readOrder;
		
		unsigned fromOp;
		unsigned toOp;
		unsigned baseIndex;
		
		tie(nJobs, nMachs, readTimes, readOrder) = readTaiJspInfo(jspTaiFormatPath);
		
		assert(nMachs > 0);
		assert(nJobs > 0);

		stageInst.procTs.resize(nJobs, vector<unsigned>(nMachs, 0));
		for(unsigned j=0; j<nJobs; j++) {
			for(unsigned u=0; u<nMachs; u++) {
				stageInst.procTs[j][readOrder[j][u]-1] = readTimes[j][u];
			}
		}
		stageInst.jobCGs.resize(nJobs, ClosedGraph<MAX_ELEMENTS>(nMachs));
		stageInst.jobSets.resize(nJobs);
		
		allStagePart = readNasPart(nasStagesDir, nJobs, nMachs);

		for(unsigned j=0; j<nJobs; j++) {
			baseIndex = 0;
		    for(unsigned bailiff=0; bailiff<allStagePart[j].size()-1; bailiff++) {
				for(unsigned prevStageI=0; prevStageI<allStagePart[j][bailiff]; prevStageI++) {
					fromOp = readOrder[j][baseIndex + prevStageI]-1;
					for(unsigned nextStageI=0; nextStageI<allStagePart[j][bailiff+1]; nextStageI++) {
						toOp = readOrder[j][baseIndex + allStagePart[j][bailiff]  + nextStageI]-1;
						//cout << "from: " << fromOp << " ; toOp: " << toOp << endl;
						stageInst.addArrow(j, Arrow(fromOp, toOp));
					}
				}
				baseIndex += allStagePart[j][bailiff];
			}
			//cout << endl;
		}

		return stageInst;
		
	}



	static vector<vector<unsigned>> readNasPart(const string & stageDir, unsigned nJobs, unsigned nMachs) {
		
		const string stageFilePath = stageDir+"/"+unsigStr(nJobs)+"_"+unsigStr(nMachs)+".txt";
		ifstream fromFile;
		vector<vector<unsigned>> nasPart(nJobs);
		unsigned partSize;

		fromFile.open(stageFilePath);
		if( ! fromFile.is_open())
			throw errorText("Could not open stage partition file: "+stageFilePath, "InstGen.hpp", "readNasPart");

		for(unsigned j=0; j<nJobs; j++) {
			fromFile >> partSize;
			assert(partSize != 0);
			nasPart[j].resize(partSize);
			for(unsigned p=0; p<partSize; p++) {
				fromFile >> nasPart[j][p];
			}
		}
		
		return nasPart;
	}



	static InstGen jspToStage(const string & filePath) {
		
		InstGen stageInst;
		vector<unsigned> stagePart;

		unsigned nJobs;
		unsigned nMachs;
		vector<vector<unsigned>> readTimes;
		vector<vector<unsigned>> readOrder;
		
		unsigned fromOp;
		unsigned toOp;
		unsigned baseIndex;
		
		tie(nJobs, nMachs, readTimes, readOrder) = readTaiJspInfo(filePath);
		
		assert(nMachs > 0);
		assert(nJobs > 0);

		stageInst.procTs.resize(nJobs, vector<unsigned>(nMachs, 0));
		for(unsigned j=0; j<nJobs; j++) {
			for(unsigned u=0; u<nMachs; u++) {
				stageInst.procTs[j][readOrder[j][u]-1] = readTimes[j][u];
			}
		}
		stageInst.jobCGs.resize(nJobs, ClosedGraph<MAX_ELEMENTS>(nMachs));
		stageInst.jobSets.resize(nJobs);
		
		for(unsigned j=0; j<nJobs; j++) {
			stagePart = genStagePartition(nMachs);

			//for(unsigned u : stagePart) cout << u << " ";
			//cout << endl << endl;
			
			baseIndex = 0;
			for(unsigned bailiff=0; bailiff<stagePart.size()-1; bailiff++) {
				for(unsigned prevStageI=0; prevStageI<stagePart[bailiff]; prevStageI++) {
					fromOp = readOrder[j][baseIndex + prevStageI]-1;
					for(unsigned nextStageI=0; nextStageI<stagePart[bailiff+1]; nextStageI++) {
						toOp = readOrder[j][baseIndex + stagePart[bailiff]  + nextStageI]-1;
						//cout << "from: " << fromOp << " ; toOp: " << toOp << endl;
						stageInst.addArrow(j, Arrow(fromOp, toOp));
					}
				}
				baseIndex += stagePart[bailiff];
			}
  
			//getchar();
		}
		
		return stageInst;
	}


	//WRONG !! REDO AND DELETE
	static InstGen fromRawPsspToMixed(const string & filePath, unsigned n0Size) {
		
		InstGen ig;

		ifstream fromFile;

		fromFile.open(filePath);
		if( ! fromFile.is_open())
			throw errorText("Could not open partitionedraw nonTaillard's instance file: "+filePath, "Parser.hpp", "fromRawPsspToMixed");

		fromFile >> ig.nJobs;
		ig.nJobs+=n0Size;
		fromFile >> ig.nMachs;

		vector<unsigned> order;
		unsigned aOrder;
		vector<unsigned> costs;
		unsigned aCost;

		ig.procTs.resize(ig.nJobs, vector<unsigned>(ig.nMachs));
		ig.jobCGs.resize(ig.nJobs, ClosedGraph<MAX_ELEMENTS>(ig.nMachs));
		ig.jobSets.resize(ig.nJobs);

		//Setting NJ
		for(unsigned j=0; j<ig.nJobs-n0Size; j++) {
			order.clear();
			costs.clear();
			for(unsigned m=0; m<ig.nMachs; m++) {
				fromFile >> aOrder;
				fromFile >> aCost;
				order.push_back(aOrder);
				costs.push_back(aCost);
			}

			for(unsigned index=0; index<ig.nMachs; index++) {
				ig.procTs[j][order[index]] = costs[index];
			}

			for(unsigned index=0; index<ig.nMachs-1; index++) {
				ig.addArrow(j, Arrow(order[index], order[index+1]));
			}
		}

		//Setting NO
		for(unsigned j=ig.nJobs-n0Size; j<ig.nJobs; j++) {
			for(unsigned index=0; index<ig.nMachs; index++) {
				ig.procTs[j][index] = rand()%98+1;//CONFIRM LIMITS
			}
		}

		fromFile.close();

		return ig;
	}
	
	

	unsigned nJobs;
	unsigned nMachs;

	vector<vector<unsigned>> procTs;
	
	vector<ClosedGraph<MAX_ELEMENTS>> jobCGs; 
	vector<set<Arrow>> jobSets;
};
//########################################################################################################################################################################################
