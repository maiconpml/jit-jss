/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#pragma once


#include "../Settings.hpp"
#include "../Utilities.hpp"

#include <boost/multi_array.hpp>


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



class PsiInst {
public:
	PsiInst() {  }
	~PsiInst() {  }



	static vector<PsiInst> readOspTaiFormat(const string & filePath) {
		vector<PsiInst> v;
		PsiInst pi;

	    unsigned J;
		unsigned M;
		multi_array<unsigned, 2> readOrder;
		multi_array<unsigned, 2> readProcs;
		string bufStr = "T";

		ifstream fromFile;

		fromFile.open(filePath);
		if( ! fromFile.is_open())
			throw errorText("Could not open osp file: "+filePath, "InstGen.hpp", "readOspTaiFormat");

		for(unsigned rep=0; rep<10; rep++) {

			while( ! bufStr.compare(":")) { assert( ! fromFile.eof());   fromFile>>bufStr; }//fast forward to :

			fromFile >> J;
			fromFile >> M;

			assert(J < 1000);
			assert(M < 1000);

			readOrder.resize(extents[J][M]);
			readProcs.resize(extents[J][M]);

			while( ! bufStr.compare(":")) { assert( ! fromFile.eof());   fromFile>>bufStr; }//fast forward to :

			for(unsigned j=0; j<J; j++) {
				for(unsigned m=0; m<M; m++) {
					fromFile >> readProcs[j][m];
					cout << readProcs[j][m] << endl;
					getchar();
					assert(readProcs[j][m] < 100);
				} // for m : M
			} // for j : J

			while( ! bufStr.compare(":")) { assert( ! fromFile.eof());   fromFile>>bufStr; }//fast forward to :

			for(unsigned j=0; j<J; j++) {
				for(unsigned m=0; m<M; m++) {
					fromFile >> readOrder[j][m];
					readOrder[j][m]--;
					assert(readOrder[j][m] < M);
				} // for m : M
			} // for j : J

			pi.procTs.resize(extents[J][M]);
			pi.J = J;
			pi.M = M;
			pi.jobSets.resize(J);

			for(unsigned j=0; j<J; j++) {
				for(unsigned m=0; m<M; m++) {
					pi.procTs[j][readOrder[j][m]] = readProcs[j][m];
				}  // for m : M
			} // for j : J

			pi.print("");
			getchar();

		} //end rep

		return v;
	}


	static void jspLawrenceParallelPsi(const string & inPath, const string & outDir, const vector<unsigned> & kV){
		ifstream fromFile;
		ofstream toFile;

		string bufString;
		string instName;

		PsiInst pi;
		multi_array<unsigned, 2> readOrder;
		multi_array<unsigned, 2> readProcs;

		fromFile.open(inPath);
		if( ! fromFile.is_open())
			throw errorText("Could not open jsp file: "+inPath, "InstGen.hpp", "sepLawrence");


		for(unsigned instIndex=1; instIndex<=15; instIndex++) {
			instName.clear();

			fromFile >> bufString;
			assert(bufString.compare("instance")==0);
			fromFile >> instName;

			fromFile >> pi.J;
			//assert(pi.J>0   &&   pi.J <= MAX_JOBS);
			fromFile >> pi.M;
			//assert(pi.M>0   &&   pi.M <= MAX_MACHS);

			pi.procTs.resize(extents[pi.J][pi.M]);
			pi.jobSets.clear();
			pi.jobSets.resize(pi.J);
			readOrder.resize(extents[pi.J][pi.M]);
			readProcs.resize(extents[pi.J][pi.M]);

			for(unsigned j=0; j<pi.J; j++) {
				for(unsigned m=0; m<pi.M; m++) {
					fromFile >> readOrder[j][m];
					assert(readOrder[j][m] < pi.M);
					fromFile >> readProcs[j][m];
				}
			}

			for(unsigned j=0; j<pi.J; j++) {
				for(unsigned bailiff=0; bailiff<pi.M-1; bailiff++) {
					for(unsigned provost=bailiff+1; provost<pi.M; provost++) {
						assert(readOrder[j][bailiff] != readOrder[j][provost]);
					}
				}
			}

			for(unsigned j=0; j<pi.J; j++) {
				for(unsigned m=0; m<pi.M; m++) {
					pi.procTs[j][readOrder[j][m]] = readProcs[j][m];
				}
			}

			for(unsigned j=0; j<pi.J; j++) {
				for(unsigned u=0; u<pi.M-1; u++) {
					pi.jobSets[j].insert(Arrow(readOrder[j][u], readOrder[j][u+1]));
				}
			}

			//pi.print(outDir+instName+".ppsi");
			for(unsigned k : kV)
				pi.printParallel(outDir+instName+"_k"+unsigStr(k)+".ppsi", k);
			//getchar();
		}
		


		fromFile.close();
	}



	static void jspLawrencePsi(const string & inPath, const string & outDir) {
		ifstream fromFile;
		ofstream toFile;

		string bufString;
		string instName;

		PsiInst pi;
		multi_array<unsigned, 2> readOrder;
		multi_array<unsigned, 2> readProcs;

		fromFile.open(inPath);
		if( ! fromFile.is_open())
			throw errorText("Could not open jsp file: "+inPath, "InstGen.hpp", "sepLawrence");


		for(unsigned instIndex=1; instIndex<=15; instIndex++) {
			instName.clear();

			fromFile >> bufString;
			assert(bufString.compare("instance")==0);
			fromFile >> instName;

			fromFile >> pi.J;
			//assert(pi.J>0   &&   pi.J <= MAX_JOBS);
			fromFile >> pi.M;
			//assert(pi.M>0   &&   pi.M <= MAX_MACHS);

			pi.procTs.resize(extents[pi.J][pi.M]);
			pi.jobSets.clear();
			pi.jobSets.resize(pi.J);
			readOrder.resize(extents[pi.J][pi.M]);
			readProcs.resize(extents[pi.J][pi.M]);

			for(unsigned j=0; j<pi.J; j++) {
				for(unsigned m=0; m<pi.M; m++) {
					fromFile >> readOrder[j][m];
					assert(readOrder[j][m] < pi.M);
					fromFile >> readProcs[j][m];
				}
			}

			for(unsigned j=0; j<pi.J; j++) {
				for(unsigned bailiff=0; bailiff<pi.M-1; bailiff++) {
					for(unsigned provost=bailiff+1; provost<pi.M; provost++) {
						assert(readOrder[j][bailiff] != readOrder[j][provost]);
					}
				}
			}

			for(unsigned j=0; j<pi.J; j++) {
				for(unsigned m=0; m<pi.M; m++) {
					pi.procTs[j][readOrder[j][m]] = readProcs[j][m];
				}
			}

			for(unsigned j=0; j<pi.J; j++) {
				for(unsigned u=0; u<pi.M-1; u++) {
					pi.jobSets[j].insert(Arrow(readOrder[j][u], readOrder[j][u+1]));
				}
			}

			pi.print(outDir+instName+".psi");
			//getchar();
		}


		fromFile.close();
	}


	//from JSP cannonic format (Ritt usual format) to JSP
	void readCannonJsp(const string & jspPath) {
		ifstream fromFile;
		multi_array<unsigned, 2> readOrder;
		multi_array<unsigned, 2> readProcs;

		fromFile.open(jspPath);
		if( ! fromFile.is_open())
			throw errorText("Could not open jsp file: "+jspPath, "InstGen.hpp", "readCannonJsp");

		fromFile >> J;
		//assert(J>0   &&   J <= MAX_JOBS);
		fromFile >> M;
		//assert(M>0   &&   M <= MAX_MACHS);

		procTs.resize(extents[J][M]);
		jobSets.resize(J);
		readOrder.resize(extents[J][M]);
		readProcs.resize(extents[J][M]);

		for(unsigned j=0; j<J; j++) {
			for(unsigned m=0; m<M; m++) {
				fromFile >> readOrder[j][m];
				assert(readOrder[j][m] < M);
				fromFile >> readProcs[j][m];
			}
		}

		for(unsigned j=0; j<J; j++) {
			for(unsigned bailiff=0; bailiff<M-1; bailiff++) {
				for(unsigned provost=bailiff+1; provost<M; provost++) {
					assert(readOrder[j][bailiff] != readOrder[j][provost]);
				}
			}
		}

		for(unsigned j=0; j<J; j++) {
			for(unsigned m=0; m<M; m++) {
				procTs[j][readOrder[j][m]] = readProcs[j][m];
			}
		}

		for(unsigned j=0; j<J; j++) {
			for(unsigned u=0; u<M-1; u++) {
				jobSets[j].insert(Arrow(readOrder[j][u], readOrder[j][u+1]));
			}
		}
	}




	//replicated jobs and machines - costs and J replicated on jobs - machines is implcit in k (kM machines total) 
	void printParallel(const string & toPath, unsigned k) const {

		assert(J > 0);
		assert(M > 0);
		assert(jobSets.size() == J);
		assert(procTs.shape()[0] == J);
		assert(procTs.shape()[1] == M);

		ofstream toFile;
		string str ="#START#\n\nTotalJobs: ";		
		str.append(unsigStr(k*J));
		str.append(" TotalMachines: ");
		str.append(unsigStr(M));
		str.append(" Replications: ");
		str.append(unsigStr(k));

		str.append("\nCosts:");
		for(unsigned j=0; j<J; j++) {
			for(unsigned rep=0; rep<k; rep++) {
				str.append("\n\t");
				for(unsigned m=0; m<M; m++) {
					str.append(unsigStr(procTs[j][m])+" ");
				}
			}
		}
		str.append("\n");
		
		for(const set<Arrow> & sa : jobSets) {
			for(unsigned rep=0; rep<k; rep++) {
				str.append("\nJob: "+unsigStr(sa.size())+"\n");
				for(const Arrow & a : sa) {
					assert(a.from < M);
					assert(a.to < M);
					str.append("\n\t"+unsigStr(a.from)+" -> "+unsigStr(a.to));
				}
				str.append("\n");
			}
		}

		if( ! toPath.empty()) {
			toFile.open(toPath);
			if( ! toFile.is_open())
				throw errorText("Could not open file to write psp instance: "+toPath, "Parser.hpp", "Parser::printPsp()");
			toFile << str; 
			toFile.close();
		} else {
			cout << str << endl;
		}
	}



	//empty toPath for cou instead of file
	void print(const string & toPath) const {
		assert(J > 0);
		assert(M > 0);
		assert(jobSets.size() == J);
		assert(procTs.shape()[0] == J);
		assert(procTs.shape()[1] == M);
		
		ofstream toFile;
		string str ="#START#\n\nTotalJobs: ";
		
		str.append(unsigStr(procTs.size()));
		str.append(" TotalMachines: ");
		str.append(unsigStr(procTs[0].size()));
		
		str.append("\nCosts:");
		for(unsigned j=0; j<J; j++) {
			str.append("\n\t");
			for(unsigned m=0; m<M; m++) {
				str.append(unsigStr(procTs[j][m])+" ");
			}
		}
		str.append("\n");
		
		for(const set<Arrow> & sa : jobSets) {
			str.append("\nJob: "+unsigStr(sa.size())+"\n");
			for(const Arrow & a : sa) {
				assert(a.from < M);
				assert(a.to < M);
				str.append("\n\t"+unsigStr(a.from)+" -> "+unsigStr(a.to));
			}
			str.append("\n");
		}
		
		if( ! toPath.empty()) {
			toFile.open(toPath);
			if( ! toFile.is_open())
				throw errorText("Could not open file to write psp instance: "+toPath, "Parser.hpp", "Parser::printPsp()");
			toFile << str; 
			toFile.close();
		} else {
			cout << str << endl;
		}
	}
	
	
	
	//creates gss instance such that any solution of it is also solution for PSS
	//reads <fromPath>.psi and prints in <toPath>.psi
	static void mapPssToGSS(const string & fromPath, const string & toPath) {
		PsiInst pssPi;
		PsiInst gssPi;
		
		pssPi.readPsiInst(fromPath);
		assert(pssPi.J > 0   &&   pssPi.M > 0);
		
		//setting basic data for gssPsi
		gssPi.J = pssPi.J;
		gssPi.M = pssPi.M;
		//gssPi.procTs = pssPi.procTs;
		vector<size_t> allShapes;
    	const size_t* shape = pssPi.procTs.shape();
    	allShapes.assign( shape, shape+pssPi.procTs.num_dimensions() );
    	gssPi.procTs.resize(allShapes);
    	gssPi.procTs = pssPi.procTs;	
		gssPi.jobSets.resize(pssPi.J);
	
		
		for(unsigned j = 0; j < pssPi.J; j++) {
			fromPartialToGssOrder(pssPi.M, pssPi.jobSets[j], gssPi.jobSets[j]);
		}
		gssPi.print(toPath);
	}
	
	
	
	static void fromPartialToGssOrder(const unsigned M, const set<Arrow> & jobSetPartial, set<Arrow> & jobSetTotal) {
		
		vector<unsigned> indeg(M, 0);
		vector<unsigned> closedSet(M, false);//closedSet[o]==true ? o is already in a stage
		unsigned numClosedSet = 0;
		
		vector<unsigned> prevStage;
		vector<unsigned> thisStage;
		
		for(const Arrow & a : jobSetPartial)
			indeg[a.to]++;
			
		while(numClosedSet < M) {
		
			//find thisStage
			thisStage.clear();
			for(unsigned o = 0; o < M; o++) {
				if(indeg[o] == 0   &&   ! closedSet[o]) {
					thisStage.push_back(o);
					closedSet[o] = true;
				}
			}

			
			//link previous stage with this
			assert( (prevStage.empty()   &&   numClosedSet==0)   || //either is first stage and no link is possible
					(! prevStage.empty()   &&   numClosedSet>0));	//or is not first stage and there must be a previous
			for(unsigned from : prevStage) {
				for(unsigned to : thisStage) {
					jobSetTotal.insert(Arrow(from, to));
				}
			}			
					
			//update indeg
			for(const Arrow & a : jobSetPartial) {
				for(unsigned o : thisStage) {
					if(a.from == o) {
						assert( ! closedSet[a.to]);
						assert(indeg[a.to] > 0);
						indeg[a.to]--;
					}
				}
			}
			
			//prep for next loop
			numClosedSet += thisStage.size();
			prevStage = thisStage;
		}
	}
	
	
	
	//creates mss instance such that any solution of it is also solution for PSS
	//reads <fromPath>.psi and prints in <toPath>.psi
	static void mapPssToMSS(const string & fromPath, const string & toPath) {
		PsiInst pssPi;
		PsiInst mssPi;
		
		pssPi.readPsiInst(fromPath);
		assert(pssPi.J > 0   &&   pssPi.M > 0);
		
		//setting basic data for mssPsi
		mssPi.J = pssPi.J;
		mssPi.M = pssPi.M;
		//mssPi.procTs = pssPi.procTs;
		vector<size_t> allShapes;
    	const size_t* shape = pssPi.procTs.shape();
    	allShapes.assign( shape, shape+pssPi.procTs.num_dimensions() );
    	mssPi.procTs.resize(allShapes);
    	mssPi.procTs = pssPi.procTs;	
		mssPi.jobSets.resize(pssPi.J);
	
		
		for(unsigned j = 0; j < pssPi.J; j++) {
			if( ! pssPi.jobSets[j].empty()) {
				//job shop job
				fromPartialToTotalOrder(pssPi.M, pssPi.jobSets[j], mssPi.jobSets[j]);
			}
			//open shop job need nothing 
		}
		//pssPi.print("");
		//mssPi.print("");
		mssPi.print(toPath);
	}
	
	
	
	//from the partial order in <jobSetsPartial> get random associated total order in jobSetsTotal
	static void fromPartialToTotalOrder(const unsigned M, const set<Arrow> & jobSetPartial, set<Arrow> & jobSetTotal) {
	
		assert( ! jobSetPartial.empty());
		
		vector<unsigned> indeg(M, 0);
		vector<unsigned> openSet;
		
		unsigned prevOp = UINT_MAX;
		unsigned thisOp;
		
		for(const Arrow & a : jobSetPartial)
			indeg[a.to]++;
		
		for(unsigned o = 0; o < M; o++) {
			if(indeg[o] == 0)
				openSet.push_back(o);
		}
		
		assert( ! openSet.empty());
		
		while( ! openSet.empty()) {
			assert(jobSetTotal.size() < M-1);
		
			random_shuffle(openSet.begin(), openSet.end());
			
			thisOp = openSet.back();
			openSet.pop_back();
			
			if(prevOp != UINT_MAX) {//there was an actual prev
				jobSetTotal.insert(Arrow(prevOp, thisOp));
			}
			
			//update indeg
			for(const Arrow & a : jobSetPartial) {
				if(a.from == thisOp) {
					assert(indeg[a.to] > 0);
					indeg[a.to]--;
					if(indeg[a.to] == 0) {
						openSet.push_back(a.to);
					}
				}
			}
			
			prevOp = thisOp;
		}
		assert(jobSetTotal.size() == M-1);
	}
	


	void readPsiInst(const string & filePath) {
		
		string bufferStr;
		unsigned bufferT;
		unsigned auxJSize;
		unsigned fstM, sndM;

		ifstream streamFromFile;
		streamFromFile.open(filePath);
		if( ! streamFromFile.is_open())
			throw errorText("Could not open instance file: "+filePath, "Instance.hpp", "Inst::parse()");

		while(true) {
			streamFromFile >> bufferStr;
			
			if(streamFromFile.eof())
				throw errorText("Expected #START# but reached eof", "Instance.hpp", "Inst::parse()");
			
			if(bufferStr.compare("#START#") == 0)
				break;
		}

		streamFromFile >> bufferStr;
		assert(bufferStr.compare("TotalJobs:") == 0);
		streamFromFile >> J;
		
		streamFromFile >> bufferStr;
		assert(bufferStr.compare("TotalMachines:") == 0);
		streamFromFile >> M;

		procTs.resize(extents[J][M]);

		//Reading proc Ts
		streamFromFile >> bufferStr;
		assert(bufferStr.compare("Costs:") == 0);
		for(unsigned j=0; j<J; j++) {
			for(unsigned m=0; m<M; m++) {
				streamFromFile >> bufferT;
				assert(bufferT != 0);
				procTs[j][m] = bufferT;
			}
		}

		jobSets.resize(J);

		//Reading precedences
		for(unsigned j=0; j<J; j++) {
			streamFromFile >> bufferStr;
			assert(bufferStr.compare("Job:") == 0);

			streamFromFile >> auxJSize;

			for(unsigned u=0; u<auxJSize; u++) {
				streamFromFile >> fstM;
				assert(fstM < M);

				streamFromFile >> bufferStr;
				assert(bufferStr.compare("->") == 0);
				
				streamFromFile >> sndM;
				assert(sndM < M);

				jobSets[j].insert(Arrow(fstM, sndM));
			}
		}
	}


	unsigned J;
	unsigned M;
    multi_array<unsigned, 2>  procTs;
	vector<set<Arrow>> jobSets;
};
