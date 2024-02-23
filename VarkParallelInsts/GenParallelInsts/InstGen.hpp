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
	Arrow(unsigned _from, unsigned _to) : from(_from), to(_to) { }
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




class InstGen {
public:
	InstGen(unsigned _J, unsigned _M) : J(_J), M(_M) {
		varK.resize(M);
		procTs.resize(extents[J][M]);
		jobOrder.resize(J, vector<unsigned>(M));
	}
	InstGen() {}
	~InstGen() {}
	
	
	
	void alloc(unsigned _J, unsigned _M) {
		J = _J;
		M = _M;
		varK.resize(M);
		procTs.resize(extents[J][M]);
		jobOrder.resize(J, vector<unsigned>(M));
	}
	
	
	
	bool isAlloced() const {
		if(J == 0  ||   J>1000) return false;
		if(M == 0  ||   M>1000) return false;
		if(varK.size() != M) return false;
		if(procTs.shape()[0] != J  ||   procTs.shape()[1] != M) return false;
		if(jobOrder.size() != J) return false;
		for(const vector<unsigned> & v : jobOrder)
			if(v.size() != M)
				return false;
		return true;
	}
	
	
	
	//prints in vark pjss format.
	//If fileName is empty uses cout instead
	void print(const string & fileName) const {
		assert(isAlloced());
		
		ofstream toFile;
		string str ="#START#\n\nTotalJobs: ";		
		str.append(unsigStr(J));
		str.append(" TotalMachines: ");
		str.append(unsigStr(M));
		str.append(" Replications: ");
		for(unsigned k : varK)
			str.append(unsigStr(k) + " ");

		str.append("\nCosts:");
		for(unsigned j=0; j<J; j++) {
			str.append("\n\t");
			for(unsigned m=0; m<M; m++) {
				assert(procTs[j][m] > 0   &&   procTs[j][m] <= 1000);
				str.append(unsigStr(procTs[j][m])+" ");
			}
		}
		str.append("\n");
		
		for(unsigned j=0; j<J; j++) {
			str.append("\nJob: "+unsigStr(M-1)+"\n");
			for(unsigned x=0; x<M-1; x++) {
				assert(jobOrder[j][x] < M);
				assert(jobOrder[j][x+1] < M);
				str.append("\n\t"+unsigStr(jobOrder[j][x])+" -> "+unsigStr(jobOrder[j][x+1]));
			}
			str.append("\n");
		}

		if( ! fileName.empty()) {
			toFile.open(fileName);
			if( ! toFile.is_open())
				throw errorText("Could not open file to write psp instance: "+fileName, "Parser.hpp", "Parser::printPsp()");
			toFile << str; 
			toFile.close();
		} else {
			cout << str << endl;
		}
	}
	
	
	
	
	//randomizes varK, procTs and jobOrder
	//expects J and M to be set
	void randomize(unsigned maxK, unsigned maxP) {
		assert(isAlloced());
		
		vector<unsigned> bufV(M);
		
		//varK
		for(unsigned m=0; m<M; m++) {
			varK[m] = 1 + (rand()%maxK);
			bufV[m] = m;//for jobOrder later
		}
		
		//procTs
		for(unsigned j=0; j<J; j++) {
			for(unsigned m=0; m<M; m++) {
				procTs[j][m] = 1 + (rand()%maxP);
			}
		}
		
		//jobOrder
		for(unsigned j=0; j<J; j++) {
			random_shuffle(bufV.begin(), bufV.end());
			jobOrder[j] = bufV;
		}
	}
	
	
	
	//generates random instace with j jobs and m machines
	//@maxK: possible number of replications of a machine [1,maxK]
	//@maxP: possible processing time for an operation is [1,maxP]
	//Prints on file given by fileName - uses cout if empty
	static void genRandInst(const string & fileName, unsigned j ,unsigned m, unsigned maxK, unsigned maxP) {
		InstGen ig(j ,m);
		
		ig.randomize(maxK, maxP);
		
		ig.print(fileName);
	}
	
	
	
	
	unsigned J;
	unsigned M;
	vector<unsigned> varK;
	multi_array<unsigned, 2>  procTs; //procTs[j][m]: processig time of oper in job j and machine m 
	vector<vector<unsigned>>  jobOrder; //jobOrder[j][x]: oper in job j, in machine jobOrder[j][x], is in position x
};



