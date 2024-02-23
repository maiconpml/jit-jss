/*
  @author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/
#pragma once

#include <boost/multi_array.hpp>
#include <climits>
#include <algorithm> 
#include <chrono>

#include "Settings.hpp"
#include "Utilities.hpp"
#include "Instance.hpp"
#include "RandomNumberGenerator.hpp"
#include "State.hpp"

using namespace boost;




#ifdef ANALYZE_MODE


class OperPairQuality {
public:
	OperPairQuality() {}
	OperPairQuality(unsigned _bailiff, unsigned _provost, unsigned _makes) : bailiff(_bailiff), provost(_provost), makes(_makes) {}
	~OperPairQuality() {}
	
	bool operator<(const OperPairQuality & oq) const {
		return this->makes  <  oq.makes;
	}
	
	string toString() const {
		string str;
		
		str.append(unsigStr(bailiff) + " "+unsigStr(provost) + " "+unsigStr(makes) + " ");
		
		return str;
	}
	
	unsigned bailiff;
	unsigned provost;
	unsigned makes;
};

//predict higher == better chance to improve
class OperPairPredict {
public:
	OperPairPredict() {}
	OperPairPredict(unsigned _bailiff, unsigned _provost, double _predict) : bailiff(_bailiff), provost(_provost), predict(_predict){}
	~OperPairPredict() {}
	
	bool operator<(const OperPairPredict & op) const {
		return this->predict  >  op.predict;
	}
	
		string toString() const {
		string str;
		
		str.append(unsigStr(bailiff) + " "+unsigStr(provost) + " "+doubleStr(predict) + " ");
		
		return str;
	}
	
	unsigned bailiff;
	unsigned provost;
	double predict;//higher == better chance to improve
};



class OperData {
public:
	OperData() : inNeigh(0), inNeighAsBailiff(0), inNeighAsProvost(0), bestInNSP(0), bestInNSPAsBailiff(0), bestInNSPAsProvost(0), improvedCurrent(0), improvedCurrentAsBailiff(0), improvedCurrentAsProvost(0) {}
	~OperData() {}
	

	string toString() const {
		
		string str = "      ";
		
		str.append(unsigStr(inNeigh) + "     ");
		str.append(unsigStr(inNeighAsBailiff) + "     ");
		str.append(unsigStr(inNeighAsProvost) + "     ");
		str.append(unsigStr(bestInNSP) + "     ");
		str.append(unsigStr(bestInNSPAsBailiff) + "     ");
		str.append(unsigStr(bestInNSPAsProvost) + "     ");
		str.append(unsigStr(improvedCurrent) + "     ");
		str.append(unsigStr(improvedCurrentAsBailiff) + "     ");
		str.append(unsigStr(improvedCurrentAsProvost) + "     ");
		
		return str;
	}
	
	
	unsigned inNeigh;
	unsigned inNeighAsBailiff;
	unsigned inNeighAsProvost;
	unsigned bestInNSP;
	unsigned bestInNSPAsBailiff;
	unsigned bestInNSPAsProvost;
	unsigned improvedCurrent;
	unsigned improvedCurrentAsBailiff;
	unsigned improvedCurrentAsProvost;
	
	static string operDataHeader;
};
string OperData::operDataHeader = "Oper inNeigh inNeighAsBailiff inNeighAsProvost bestInNSP bestInNSPAsBailiff bestInNSPAsProvost improvedCurrentAsBailiff inNeighAsProvost";










class Analyzer {
public:
	Analyzer()  {}
	~Analyzer() {}
	
	
	
	void set() {
		alloc();
		setOpersByJobTail();
	}
	
	
	
	void alloc() {
		allOperData.resize(inst.O);
		opersByJobTail.resize(inst.O);
	}
	
	
	//prioActual[x]: operation x was in order orderActual[x] of quality
	//prioPredict same
	//UINFTY if not present in critical path
	// rank is from 1 - n --- no repeats
	static double simpleKendallW(const vector<unsigned> & prioActual, const vector<unsigned> & prioPredict) {
		assert(prioActual.size() == inst.O);
		assert(prioPredict.size() == inst.O);
		double w;	
		vector<unsigned> R(inst.O, 0);
		double meanR = 0.0;
		double S = 0.0;
		unsigned n = 0;
		
		//computing Ri
		for(unsigned o=1; o<inst.O; o++) {
			if(prioActual[o] < UINFTY) {
				assert(prioPredict[o] < UINFTY);
				R[o] = prioActual[o] + prioPredict[o];
				n++;
				meanR += R[o];
			}
		}
		
		//meanR
		assert(n > 0);
		meanR = meanR/(double)(n);

		
		for(unsigned o=1; o<inst.O; o++) {
			if(R[o] > 0)
				S += pow((double)(R[o]) - meanR, 2.0);
		}
		
		w = (3 * S)/(double)(pow(n, 3.0) - n);		

		assert(w>=0    &&    w<=1.0);
		
		return w;
	} 
	
	
	//prioActual[bailiff][provost]
	static double pairKendallW(const vector<vector<unsigned>> & prioActual, const vector<vector<unsigned>> & prioPredict,  vector<vector<unsigned>> & R, const vector<OperPairQuality> & qualityV) {
		assert(R.size() == inst.O);
		double w;	
		double meanR = 0.0;
		double S = 0.0;
		unsigned n = qualityV.size();
		
		if(n==1)
			return 1.0;
		
		
		for(const OperPairQuality & opq : qualityV) {
			assert(prioActual[opq.bailiff][opq.provost] > 0);
			assert(prioActual[opq.bailiff][opq.provost] <=  n);
			assert(prioPredict[opq.bailiff][opq.provost] > 0);
			assert(prioPredict[opq.bailiff][opq.provost] <=  n);
			R[opq.bailiff][opq.provost] = prioActual[opq.bailiff][opq.provost] + prioPredict[opq.bailiff][opq.provost];
			meanR += R[opq.bailiff][opq.provost];
		}
				
		//meanR
		assert(n > 0);
		meanR = meanR/(double)(n);

		
		for(const OperPairQuality & opq : qualityV) {
			assert(R[opq.bailiff][opq.provost] >= 2);
			assert(R[opq.bailiff][opq.provost] <=  2*n);
			S += pow((double)(R[opq.bailiff][opq.provost]) - meanR, 2.0);
		}
		
		w = (3 * S)/(double)(pow(n, 3.0) - n);
		
		//cout << n << endl;
		//cout << "\nw: " << w << endl;  

		assert(w>=0    &&    w<=1.0);
		
		return w;
	}
	
	
	
	void setOpersByJobTail() {
#ifndef NDEBUG
		assert(inst.O > 0    &&     inst.O<100000);
		vector<bool> closedSet(inst.O, false);
#endif //NDEBUG
		
		vector<pair<unsigned, unsigned>> v; 
		vector<unsigned> operJobTail(inst.O);
		unsigned jobSum;
		unsigned curOp;
		
		//setting job tails
		for(unsigned aTail : inst.leafs) {
#ifndef NDEBUG
			assert( ! closedSet[aTail]);
			closedSet[aTail] = true;
#endif //NDEBUG
			jobSum = 0;
			curOp = aTail;

			for(unsigned u=0; u<inst.M; u++) {
				assert(curOp>0    &&    curOp<inst.O);
				operJobTail[curOp] = jobSum;
				jobSum += inst.P[curOp];
				
				if(u!=inst.M-1) {
					assert(inst.prev[curOp].size() == 1);
					curOp = inst.prev[curOp][0];
				}
			}
		}
		
		for(unsigned o=1; o<inst.O; o++)
			v.push_back(pair<unsigned, unsigned>(operJobTail[o], o));
		
		sort(v.begin(), v.end());
		
		//for(unsigned )
	}
	
	
	//quality v is only useeful to find cands - not using cands because NSP deletes the selected one
	//the predition strength is % of improve best of bailif  as bailiff and of provost as provost
	void predictInPosImproveBest(vector<vector<unsigned>> & prioPredict, const vector<pair<unsigned, unsigned>> & cands) const {
#ifndef NDEBUG
		assert( ! cands.empty());
		assert(prioPredict.size() == inst.O);
		for(const vector<unsigned> & v : prioPredict)
			assert(v.size() == inst.O);
#endif //NDEBUG
		vector<OperPairPredict> predictV;
		double bailiffPredict;
		double provostPredict;
		unsigned bailiff;
		unsigned provost;
		
		
		for(const pair<unsigned, unsigned> & aCand : cands) {
			bailiff = aCand.first;
			provost = aCand.second;
			if(allOperData[bailiff].inNeighAsBailiff == 0) {
				bailiffPredict = 0.0;
			} else {
				bailiffPredict = (allOperData[bailiff].improvedCurrentAsBailiff)/(double)(allOperData[bailiff].inNeighAsBailiff);
			}
			
			if(allOperData[provost].inNeighAsProvost == 0) {
				provostPredict = 0.0;
			} else {
				provostPredict = (allOperData[provost].improvedCurrentAsProvost)/(double)(allOperData[provost].inNeighAsProvost);
			}
		
			predictV.push_back(OperPairPredict(bailiff, provost, bailiffPredict+provostPredict));
		}
		
		sort(predictV.begin(), predictV.end());
		
		//for(const OperPairPredict & opp : predictV)
		//	cout << opp.toString() << endl;
		//getchar();
		
		for(unsigned u=0; u<predictV.size(); u++) {
			prioPredict[predictV[u].bailiff][predictV[u].provost] = u+1;
		}
		

	}
	
	
	
	
	//the predition strength is % of improve best of bailif  as bailiff and of provost as provost
	void predictInPosBestNSP(vector<vector<unsigned>> & prioPredict, const vector<pair<unsigned, unsigned>> & cands) const {
#ifndef NDEBUG
		assert( ! cands.empty());
		assert(prioPredict.size() == inst.O);
		for(const vector<unsigned> & v : prioPredict)
			assert(v.size() == inst.O);
#endif //NDEBUG
		vector<OperPairPredict> predictV;
		double bailiffPredict;
		double provostPredict;
		unsigned bailiff;
		unsigned provost;
		
		
		for(const pair<unsigned, unsigned> & aCand : cands) {
			bailiff = aCand.first;
			provost = aCand.second;
			if(allOperData[bailiff].inNeighAsBailiff == 0) {
				bailiffPredict = 0.0;
			} else {
				bailiffPredict = (allOperData[bailiff].bestInNSPAsBailiff)/(double)(allOperData[bailiff].inNeighAsBailiff);
			}
			
			if(allOperData[provost].inNeighAsProvost == 0) {
				provostPredict = 0.0;
			} else {
				provostPredict = (allOperData[provost].bestInNSPAsProvost)/(double)(allOperData[provost].inNeighAsProvost);
			}
		
			predictV.push_back(OperPairPredict(bailiff, provost, bailiffPredict+provostPredict));
		}
		
		sort(predictV.begin(), predictV.end());
		
		for(unsigned u=0; u<predictV.size(); u++) {
			prioPredict[predictV[u].bailiff][predictV[u].provost] = u+1;
		}
		

	}
	
	
	
	void sortCandsByPredict(vector<pair<unsigned,unsigned>> & cands) const {
		assert( ! cands.empty());
		vector<OperPairPredict> predictV;
		double bailiffPredict;
		double provostPredict;
		unsigned bailiff;
		unsigned provost;
		
		
		for(const pair<unsigned, unsigned> & aCand : cands) {
			bailiff = aCand.first;
			provost = aCand.second;
			if(allOperData[bailiff].inNeighAsBailiff == 0) {
				bailiffPredict = 0.0;
			} else {
				bailiffPredict = (allOperData[bailiff].bestInNSPAsBailiff)/(double)(allOperData[bailiff].inNeighAsBailiff);
			}
			
			if(allOperData[provost].inNeighAsProvost == 0) {
				provostPredict = 0.0;
			} else {
				provostPredict = (allOperData[provost].bestInNSPAsProvost)/(double)(allOperData[provost].inNeighAsProvost);
			}
		
			predictV.push_back(OperPairPredict(bailiff, provost, bailiffPredict+provostPredict));
		}
		
		sort(predictV.begin(), predictV.end());
		
		cands.clear();
		for(const OperPairPredict & opp : predictV) {
			//cout << "\t" << opp.toString() << endl;
			cands.push_back(pair<unsigned, unsigned>(opp.bailiff, opp.provost));
		}
	}
	
	
	
	
	
	string toString() const {
		string str ="";
	
		str.append(OperData::operDataHeader );
		for(unsigned u=1; u<inst.O; u++) {
			str.append("\n" + unsigStr(u) + " " + allOperData[u].toString());
		}
	
		return str;
	}
	
	
	//static stuff
	vector<unsigned> opersByJobTail;//all operations ordered by the INCREASING size of its static job tail 
	vector<OperData> allOperData;
};


class Observer {
public:
	Observer() : sumKW(0), kwTimes(0), setMetaTimes(0) { }
	~Observer() { }
	
	
	void set(unsigned _minLearnPeriod, double _minLearnKW, double _trimmSize) {
		analyzer.set();
		minLearnPeriod = _minLearnPeriod;
		minLearnKW = _minLearnKW;
		trimmSize = _trimmSize;
	}
	
	//returns is you should start trimming
	bool addKendallW(double aKendalW) {
		sumKW += aKendalW;
		kwTimes++;
		//buffKW = sumKW/(double)(kwTimes);
		//if(buffKW > maxKW)
			//maxKW = buffKW;
		
		//XXX	
		//cout << "\t\tkwTimes: " << kwTimes << "   avrgKW: " << (sumKW/(double)(kwTimes)) << endl;
		//cout << "\t\tminLearnPeriod: " << minLearnPeriod <<"  minLearnKW:" << minLearnKW << endl;
		return minLearnPeriod<kwTimes   &&   minLearnKW<(sumKW/(double)(kwTimes));
	}
	
	
	
	
	void forget() {
		sumKW=0.0;
		kwTimes = 0;
	}
	
	
	
	unsigned minLearnPeriod;
	double minLearnKW;
	double trimmSize;

	Analyzer analyzer;
	
	//kendall 2 fields
	double sumKW;  //sum of kenda w average kendal w is sumKW/(double)kwTimes
	unsigned kwTimes; //how many times kendall w was computed
	//double buffKW;
	//double maxKW;
	
	unsigned setMetaTimes; //times setMeta was called
};
Observer observer;

#endif //ANALYZE_MODE

