/*
  @author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/
#pragma once

#include "Settings.hpp"
#include "Utilities.hpp"
#include "Instance.hpp"
//#include "Analyzer.hpp"
#include "Schedule.hpp"


#ifdef VAR_PARALLEL_MACHS


class Cand {
public:
	Cand() : bailiff(0), provost(0), direction(0) { }
	Cand(unsigned _bailiff, unsigned _provost, unsigned _direction) : bailiff(_bailiff), provost(_provost), direction(_direction) { assert(direction == BACK   ||   direction == FORTH); }
	~Cand() { }


	
	string toString() const {
		assert(direction == BACK   ||   direction == FORTH);
		return unsigStr(bailiff) + " " + unsigStr(provost) + (direction==BACK ? " B" : " F");
	}



	bool operator==(const Cand & cand) const {
		if(bailiff != cand.bailiff) return false;
		if(provost != cand.provost) return false;
		if(direction != cand.direction) return false;
		return true;
	}


	

	bool isValid() const {
		if(bailiff==0   ||   bailiff>=inst.O) return false;
		if(provost==0   ||   provost>=inst.O) return false;
	      if(direction!=BACK   &&   direction!=FORTH) return false;
		return true;
	}


	
	unsigned bailiff;
	unsigned provost;
	unsigned direction;//BACK or FORTH
};








//information abou the machines (buckets) of all stages (machines)
class Stages {
public:
	Stages() { }
	~Stages() { }
	
	
	
	bool isAlloc() const {
		assert(inst.O > 0);
		if(bucketLimits.shape()[0] != inst.M)
			return false;
		if(bucketLimits.shape()[1] != inst.maxK)
			return false;
		if(bucketTops.shape()[0] != inst.M)
			return false;
		if(bucketTops.shape()[1] != inst.maxK)
			return false;
		return true;
	}
	
	
	
	void clear() {
		fill(bucketLimits.data(),bucketLimits.data()+bucketLimits.num_elements(), 0);
		fill(bucketTops.data(),bucketTops.data()+bucketTops.num_elements(), 0);
	}
	
	
	
	void alloc() {
		bucketLimits.resize(boost::extents[inst.M][inst.maxK]);
		bucketTops.resize(boost::extents[inst.M][inst.maxK]);
	}
	
	
	
	multi_array<unsigned, 2> bucketLimits;
	multi_array<unsigned, 2> bucketTops;
};








//data structure to do topological walks
class Utility {
public:
	Utility() { }
	~Utility() { }
	
	
	
	void alloc() {
		assert(inst.O > 0);
		Q.resize(inst.O);
		indeg.resize(inst.O);
	}
	
	
	
	bool isAlloc() const {
		assert(inst.O > 0);
		
		if(Q.size() != inst.O)
			return false;
		if(indeg.size() != inst.O) 
			return false;
		
		return true;
	}	
	
	
	
	void reset() {
		fill(Q.begin(), Q.end(), 0);
		fill(indeg.begin(), indeg.end(), 0);
	}



	vector<unsigned> Q;
	vector<unsigned> indeg;
};








//complete selection
class Select {
public:
	Select() { }
	~Select() { }
	

	
	void alloc() {
		mach.resize(inst.O);
		_mach.resize(inst.O);
	}



	bool isAlloc() {
		assert(mach.size() == 0   ||   mach.size() == inst.O);
		assert(mach.size() == _mach.size());
		return mach.size() == inst.O;
	}



	bool operator==(const Select & _s) const {
		assert(inst.O > 0);
		assert(_s._mach.size() == inst.O);
		assert(_s.mach.size() == inst.O);
		assert(_mach.size() == inst.O);
		assert(mach.size() == inst.O);
		
		for(unsigned u=0; u<inst.O; u++) {
			if(_s._mach[u] != _mach[u]) return false;
			if(_s.mach[u] != mach[u]) return false;
		}
		return true;
	}

	
	
	string toString() const {
		string str = "";
		str.append("\nmach:");
		for(unsigned o=1; o<mach.size(); o++) {
			str.append("  "+unsigStr(o)+"-"+unsigStr(mach[o]));
			assert(mach[o]==0    ||  _mach[mach[o]] == o);
		}
		return str;
	}

	
	
	void reset() {
		fill(mach.begin(), mach.end(), 0);
		fill(_mach.begin(), _mach.end(), 0);
	}



	void transpose(unsigned o1, unsigned o2) {
		assert(o1 != 0);
		assert(o1 < inst.O);
		assert(o2 != 0);
		assert(o2 < inst.O);
		assert(o1 != o2);
		assert(inst.operToM[o1] == inst.operToM[o2]);
		assert(mach[o1]==o2);
		assert(_mach[o2]==o1);

		unsigned prev = _mach[o1];
		unsigned post = mach[o2];

		//forward
		if(prev != 0)
			mach[prev] = o2;
		mach[o1] = post;
		mach[o2] = o1;
		//backward
		_mach[o1] = o2;
		_mach[o2] = prev;
		if(post != 0)
			_mach[post] = o1;
	}
	


	void shift(const Cand & cand)  { shift(cand.bailiff, cand.provost, cand.direction); }
	
	

	void shift(unsigned bailiff, unsigned provost, unsigned direction) {
		assert(bailiff != 0);
		assert(provost != 0);
		assert(provost != bailiff);
		assert(inst.operToM[bailiff] = inst.operToM[provost]);
		assert(direction == FORTH   ||   direction == BACK);
		
		unsigned x;//adjacent to provost while it moves

		if(direction == FORTH) {//o1==provost         o2==bailiff	
			x = mach[provost];
			while(x != bailiff) {
				assert(x != 0);
				swap(provost, x);
				assert(x != mach[provost]);
				x = mach[provost];
			}
			transpose(provost, bailiff);
		} else {//direction == BACK      o1==bailiff         o2 == provost
			x = _mach[provost];
			while(x != bailiff) {
				assert(x != 0);
				swap(x, provost);
				assert(x != _mach[provost]);
				x = _mach[provost];
			}
			transpose(bailiff, provost);
		}
	}


	
	//fast find makes - do not change meta - UINFTY if cycle
	unsigned fastMakes(Utility & utility, Stages & stages, vector<unsigned> & auxDists) const {
#ifndef NDEBUG
		assert(auxDists.size() == inst.O);
		assert(utility.isAlloc());
		assert(stages.isAlloc());
#endif //NDEBUG
		unsigned qInsert = 0;
		unsigned qAccess = 0;
		unsigned curOp;
		unsigned curMach;
		unsigned newDist;
		unsigned chosenBucket;
		unsigned bucketHead;
		unsigned jobHead;
		bool foundNonDelayBucket;
		unsigned prevJ;
		unsigned nextJ;
		unsigned nextM;
		
		unsigned foundMakes = 0;

		//clearing aux
		utility.reset();
		stages.clear();
		fill(auxDists.begin(), auxDists.end(), 0);
		
		//intialize indeg and put minimals in queue
		for(unsigned o=1; o<inst.O; o++) {
			if( ! inst.prev[o].empty())
				utility.indeg[o]++;
			if(_mach[o] != 0)
				utility.indeg[o]++;
			if(utility.indeg[o] == 0)
				utility.Q[qInsert++] = o;
		}
		
		//forwards topological walk main loop
		while(qAccess < qInsert) {
			//GET NEW CUR OP
			assert(qAccess<utility.Q.size());
			curOp = utility.Q[qAccess++];
			assert(utility.indeg[curOp] == 0);
			
			//AUX FOR CUR OP
			curMach = inst.operToM[curOp];
			if( ! inst.prev[curOp].empty()) {
				assert(inst.prev[curOp].size() == 1   &&   inst.prev[curOp][0] > 0   &&   inst.prev[curOp][0] < inst.O);
				prevJ = inst.prev[curOp][0];
				jobHead = auxDists[prevJ] + inst.P[prevJ];
			} else {
				prevJ = 0;
				jobHead = 0;
			}
			if( ! inst.next[curOp].empty()) {
				assert(inst.next[curOp].size() == 1   &&   inst.next[curOp][0] > 0   &&   inst.next[curOp][0] < inst.O);
				nextJ = inst.next[curOp][0];
			} else {
				nextJ = 0;
			}
			nextM = mach[curOp];
			
			//RELEASE NEXT OPS
			if( nextJ ) {
				assert(utility.indeg[nextJ] > 0);
				utility.indeg[nextJ]--;
				if(utility.indeg[nextJ] == 0) {
					assert(qInsert<utility.Q.size()-1);
					utility.Q[qInsert++] = nextJ;
				}
			}
			if( nextM ) {
				assert(utility.indeg[nextM] > 0);
				utility.indeg[nextM]--;
				if(utility.indeg[nextM] == 0) {
					assert(qInsert<utility.Q.size()-1);
					utility.Q[qInsert++] = nextM;
				}
			}
			
			//FIND PROPPER BUCKET
			//getting biggest bucket that doesnt delay curOp
			foundNonDelayBucket = false;
			for(unsigned bin=0; bin<inst.varK[curMach]; bin++) {
				//if(bucketLimits[curMach][bin] <= jobHead) {
				if(stages.bucketLimits[curMach][bin] < jobHead) {
					if( ! foundNonDelayBucket) {//first non-delay bucket found
						chosenBucket = bin;
						foundNonDelayBucket = true;
					} else {//other non delay bucket previously found
						if(stages.bucketLimits[curMach][bin] > stages.bucketLimits[curMach][chosenBucket])
							chosenBucket = bin; //bin's bucket still smaller the jobHead but bigger then prev chosen
					}
				}
			}
			// all buckets delay curOp? -> get smallest one then
			if( ! foundNonDelayBucket) {
				chosenBucket = 0;
				for(unsigned bin=1; bin<inst.varK[curMach]; bin++) {
					if(stages.bucketLimits[curMach][bin] < stages.bucketLimits[curMach][chosenBucket])
						chosenBucket = bin;
				}
			}
			bucketHead = stages.bucketLimits[curMach][chosenBucket];
			
			//UPDATE dists  - foundMakes
			auxDists[curOp] = max(jobHead, bucketHead);
			newDist = auxDists[curOp] + inst.P[curOp];
			if(foundMakes < newDist) {
				foundMakes = newDist;
			}
			
			//UPDATE BUCKET
			stages.bucketLimits[curMach][chosenBucket] = newDist;
			stages.bucketTops[curMach][chosenBucket] = curOp;
		}//while(qAccess < qInsert) -- forwards topological walk main loop
		
		if(qAccess<inst.O-1)
			return UINFTY;
		else
			return foundMakes;
	}


	
	vector<unsigned> mach;//machine successor
	vector<unsigned> _mach;//machine predecessor
};








//information about tight arcs that can tell us the critical operations to change
class DG {
public:
#ifndef NDEBUG
	DG() : dirty(true) { }
#else
	DG() { }
#endif //NDEBUG
	~DG() { }
	
	
	
	bool operator==(const DG & dgvk) const {
		assert(this->isAlloc());
		assert(dgvk.isAlloc());
		assert(dgvk.select.mach.size() == this->select.mach.size());
		assert(dgvk.select._mach.size() == this->select._mach.size());
		assert(dgvk.select.mach.size() == this->select._mach.size());
		
		for(unsigned u=0; u<dgvk.select.mach.size(); u++) {
			if(dgvk.select.mach[u] != this->select.mach[u])
				return false;
		}
		return true;
	}
	
	
	
	string toString() const {
		assert(isAlloc());
		string str = "";
		
		str.append("millis: "+ unsigStr(millisecs));
		str.append("   makes: "+ unsigStr(makes));
		
		str.append("\nheads:");
		for(unsigned x=0; x<heads.size(); x++)
	 		str.append(" "+unsigStr(x)+ "(" +unsigStr(heads[x])+")");

	 	str.append(select.toString());
	
	 	str.append("\nprevJobTightArc:");
	 	for(unsigned x=0; x<prevJobTightArc.size(); x++)
	 		str.append(" "+unsigStr(x)+"("+unsigStr(prevJobTightArc[x])+")");
	 		
	 	str.append("\nnextJobTightArc:");
	 	for(unsigned x=0; x<nextJobTightArc.size(); x++)
	 		str.append(" "+unsigStr(x)+"("+unsigStr(nextJobTightArc[x])+")");
	 	
	 	str.append("\nmachHolds:");
	 	for(unsigned x=0; x<machHolds.size(); x++) {
	 		str.append("\n\t"+unsigStr(x)+": ");
	 		for(unsigned o : machHolds[x]) {
	 			str.append(" " + unsigStr(o));
	 		}
	 	}
	 	
	 	str.append("\nmachHeldBy:");
	 	for(unsigned x=0; x<machHeldBy.size(); x++) {
	 		str.append("\n\t"+unsigStr(x)+": ");
	 		for(unsigned o : machHeldBy[x]) {
	 			str.append(" " + unsigStr(o));
	 		}
	 	}

		return str;
	}
	
	
	
	void alloc() {
		assert(inst.O > 0);
		assert(heads.empty());
		assert(prevJobTightArc.empty());
		assert(nextJobTightArc.empty());
		assert(machHolds.empty());
		assert(machHeldBy.empty());
		
		millisecs = UINFTY;
		makes = 0;
		
		select.alloc();
		
		heads.resize(inst.O);
		
		prevJobTightArc.resize(inst.O);
		nextJobTightArc.resize(inst.O);
		
		machHolds.resize(inst.O);
		for(vector<unsigned> & v : machHolds)
			v.reserve(inst.maxK);
		machHeldBy.resize(inst.O);
		
		for(vector<unsigned> & v : machHeldBy)
			v.reserve(inst.maxK);

		clear();
	}
	
	
	
	bool isAlloc() const {
		assert(inst.O > 0);
		
		if(heads.size() != inst.O)
			return false;
		if(prevJobTightArc.size() != inst.O)
			return false;
		if(nextJobTightArc.size() != inst.O) 
			return false;
		if(machHolds.size() != inst.O) 
			return false;
		if(machHeldBy.size() != inst.O) 
			return false;
		
		return true;
	}
	
	
	
	void clear() {
		select.reset();
		clearMeta();
#ifndef NDEBUG
		dirty = true;
#endif //NDEBUG
	}
	
	
	
	void clearMeta() {
		makes = 0;
		fill(heads.begin(), heads.end(), 0);
		fill(prevJobTightArc.begin(), prevJobTightArc.end(), 0);
		fill(nextJobTightArc.begin(), nextJobTightArc.end(), 0);
		for(vector<unsigned> & v : machHolds)
			v.clear();
		for(vector<unsigned> & v : machHeldBy)
			v.clear();
	}
	
	
	
	//calls job tail heuristic from Schedule and converts the result into a DG
	//DOES NOT SET META
	void jobTailHeuristic() {
		assert(isAlloc());
	
		Schedule sch;
		vector<vector<unsigned>> allMachOrder(inst.M);
			
		sch.alloc();
		sch.dispatchLongestJobTail(SMALLER_GAP);
		sch.getMachOrders(allMachOrder);

		//XXX REMOVE
		cout << sch.toString() << endl;
		
		importFromMachOrder(allMachOrder);
	}
	
	
	
	//receives the order of the machines and sets the state accordingly
	//designed to be used in conjunction with Schedule::getMachOrders()
	//Does not setMeta
	void importFromMachOrder(const vector<vector<unsigned>> & allMachOrder) {
		assert(isAlloc());
		assert(allMachOrder.size() == inst.M);
		
		unsigned curOp;
		
		//set machines
		for(unsigned m=0; m<inst.M; m++) {
			assert(allMachOrder[m].size() == inst.J);
			
			for(unsigned index=0; index<inst.J; index++) {
				curOp = allMachOrder[m][index];
				assert(curOp>0    &&    curOp<inst.O);
				
				if(index != 0)
					select._mach[curOp] = allMachOrder[m][index-1];
				else
					select._mach[curOp] = 0;
					
				if(index != inst.J-1)
					select.mach[curOp] = allMachOrder[m][index+1];
				else
					select.mach[curOp] = 0;
			}
		}
#ifndef NDEBUG
		dirty = true;
#endif //NDEBUG
	}
	
	
	
	void importFromSelect(const Select & _select) {
		select = _select;
#ifndef NDEBUG
		dirty = true;
#endif //NDEBUG
	}
	
	

	//sets all fields based on mach and _mach
	//@return: valid?
	//@utility: vector<unsigned> Q   vector<unsigned> indeg  
	//@stages:   multi_array<unsigned, 2> bucketLimits    multi_array<unsigned, 2> bucketTops
	bool setMeta(vector<Cand> & cands, Utility & utility, Stages & stages) {
#ifndef NDEBUG
		assert(isAlloc());
		assert(dirty);
		assert(utility.isAlloc());
		assert(stages.isAlloc());
#endif //NDEBUG
		unsigned qInsert = 0;
		unsigned qAccess = 0;
		unsigned curOp;
		unsigned curMach;
		unsigned newDist;
		unsigned chosenBucket;
		unsigned bucketHead;
		unsigned jobHead;
		bool foundNonDelayBucket;
		unsigned prevJ;
		unsigned nextJ;
		unsigned prevM;
		unsigned nextM;
	
		cands.clear();
		
		//minimal/maximal elements of DG that are part of a critical path
		//vector<unsigned> roots;
		vector<unsigned> critLeafs;
		vector<bool> critOpers(inst.O, false);//is in some critical path?
		
		makes = 0;	
		
		//clearing aux
		clearMeta();
		utility.reset();
		stages.clear();
		
		//intialize indeg and put minimals in queue
		for(unsigned o=1; o<inst.O; o++) {
			if( ! inst.prev[o].empty())
				utility.indeg[o]++;
			if(select._mach[o] != 0)
				utility.indeg[o]++;
			if(utility.indeg[o] == 0)
				utility.Q[qInsert++] = o;
		}
		
		//forwards topological walk main loop
		while(qAccess < qInsert) {
			//GET NEW CUR OP
			assert(qAccess<utility.Q.size());
			curOp = utility.Q[qAccess++];
			assert(utility.indeg[curOp] == 0);
			
			//AUX FOR CUR OP
			curMach = inst.operToM[curOp];
			if( ! inst.prev[curOp].empty()) {
				assert(inst.prev[curOp].size() == 1   &&   inst.prev[curOp][0] > 0   &&   inst.prev[curOp][0] < inst.O);
				prevJ = inst.prev[curOp][0];
				jobHead = heads[prevJ] + inst.P[prevJ];
			} else {
				prevJ = 0;
				jobHead = 0;
			}
			if( ! inst.next[curOp].empty()) {
				assert(inst.next[curOp].size() == 1   &&   inst.next[curOp][0] > 0   &&   inst.next[curOp][0] < inst.O);
				nextJ = inst.next[curOp][0];
			} else {
				nextJ = 0;
			}
			prevM = select._mach[curOp];
			nextM = select.mach[curOp];
			
			//RELEASE NEXT OPS
			if( nextJ ) {
				assert(utility.indeg[nextJ] > 0);
				utility.indeg[nextJ]--;
				if(utility.indeg[nextJ] == 0) {
					assert(qInsert<utility.Q.size()-1);
					utility.Q[qInsert++] = nextJ;
				}
			}
			if( nextM ) {
				assert(utility.indeg[nextM] > 0);
				utility.indeg[nextM]--;
				if(utility.indeg[nextM] == 0) {
					assert(qInsert<utility.Q.size()-1);
					utility.Q[qInsert++] = nextM;
				}
			}
			
			//FIND PROPPER BUCKET
			//getting biggest bucket that doesnt delay curOp
			foundNonDelayBucket = false;
			for(unsigned bin=0; bin<inst.varK[curMach]; bin++) {
				//if(bucketLimits[curMach][bin] <= jobHead) {
				if(stages.bucketLimits[curMach][bin] < jobHead) {
					if( ! foundNonDelayBucket) {//first non-delay bucket found
						chosenBucket = bin;
						foundNonDelayBucket = true;
					} else {//other non delay bucket previously found
						if(stages.bucketLimits[curMach][bin] > stages.bucketLimits[curMach][chosenBucket])
							chosenBucket = bin; //bin's bucket still smaller the jobHead but bigger then prev chosen
					}
				}
			}
			// all buckets delay curOp? -> get smallest one then
			if( ! foundNonDelayBucket) {
				chosenBucket = 0;
				for(unsigned bin=1; bin<inst.varK[curMach]; bin++) {
					if(stages.bucketLimits[curMach][bin] < stages.bucketLimits[curMach][chosenBucket])
						chosenBucket = bin;
				}
			}
			bucketHead = stages.bucketLimits[curMach][chosenBucket];
			
			//UPDATE dists  - makes - leafs
			heads[curOp] = max(jobHead, bucketHead);
			newDist = heads[curOp] + inst.P[curOp];
			if(makes < newDist) {
				makes = newDist;
				critLeafs.clear();
				critLeafs.push_back(curOp);
			} else if (makes == newDist) {
				assert( ! critLeafs.empty());
				critLeafs.push_back(curOp);
			}
			
			//update prevJobTightArc & nextJobTightArc- find all job tight arcs
			if( prevJ ) {
				assert(jobHead <= heads[curOp]);
				assert(prevJobTightArc[curOp] == 0);
				if(jobHead == heads[curOp]) {
					prevJobTightArc[curOp] = prevJ;
					nextJobTightArc[prevJ] = curOp;
				}
			}
			
			//update machHolds & machHeldBy- find all machine tight arcs
			if(bucketHead > 0) {
				assert(prevM);
				for(unsigned bin=0; bin<inst.varK[curMach]; bin++) {
					//assert(stages.bucketLimits[curMach][bin] <= heads[curOp]);
					if(stages.bucketLimits[curMach][bin] >= jobHead) {
						machHolds[stages.bucketTops[curMach][bin]].push_back(curOp);
						machHeldBy[curOp].push_back(stages.bucketTops[curMach][bin]);
					}
				}
			}
			
			//UPDATE BUCKET
			stages.bucketLimits[curMach][chosenBucket] = newDist;
			stages.bucketTops[curMach][chosenBucket] = curOp;
		}//while(qAccess < qInsert) -- forwards topological walk main loop
		
#ifndef NDEBUG
		Utility  auxUtility;
		auxUtility.alloc();
		Stages  auxStages;
		auxStages.alloc();
		vector<unsigned>  auxAuxDists(inst.O, 0);
		unsigned auxMakes = select.fastMakes(auxUtility, auxStages, auxAuxDists);
		if(qAccess==inst.O-1)
			assert(makes == auxMakes);
		else
			assert(UINFTY == auxMakes);
#endif //NDEBUG
		
		if(qAccess<inst.O-1)
			return false;
		
		utility.reset();
		stages.clear();
		qInsert = 0;
		qAccess = 0;
		
		//intialize critOpers
		for(unsigned o : critLeafs)
			critOpers[o] = true;
		
		//intialize indeg and put minimals in queue
		for(unsigned o=1; o<inst.O; o++) {
			if( ! inst.next[o].empty())
				utility.indeg[o]++;
			if(select.mach[o] != 0)
				utility.indeg[o]++;
			if(utility.indeg[o] == 0)
				utility.Q[qInsert++] = o;
		}
		
		//backward topological walk main loop
		while(qAccess < qInsert) {
			//GET NEW CUR OP
			assert(qAccess<utility.Q.size());
			curOp = utility.Q[qAccess++];
			assert(utility.indeg[curOp] == 0);
			
			//AUX FOR CUR OP
			curMach = inst.operToM[curOp];
			if( ! inst.prev[curOp].empty()) {
				assert(inst.prev[curOp].size() == 1   &&   inst.prev[curOp][0] > 0   &&   inst.prev[curOp][0] < inst.O);
				prevJ = inst.prev[curOp][0];
			} else {
				prevJ = 0;
			}
			if( ! inst.next[curOp].empty()) {
				assert(inst.next[curOp].size() == 1   &&   inst.next[curOp][0] > 0   &&   inst.next[curOp][0] < inst.O);
				nextJ = inst.next[curOp][0];
			} else {
				nextJ = 0;
			}
			prevM = select._mach[curOp];
			nextM = select.mach[curOp];
			
			//RELEASE PREV OPS
			if( prevJ ) {
				assert(utility.indeg[prevJ] > 0);
				utility.indeg[prevJ]--;
				if(utility.indeg[prevJ] == 0) {
					assert(qInsert<utility.Q.size()-1);
					utility.Q[qInsert++] = prevJ;
				}
			}
			if( prevM ) {
				assert(utility.indeg[prevM] > 0);
				utility.indeg[prevM]--;
				if(utility.indeg[prevM] == 0) {
					assert(qInsert<utility.Q.size()-1);
					utility.Q[qInsert++] = prevM;
				}
			}
			
			if( critOpers[curOp] ) {
				if( prevJobTightArc[curOp] ) {
					critOpers[prevJobTightArc[curOp]] = true;
					//Add to neighbourhood in here for partial shop or open shop
				}
				for(unsigned o : machHeldBy[curOp]) {
					//assert( ! critOpers[o] );
					critOpers[o] = true;
					//Add to neighborhood:     o <-> curOp
					//case LEP  -  (e,i): tight job arc (Pj(e), e)
					//     o==e   curOp==i    -> tight job arc (Pj(o), o)
					if(prevJobTightArc[o]) {
						assert(inst.prev[o].size() == 1);
						assert(nextJobTightArc[inst.prev[o][0]]);
						//bailiff == e == o      provost == i == curOp
						cands.push_back(Cand(o, curOp, BACK));
						//case REP  -  (i,e): tight job arc (e, Sj(e))
						//     o==i   curOp==e    -> tight job arc (curOp, Sj(curOp))
					} else if(nextJobTightArc[curOp]) {
						assert(inst.next[o].size() == 1);
						//assert(prevJobTightArc[inst.next[o][0]]);//NOT TURE FOR PARALLEL !
						//bailiff == e == curOp      provost == i == o
						cands.push_back(Cand(curOp, o, FORTH));
					}
					
				}
			}
			
		}//while(qAccess < qInsert) -- backward topological walk main loop
		
		assert(qAccess==inst.O-1);//after backward while or else would already have finished in first while
		
#ifndef NDEBUG
		dirty = false;
#endif //NDEBUG
		
		return true;
	}
	
	
	
	//expects o1 to be before o2
	//DOES NOT update meta
	void shift(unsigned bailiff, unsigned provost, unsigned direction) {
		select.shift(bailiff, provost, direction);
#ifndef NDEBUG
		dirty = true;
#endif //NDEBUG
	}

	

	unsigned millisecs;
	
	Select select;

	unsigned makes;//UINFTY if invalid, 0 if none defined
	
	vector<unsigned> heads; //start times as well
	//vector<unsigned> tails; //XXX remove tails if trimming is not used

	
	//0 if there is no tight arcs for job
	vector<unsigned> prevJobTightArc;
	vector<unsigned> nextJobTightArc;

	//machHolds[x]: x held back the elements in machHolds[x]
	vector<vector<unsigned>> machHolds;// may hold any number of other perations back
	//whenever a new operation found only delay buckets all fortier operations of the buckets add the new op to machHolds
	
	//machHeldBy[x]: operations that held x back when scheduled in the machine
	vector<vector<unsigned>> machHeldBy; //from 0 to varK[inst.operToM[x]] opers in vector
	//whenever a new operation found only delay buckets the new op adds all fortier operations of the buckets add 
	//sum of sizes of vecs in machHeldBy must be == to the sum of machHeld
	
	//vector<unsigned> bestBucketLimit;//XXX USe for trimming - when o was scheduled what was best H_M(o)? 
	
#ifndef NDEBUG
	bool dirty;
#endif //NDEBUG
};


#endif //VAR_PARALLEL_MACHS
