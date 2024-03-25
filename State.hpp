/*
  @author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/
#pragma once

#include <boost/multi_array.hpp>
#include <climits>
#include <algorithm> 
#include <chrono>
#include <set>
#include <queue>

#include "Settings.hpp"
#include "Utilities.hpp"
#include "Instance.hpp"
#include "RandomNumberGenerator.hpp"
#include "Schedule.hpp"

#include "Analyzer.hpp"

using namespace boost;


class State {
public:
	State() {  }
	~State() {  }


	bool operator==(const State & other) const {
		assert(verify());
		assert(other.verify());

		for(unsigned pos=0; pos<inst.O; pos++) {
			if(job[pos] != other.job[pos]) return false;
			if(_job[pos] != other._job[pos]) return false;
			if(mach[pos] != other.mach[pos]) return false;
			if(_mach[pos] != other._mach[pos]) return false;
		}

		return true;
	}



	void alloc() {
		assert(inst.O != 0);
		assert( ! isAlloced());
		job.resize(inst.O);
		_job.resize(inst.O);
		mach.resize(inst.O);
		_mach.resize(inst.O);
	}



	void clear() {
		assert(inst.O != 0);
		assert(mach.size() == inst.O);
		assert(job.size() == inst.O);

		//roots.clear();
		//leafs.clear();
		fill(job.begin(), job.end(), 0);
		fill(_job.begin(), _job.end(), 0);
		fill(mach.begin(), mach.end(), 0);
		fill(_mach.begin(), _mach.end(), 0);
	}



	string toString() const {
		assert(job[0] == 0);
		assert(_job[0] == 0);
		assert(mach[0] == 0);
		assert(_mach[0] == 0);

		string str = "makes: " + unsigStr(makes);

		str.append("\njob: ");
		for(unsigned o=1; o<job.size(); o++) {
			str.append("  "+unsigStr(o)+"->"+unsigStr(job[o]));
			assert(job[o]==0    ||   _job[job[o]] == o);
		}

		str.append("\nmach:");
	    for(unsigned o=1; o<mach.size(); o++) {
		    str.append("  "+unsigStr(o)+"->"+unsigStr(mach[o]));
			assert(mach[o]==0    ||  _mach[mach[o]] == o);
		}

		return str;
	}
	//expects fully defined state
	string easyString() const {
		assert(job[0] == 0);
		assert(_job[0] == 0);
		assert(mach[0] == 0);
		assert(_mach[0] == 0);

		unsigned root;
		unsigned curOp;

		string str = "makes: " + unsigStr(makes);

		str.append("\njob:");
		for(unsigned j=0; j<inst.J; j++) {
			root = 0;
			for(unsigned o : inst.jobOpers[j]) {
				if(_job[o] == 0) {
					assert(root == 0);
					root = o;
				}
			}
			curOp = root;
			str.append("\n\t");
			for(unsigned u=0; u<inst.jSize[j]; u++) {
				assert(curOp != 0);
				str.append(unsigStr(curOp)+" ");
				curOp = job[curOp];
			}
		}

		str.append("\nmach:");
		for(unsigned m=0; m<inst.M; m++) {
			root = 0;
			for(unsigned o : inst.machOpers[m]) {
				if(_mach[o] == 0) {
					assert(root == 0);
					root = o;
				}
			}
			curOp = root;
			str.append("\n\t");
			for(unsigned u=0; u<inst.mSize[m]; u++) {
				assert(curOp != 0);
				str.append(unsigStr(curOp)+" ");
				curOp = mach[curOp];
			}
		}

		return str;
	}



	//x and y in same job 
	//max(head(Pm[y]), head(Pj[x]))+procT[X]+procT[Y]+max(tail(Sj]y]), tail(Sm[x]))
	//max(head(Pj[x], head(Pm[y])) +procT[y] + tail(Sm[x])
	//head(Pm[x]) + procT[x] +max(tail(Sj[y], tail(Sm[y])
	unsigned swapLowerBound(const unsigned x, const unsigned y, const vector<unsigned> & heads, const vector<unsigned> & tails) const {
		assert(heads.size() == inst.O);
		assert(tails.size() == inst.O);
		assert(x!=0);
		assert(y!=0);
		assert(inst.operToJ[x]==inst.operToJ[y]   ||   inst.operToM[x]==inst.operToM[y]);
		assert( (job[x] == y   && _job[y]==x)   ||   (mach[x] == y   && _mach[y]==x) );

		unsigned hA1;
		unsigned hA2;
		unsigned tB1;
		unsigned tB2;
		unsigned hAlpha;
		unsigned tBetha;
		unsigned hGamma;
		unsigned tDelta;
		unsigned lb;

		if(inst.operToJ[x] == inst.operToJ[y])  {
			hA1 = (_job[x]!=0   ?   heads[_job[x]] + inst.P[_job[x]]   :   0);
			hA2 = (_mach[y]!=0   ?   heads[_mach[y]] + inst.P[_mach[y]]   :   0);
			tB1 =  (job[y]!=0   ?   tails[job[y]] + inst.P[job[y]]   :   0);
			tB2 =  (mach[x]!=0   ?   tails[mach[x]] + inst.P[mach[x]]   :   0);
			hGamma = (_mach[x]!=0   ?   heads[_mach[x]] + inst.P[_mach[x]]   :   0);
			tDelta = (mach[y]!=0   ?   tails[mach[y]] + inst.P[mach[y]]   :   0);
		} else {//same mach
			hA1 = (_mach[x]!=0   ?   heads[_mach[x]] + inst.P[_mach[x]]   :   0);
			hA2 = (_job[y]!=0   ?   heads[_job[y]] + inst.P[_job[y]]   :   0);
			tB1 =  (mach[y]!=0   ?   tails[mach[y]] + inst.P[mach[y]]   :   0);
			tB2 =  (job[x]!=0   ?   tails[job[x]] + inst.P[job[x]]   :   0);
			hGamma = (_job[x]!=0   ?   heads[_job[x]] + inst.P[_job[x]]   :   0);
			tDelta = (job[y]!=0   ?   tails[job[y]] + inst.P[job[y]]   :   0);
		}
		hAlpha = max(hA1, hA2);
		tBetha = max(tB1, tB2);

		lb = hAlpha + inst.P[y] + inst.P[x] + tBetha;
		lb = max(lb, hAlpha +inst.P[y]+ tDelta);
		lb = max(lb, hGamma + inst.P[x] + tBetha);

		return lb;
	}



	void swap(unsigned o1, unsigned o2) {
		assert(o1 != 0);
		assert(o2 != 0);
		assert(o1 != o2);

		unsigned prev;
		unsigned post;
#ifdef EXTENDED_JSS
		assert(inst.operToM[o1] == inst.operToM[o2]);
		//MACH swap
		assert(mach[o1]==o2);
		assert(_mach[o2]==o1);
		prev = _mach[o1];
		post = mach[o2];
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
#endif //EXTENDED_JSS
#ifdef VANILLA_PSS
		if(inst.operToJ[o1] == inst.operToJ[o2]) {
			//JOB swap
			assert(job[o1]==o2);
			assert(_job[o2]==o1);
			prev = _job[o1];
			post = job[o2];
			//forward
			if(prev != 0)
				job[prev] = o2;
			job[o1] = post;
			job[o2] = o1;
			//backward
			_job[o1] = o2;
			_job[o2] = prev;
			if(post != 0)
				_job[post] = o1;
		} else {
			assert(inst.operToM[o1] == inst.operToM[o2]);
			//MACH swap
			assert(mach[o1]==o2);
			assert(_mach[o2]==o1);
			prev = _mach[o1];
			post = mach[o2];
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
#endif //VANILLA_PSS
	}



	//moves o1 next to o2 passing through it
	//o2 stays still while o1 is moved until it goes through o2
	//expects that moving fomr o1 in the <direction> will hit oBailiff eventually
	void shift(unsigned type, unsigned direction, unsigned oBailiff, unsigned oProvost) {
#ifndef NDEBUG
		unsigned iters = 0;
		assert(oProvost != 0);
		assert(oBailiff != 0);
		assert(direction==FORTH   ||   direction==BACK);
		assert((type==JOB && inst.operToJ[oProvost] == inst.operToJ[oBailiff])
			   ||  (type==MACH && inst.operToM[oProvost] == inst.operToM[oBailiff]));
#endif


		//for each swap: BEFORE alpha -> x -> y -> omega
		//                       AFTER   alpha -> y -> x -> omega
		unsigned alpha;
		unsigned x;
		unsigned y;
		unsigned omega;

		if(direction == FORTH) {
			x = oProvost;
			if(type == JOB)
				y = job[x];
			else
				y = mach[x];
		} else {//BACK
			y = oProvost;
			if(type == JOB)
				x = _job[y];
			else
				x = _mach[y];
		}

		//job shift
		if(type == JOB) {
		    while(true) {
#ifndef NDEBUG
				iters++;
				assert(x!=0);
				assert(y!=0);
				assert(iters<inst.O);
#endif
				alpha = _job[x];
				omega = job[y];
				//normal arcs
				if(alpha != 0)
					job[alpha] = y;
				job[y] = x;
				job[x] = omega;
				//inverted arcs
				if(omega != 0)
					_job[omega] = x;
				_job[x] = y;
				_job[y] = alpha;

				if(direction == FORTH) {
					if(y==oBailiff)
						break;
					y = job[x];
				} else {
					if(x==oBailiff)
						break;
					x = _job[y];
				}
			}//while true for shifting until swap oProvost and oBailiff
		} else { // mach shift
			while(true) {
#ifndef NDEBUG
				iters++;
				assert(x!=0);
				assert(y!=0);
				assert(iters<inst.O);
#endif
				alpha = _mach[x];
				omega = mach[y];
				//normal arcs
				if(alpha != 0)
					mach[alpha] = y;
				mach[y] = x;
				mach[x] = omega;
				//inverted arcs
				if(omega != 0)
					_mach[omega] = x;
				_mach[x] = y;
				_mach[y] = alpha;

				if(direction == FORTH) {
					if(y==oBailiff)
						break;
					y = mach[x];
				} else {
					if(x==oBailiff)
						break;
					x = _mach[y];
				}
			}//while true for shifting until swap oProvost and oBailiff
		}
	}	



	void makeRand(vector<unsigned> & indeg, vector<unsigned> & Q) {
#ifndef NDEBUG
		assert(isAlloced());
		assert(isClear());
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);
		fill(indeg.begin(), indeg.end(), 0);
#endif

		unsigned qInsert;
		unsigned qAccess;
		unsigned curOp;
		unsigned lastOp;
		unsigned rootOp;

		unsigned swapCount;

		vector<unsigned> opList;
		opList.reserve(inst.O);

		vector<unsigned> mList(inst.M);

		//insert JOBS in topo order
		for(unsigned j=0; j<inst.J; j++) {
			qInsert = 0;
			qAccess = 0;

			//set indeg
			opList = inst.jobOpers[j];
			random_shuffle(opList.begin(), opList.end());
			for(unsigned o : opList) {
				assert(inst.operToJ[o] == j);
				assert(indeg[o] == 0);
				indeg[o] = inst.prev[o].size();
				if(indeg[o] == 0)
					Q[qInsert++] = o;
			}
			assert(qInsert>0);

			curOp = 0;

			//topological walk
			while(qAccess < qInsert) {
				assert(qAccess<inst.M);
				assert(qAccess<Q.size());

				lastOp = curOp;
				curOp = Q[qAccess++];

				assert(curOp != 0);
				assert(curOp < inst.O);
				assert(indeg[curOp] == 0);
				assert(inst.operToJ[curOp] == j);

				if(lastOp != 0) job[lastOp] = curOp;
				_job[curOp] = lastOp;

				//job precs
				opList = inst.next[curOp];
				random_shuffle(opList.begin(), opList.end());
				for(unsigned n : opList) {
					assert(indeg[n]>0);
					indeg[n]--;
					if(indeg[n] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = n;
					}
				}
			}
			assert(qAccess == inst.jSize[j]);
			assert( ! hasCycle(indeg, Q));
		}
		
		//insert machines wthout making a cycle
		for(unsigned m=0; m<inst.M; m++) mList[m] = m;
		random_shuffle(mList.begin(), mList.end());
		for(unsigned m : mList) {
			opList =  inst.machOpers[m];
			random_shuffle(opList.begin(), opList.end());
			rootOp = 0;
			
			for(unsigned o : opList) {
				assert(o != 0);
				assert( ! hasCycle(indeg, Q));
				assert(_mach[rootOp] == 0);
				assert(_mach[o] == 0);

				//insert o in first postion

				mach[o] = rootOp;
				if(rootOp != 0) _mach[rootOp] = o;

				swapCount = 0;

				//swap until no cycle
				while(hasCycle(indeg, Q)) {
					assert(mach[o] != 0); //must have a valid place
					swap(o, mach[o]);
					swapCount++;
				}

				if(swapCount == 0) rootOp = o;
			}
		}	   
	}



	//partial state has cycles?
	bool hasCycle(vector<unsigned> & indeg, vector<unsigned> & Q) const {

		assert(isAlloced());
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);

		unsigned qInsert = 0;
		unsigned qAccess = 0;

		unsigned curOp;
		unsigned newOp;

		fill(indeg.begin(), indeg.end(), 0);

		for(unsigned o=1; o<inst.O; o++) {
			indeg[o] = inst.prev[o].size();
			if(_job[o] != 0)
				indeg[o]++;
			if(_mach[o] != 0)
				indeg[o]++;
			if(indeg[o] == 0)
				Q[qInsert++] = o;
		}


		while(qAccess < qInsert) {
			assert(qAccess<Q.size());
			curOp = Q[qAccess++];
			assert(indeg[curOp] == 0);

			//job order
			newOp = job[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
			}
			//mach order
			newOp = mach[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
			}
			//job precs
			for(unsigned n : inst.next[curOp]) {
				assert(indeg[n]>0);
				indeg[n]--;
				if(indeg[n] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = n;
				}
			}
		}
		assert(qAccess<=inst.O-1);
		return qAccess<inst.O-1;
	}

	//expects compete state
	//@return: cycle?
	bool findTopo(vector<unsigned> & topo, vector<unsigned> & indeg) {
		assert(isAlloced());
		assert(indeg.size() == inst.O);
		assert(topo.size() == inst.O);

		unsigned qInsert = 0;
		unsigned qAccess = 0;

		unsigned curOp;
		unsigned newOp;

		fill(indeg.begin(), indeg.end(), 0);

		for(unsigned o=1; o<inst.O; o++) {
			if(_job[o] != 0)
				indeg[o]++;
			if(_mach[o] != 0)
				indeg[o]++;
			if(indeg[o] == 0)
				topo[qInsert++] = o;
		}	

		while(qAccess < qInsert) {
			assert(qAccess<topo.size());
			curOp = topo[qAccess++];
			assert(indeg[curOp] == 0);

			//job order
			newOp = job[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<topo.size()-1);
					topo[qInsert++] = newOp;
				}
			}
			//mach order
			newOp = mach[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<topo.size()-1);
					topo[qInsert++] = newOp;
				}
			}
		}
		assert(qAccess<=inst.O-1);
		return qAccess<inst.O-1;
	}


	//says if two opers from same item can reach one another -  ignores order of item but not precs - if not enforce then guaranteed obey order and do not form cycle;
	//make an "enforce" for job and one for mach with dif sizes - type JOB then size is M (contrary)
	//@enforce[x][y] => x->y  ---  x and y are not op index but the index of the machine or job
	//@operPos: position in topological order of each operation of item
	void mustEnforce(vector<vector<unsigned>> & enforce, vector<unsigned> & enfIndeg, unsigned type, unsigned index,  vector<unsigned> & indeg, vector<unsigned> & Q, vector<unsigned> & operPos, vector<bool> & reach) const {
#ifndef NDEBUG
		assert(enforce.size() == inst.O);
		for(const vector<unsigned> & v : enforce) assert(v.capacity() == max(inst.J, inst.M));
		assert(enfIndeg.size() == inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);
		assert(operPos.size() == inst.O);
		assert(reach.size() == inst.O);
#endif


		unsigned qInsert = 0;
		unsigned qAccess = 0;

		unsigned curOp;
		unsigned newOp;

		unsigned maxPos = 0; //last position of an oper 

		fill(indeg.begin(), indeg.end(), 0);
		//reseting enforce and enfIndeg
		for(unsigned o : (type==JOB ? inst.jobOpers[index] : inst.machOpers[index])) {
			enforce[o].clear();
			enfIndeg[o] = 0;
		}

		for(unsigned o=1; o<inst.O; o++) {
			if(type!=JOB   ||    inst.operToJ[o]!=index) {
				if(_job[o] != 0)
					indeg[o]++;
			} else {//type==JOB   &&    inst.operToJ[curOp]==index
				indeg[o] += inst.prev[o].size();
			}
			if(_mach[o]!=0    &&   (type!=MACH    ||    inst.operToM[o]!=index))
				indeg[o]++;
			if(indeg[o] == 0)
				Q[qInsert++] = o;
		}

		assert(qInsert>0);

		//preparing topological order
		while(qAccess < qInsert) {
			assert(qAccess<Q.size());
			curOp = Q[qAccess];
			assert(indeg[curOp] == 0);

			if(index   ==   (type==JOB ? inst.operToJ[curOp] : inst.operToM[curOp])) {
				operPos[curOp] = qAccess;
				maxPos = max(maxPos, qAccess);
			}

			qAccess++;

			//from JOB
			if(type!=JOB   ||    inst.operToJ[curOp]!=index) {
				newOp = job[curOp];
			    if(newOp != 0) {
					assert(indeg[newOp] > 0);
					indeg[newOp]--;
					if(indeg[newOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = newOp;
					}
				}
			} else {//type==JOB   &&    inst.operToJ[curOp]==index
				for(unsigned nOp : inst.next[curOp]) {
					assert(indeg[nOp] > 0);
					indeg[nOp]--;
					if(indeg[nOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = nOp;
					}
				}
			}

			//from MACH
			newOp = mach[curOp];
			if(newOp != 0    &&   (type!=MACH   ||    inst.operToM[curOp]!=index)) {
				assert(indeg[newOp] > 0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
			}
		}
		assert(qAccess == inst.O-1);
		//Here Q is a topological order, operPos has the positions of the item's opers, and maxPos is the biggest

		//propagating forward from each element in item
		for(unsigned fromOp : (type==JOB ? inst.jobOpers[index] : inst.machOpers[index])) {
			assert(Q[operPos[fromOp]] = fromOp);			
			fill(reach.begin(), reach.end(), false);
			reach[fromOp] = true;
			for(unsigned aPos=operPos[fromOp]; aPos<=maxPos; aPos++) {
				curOp = Q[aPos];
				if(reach[curOp]) {
					if(index == (type==JOB ? inst.operToJ[curOp] : inst.operToM[curOp])) {
						if(curOp != fromOp) {
							enforce[fromOp].push_back(curOp);
							enfIndeg[curOp]++;
						}
					}

					//from JOB
					if(type!=JOB   ||    inst.operToJ[curOp]!=index) {
						newOp = job[curOp];
						if(newOp != 0)
							reach[newOp] = true;	
					} else {//type==JOB   &&    inst.operToJ[curOp]==index
						for(unsigned nOp : inst.next[curOp])
							reach[nOp] = true;	
					}

					//from MACH
					newOp = mach[curOp];
					if(newOp != 0    &&   (type!=MACH   ||    inst.operToM[curOp]!=index)) {
						reach[newOp] = true;
					}
				}
			}
		}
	}
	


	void bubblePerturb(double p, unsigned perturbSize, vector<vector<unsigned>> & enforce, vector<unsigned> & enfIndeg, vector<unsigned> & indeg, vector<unsigned> & Q, vector<unsigned> & operPos, vector<bool> & reach, const ElementProb & ep) {
#ifndef NDEBUG
		assert(enforce.size() == inst.O);
		for(const vector<unsigned> & v : enforce) assert(v.capacity() == max(inst.J, inst.M));
		assert(enfIndeg.size() == inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);
		assert(operPos.size() == inst.O);
		assert(reach.size() == inst.O);
		assert(p>0.0 && p<1.0);
		assert(perturbSize <= inst.O);
#endif

		unsigned root = 0;
		unsigned curOp;
		unsigned type;
		vector<unsigned> next;
		vector<unsigned> prev;
		unsigned lastInsertedOp;
		unsigned bubbleTries = 0;

		makes = 0;

		unsigned index;
		for(unsigned iter=0; iter<perturbSize; iter++) {
			//selecting element to bubble

			tie(type, index) = ep.randElement();

			mustEnforce(enforce, enfIndeg, type, index, indeg, Q, operPos, reach);
 			//enforce is poset so no cycles and precs are ok, enfIndeg is indeg of such poset

			//finding root
#ifdef VANILLA_PSS
			for(unsigned u : (type==JOB ? inst.posets[index].roots  :  inst.machOpers[index])) {
				curOp = (type==JOB ? inst.jmToIndex[index][u] : u);
				if((type==JOB ? _job[curOp]  :  _mach[curOp]) == 0) {
					root = curOp;
					break;
				}
			}
#endif //VANILLA_PSS
#ifdef EXTENDED_JSS
			for(unsigned u : inst.machOpers[index]) {
				assert(type==MACH);
				curOp = u;
				if((type==JOB ? _job[curOp]  :  _mach[curOp]) == 0) {
					root = curOp;
					break;
				}
			}
#endif //EXTENDED_JSS
			assert(root!=0);
			assert(enfIndeg[root] == 0);

#ifndef NDEBUG
			//confirm that order is feasible accroing to must enforce
			unsigned a = root;
			unsigned b;
			for(unsigned bailiff=0; bailiff<(type==JOB  ?  inst.M : inst.J)-1; bailiff++) {
				assert(a != 0);
				b = (type==JOB  ?  job[a] : mach[a]);
				for(unsigned provost=bailiff+1; provost<(type==JOB  ?  inst.M : inst.J); provost++) {
					assert(b != 0);
					for(unsigned e : enforce[b])
						assert(e != a);
					b = (type==JOB  ?  job[b] : mach[b]);
				}
				a = (type==JOB  ?  job[a] : mach[a]);
			}
#endif
			
			//Bubbling
			next = (type==JOB ? job  :  mach);
			prev = (type==JOB ? _job  :  _mach);
			lastInsertedOp = 0;
			while(root != 0) {
				curOp = root;
				while(curOp != 0) {
					bubbleTries++;
					if(bubbleTries > 10000) throw errorText("bubble stuck!", "bubblePerturb", "State.hpp");
					if(enfIndeg[curOp]==0   &&   random_number<double>(0.0, 1.0)<=p) {
						//update indeg
						for(unsigned o : enforce[curOp]) {
							assert(enfIndeg[o] > 0);
							enfIndeg[o]--;
						}
						//insert in actual order
						if(lastInsertedOp != 0) 
							(type==JOB ? job[lastInsertedOp]=curOp  :  mach[lastInsertedOp]=curOp);
						(type==JOB ? _job[curOp]=lastInsertedOp  :  _mach[curOp]=lastInsertedOp);
						//removing inserted oper from aux (prev, next)
						next[prev[curOp]] = next[curOp];
						prev[next[curOp]] = prev[curOp];
						//updateing root 
						if(curOp == root)
							root = next[curOp];
						lastInsertedOp = curOp;
					}
					curOp = next[curOp];
				}
			}
			(type==JOB ? job[lastInsertedOp]=0  :  mach[lastInsertedOp]=0);
		}
	}
	




	static void removeOperFromOrder(unsigned targOp, vector<unsigned> & next, vector<unsigned> & prev) {
		assert(targOp!=0);

		const unsigned prevOp = prev[targOp];
		assert(prevOp==0   ||   next[prevOp] == targOp);
		const unsigned postOp = next[targOp];
		assert(postOp==0   ||   prev[postOp] == targOp);

		next[targOp] = 0;
		prev[targOp] = 0;
		if(prevOp != 0) prev[prevOp] = postOp;
		if(postOp != 0) prev[postOp] = prevOp;
	}



	void removeOper(unsigned targOp) {
		assert(targOp!=0);
		
		const unsigned prevJobOp = _job[targOp];
		assert(prevJobOp==0   ||   job[prevJobOp] == targOp);
		const unsigned postJobOp = job[targOp];
		assert(postJobOp==0   ||   _job[postJobOp] == targOp);
		const unsigned prevMachOp = _mach[targOp];
		assert(prevMachOp==0   ||   mach[prevMachOp] == targOp);
		const unsigned postMachOp = mach[targOp];
		assert(postMachOp==0   ||   _mach[postMachOp] == targOp);

		job[targOp] = 0;
		_job[targOp] = 0;
		mach[targOp] = 0;
		_mach[targOp] = 0;
		if(prevJobOp != 0) job[prevJobOp] = postJobOp;
		if(postJobOp != 0) _job[postJobOp] = prevJobOp;
		if(prevMachOp != 0) mach[prevMachOp] = postMachOp;
		if(postMachOp != 0) _mach[postMachOp] = prevMachOp;
	}



	void removeOperFromMach(unsigned targOp) {
		assert(targOp!=0);

		const unsigned prevMachOp = _mach[targOp];
		assert(prevMachOp==0   ||   mach[prevMachOp] == targOp);
		const unsigned postMachOp = mach[targOp];
		assert(postMachOp==0   ||   _mach[postMachOp] == targOp);

		mach[targOp] = 0;
		_mach[targOp] = 0;
		if(prevMachOp != 0) mach[prevMachOp] = postMachOp;
		if(postMachOp != 0) _mach[postMachOp] = prevMachOp;
	}



	//insert newOp right after prevOp and before postOp --- needs both in case prev is 0 needs to know cur root
	//care if using roots and leafs - must adjust outside
	void insertOper(unsigned newOp, unsigned prevJobOp, unsigned postJobOp, unsigned prevMachOp, unsigned postMachOp) {
		insertJobOper(newOp, prevJobOp, postJobOp);
		insertMachOper(newOp, prevMachOp, postMachOp);
	}
	void insertJobOper(unsigned newOp, unsigned prevJobOp, unsigned postJobOp) {
		assert(newOp != 0);
		assert(prevJobOp==0   ||   job[prevJobOp]==postJobOp);
		assert(postJobOp==0   ||   _job[postJobOp]==prevJobOp);

		job[newOp] = postJobOp;
		if(postJobOp!=0) _job[postJobOp] = newOp;
		_job[newOp] = prevJobOp;
		if(prevJobOp!=0) job[prevJobOp] = newOp;
	}
	void insertMachOper(unsigned newOp, unsigned prevMachOp, unsigned postMachOp) {
		assert(newOp != 0);
		assert(prevMachOp==0   ||   mach[prevMachOp]==postMachOp);
		assert(postMachOp==0   ||   _mach[postMachOp]==prevMachOp);

		mach[newOp] = postMachOp;
		if(postMachOp!=0) _mach[postMachOp] = newOp;
		_mach[newOp] = prevMachOp;
		if(prevMachOp!=0) mach[prevMachOp] = newOp;
	}



	bool compHeadsComplete(vector<unsigned> & heads, vector<unsigned> & indeg, vector<unsigned> & Q) const {
		assert(isAlloced());
		assert(heads.size() == inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);

		unsigned qInsert = 0;
		unsigned qAccess = 0;

		unsigned curOp;
		unsigned newOp;

		unsigned newMax;

		fill(heads.begin(), heads.end(), 0);
		fill(indeg.begin(), indeg.end(), 0);

		for(unsigned o=1; o<inst.O; o++) {
			if(_job[o] != 0)
				indeg[o]++;
			if(_mach[o] != 0)
				indeg[o]++;
			if(indeg[o] == 0)
				Q[qInsert++] = o;
		}

		assert(qInsert>0);

		while(qAccess < qInsert) {
			assert(qAccess<Q.size());
			curOp = Q[qAccess++];
			assert(indeg[curOp] == 0);
			newMax = heads[curOp] + inst.P[curOp];

			//job order
			newOp = job[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
				heads[newOp] = max(heads[newOp], newMax);
			}
			//mach order
			newOp = mach[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
				heads[newOp] = max(heads[newOp], newMax);
			}
		}
	
		assert(qAccess<=inst.O-1);
		return qAccess<inst.O-1;
	}

	
	

	bool compHeads(vector<unsigned> & heads, vector<unsigned> & indeg, vector<unsigned> & Q) const {
		assert(isAlloced());
		assert(heads.size() == inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);

		unsigned qInsert = 0;
		unsigned qAccess = 0;

		unsigned curOp;
		unsigned newOp;

		unsigned newMax;

		fill(heads.begin(), heads.end(), 0);

		for(unsigned o=1; o<inst.O; o++) {
			indeg[o] = inst.prev[o].size();
			if(_job[o] != 0)
				indeg[o]++;
			if(_mach[o] != 0)
				indeg[o]++;
			if(indeg[o] == 0)
				Q[qInsert++] = o;
		}

		//assert(qInsert>0);//cycle may reate no first op in topo

		while(qAccess < qInsert) {
			assert(qAccess<Q.size());
			curOp = Q[qAccess++];
			assert(indeg[curOp] == 0);
			newMax = heads[curOp] + inst.P[curOp];

			//job order
			newOp = job[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
				heads[newOp] = max(heads[newOp], newMax);
			}
			//mach order
			newOp = mach[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
				heads[newOp] = max(heads[newOp], newMax);
			}
			//job precs
			for(unsigned n : inst.next[curOp]) {
				assert(indeg[n]>0);
				indeg[n]--;
				if(indeg[n] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = n;
				}
				heads[n] = max(heads[n], newMax);
			}
		}
	
		assert(qAccess<=inst.O-1);
		return qAccess<inst.O-1;
	}



	bool compTailsComplete(vector<unsigned> & tails, vector<unsigned> & indeg, vector<unsigned> & Q) const {
		assert(isAlloced());
		assert(tails.size() == inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);
		unsigned qInsert = 0;
		unsigned qAccess = 0;

		unsigned curOp;
		unsigned newOp;

		unsigned newMax;

		fill(tails.begin(), tails.end(), 0);
		fill(indeg.begin(), indeg.end(), 0);

		for(unsigned o=1; o<inst.O; o++) {
			if(job[o] != 0)
				indeg[o]++;
			if(mach[o] != 0)
				indeg[o]++;
			if(indeg[o] == 0)
				Q[qInsert++] = o;
		}

		assert(qInsert>0);

		while(qAccess < qInsert) {
			assert(qAccess<Q.size());
			curOp = Q[qAccess++];
			assert(indeg[curOp] == 0);
			newMax = tails[curOp] + inst.P[curOp];

			//job order
			newOp = _job[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
				tails[newOp] = max(tails[newOp], newMax);
			}
			//mach order
			newOp = _mach[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
				tails[newOp] = max(tails[newOp], newMax);
			}
		}
	
		assert(qAccess<=inst.O-1);
		return qAccess<inst.O-1;
	}




	bool compTails(vector<unsigned> & tails, vector<unsigned> & indeg, vector<unsigned> & Q) const {
		assert(isAlloced());
		assert(tails.size() == inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);
		unsigned qInsert = 0;
		unsigned qAccess = 0;

		unsigned curOp;
		unsigned newOp;

		unsigned newMax;

		fill(tails.begin(), tails.end(), 0);

		for(unsigned o=1; o<inst.O; o++) {
			indeg[o] = inst.next[o].size();
			if(job[o] != 0)
				indeg[o]++;
			if(mach[o] != 0)
				indeg[o]++;
			if(indeg[o] == 0)
				Q[qInsert++] = o;
		}

		//assert(qInsert>0);

		while(qAccess < qInsert) {
			assert(qAccess<Q.size());
			curOp = Q[qAccess++];
			assert(indeg[curOp] == 0);
			newMax = tails[curOp] + inst.P[curOp];

			//job order
			newOp = _job[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
				tails[newOp] = max(tails[newOp], newMax);
			}
			//mach order
			newOp = _mach[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
				tails[newOp] = max(tails[newOp], newMax);
			}
			//job precs
			for(unsigned n : inst.prev[curOp]) {
				assert(indeg[n]>0);
				indeg[n]--;
				if(indeg[n] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = n;
				}
				tails[n] = max(tails[n], newMax);
			}
		}
	
		assert(qAccess<=inst.O-1);
		return qAccess<inst.O-1;
	}



	//can reach fromIndex toIndex in Q?
	bool propagateReach(const vector<unsigned> & Q, vector<bool> & reach, unsigned fromIndexA, unsigned fromIndexB, unsigned toIndexA, unsigned toIndexB) const {
		unsigned curOp;
		unsigned newOp;
		unsigned minIndex = min(fromIndexA, fromIndexB);
		unsigned maxIndex = max(toIndexA, toIndexB);

		unsigned targA = Q[toIndexA];
		unsigned targB = Q[toIndexB];

		//cout << "\t\t\t\ttargA: " << targA << "  ;  targB: " << targB << endl;

		fill(reach.begin(), reach.end(), false);
		if(Q[fromIndexA] != 0)
			reach[Q[fromIndexA]] = true;
		if(Q[fromIndexB] != 0)
			reach[Q[fromIndexB]] = true;

		for(unsigned topoIndex=minIndex; topoIndex<maxIndex; topoIndex++) {
			curOp = Q[topoIndex];
			//assert(curOp != targA   && curOp != targB);
			if(reach[curOp]) {
				//job order
				newOp = job[curOp];
				if(newOp != 0) {
					if(newOp == targA   ||   newOp == targB) {
						//tcout << "\t\t\t\t\tX - newOp" << newOp << " ; curOp: " << curOp  << endl;
						return true;
					}
					reach[newOp] = true;
				}
				//mach order
				newOp = mach[curOp];
				if(newOp != 0) {
					if(newOp == targA   ||   newOp == targB) {
						//cout << "\t\t\t\t\tY - newOp"  << newOp << " ; curOp: " << curOp  << endl;
						return true;
					}
					reach[newOp] = true;
				}
				//job precs
				for(unsigned n : inst.next[curOp]) {
				    assert(n != 0);
					if(n == targA   ||   n == targB) {
						//cout << "\t\t\t\t\tZ - newOp"  << newOp << " ; curOp: " << curOp  << endl;
						return true;
					}
					reach[n] = true;
				}
			}
		}

		return false;
	}

#ifdef EXTENDED_JSS
	void insertInsaPos(unsigned curOp, vector<unsigned> & indeg, vector<unsigned> & Q, vector<bool> & inOpers, vector<unsigned> & jobRoot, vector<unsigned> & jobLeaf, vector<unsigned> & machRoot, vector<unsigned> & machLeaf, vector<unsigned> & heads, vector<unsigned> & tails, vector<unsigned> & invertedQ, vector<bool> & reach) {
#ifndef NDEBUG
		bool cycle;
		vector<unsigned> auxQ(inst.O);
		vector<unsigned> auxIndeg(inst.O);
		assert(curOp != 0   && curOp<inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);
		assert(inOpers.size() == inst.O);
		assert( ! inOpers[curOp]);
		assert(machRoot.size() == inst.M);
		assert(machLeaf.size() == inst.M);
		for(unsigned mr : machRoot) {
			assert(mr==0   ||   inOpers[mr]);
			assert(_mach[mr]==0);
		}
		for(unsigned ml : machLeaf) {
			assert(ml==0   ||   inOpers[ml]);
			assert(mach[ml]==0);
		}
		bool found;
		for(unsigned m=0; m<inst.M; m++) {
			assert((machRoot[m] && machLeaf[m])   ||   ( ! machRoot[m] && ! machLeaf[m]));
			found=false;
			for(unsigned o : inst.machOpers[m]) found |= inOpers[o];
			if(found) {
				assert(machRoot[m] != 0);
				assert(machLeaf[m] != 0);
			} else {
				assert(machRoot[m] == 0);
				assert(machLeaf[m] == 0);
			}
		}
#endif

		//XXX
		cout << "\n" << toString() << endl;


		
		
		unsigned curCost;
		const unsigned curJob = inst.operToJ[curOp];
		const unsigned curMach = inst.operToM[curOp];


		//XXX
		cout << "curOp: " << curOp << " ; curJob: " << curJob << " ; curMach: " << curMach << endl;

		unsigned bestCost = UINT_MAX;
		//unsigned bufOp;
		unsigned jobHead;
		unsigned jobTail;

		unsigned bestJobPrev = 0;
		unsigned bestJobNext = 0;
		unsigned bestMachPrev = 0;
		unsigned bestMachNext = 0;

		assert( ! inOpers[curOp]);
		inOpers[curOp] = true;
	 
		
#ifndef NDEBUG
		cycle = 
#endif
			compTails(tails, indeg, Q);
		assert( ! cycle);
#ifndef NDEBUG
		cycle = 
#endif
			compHeads(heads, indeg, Q);//important that heads is after tails - Q contains a topo sort
		assert( ! cycle);

		//finding inverted table of Q
		for(unsigned u=0; u<inst.O; u++)
			invertedQ[Q[u]] = u;

		//finding job heads and tails --- either from order or from precs
		jobHead = heads[_job[curOp]]+inst.P[_job[curOp]];
		for(unsigned p : inst.prev[curOp])
			jobHead = max(jobHead, heads[p]+inst.P[p]);			
		jobTail = tails[job[curOp]] + inst.P[job[curOp]];
		for(unsigned n : inst.next[curOp])
			jobTail = max(jobTail, tails[n] + inst.P[n]);


		insertMachOper(curOp, 0, machRoot[curMach]);//insert in 1st mach position
		while(true) {//mach runthrough
			curCost = jobHead + jobTail;
			curCost = max(curCost,  jobHead + tails[mach[curOp]] +inst.P[mach[curOp]]);
			curCost = max(curCost,   heads[_mach[curOp]] + inst.P[_mach[curOp]] + jobTail);
			curCost = max(curCost,   heads[_mach[curOp]] + inst.P[_mach[curOp]] + tails[mach[curOp]] +inst.P[mach[curOp]]);

			if(bestCost > curCost) {
				assert(Q[invertedQ[job[curOp]]] == job[curOp]);
				assert(Q[invertedQ[_job[curOp]]] == _job[curOp]);
				assert(Q[invertedQ[mach[curOp]]] == mach[curOp]);
				assert(Q[invertedQ[_mach[curOp]]] == _mach[curOp]);
				if( ! propagateReach(Q, reach, invertedQ[mach[curOp]], invertedQ[job[curOp]], invertedQ[_job[curOp]],  invertedQ[_mach[curOp]])) {
					bestCost = curCost;
					bestJobPrev = _job[curOp];
					bestJobNext = job[curOp];
					bestMachPrev = _mach[curOp];
					bestMachNext = mach[curOp];					
					assert( ! hasCycle(auxIndeg, auxQ));
				} else {	
					assert(hasCycle(auxIndeg, auxQ));
				}
			}

			if(mach[curOp] == 0)
				break;
			//cout << "\t\tswapping:  " << curOp << " " << mach[curOp];
			swap(curOp, mach[curOp]);
			//cout << "   done !!!!" << endl;
		}//mach runthrough
		removeOperFromMach(curOp);

		assert(bestCost < UINT_MAX);

		insertMachOper(curOp, bestMachPrev, bestMachNext);//insert in 1st mach position

		//remove from last postion and reinsert in best place
		//removeOper(curOp);
		//insertOper(curOp, bestJobPrev, bestJobNext, bestMachPrev, bestMachNext);

		assert(partialVerify(inOpers));

		//adjust roots and leafs
		if(bestJobNext == jobRoot[curJob])
			jobRoot[curJob] = curOp;
		if(bestJobPrev == jobLeaf[curJob])
			jobLeaf[curJob] = curOp;
		if(bestMachNext == machRoot[curMach])
			machRoot[curMach] = curOp;
		if(bestMachPrev == machLeaf[curMach])
			machLeaf[curMach] = curOp;
	}
#endif //EXTENDED_JSS
#ifdef VANILLA_PSS
	void insertInsaPos(unsigned curOp, vector<unsigned> & indeg, vector<unsigned> & Q, vector<bool> & inOpers, vector<unsigned> & jobRoot, vector<unsigned> & jobLeaf, vector<unsigned> & machRoot, vector<unsigned> & machLeaf, vector<unsigned> & heads, vector<unsigned> & tails, vector<unsigned> & invertedQ, vector<bool> & reach) {
#ifndef NDEBUG
		bool cycle;
		vector<unsigned> auxQ(inst.O);
		vector<unsigned> auxIndeg(inst.O);
		assert(curOp != 0   && curOp<inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);
		assert(inOpers.size() == inst.O);
		assert( ! inOpers[curOp]);
		assert(jobRoot.size() == inst.J);
		assert(jobLeaf.size() == inst.J);
		assert(machRoot.size() == inst.M);
		assert(machLeaf.size() == inst.M);
		for(unsigned jr : jobRoot) {
			assert(jr==0   ||   inOpers[jr]);
			assert(_job[jr]==0);
		}
		for(unsigned jl : jobLeaf) {
			assert(jl==0   ||   inOpers[jl]);
			assert(job[jl]==0);
		}
		for(unsigned mr : machRoot) {
			assert(mr==0   ||   inOpers[mr]);
			assert(_mach[mr]==0);
		}
		for(unsigned ml : machLeaf) {
			assert(ml==0   ||   inOpers[ml]);
			assert(mach[ml]==0);
		}
		bool found;
		for(unsigned j=0; j<inst.J; j++) {
			assert((jobRoot[j] && jobLeaf[j])   ||   ( ! jobRoot[j] && ! jobLeaf[j]));
			found=false;
			for(unsigned o : inst.jobOpers[j]) found |= inOpers[o];
			if(found) {
				assert(jobRoot[j] != 0);
				assert(jobLeaf[j] != 0);
			} else {
				assert(jobRoot[j] == 0);
				assert(jobLeaf[j] == 0);
			}
		}
		for(unsigned m=0; m<inst.M; m++) {
			assert((machRoot[m] && machLeaf[m])   ||   ( ! machRoot[m] && ! machLeaf[m]));
			found=false;
			for(unsigned o : inst.machOpers[m]) found |= inOpers[o];
			if(found) {
				assert(machRoot[m] != 0);
				assert(machLeaf[m] != 0);
			} else {
				assert(machRoot[m] == 0);
				assert(machLeaf[m] == 0);
			}
		}
#endif
		
		unsigned curCost;
		const unsigned curJob = inst.operToJ[curOp];
		const unsigned curMach = inst.operToM[curOp];

		unsigned bestCost = UINT_MAX;
		unsigned bufOp;
		unsigned jobHead;
		unsigned jobTail;

		//opers that delimit where we can insert curOp in job
		unsigned fromHere = 0;
		unsigned toHere = 0;
		bool leftPrecFound = false;
		unsigned nextOper;

		unsigned bestJobPrev = 0;
		unsigned bestJobNext = 0;
		unsigned bestMachPrev = 0;
		unsigned bestMachNext = 0;

		assert( ! inOpers[curOp]);
		inOpers[curOp] = true;

		//checking limits in job for curOp
		bufOp = jobRoot[curJob];
		while(bufOp != 0) {
			if(inst.hasPrec(bufOp, curOp)) {
				fromHere = bufOp;
				leftPrecFound = true;
			}
			if(inst.hasPrec(curOp, bufOp)) {
				toHere = bufOp;
				break;
			}
			bufOp = job[bufOp];
		}
		if(! leftPrecFound) {//no precs -> fromHere is null and next is root (may be also null if empty)
			assert(fromHere == 0);
			nextOper = jobRoot[curJob];
		} else {// precs found
			assert(jobRoot[curJob] != 0);//not empty
			assert(fromHere != 0);//something is actually prec to curOp
			nextOper = job[fromHere];
		}
		
#ifndef NDEBUG
		cycle = 
#endif
			compTails(tails, indeg, Q);
		assert( ! cycle);
#ifndef NDEBUG
		cycle = 
#endif
			compHeads(heads, indeg, Q);//important that heads is after tails - Q contains a topo sort
		assert( ! cycle);

		//finding inverted table of Q
		for(unsigned u=0; u<inst.O; u++)
			invertedQ[Q[u]] = u;
		//evaluating position pairs
		insertJobOper(curOp, fromHere, nextOper);//insert in 1st valid job position
		while(true) {//job runthrough
			//finding job heads and tails --- either from order or from precs
			jobHead = heads[_job[curOp]]+inst.P[_job[curOp]];
			for(unsigned p : inst.prev[curOp])
				jobHead = max(jobHead, heads[p]+inst.P[p]);			
			jobTail = tails[job[curOp]] + inst.P[job[curOp]];
			for(unsigned n : inst.next[curOp])
				jobTail = max(jobTail, tails[n] + inst.P[n]);		

			insertMachOper(curOp, 0, machRoot[curMach]);//insert in 1st mach position
			while(true) {//mach runthrough
				curCost = jobHead + jobTail;
				curCost = max(curCost,  jobHead + tails[mach[curOp]] +inst.P[mach[curOp]]);
				curCost = max(curCost,   heads[_mach[curOp]] + inst.P[_mach[curOp]] + jobTail);
				curCost = max(curCost,   heads[_mach[curOp]] + inst.P[_mach[curOp]] + tails[mach[curOp]] +inst.P[mach[curOp]]);

				if(bestCost > curCost) {
					assert(Q[invertedQ[job[curOp]]] == job[curOp]);
					assert(Q[invertedQ[_job[curOp]]] == _job[curOp]);
					assert(Q[invertedQ[mach[curOp]]] == mach[curOp]);
					assert(Q[invertedQ[_mach[curOp]]] == _mach[curOp]);
					if( ! propagateReach(Q, reach, invertedQ[mach[curOp]], invertedQ[job[curOp]], invertedQ[_job[curOp]],  invertedQ[_mach[curOp]])) {
						bestCost = curCost;
						bestJobPrev = _job[curOp];
						bestJobNext = job[curOp];
						bestMachPrev = _mach[curOp];
						bestMachNext = mach[curOp];						
						assert( ! hasCycle(auxIndeg, auxQ));
					} else {
						assert(hasCycle(auxIndeg, auxQ));
					}
				}

				if(mach[curOp] == 0)
					break;
				swap(curOp, mach[curOp]);
			}//mach runthrough
			removeOperFromMach(curOp);

			if(job[curOp] == toHere)
				break;
			assert(job[curOp] != 0);
			swap(curOp, job[curOp]);
		}//job runthrough
		assert(bestCost < UINT_MAX);

		//remove from last postion and reinsert in best place
		removeOper(curOp);
		insertOper(curOp, bestJobPrev, bestJobNext, bestMachPrev, bestMachNext);

		assert(partialVerify(inOpers));

		//adjust roots and leafs
		if(bestJobNext == jobRoot[curJob])
			jobRoot[curJob] = curOp;
		if(bestJobPrev == jobLeaf[curJob])
			jobLeaf[curJob] = curOp;
		if(bestMachNext == machRoot[curMach])
			machRoot[curMach] = curOp;
		if(bestMachPrev == machLeaf[curMach])
			machLeaf[curMach] = curOp;
	}
#endif //VANILLA_PSS


	//never removes the last oper of a job or mach
	//exects inOpers to be all true - computes roots and leafs
	void insaPerturb(unsigned perturbSize, vector<unsigned> & indeg, vector<unsigned> & Q, vector<bool> & inOpers, vector<unsigned> & jobRoot, vector<unsigned> & jobLeaf, vector<unsigned> & machRoot, vector<unsigned> & machLeaf, vector<unsigned> & heads, vector<unsigned> & tails, vector<unsigned> & invertedQ, vector<bool> & reach) {
		unsigned targOp;
		vector<unsigned> chosenOps;
		chosenOps.reserve(perturbSize);
		unsigned targJ;
		unsigned targM;
		unsigned iters = 0;

		if(perturbSize > inst.O-inst.J-inst.M)
			perturbSize = inst.O-inst.J-inst.M;

		//setting roots & leafs
		for(unsigned o=1; o<inst.O; o++) {
			assert(inOpers[o]);
			if(_job[o]==0) jobRoot[inst.operToJ[o]] = o;
			if(job[o]==0) jobLeaf[inst.operToJ[o]] = o;
			if(_mach[o]==0) machRoot[inst.operToM[o]] = o;
			if(mach[o]==0) machLeaf[inst.operToM[o]] = o;
		}

		//choosing opers to perturb
		assert(perturbSize<inst.O-1);
		while(chosenOps.size()<perturbSize) {
			iters++;
			if(iters>10000) throw errorText("stuck in insa perturb", "State.hpp", "insaPerturb()");
			targOp = rand()%(inst.O-1)+1;
			targJ = inst.operToJ[targOp];
			targM = inst.operToM[targOp];
			//not removed still and not last of a Job or Machine
			if(inOpers[targOp]   &&   jobRoot[targJ]!=jobLeaf[targJ]   &&   machRoot[targM]!=machLeaf[targM]) {
				chosenOps.push_back(targOp);
				inOpers[targOp] = false;
				if(jobRoot[targJ] == targOp) {
					assert(job[targOp] != 0);
					jobRoot[targJ] = job[targOp];
				}
				if(jobLeaf[targJ] == targOp) {
					assert(_job[targOp] != 0);
					jobLeaf[targJ] = _job[targOp];
				}
				if(machRoot[targM] == targOp) {
					assert(mach[targOp] != 0);
					machRoot[targM] = mach[targOp];
				}
				if(machLeaf[targM] == targOp) {
					assert(_mach[targOp] != 0);
					machLeaf[targM] = _mach[targOp];
				}
				removeOper(targOp);
			}
		}

		//reinserting them insa style
		for(unsigned o : chosenOps)
			insertInsaPos(o, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails, invertedQ, reach);
	}




	void insaPsp(bool bigJobFirst, vector<unsigned> & indeg, vector<unsigned> & Q, vector<unsigned> & heads, vector<unsigned> & tails, vector<unsigned> & invertedQ, vector<bool> & reach) {
		assert(isAlloced());
		assert(isClear());
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);

		unsigned curOp;

		unsigned curCost;
		unsigned bigJob = UINT_MAX;
		unsigned bigJobCost = 0;

		vector<bool> inOpers(inst.O, false);
		vector<unsigned> jobRoot(inst.J, 0);
		vector<unsigned> jobLeaf(inst.J, 0);
		vector<unsigned> machRoot(inst.M, 0);
		vector<unsigned> machLeaf(inst.M, 0);
		
		vector<pair<unsigned, unsigned>> costVoper;
		costVoper.reserve(inst.O);

#ifdef EXTENDED_JSS
		//setting joborder that will not be changed
		for(unsigned o=1; o<inst.O; o++) {
			if( ! inst.next[o].empty()) {
				assert(inst.next[o].size() == 1);
				job[o]  = inst.next[o][0];
			} else {
				job[o] = 0;
			}

			if( ! inst.prev[o].empty()) {
				assert(inst.prev[o].size() == 1);
				_job[o]  = inst.prev[o][0];
			} else {
				_job[o] = 0;
			}
		}
#endif //EXTENDED_JSS

		//inserting big job first
		if(bigJobFirst) {
			for(unsigned j=0; j<inst.J; j++) {
				curCost = 0;
				for(unsigned o : inst.jobOpers[j])
					curCost += inst.P[o];
				if(curCost > bigJobCost) {
					bigJobCost = curCost;
					bigJob = j;
				}
			}
			assert(bigJob < inst.J);
			for(unsigned o : inst.jobOpers[bigJob]) {
				insertInsaPos(o, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails, invertedQ, reach);
			}
		}

		//ordering the rest of the opers
		for(unsigned o=1; o<inst.O; o++) {
			if( ! inOpers[o])
				costVoper.push_back(pair<unsigned, unsigned>(inst.P[o], o));
		}
		sort(costVoper.rbegin(), costVoper.rend());

		assert( (bigJobFirst && costVoper.size()==(inst.O-1-inst.jSize[bigJob])) 
				||   ( ! bigJobFirst && costVoper.size()==(inst.O-1))); 

		//inserting missing ops
		for(unsigned pos=0; pos<costVoper.size(); pos++) {
			curOp = costVoper[pos].second;
			insertInsaPos(curOp, indeg, Q, inOpers, jobRoot, jobLeaf, machRoot, machLeaf, heads, tails,invertedQ, reach);
		}
	}

	void gifflerThompson() {

		set<unsigned> ready;
		vector<unsigned> ready0;
		vector<unsigned> ready1;
		vector<unsigned> jobStartTimes(inst.J, 0);
		vector<unsigned> machStartTimes(inst.M, 0);
		vector<unsigned> jobLeafs(inst.J, 0);
		vector<unsigned> jobLeafs2(inst.J, 0);
 		vector<unsigned> machLeafs(inst.M, 0);
		vector<unsigned> machLeafs2(inst.M, 0);

		unsigned completionTime;
		unsigned earlCompletion;
		unsigned mach;
		unsigned auxStartTime;

		for (unsigned op : inst.roots) {
			ready.insert(op);
		}

		while (!ready.empty()) {

			earlCompletion = UINT_MAX;

			ready0.clear();
			ready1.clear();

			for (unsigned op : ready) {
				assert(op < inst.O);
				completionTime = max(jobStartTimes[inst.operToJ[op]], machStartTimes[inst.operToM[op]]) + inst.P[op];
				if (completionTime < earlCompletion) {
					earlCompletion = completionTime;
					mach = inst.operToM[op];
				}
			}

			for (unsigned op : ready) {
				if (inst.operToM[op] == mach) {
					ready0.push_back(op);
				}
			}

			for (unsigned op : ready0) {
				auxStartTime = max(jobStartTimes[inst.operToJ[op]], machStartTimes[inst.operToM[op]]);
				if (auxStartTime < earlCompletion) {
					ready1.push_back(op);
				}
			}

			assert(ready1.size() > 0);
			unsigned op = ready1[0];

			for (unsigned i = 1; i < ready1.size(); ++i) {
				if (inst.deadlines[op] > inst.deadlines[ready1[i]]) {
					op = ready1[i];
				}
			}

			if (jobLeafs[inst.operToJ[op]]) {
				insertJobOper(jobLeafs[inst.operToJ[op]], jobLeafs2[inst.operToJ[op]], op);
				jobLeafs2[inst.operToJ[op]] = jobLeafs[inst.operToJ[op]];
			}

			jobLeafs[inst.operToJ[op]] = op;

			if (machLeafs[inst.operToM[op]]) {
				insertMachOper(machLeafs[inst.operToM[op]], machLeafs2[inst.operToM[op]], op);
				machLeafs2[inst.operToM[op]] = machLeafs[inst.operToM[op]];
			}

			machLeafs[inst.operToM[op]] = op;
			
			ready.erase(op);

			if (inst.next[op].size()) {
				ready.insert(inst.next[op][0]);
			}

			auxStartTime = max(jobStartTimes[inst.operToJ[op]], machStartTimes[inst.operToM[op]]);
			jobStartTimes[inst.operToJ[op]] = auxStartTime + inst.P[op];
			machStartTimes[inst.operToM[op]] = auxStartTime + inst.P[op];
		}
	}


	static void computeCritic(vector<unsigned> & critic, unsigned lastOp, const vector<unsigned> & prev) {
		assert(lastOp != 0);
		assert(lastOp < inst.O);
		assert(critic.capacity() == inst.O);

		critic.clear();

		while(lastOp != 0) {

			assert(critic.size() < inst.O);
			critic.push_back(lastOp);
			lastOp = prev[lastOp];
		}
		reverse(critic.begin(), critic.end());
	}



	static void blockBegins(vector<unsigned> & jobBb, vector<unsigned> & machBb, const vector<unsigned> & critic) {
		assert( ! critic.empty());
		assert(critic.size() < inst.O);
		assert(jobBb.capacity() == inst.O);
		assert(machBb.capacity() == inst.O);

		jobBb.clear();
		machBb.clear();
		jobBb.push_back(0);
		machBb.push_back(0);

		for(unsigned pos=0; pos<critic.size()-1; pos++) {
			assert(critic[pos] != 0);
			assert(inst.operToJ[critic[pos]] == inst.operToJ[critic[pos+1]]
				   || inst.operToM[critic[pos]] == inst.operToM[critic[pos+1]]);
#ifdef VANILLA_PSS
			//job BBs
			if(inst.operToJ[critic[pos]] != inst.operToJ[critic[pos+1]])
				jobBb.push_back(pos+1);
#endif //VANILLA_PSS 
			//mach BBs
			if(inst.operToM[critic[pos]] != inst.operToM[critic[pos+1]])
				machBb.push_back(pos+1);
		}
	}

	static void fillCandidatesTest1(vector<pair<unsigned, unsigned>> & cands, vector<unsigned> & mach, unsigned lastOp) {

		assert(mach.size() == inst.O);
		assert(cands.capacity() == inst.O);

		for (unsigned currentOp = 1; currentOp < inst.O; ++currentOp) {

			if (inst.operToJ[currentOp] != inst.operToJ[mach[currentOp]] && mach[currentOp]) {
				cands.push_back(pair<unsigned, unsigned>(currentOp, mach[currentOp]));
			}
		}
	}

	static void fillCandidatesN5(vector<pair<unsigned, unsigned>> & cand, vector<unsigned> & jobBb, vector<unsigned> & machBb, const vector<unsigned> & critic) {
		assert(cand.capacity() == inst.O);
		assert(jobBb.capacity() == inst.O);
		assert(machBb.capacity() == inst.O);
		assert( ! critic.empty());

		unsigned lOp;
		unsigned rOp;

		cand.clear();

		blockBegins(jobBb, machBb, critic);

#ifdef PARALLEL_MACHS
		if(jobBb.size() <= 1    ||    machBb.size() <= 1)
			return;
#endif //PARALLEL_MACHS
#ifdef VAR_PARALLEL_MACHS
		if(jobBb.size() <= 1    ||    machBb.size() <= 1)
			return;
#endif //VAR_PARALLEL_MACHS


#ifdef VANILLA_PSS
		assert(jobBb.size() > 1); //or else in job lower bound. why still calling neighbourhhod still?
		assert(jobBb[0] == 0);
		assert(jobBb[1] > 0);
		assert(jobBb[1] < critic.size());
		//JOB BLOCKS
		//first job block
		if(jobBb[0]+2 <= jobBb[1]) {//has 2 or more opers
			lOp = critic[jobBb[1]-2];
			rOp = critic[jobBb[1]-1];
			assert(inst.operToJ[lOp] == inst.operToJ[rOp]);
			if( ! inst.hasPrec(lOp, rOp))
				cand.push_back(pair<unsigned, unsigned>(lOp, rOp));//internal edge
		}
		//middle job blocks
		for(unsigned bbBaseInd=1;  bbBaseInd<jobBb.size()-1; bbBaseInd++) {
			assert(jobBb[bbBaseInd] < jobBb[bbBaseInd+1]);
			assert(jobBb[bbBaseInd+1]<critic.size());
			if(jobBb[bbBaseInd]+2 <= jobBb[bbBaseInd+1]) { //block has size at least 2
				//left edge
				lOp = critic[jobBb[bbBaseInd]];
				rOp = critic[jobBb[bbBaseInd]+1];
				if( ! inst.hasPrec(lOp, rOp))
					cand.push_back(pair<unsigned, unsigned>(lOp, rOp));

				if(jobBb[bbBaseInd]+3 <= jobBb[bbBaseInd+1]) {//block has size at least 3
					//right edge
					lOp = critic[jobBb[bbBaseInd+1]-2];
					rOp = critic[jobBb[bbBaseInd+1]-1];
					if( ! inst.hasPrec(lOp, rOp))
						cand.push_back(pair<unsigned, unsigned>(lOp, rOp));
				}
			}
		}
		//last job block
		if(jobBb.back()+2<=critic.size()) {//has 2 or more opers
			lOp = critic[jobBb.back()];
			rOp = critic[jobBb.back()+1];
			if( ! inst.hasPrec(lOp, rOp))
				cand.push_back(pair<unsigned, unsigned>(lOp, rOp));
		}
#endif //VANILLA_PSS
		//MACH BLOCKS
		assert(machBb.size() > 1);//or else in machine lower bound. why still calling neighbourhhod still?
		assert(machBb[0] == 0);
		assert(machBb[1] > 0);
		assert(machBb[1] < critic.size());
		
		//first mach block
		if(machBb[0]+2 <= machBb[1]) {//has 2 or more opers
			lOp = critic[machBb[1]-2];
			rOp = critic[machBb[1]-1];
			assert(inst.operToM[lOp] == inst.operToM[rOp]);
			cand.push_back(pair<unsigned, unsigned>(lOp, rOp));
		}
		//middle mach blocks
		for(unsigned bbBaseInd=1;  bbBaseInd<machBb.size()-1; bbBaseInd++) {
			assert(machBb[bbBaseInd] < machBb[bbBaseInd+1]);
			assert(machBb[bbBaseInd+1]<critic.size());
			if(machBb[bbBaseInd]+2 <= machBb[bbBaseInd+1]) { //block has size at least 2
				//left edge
				lOp = critic[machBb[bbBaseInd]];
				rOp = critic[machBb[bbBaseInd]+1];
				cand.push_back(pair<unsigned, unsigned>(lOp, rOp));

				if(machBb[bbBaseInd]+3 <= machBb[bbBaseInd+1]) {//block has size at least 3
					//right edge
					lOp = critic[machBb[bbBaseInd+1]-2];
					rOp = critic[machBb[bbBaseInd+1]-1];
					cand.push_back(pair<unsigned, unsigned>(lOp, rOp));
				}
			}
		}
		//last mach block
		if(machBb.back()+2<=critic.size()) {//has 2 or more opers
			lOp = critic[machBb.back()];
			rOp = critic[machBb.back()+1];
			cand.push_back(pair<unsigned, unsigned>(lOp, rOp));//internal edge
		}


		assert( ! cand.empty());
		assert(cand.size() < inst.O);
	}



	//@type: shared element: MACH or JOB
	//@bailiff and provost: static and moving ops
	//@return: guarateed not to create cycle?
	bool nonCycleN6Criterium(unsigned direction, unsigned bailiff, unsigned provost, unsigned type, const vector<unsigned> & heads, const vector<unsigned> & tails) const {

		//cout << "\t\tTesting: "  << (type==MACH ? "  MACH" : "  JOB") << (direction==FORTH ? ",FORTH,b" : ",BACK,b") << bailiff << ",p" << provost << endl;

		
		if(direction == FORTH) { //u is provost - v is bailiff
			//from Zhang(2007): u - v in crit: Moving u after v is acyclic if L(v, *) >= L(js[u], *) ---- L(u,v) is sum in interval [u,v)
			if(type == JOB) {
				assert(inst.operToJ[bailiff] == inst.operToJ[provost]);
				return tails[bailiff] + inst.P[bailiff] >= tails[mach[provost]] + inst.P[mach[provost]];
			} else {
				assert(type == MACH);
				assert(inst.operToM[bailiff] == inst.operToM[provost]);
				return tails[bailiff] + inst.P[bailiff] >= tails[job[provost]] + inst.P[job[provost]];
			}
		} else { //BACK: u is bailiff - v is provost
			assert(direction == BACK);
			//from Zhang(2007): u - v in crit: Moving v before u is acyclic if L(0,u) + p_u >= L(0, jp[v]) + p_jp[v] ---- L(u,v) is sum in interval [u,v)			
			if(type == JOB) {
				assert(inst.operToJ[bailiff] == inst.operToJ[provost]);
				return heads[bailiff] + inst.P[bailiff] >= heads[_mach[provost]] + inst.P[_mach[provost]];
			} else {
				assert(type == MACH);
				assert(inst.operToM[bailiff] == inst.operToM[provost]);
				return heads[bailiff] + inst.P[bailiff] >= heads[_job[provost]] + inst.P[_job[provost]];
			}
		}
	}




	//Luv(0, *) = max {Luv(0,w) + Luv(w, *)  :  w in Q} - Q is {u,...., V}
	//
/*	unsigned guidepost(const vector<unsigned> & heads, const vector<unsigned> & tails) const {


	}
	*/
	
	

	//Balas Vazacopoulos 1998
	//cand: type - direction - bailiff - provost
	void fillCandidatesN6(vector<tuple<unsigned, unsigned, unsigned, unsigned>> & cands, vector<unsigned> & jobBb, vector<unsigned> & machBb, const vector<unsigned> & critic, const vector<unsigned> & heads, const vector<unsigned> & tails) const {
		assert(cands.capacity() == ((inst.O)^2));
		assert(jobBb.capacity() == inst.O);
		assert(machBb.capacity() == inst.O);
		assert( ! critic.empty());

		cands.clear();

		blockBegins(jobBb, machBb, critic);

		/*
		cout << "\n\n\nCalled N6";
		cout << toString() << endl;
		cout << "Crit: ";
		for(unsigned o : critic) cout << o << "(" << inst.operToJ[o] << "," << inst.operToM[o] << ")  ";
		cout << endl;
		*/

		unsigned leftEdgePos;
		unsigned rightEdgePos;
#ifdef VANILLA_PSS
		bool canShift;
#endif
		assert(jobBb.size() > 1); //or else in job lower bound. why still calling neighbourhhod still?
		assert(machBb.size() > 1);
		assert(jobBb[0] == 0);
		assert(machBb[0] == 0);
		assert(jobBb[1] > 0);
		assert(machBb[1] > 0);
		assert(jobBb[1] < critic.size());
		assert(machBb[1] < critic.size());

		//adding dummy final block for symmery in code
		jobBb.push_back(critic.size());
		machBb.push_back(critic.size());


		//cout << "jobBb:";
		//for(unsigned edge : jobBb) cout << " "<< edge;
		//cout << endl;

		//cout << "machBb:";
		//for(unsigned edge : machBb) cout << " "<< edge;
		//cout << endl;
		

#ifdef VANILLA_PSS
		//JOB BLOCKS
		for(unsigned bbBaseInd=0;  bbBaseInd<jobBb.size()-1; bbBaseInd++) {
			leftEdgePos = jobBb[bbBaseInd];
			rightEdgePos = jobBb[bbBaseInd+1]-1;
			assert(leftEdgePos <= rightEdgePos);
			assert(rightEdgePos <= critic.size());


			// 1 - foward to right edge: aPos -> rightEdgePos
		    for(unsigned aPos=leftEdgePos; aPos+1<=rightEdgePos; aPos++) {

				//cout << "\taPos -> rightEdgePos" << endl;
				//cout << "\t\t: " << critic[aPos] << "->" << critic[rightEdgePos] << endl;
				
				assert(leftEdgePos < rightEdgePos); //block is longer than 1
				canShift = true;
				for(unsigned otherPos=leftEdgePos+1; otherPos<=rightEdgePos; otherPos++) {
					if(inst.hasPrec(critic[aPos] , critic[otherPos])) {
						canShift = false;
						break;
					}
				}

				if(canShift)
					if(nonCycleN6Criterium(FORTH, critic[rightEdgePos], critic[aPos], JOB, heads, tails))
						cands.push_back(make_tuple(JOB, FORTH, critic[rightEdgePos], critic[aPos]));
			}

			//cands:  type - direction - bailiff - provost
			// 2 - sPos backward to left edge: leftEdgePos <- aPos
			if(rightEdgePos-leftEdgePos >= 2) { //size of block larger than 3 - else FORT shift is same
				for(unsigned aPos=leftEdgePos+1; aPos<=rightEdgePos; aPos++) {


					//cout << "\tleftEdgePos <- aPos" << endl;
					//cout << "\t\t: " << critic[leftEdgePos] << "<-" << critic[aPos] << endl;
					
					canShift = true;
					for(unsigned otherPos=leftEdgePos; otherPos<rightEdgePos; otherPos++) {
						if(inst.hasPrec(critic[otherPos] , critic[aPos])) {
							canShift = false;
							break;
						}
					}

					if(canShift)
						if(nonCycleN6Criterium(BACK, critic[leftEdgePos], critic[aPos], JOB, heads, tails))
							cands.push_back(make_tuple(JOB, BACK, critic[leftEdgePos], critic[aPos]));
				}//for(unsigned aPos=leftEdgePos+1; aPos<=rightEdgePos; aPos++)
			}//if(rightEdgePos-leftEdgePos >= 2)
		}//for(unsigned bbBaseInd=0;  bbBaseInd<jobBb.size()-1; bbBaseInd++)
#endif //VANILLA_PSS

		//MACHINE BLOCKS
		for(unsigned bbBaseInd=0;  bbBaseInd<machBb.size()-1; bbBaseInd++) {
			leftEdgePos = machBb[bbBaseInd];
			rightEdgePos = machBb[bbBaseInd+1]-1;
			assert(leftEdgePos <= rightEdgePos);
			assert(rightEdgePos <= critic.size());


			// 1 - foward to right edge: aPos -> rightEdgePos
		    for(unsigned aPos=leftEdgePos; aPos+1<=rightEdgePos; aPos++) {
				
				//cout << "\taPos -> rightEdgePos" << endl;
				//cout << "\t\t: " << critic[aPos] << "->" << critic[rightEdgePos] << endl;
				
				assert(leftEdgePos < rightEdgePos); //block is longer than 1
				if(nonCycleN6Criterium(FORTH, critic[rightEdgePos], critic[aPos], MACH, heads, tails))
					cands.push_back(make_tuple(MACH, FORTH, critic[rightEdgePos], critic[aPos]));
			}

			//cands:  direction - bailiff - provost
			// 2 - sPos backward to left edge: leftEdgePos <- aMidPos
			if(rightEdgePos-leftEdgePos >= 2) { //size of block larger than 3 - else FORT shift is same
				for(unsigned aPos=leftEdgePos+1; aPos<=rightEdgePos; aPos++) {

					//cout << "\tleftEdgePos <- aPos" << endl;
					//cout << "\t\t: " << critic[leftEdgePos] << "<-" << critic[aPos] << endl;
					
					if(nonCycleN6Criterium(BACK, critic[leftEdgePos], critic[aPos], MACH, heads, tails))
						cands.push_back(make_tuple(MACH, BACK, critic[leftEdgePos], critic[aPos]));
				}
			}
		}
		/*
		cout << "cands:";
		for(tuple<unsigned, unsigned, unsigned, unsigned> & t : cands) cout << (get<0>(t)==MACH ? "  MACH" : "  JOB") << (get<1>(t)==FORTH ? ",FORTH,b" : ",BACK,b") << get<2>(t) << ",p" << get<3>(t);
		cout << endl;
		getchar();
		*/
	}



	

	//Zhang 2007
	//cand: direction - bailiff - provost
	///XXX CHECK WHEN CANDS REPEAT -
	// -- instance block ize two all the same
	// -- other situations???
	//TODO SPECIAL CASES WHERE SIZE IS 2 or 3 WLL REPEAT CANDIDATES (FORTH AND BACK IN ADJ ARE SAME)
	void fillCandidatesN7(vector<tuple<unsigned, unsigned, unsigned>> & cands, vector<unsigned> & jobBb, vector<unsigned> & machBb, const vector<unsigned> & critic, const vector<unsigned> & heads, const vector<unsigned> & tails) const {
		assert(cands.capacity() == (inst.O)*(inst.O));
		assert(jobBb.capacity() == inst.O);
		assert(machBb.capacity() == inst.O);
		assert( ! critic.empty());

		cands.clear();

		blockBegins(jobBb, machBb, critic);

		unsigned leftEdgePos;
		unsigned rightEdgePos;
		unsigned leftMiddlePos;
		unsigned rightMiddlePos;
		unsigned blockSize;

		bool canShift;



		assert(jobBb.size() > 1); //or else in job lower bound. why still calling neighbourhhod still?
		assert(machBb.size() > 1);
		assert(jobBb[0] == 0);
		assert(machBb[0] == 0);
		assert(jobBb[1] > 0);
		assert(machBb[1] > 0);
		assert(jobBb[1] < critic.size());
		assert(machBb[1] < critic.size());

		//adding dummy final block for symmery in code
		jobBb.push_back(critic.size());
		machBb.push_back(critic.size());
#ifdef VANILLA_PSS
		//XXX TODO
		//JOB BLOCKS
		for(unsigned bbBaseInd=0;  bbBaseInd<jobBb.size()-1; bbBaseInd++) {
			leftEdgePos = machBb[bbBaseInd]-1;
			rightEdgePos = machBb[bbBaseInd+1];
			leftMiddlePos = leftEdgePos + 1;
			assert(rightEdgePos>0);
			rightMiddlePos = rightEdgePos -1;

			//cands:  direction - bailiff - provost
			// 1- left edge forward to middle: leftEdgePos -> aMidPos
			for(unsigned aMidPos=leftMiddlePos; aMidPos<=rightMiddlePos; aMidPos++) {
				if(inst.hasPrec(critic[leftEdgePos], critic[aMidPos]))
					break;

				if(nonCycleN6Criterium(FORTH, critic[aMidPos], critic[leftEdgePos], JOB, heads, tails))
					cands.push_back(make_tuple(FORTH, critic[aMidPos], critic[leftEdgePos]));
			}

			// 2 - right edge backward to middle: aMidPos <- rightEdgePos
			for(unsigned aMidPos=rightMiddlePos; aMidPos>=leftMiddlePos; aMidPos--) {			
				if(inst.hasPrec(critic[aMidPos], critic[rightEdgePos]))
					break;
				if(nonCycleN6Criterium(BACK, critic[aMidPos], critic[rightEdgePos], JOB, heads, tails))
					cands.push_back(make_tuple(BACK, critic[aMidPos], critic[rightEdgePos]));
			}

			// 3 - midle backward to left edge: leftEdgePos <- aMidPos
			for(unsigned aMidPos=leftMiddlePos; aMidPos<=rightMiddlePos; aMidPos++) {
				canShift = true;
				for(unsigned o=aMidPos-1; o>=leftEdgePos; o--) {
					if(inst.hasPrec(critic[o] , critic[aMidPos])) {
						canShift = false;
						break;
					}
				}
					
				if(canShift)
					if(nonCycleN6Criterium(BACK, critic[leftEdgePos], critic[aMidPos], JOB, heads, tails))
						cands.push_back(make_tuple(BACK, critic[leftEdgePos], critic[aMidPos]));
			}

			// 4 - midle foward to right edge: aMidPos -> rightEdgePos
		    for(unsigned aMidPos=leftMiddlePos; aMidPos<=rightMiddlePos; aMidPos++) {
				canShift = true;
				for(unsigned o=aMidPos+1; o<=rightEdgePos; o++) {
					if(inst.hasPrec(critic[aMidPos] , critic[o])) {
						canShift = false;
						break;
					}
				}

				if(canShift)
					if(nonCycleN6Criterium(FORTH, rightEdgePos, aMidPos, JOB, heads, tails))
						cands.push_back(make_tuple(FORTH, rightEdgePos, aMidPos));
			}
		}
#endif //VANILLA_PSS

		//MACHINE BLOCKS
		for(unsigned bbBaseInd=0;  bbBaseInd<machBb.size()-1; bbBaseInd++) {
			leftEdgePos = machBb[bbBaseInd];
			rightEdgePos = machBb[bbBaseInd+1]-1;
			leftMiddlePos = leftEdgePos + 1;
			assert(rightEdgePos>0);
			rightMiddlePos = rightEdgePos -1;

			assert(rightEdgePos >= leftEdgePos);
			blockSize = rightEdgePos - leftEdgePos + 1;

			if(blockSize == 2) {
				assert(rightEdgePos == leftEdgePos+1);
				if(nonCycleN6Criterium(FORTH, critic[rightEdgePos], critic[leftEdgePos], MACH, heads, tails))
					cands.push_back(make_tuple(FORTH, critic[rightEdgePos], critic[leftEdgePos]));
			}

			//cands:  direction - bailiff - provost
			for(unsigned aMidPos=leftMiddlePos; aMidPos<=rightMiddlePos; aMidPos++) {
				assert(blockSize >=3);
				//left edge forward to middle: leftEdgePos -> aMidPos
				if(nonCycleN6Criterium(FORTH, critic[aMidPos], critic[leftEdgePos], MACH, heads, tails))
					cands.push_back(make_tuple(FORTH, critic[aMidPos], critic[leftEdgePos]));

				//right edge backward to middle: aMidPos <- rightEdgePos
				if(nonCycleN6Criterium(BACK, critic[aMidPos], critic[rightEdgePos], MACH, heads, tails))
					cands.push_back(make_tuple(BACK, critic[aMidPos], critic[rightEdgePos]));

				//midle backward to left edge: leftEdgePos <- aMidPos
				if(nonCycleN6Criterium(BACK, critic[leftEdgePos], critic[aMidPos], MACH, heads, tails))
					cands.push_back(make_tuple(BACK, critic[leftEdgePos], critic[aMidPos]));

				//midle foward to right edge: aMidPos -> rightEdgePos
				if(nonCycleN6Criterium(FORTH, critic[rightEdgePos], critic[aMidPos], MACH, heads, tails))
					cands.push_back(make_tuple(FORTH, critic[rightEdgePos], critic[aMidPos]));
			}
		}
	}
	


	//N5 neighbourhood first improvement - randomizes
	void localSearch( vector<unsigned> & dists, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, vector<unsigned> & jobBb, vector<unsigned> & machBb, vector<pair<unsigned, unsigned>> & cands, vector<unsigned> & critic, unsigned lowerBound) {
#ifndef NDEBUG
		bool cycle;
		assert(dists.size() == inst.O);
		assert(prev.size() == inst.O);
		assert(prev.size() == inst.O);
		assert(Q.size() == inst.O);
		assert(cands.capacity() == inst.O);
		assert(critic.capacity() == inst.O);
		assert(jobBb.capacity() == inst.O);
		assert(machBb.capacity() == inst.O);
#endif
		unsigned lastOp;
		unsigned thisMakes;

		setMeta(dists, lastOp, prev, indeg, Q);

		while(true) {
			assert(makes >= lowerBound);
			if(makes == lowerBound)
				break;
			thisMakes = makes;
			computeCritic(critic, lastOp, prev);
			assert( ! critic.empty());
			assert(critic.size() < inst.O);
			fillCandidatesN5(cands, jobBb, machBb, critic);
			random_shuffle(cands.begin(), cands.end());

			for(const pair<unsigned, unsigned> & p : cands) {
				swap(p.first, p.second);
#ifndef NDEBUG
				cycle = 
#endif
					setMeta(dists, lastOp, prev, indeg, Q);
				assert( ! cycle);

				if(thisMakes > makes) {
					break; //first improvement
				} else {
					swap(p.second, p.first); //undoing swap
				}
			}
			if(thisMakes <= makes)
				break;
		}
		setMeta(dists, lastOp, prev, indeg, Q);
		assert(makes == thisMakes);
	}




	//maxLen = maxD * maxC
	static bool detectRepeat(vector<int> & oldValues, unsigned & posCurMakes, unsigned & cycleLastPos, unsigned & cycleL, unsigned newMakes, bool isNewBest, unsigned maxD, unsigned maxLen) {
		assert(oldValues.size() == maxD);
		assert(maxD != 0);


		// new best solution: reset cycle detector
		if (isNewBest) { 
			posCurMakes = 0;
			cycleLastPos = 0;
			cycleL = 0;
			for (unsigned u=0; u<maxD; u++)
				oldValues[u]=-1; 
		}
		// store current makespan in ring buffer

		++posCurMakes;
		posCurMakes = posCurMakes % maxD; 
		oldValues[posCurMakes]=newMakes; //posCurMakes is latest makes pos
		//fprintf(stderr, "cur makes: %d.\n", newMakes);
		// if current makespan occurred earlier
		if (cycleLastPos) {
			// if current make is also the same
			if (!(oldValues[(posCurMakes+cycleLastPos)%maxD]-newMakes)) {//maxD is lenght of ring buffer
				// if this has been the case for maxc periods: we have cycle
				if (++cycleL >= maxLen) {//cycleL is number of times it was already found to be equal
					//printf(    "oldValues d [%3d] (%8ld)   %6.2f'\n",cycleLastPos,iter,showtime(timep));
					//getchar();
					return true;
				}
				// otherwise: keep watching
				return false;// it is false but is finidng a cycle just not big enough yet
			}
		}
		// not detecting a cycle, or current cycle is not a cycle: reset, find last occurrence of current makespan
		cycleL=cycleLastPos=0; 
		for (unsigned u=1; u<maxD; u++) {
			if (!(oldValues[(posCurMakes+u)%maxD]-newMakes)) 
				cycleLastPos=u; //cycleLastPos is position of the of last makespan - if cycleLastPos is 0 then it wasnt in list
		}
		return false; 
	}

	

	

	//does not set critic but do set makes
	//@return: cycle?
	bool setMeta(vector<unsigned> & dists, unsigned & lastOp, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q)  {
		assert(isAlloced());
		assert(dists.size() == inst.O);
		assert(prev.size() == inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);

		unsigned qInsert = 0;
		unsigned qAccess = 0;

		unsigned curOp;
		unsigned newOp;

		unsigned newMax;

		makes = 0;

		fill(dists.begin(), dists.end(), 0);
		fill(indeg.begin(), indeg.end(), 0);

		for(unsigned o=1; o<inst.O; o++) {
			if(_job[o] != 0)
				indeg[o]++;
			if(_mach[o] != 0)
				indeg[o]++;
			if(indeg[o] == 0) {
				prev[o] = 0;
				Q[qInsert++] = o;
			}
		}

		//assert(qInsert>0);

		while(qAccess < qInsert) {
			assert(qAccess<Q.size());
			curOp = Q[qAccess++];

			assert(indeg[curOp] == 0);

			newMax = dists[curOp] + inst.P[curOp];
			if(makes < newMax) {
				makes = newMax;
				lastOp = curOp;
			}

			//from JOB
			newOp = job[curOp];
			if(newOp != 0) {
				assert(indeg[newOp] > 0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
				if(dists[newOp] < newMax) {
					dists[newOp] = newMax;	
					prev[newOp] = curOp;
				}			
			}

			//from MACH
			newOp = mach[curOp];
			if(newOp != 0) {
				assert(indeg[newOp] > 0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
				if(dists[newOp] < newMax) {
					dists[newOp] = newMax; 
					prev[newOp] = curOp;
				}			
			}
		}
		inst.calcPenalties(dists, ePenalty,lPenalty);
		penalties = ePenalty + lPenalty;
		return qAccess<inst.O-1;
	}

    //##############################################################
	//HALF MAKES functions vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
	//##############################################################
#ifdef HALF_MAKES
	//problematic - order in queue is not ok with two swaps
	//CONST DISTS - put dirty dists
	unsigned halfMakes(const vector<unsigned> & heads, vector<unsigned> & dists, const vector<unsigned> & invertTable, vector<unsigned> & Q, unsigned o1, unsigned o2) const {
		assert(isAlloced());
		assert(dists.size() == inst.O);
		assert(invertTable.size() == inst.O);
		assert(Q.size() == inst.O);
		assert(invertTable[o1] < invertTable[o2]);
		assert(o1!=0);
		assert(o2!=0);

		assert(Q[invertTable[o1]] == o1);
		assert(Q[invertTable[o2]] == o2);


		unsigned curOp;
		unsigned fromJob;
		unsigned fromMach;
		unsigned foundMakes = 0;

		vector<unsigned> fixOps = {o2, o1};

		copy(heads.begin(), heads.end(), dists.begin());

		for(unsigned o : fixOps) {
			if(_job[o] != 0)
				fromJob = dists[_job[o]] + inst.P[_job[o]];
			else
				fromJob = 0;
			if(_mach[o] != 0)
				fromMach = dists[_mach[o]] + inst.P[_mach[o]];
			else
				fromMach = 0;
			dists[o] = max(fromJob, fromMach);
		}

		Q[invertTable[o1]] = o2;
		Q[invertTable[o2]] = o1;

		for(unsigned index=invertTable[o1]+1; index<inst.O; index++) {
			curOp = Q[index];

			//cout << "\t\tcurOp: " << curOp << endl;

			if(_job[curOp] != 0)
				fromJob = dists[_job[curOp]] + inst.P[_job[curOp]];
			else
				fromJob = 0;

			if(_mach[curOp] != 0)
				fromMach = dists[_mach[curOp]] + inst.P[_mach[curOp]];
			else
				fromMach = 0;

			dists[curOp] = max(fromJob, fromMach);

			//cout << "\t\t\tdists: " << dists[curOp] << endl;
		}

		//update makes from leafs
		for(unsigned l : inst.leafs)
			foundMakes  = max(foundMakes, dists[l] + inst.P[l]);

		Q[invertTable[o1]] = o1;
		Q[invertTable[o2]] = o2;


#ifndef NDEBUG
		State dummyState = *this;
		vector<unsigned> dummyDists(inst.O);
		unsigned dummyLastOp(inst.O);
		vector<unsigned> dummyPrev(inst.O);
		vector<unsigned> dummyIdeg(inst.O);
		vector<unsigned> dummyQ(inst.O);
		dummyState.setMeta(dummyDists, dummyLastOp, dummyPrev, dummyIdeg, dummyQ);
	    for(unsigned u=0; u<inst.O; u++) {
			assert(dummyDists[u] == dists[u]);
		}
		assert(dummyState.makes == foundMakes);
#endif //NDEBUG

		return foundMakes;
	}
	


	//normal computing of heads but stores invert table for Q
	bool halfHeadsComplete(vector<unsigned> & heads, vector<unsigned> & invertTable, vector<unsigned> & indeg, vector<unsigned> & Q) const {
		assert(isAlloced());
		assert(heads.size() == inst.O);
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);

		unsigned qInsert = 0;
		unsigned qAccess = 0;

		unsigned curOp;
		unsigned newOp;

		unsigned newMax;

		fill(heads.begin(), heads.end(), 0);
		fill(indeg.begin(), indeg.end(), 0);

		for(unsigned o=1; o<inst.O; o++) {
			if(_job[o] != 0)
				indeg[o]++;
			if(_mach[o] != 0)
				indeg[o]++;
			if(indeg[o] == 0)
				Q[qInsert++] = o;
		}

		assert(qInsert>0);

		while(qAccess < qInsert) {
			assert(qAccess<Q.size());
			curOp = Q[qAccess++];
			invertTable[curOp] = qAccess-1;

			assert(indeg[curOp] == 0);
			newMax = heads[curOp] + inst.P[curOp];

			//job order
			newOp = job[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
				heads[newOp] = max(heads[newOp], newMax);
			}
			//mach order
			newOp = mach[curOp];
			if(newOp != 0) {
				assert(indeg[newOp]>0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
				heads[newOp] = max(heads[newOp], newMax);
			}
		}
	
		assert(qAccess<=inst.O-1);

#ifndef NDEBUG
		State dummyState = *this;
		vector<unsigned> dummyDists(inst.O);
		unsigned dummyLastOp(inst.O);
		vector<unsigned> dummyPrev(inst.O);
		vector<unsigned> dummyIdeg(inst.O);
		vector<unsigned> dummyQ(inst.O);
		dummyState.setMeta(dummyDists, dummyLastOp, dummyPrev, dummyIdeg, dummyQ);
		for(unsigned u=0; u<inst.O; u++)
			assert(dummyDists[u] == heads[u]);
#endif //NDEBUG

		return qAccess<inst.O-1;
	}
#endif //HALF_MAKES
	//##############################################################
	//HALF MAKES functions ^^^^^^^^^^^^^^^^^^^^^^^
	//##############################################################


	//computes the new makes after swap, expects dists the be the value before the swap and the swap to be performed in DG
	unsigned swapMakes(const vector<unsigned> & dists, vector<unsigned> & dirtyDists, vector<bool> & dirty, unsigned o1, unsigned o2) const {
		assert(isAlloced());
		assert(dists.size() == inst.O);
		assert(dirty.size() == inst.O);
		assert((mach[o2] == o1 && _mach[o1] == o2)    ||    (job[o2] == o1 && _job[o1] == o2));

		unsigned curOp;

		queue<unsigned> Q;

		unsigned newDist;

		bool updated;

		unsigned distFromJ;
		unsigned distFromM;

		//XXX
		//cout << "\n\no1: " << o1 << "  ;  o2: " << o2 << endl;
		//cout << "start dists: ";
		//for(unsigned o=1; o<inst.O; o++)
		//	 cout << o << ":" << dists[o] << "[" << inst.P[o] << "]   ";
		//cout << endl;
		//XXX

		fill(dirty.begin(), dirty.end(), false);

		//calculating new dist for  o2
		if(_job[o2] != 0)
			distFromJ = dists[_job[o2]]+inst.P[_job[o2]];
		else
			distFromJ = 0;
		if(_mach[o2] != 0)
			distFromM = dists[_mach[o2]]+inst.P[_mach[o2]];
		else
			distFromM = 0;
		newDist = max(distFromJ, distFromM);
		updated = false;
		//update o2
		if(newDist != dists[o2]) {
			dirtyDists[o2] = newDist;
			dirty[o2] = true;
		}
		//update o1 
		if( ! updated) {
			//calculating new dist for  o1
			if(_job[o1] != 0)
				distFromJ =  (dirty[_job[o1]] ? dirtyDists[_job[o1]] : dists[_job[o1]])+inst.P[_job[o1]];
			else
				distFromJ = 0;
			if(_mach[o1] != 0)
				distFromM =   (dirty[_mach[o1]] ? dirtyDists[_mach[o1]] : dists[_mach[o1]])+inst.P[_mach[o1]];
			else
				distFromM = 0;
			newDist = max(distFromJ, distFromM);
			//update o1
			if(newDist != dists[o1]) {
				dirtyDists[o1] = newDist;
				dirty[o1] = true;
			}
		}

		//inserting possible propagations of change
		if(job[o2] != 0   &&   job[o2] != o1)
			Q.push(job[o2]);
		if(mach[o2] != 0    &&   mach[o2] != o1)
			Q.push(mach[o2]);
		if(job[o1] != 0)
			Q.push(job[o1]);
		if(mach[o1] != 0)
			Q.push(mach[o1]);

		 

		while( ! Q.empty()) {
			assert(Q.size() < (inst.O*inst.O));
			curOp = Q.front();
			Q.pop();
			updated = false;

			//Computing dist for curOp
			if(_job[curOp] != 0)
				distFromJ = (dirty[_job[curOp]] ? dirtyDists[_job[curOp]] : dists[_job[curOp]])+inst.P[_job[curOp]];
			else
				distFromJ = 0;
			if(_mach[curOp] != 0)
				distFromM = (dirty[_mach[curOp]] ? dirtyDists[_mach[curOp]] : dists[_mach[curOp]])+inst.P[_mach[curOp]];
			else
				distFromM = 0;
			newDist = max(distFromJ, distFromM);

			//cout << "\t\tdistFromJ: " << distFromJ << "  ;  distFromM: " << distFromM <<"  ;  newDist: " << newDist << endl; 


			//is curOp dirty?
			if(! dirty[curOp]) {
				if(newDist != dists[curOp]) {
					dirtyDists[curOp] = newDist;
					dirty[curOp] = true;
					updated = true;
				}
			} else {
				if(newDist != dirtyDists[curOp]) {
					dirtyDists[curOp] = newDist;
					updated = true;
				}
			}
			 

			if(updated) {
				if(job[curOp] != 0)
					Q.push(job[curOp]);
				if(mach[curOp] != 0)
					Q.push(mach[curOp]);
			}
		}
			 
		//update makes from leafs
		for(unsigned l : inst.leafs)
			newDist = max(newDist, (dirty[l] ? dirtyDists[l]+inst.P[l] : dists[l]+inst.P[l]));

		// cout << "FINISHE DELTA SWAP" << endl;
		 
#ifndef NDEBUG
		State dummyState = *this;
		vector<unsigned> dummyDists(inst.O);
		unsigned lastOp;
		vector<unsigned> prev(inst.O);
		vector<unsigned> indeg(inst.O);
		vector<unsigned> dummyQ(inst.O);
		bool dummyCycle = dummyState.setMeta(dummyDists, lastOp, prev, indeg, dummyQ);
		assert( ! dummyCycle);
		//XXX
		if(newDist != dummyState.makes) {
			cout << "newDist: " << newDist << "  ;  dummyState.makes: " << dummyState.makes << endl;
			cout << easyString() << endl;
			cout << "basic: ";
			for(unsigned o=1; o<inst.O; o++)
				cout << o << ":" << dummyDists[o] << "[" << inst.P[o] << "]   ";
			cout << endl;
			cout << "delta: ";
			for(unsigned o=1; o<inst.O; o++)
				cout << o << ":" << (dirty[o] ? dirtyDists[o] : dists[o]) << "[" << inst.P[o] << "]   ";
			cout << endl;
		}
		//XXX
		assert(newDist == dummyState.makes);
#endif

		//cout << "FINISHED ASSERTING" << endl;

		return newDist;
	}




		//@return: starts
	vector<unsigned> genSchedule()  {

		vector<unsigned> indeg(inst.O, 0);
		vector<unsigned> Q(inst.O, 0);
		vector<unsigned> starts(inst.O, 0);
		vector<unsigned> prev(inst.O, 0);
		unsigned lastOp;

#ifndef NDEBUG
		bool cycle = 
#endif
			setMeta(starts, lastOp, prev, indeg, Q);
		assert(! cycle);

		return starts;
	}



	bool verifySchedule() {

		unsigned testMakes = makes;//genSchedule may change makes - store here to be sure
		vector<unsigned> schedule = genSchedule();

		return inst.verifySchedule(schedule, testMakes);
	}

	bool printPenaltys(){

		vector<unsigned> schedule = genSchedule();
		cout << makes << " ";

		inst.printPenaltys(schedule,makes);
		
	}

	


#ifndef NDEBUG
		bool isAlloced() const {
		assert(inst.O != 0);
		if(mach.size() != inst.O) return false;
		if(job.size() != inst.O) return false;

		return true;
	}



		bool isClear() const {
		assert(inst.O != 0);
		assert(isAlloced());
		//if( ! roots.empty()) return false;
		for(unsigned o : job)
			if(o != 0) return false;

		for(unsigned o : _job)
			if(o != 0) return false;

		for(unsigned o : mach) 
			if(o != 0) return false;

		for(unsigned o : _mach) 
			if(o != 0) return false;

		return true;
	}



	//not all opers in state yet
	bool partialVerify(const vector<bool> & inOpers) const {
		assert(isAlloced());
		assert( ! inOpers[0]);
		vector<unsigned> jobLeafs(inst.J, 0);
		vector<unsigned> machLeafs(inst.M, 0);
		vector<unsigned> jobRoots(inst.J, 0);
		vector<unsigned> machRoots(inst.M, 0);
		unsigned j;
		unsigned m;

		if(job[0] != 0    ||   _job[0] !=0) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "job[0] != 0    ||   _job[0] !=0" << endl;
			return false;
		}
		if(mach[0] != 0    ||   _mach[0] !=0) {
			cout << "\n\t\tERROR !\n" << toString() << endl << "mach[0] != 0    ||   _mach[0] !=0" << endl;
			return false;
		}

		for(unsigned o=1; o<inst.O; o++) {
			j = inst.operToJ[o];
			m = inst.operToM[o];
#ifdef VANILLA_PSS
			//not in must have no next or prev
			if( ! inOpers[o]) {
				if(job[o] != 0   ||    _job[o] != 0) {
					cout << "\n\t\tERROR !\n" << toString() << endl << "job[o] != 0   ||    _job[o] != 0" << " ; o:" << o << endl;
					return false;
				}

				if(mach[o] != 0   ||    _mach[o] != 0) {
					cout << "\n\t\tERROR !\n" << toString() << endl << "mach[o] != 0   ||    _mach[o] != 0" << " ; o:" << o << endl;
					return false;
				}
			}
#endif //#ifdef VANILLA_PSS
			//one leaf per job
			if(inOpers[o]   &&   job[o]==0) {
				jobLeafs[j]++;
				if(jobLeafs[j] > 1) {
					cout << "\n\t\tERROR !\n" << toString() << endl << "jobLeafs[inst.operToJ[o]] > 1" << " ; o:" << o << " ; j:" << j << endl;
					return false;
				}
			}
			//one root per job
			if(inOpers[o]   &&   _job[o]==0) {
				jobRoots[j]++;
				if(jobRoots[j] > 1) {
					cout << "\n\t\tERROR !\n" << toString() << endl << "jobRoots[inst.operToJ[o]] > 1" << " ; o:" << o << " ; j:" << j << endl;
					return false;
				}
			}
			//one leaf per mach
			if(inOpers[o]   &&   mach[o]==0) {
				machLeafs[m]++;
				if(machLeafs[m] > 1) {
					cout << "\n\t\tERROR !\n" << toString() << endl << "machLeafs[m] > 1" << " ; o:" << o << " ; m:" << m << endl;
					return false;
				}
			}
			//one root per mach
			if(inOpers[o]   &&   _mach[o]==0) {
				machRoots[m]++;
				if(machRoots[m] > 1) {
					cout << "\n\t\tERROR !\n" << toString() << endl << "machRoots[m] > 1" << " ; o:" << o << " ; m:" << m << endl;
					return false;
				}
			}

			//opers in range
			if(job[o] >= inst.O) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "job[o] >= inst.O" << " ; o:" << o << endl;
				return false;
			}
			if(_job[o]  >= inst.O) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "_job[o]  >= inst.O" << " ; o:" << o << endl;
				return false;
			}
			if(mach[o] >= inst.O) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "mach[o] >= inst.O" << " ; o:" << o << endl;
				return false;
			}
			if(_mach[o] >= inst.O) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "_mach[o] >= inst.O" << " ; o:" << o << endl;
				return false;
			}

			//same job
			if(job[o]!=0   &&   inst.operToJ[o] != inst.operToJ[job[o]]) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "job[o]!=0   &&   inst.operToJ[o] != inst.operToJ[job[o]]" << " ; o:" << o << endl;
				return false;
			}
			if(_job[o]!=0   &&   inst.operToJ[o] != inst.operToJ[_job[o]]) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "_job[o]!=0   &&   inst.operToJ[o] != inst.operToJ[_job[o]]" << " ; o:" << o << endl;
				return false;
			}

			//same mach
			if(mach[o]!=0   &&   inst.operToM[o] != inst.operToM[mach[o]]) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "mach[o]!=0   &&   inst.operToM[o] != inst.operToM[mach[o]]" << " ; o:" << o << endl;
				return false;
			}
			if(_mach[o]!=0   &&   inst.operToM[o] != inst.operToM[_mach[o]]) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "_mach[o]!=0   &&   inst.operToM[o] != inst.operToM[_mach[o]]" << " ; o:" << o << endl;
				return false;
			}

			//check precs
			if(job[o]!=0   &&    inst.hasPrec(job[o], o)) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "job[o]!=0   &&    inst.hasPrec(job[o], o)" << " ; o:" << o << endl;
				return false;
			}
			if(_job[o]!=0   &&    inst.hasPrec(o, _job[o])) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "_job[o]!=0   &&    inst.hasPrec(o, _job[o])" << " ; o:" << o << endl;
				return false;
			}

			//invert compatible
			if(_job[o]!=0   &&   job[_job[o]] != o) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "_job[o]!=0   &&   ! job[_job[o]] == o" << " ; o:" << o << endl;
				return false;
			}
			if(job[o]!=0   &&   _job[job[o]] != o) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "job[o]!=0   &&   ! _job[job[o]] == o" << " ; o:" << o << endl;
				return false;
			}
			if(_mach[o]!=0   &&   mach[_mach[o]] != o) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "_mach[o]!=0   &&   ! mach[_mach[o]] == o" << " ; o:" << o << endl;
				return false;
			}
			if(mach[o]!=0   &&   _mach[mach[o]] != o) {
				cout << "\n\t\tERROR !\n" << toString() << endl << "mach[o]!=0   &&   ! _mach[mach[o]] == o" << " ; o:" << o << endl;
				return false;
			}

			for(unsigned provost=o+1; provost<inst.O; provost++) {
				//no repeat of opers
				if(job[o]!=0   &&   job[o] == job[provost]) {
					cout << "\n\t\tERROR !\n" << toString() << endl << "job[o]!=0   &&   job[o] == job[provost]" << " ; o:" << o << " ; provost:" << provost << endl;
					return false;
				}
				if(_job[o]!=0   &&   _job[o] == _job[provost]) {
					cout << "\n\t\tERROR !\n" << toString() << endl << "_job[o]!=0   &&   _job[o] == _job[provost]" << " ; o:" << o << " ; provost:" << provost << endl;
					return false;
				}
				if(mach[o]!=0   &&   mach[o] == mach[provost]) {
					cout << "\n\t\tERROR !\n" << toString() << endl << "mach[o]!=0   &&   mach[o] == mach[provost]" << " ; o:" << o << " ; provost:" << provost << endl;
					return false;
				}
				if(_mach[o]!=0   &&   _mach[o] == _mach[provost]) {
					cout << "\n\t\tERROR !\n" << toString() << endl << "_mach[o]!=0   &&   _mach[o] == _mach[provost]" << " ; o:" << o << " ; provost:" << provost << endl;
					return false;
				}
			}
		}

		return true;
	}



		bool verify() const {
		vector<bool> inOpers(inst.O, true);
		inOpers[0] = false;
		return partialVerify(inOpers);
	}
#endif //NDEBUG


#ifdef PARALLEL_MACHS

	//expects o1 to be before o2
	void shiftForth(unsigned o1, unsigned o2) {
		assert(inst.operToM[o1] == inst.operToM[o2]);

		unsigned next = mach[o1];
		
		while(next != o2) {
			assert(next != 0);
			swap(o1, next);
			assert(next != mach[o1]);
			next = mach[o1];
		}
		
		swap(o1, o2);
	}




		//expects o1 to be before o2
		void shiftBack(unsigned o1, unsigned o2) {
			assert(o1 != o2);
			assert(o1 != 0);
			assert(o2 != 0);
			assert(inst.operToM[o1] == inst.operToM[o2]);

			unsigned prev = _mach[o2];
		
			while(prev != o1) {
				assert(prev != 0);
				swap(prev, o2);
				assert(prev != _mach[o2]);
				prev = _mach[o2];
			}
		
			swap(o1, o2);
		}


	


		//@return: cycle?
		//bucketUsage[x].first: oper x start --- bucketUsage[x].second: operx used this bucket
		bool findBucketUsage(vector<pair<unsigned, unsigned>> & bucketUsage, vector<unsigned> & dists, unsigned & lastOp, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & buckets)  {
			assert(isAlloced());
			assert(dists.size() == inst.O);
			assert(prev.size() == inst.O);
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(buckets.shape()[0] == inst.M);
			assert(buckets.shape()[1] == inst.K);


			unsigned qInsert = 0;
			unsigned qAccess = 0;

			unsigned curOp;
			unsigned newOp;

			unsigned newMax;

			unsigned chosenBucket;
			unsigned bucketHead;
			unsigned jobHead;
			bool foundNonDelayBucket;

			makes = 0;

			fill(dists.begin(), dists.end(), 0);
			fill(indeg.begin(), indeg.end(), 0);
			fill(buckets.data(),buckets.data()+buckets.num_elements(), 0);

			for(unsigned o=1; o<inst.O; o++) {
				if(_job[o] != 0)
					indeg[o]++;
				if(_mach[o] != 0)
					indeg[o]++;
				if(indeg[o] == 0) {
					prev[o] = 0;
					Q[qInsert++] = o;
				}
			}

			assert(qInsert>0);

			while(qAccess < qInsert) {
				assert(qAccess<Q.size());
				curOp = Q[qAccess++];

				cout << "\t\tcurOp: " << curOp << endl;

				assert(indeg[curOp] == 0);

				newMax = dists[curOp] + inst.P[curOp];
				if(makes < newMax) {
					makes = newMax;
					lastOp = curOp;
				}

				//from JOB
				newOp = job[curOp];
				if(newOp != 0) {
					assert(indeg[newOp] > 0);
					indeg[newOp]--;
					if(indeg[newOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = newOp;
					}
					if(dists[newOp] < newMax) {
						dists[newOp] = newMax;	
						prev[newOp] = curOp;
					}			
				}

				//put it in biggest bin that still generates the lowest machine head:
				//if exists at least one bin that is lower the the jobHead of curOp put it in the biggest bin ammong it
				//if all bins are higher the job head put it in the lowest one
				jobHead = dists[_job[curOp]] + inst.P[_job[curOp]];
				assert(jobHead <= dists[curOp]);

				//getting biggest bucket that doesnt delay curOp
				foundNonDelayBucket = false;
				for(unsigned bin=0; bin<inst.K; bin++) {
					if(buckets[inst.operToM[curOp]][bin] < jobHead) {
						if( ! foundNonDelayBucket) {
							chosenBucket = bin;
							foundNonDelayBucket = true;
						} else {
							if(buckets[inst.operToM[curOp]][bin] > buckets[inst.operToM[curOp]][chosenBucket])
								chosenBucket = bin; //bins bucket still smaller the jobHead but bigger then prev chosen
						}
					}
				}

				// all buckets delay curOp? -> get smallest one then
				if( ! foundNonDelayBucket) {
					chosenBucket = 0;
					for(unsigned bin=1; bin<inst.K; bin++) {
						if(buckets[inst.operToM[curOp]][bin] < buckets[inst.operToM[curOp]][chosenBucket])
							chosenBucket = bin;
					}
				}

				//filling chosen bucket
				buckets[inst.operToM[curOp]][chosenBucket] = newMax;
				bucketUsage[curOp].first = dists[curOp];
				bucketUsage[curOp].second = chosenBucket;

				//from MACH
				newOp = mach[curOp];
				if(newOp != 0) {
					assert(indeg[newOp] > 0);
					indeg[newOp]--;
					if(indeg[newOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = newOp;
					}
					/*
					//put it in biggest bin that still generates the lowest machine head:
					//if exists at least one bin that is lower the the jobHead of curOp put it in the biggest bin ammong it
					//if all bins are higher the job head put it in the lowest one
					jobHead = dists[_job[curOp]] + inst.P[_job[curOp]];
					assert(jobHead <= dists[curOp]);

					//getting biggest bucket that doesnt delay curOp
					foundNonDelayBucket = false;
					for(unsigned bin=0; bin<inst.K; bin++) {
					if(buckets[inst.operToM[curOp]][bin] < jobHead) {
					if( ! foundNonDelayBucket) {
					chosenBucket = bin;
					foundNonDelayBucket = true;
					} else {
					if(buckets[inst.operToM[curOp]][bin] > buckets[inst.operToM[curOp]][chosenBucket])
					chosenBucket = bin; //bins bucket still smaller the jobHead but bigger then prev chosen
					}
					}
					}

					// all buckets delay curOp? -> get smallest one then
					if( ! foundNonDelayBucket) {
					chosenBucket = 0;
					for(unsigned bin=1; bin<inst.K; bin++) {
					if(buckets[inst.operToM[curOp]][bin] < buckets[inst.operToM[curOp]][chosenBucket])
					chosenBucket = bin;
					}
					}

					//filling chosen bucket
					buckets[inst.operToM[curOp]][chosenBucket] = newMax;
					bucketUsage[curOp].first = dists[curOp];
					bucketUsage[curOp].second = chosenBucket;
					*/

					//setting prev of next oper
					if(dists[newOp] < newMax) {
						//dists[newOp] = newMax; 
						prev[newOp] = curOp;
					}


					//setting heads of next in machine order
					bucketHead = buckets[inst.operToM[curOp]][0];
					for(unsigned bin=1; bin<inst.K; bin++)
						bucketHead = min(bucketHead, buckets[inst.operToM[curOp]][bin]);
					if(dists[newOp] < bucketHead) {
						dists[newOp] = bucketHead; 
						//prev[newOp] = curOp;
					}

				}
			}

			return qAccess<inst.O-1;
		}






		//commputes dists of each oper compressing in the first gap found
		//gap x stored in gapL[m][k][x] - first start - second finish
		//@return: found makes - UINT_MAX if cycle
		unsigned headParaCompress(vector<unsigned> & dists, vector<vector<vector<pair<unsigned, unsigned>>>> & allGaps, unsigned & lastOp, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q) const {
			assert(isAlloced());
			assert(dists.size() == inst.O);
			assert(prev.size() == inst.O);
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(allGaps.size() == inst.M);

			unsigned qInsert = 0;
			unsigned qAccess = 0;

			unsigned curOp;
			unsigned jobHead;

			unsigned foundMakes = 0;

			for(unsigned m=0; m<inst.M; m++) {
				for(unsigned k=0; k<inst.K; k++) {
					allGaps[m][k] = {pair<unsigned, unsigned>(0, UINT_MAX-1)};
				}
			}
			fill(dists.begin(), dists.end(), 0);
			fill(indeg.begin(), indeg.end(), 0);


			for(unsigned o=1; o<inst.O; o++) {	
				if(_job[o] != 0)
					indeg[o]++;
				if(_mach[o] != 0)
					indeg[o]++;
				if(indeg[o] == 0) {
					prev[o] = 0;
					Q[qInsert++] = o;	
				}
			}




			while(qAccess < qInsert) {
				assert(qAccess<Q.size());
				curOp = Q[qAccess++];
				assert(curOp!=0);
				assert(indeg[curOp] == 0);

				//free next ops
				if(job[curOp] != 0) {
					assert(indeg[job[curOp]] > 0);
					indeg[job[curOp]]--;
					if(indeg[job[curOp]] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = job[curOp];
					}
				}
				if(mach[curOp] != 0) {
					assert(indeg[mach[curOp]] > 0);
					indeg[mach[curOp]]--;
					if(indeg[mach[curOp]] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = mach[curOp];
					}
				}

				jobHead = dists[_job[curOp]]+inst.P[_job[curOp]];

				//vector<vector<pair<unsigned, unsigned>>> & machGaps, unsigned jobHead, unsigned oper)
				dists[curOp] = insertInGap(allGaps[inst.operToM[curOp]], jobHead, curOp);

				assert(jobHead <= dists[curOp]);
				if(jobHead == dists[curOp])
					prev[curOp] = _job[curOp];
				else
					prev[curOp] = _mach[curOp];

				if(foundMakes < dists[curOp]+inst.P[curOp]) {
					foundMakes = dists[curOp]+inst.P[curOp];
					lastOp = curOp;
				}
			}

			assert(qAccess<=inst.O-1);
			if(qAccess<inst.O-1)
				return UINT_MAX;
			else 
				return foundMakes;
		}



		static string printBucket(const vector<vector<pair<unsigned, unsigned>>> & machGaps) {
			string str = "";

			for(unsigned k=0; k<inst.K; k++) {
				str.append("\tk: " + unsigStr(k)+" - ");
				for(unsigned x=0; x<machGaps[k].size(); x++) {
					str.append("["+unsigStr(machGaps[k][x].first)+","+unsigStr(machGaps[k][x].second)+"]    ");
				}
				str.append("\n");
			}

			return str;
		}




		//@return: head[oper]
		//updates gaps on machGaps
		static unsigned insertInGap(vector<vector<pair<unsigned, unsigned>>> & machGaps, unsigned jobHead, unsigned oper)  {
			vector<pair<unsigned, unsigned>>::iterator chosenIter = machGaps[0].end();
			vector<pair<unsigned, unsigned>>::iterator iter;
			unsigned chosenBucket = UINT_MAX;
			bool gapFound = false;
			bool nonDelayGapFound = false;
			unsigned operDist;
			unsigned prevFill;

			for(unsigned k=0; k<inst.K; k++) {
				assert(machGaps[k].back().second == UINT_MAX-1);

				//finding first possible bucket to insert
				iter = upper_bound (machGaps[k].begin(), machGaps[k].end(), pair<unsigned, unsigned>(jobHead, UINT_MAX));
				if(iter != machGaps[k].begin()) iter--;


				assert(iter==machGaps[k].begin()   ||   iter->first <= jobHead);
				assert(iter==machGaps[k].begin()   ||   (iter-1)->first < jobHead);
				while(iter != machGaps[k].end()) {


					if(nonDelayGapFound   &&   iter->first > jobHead) {
						assert(gapFound);


						break; //can stop search - got a non-delay and all next isnertions will be delay
					}
					if(gapFound    &&    iter->first > jobHead     &&   iter->first > chosenIter->first) {

						break; //can stop search - got a gap   &&   next will all be delay   &&   it will start after the chosen gap
					}
					//is end of gap enough if we put oper in here?
					if(iter->second >= max(iter->first, jobHead) + inst.P[oper]) { //found viable gap


						if(iter->first <= jobHead) { //is non delay?


							if( (! nonDelayGapFound)       ||    //first non delay? -> go 
								(chosenIter->second > iter->second)      ||    //there was another non delay -- choose lowest end --  (will preffer closed gaps)
								(chosenIter->second == iter->second   &&   chosenIter->first < iter->first)  ) {   //there was another non delay - it had the same end -> get the one with highest begin
								//iter is now chosen iter
								//cout << "\t\t\t\t\tChose it" << endl;
								chosenIter = iter;
								chosenBucket = k;
							}
							nonDelayGapFound = true;
						} else { //iter is a delay gap
							if( (! gapFound)   ||//first gap whatsoever always chosen
								((! nonDelayGapFound)    &&    (chosenIter->first > iter->first)) ) { //there was no nondelay then choose lowest begin
								chosenIter = iter;
								chosenBucket = k;
							}
						} 
						gapFound = true;
					}  //end if found viable gap
					iter++;
				} //end while iter != machGaps[k].end()
				assert(gapFound);
				assert(chosenBucket < inst.K);
				assert(chosenIter != machGaps[0].end());
			}// end for k<inst.K

			//apply found insertion
			if(nonDelayGapFound) {
				assert(chosenIter->first <= jobHead);
				if(chosenIter->second == UINT_MAX-1) { //open non-delay gap
					operDist = jobHead;
					prevFill = chosenIter->first;
					//new open gap begins from end of curOp
					chosenIter->first = jobHead + inst.P[oper];
					//maybe there is another before
					if(jobHead > prevFill)
						machGaps[chosenBucket].insert(chosenIter, pair<unsigned, unsigned>(prevFill, jobHead));
				} else { //closed non-delay gap
					operDist = chosenIter->second - inst.P[oper]; //glue to the top of the gap
					assert(operDist >= chosenIter->first);
					if(operDist > chosenIter->first) { //gap was not completelly filled
						chosenIter->second = operDist;
					} else { //gap was filled exacly - close it
						machGaps[chosenBucket].erase(chosenIter);
					}
				}
			} else { //only delay gap found
				assert(chosenIter->first > jobHead);
				operDist = chosenIter->first;
				if(operDist +inst.P[oper] < chosenIter->second) {//gap not filled completely
					chosenIter->first = operDist +inst.P[oper];
				} else {//gap filled completelly - remove it
					machGaps[chosenBucket].erase(chosenIter);
				}
			}
		
			assert(operDist >= jobHead);
			return operDist;
		}
	


		//@return: cycle?
		//XXX REMOVE THIS Sub with completeParaMeta and compMakesPara
		bool setMetaParallelBU(vector<unsigned> & dists, unsigned & lastOp, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & buckets)  {
			assert(isAlloced());
			assert(dists.size() == inst.O);
			assert(prev.size() == inst.O);
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(buckets.shape()[0] == inst.M);
			assert(buckets.shape()[1] == inst.K);


			unsigned qInsert = 0;
			unsigned qAccess = 0;

			unsigned curOp;
			unsigned newOp;

			unsigned newMax;

			unsigned chosenBucket;
			unsigned bucketHead;
			unsigned jobHead;
			bool foundNonDelayBucket;

			makes = 0;

			fill(dists.begin(), dists.end(), 0);
			fill(indeg.begin(), indeg.end(), 0);
			fill(buckets.data(),buckets.data()+buckets.num_elements(), 0);

			for(unsigned o=1; o<inst.O; o++) {
				if(_job[o] != 0)
					indeg[o]++;
				if(_mach[o] != 0)
					indeg[o]++;
				if(indeg[o] == 0) {
					prev[o] = 0;
					Q[qInsert++] = o;
				}
			}

			//assert(qInsert>0);

			while(qAccess < qInsert) {
				assert(qAccess<Q.size());
				curOp = Q[qAccess++];

				assert(indeg[curOp] == 0);

				newMax = dists[curOp] + inst.P[curOp];
				if(makes < newMax) {
					makes = newMax;
					lastOp = curOp;
				}

				//from JOB
				newOp = job[curOp];
				if(newOp != 0) {
					assert(indeg[newOp] > 0);
					indeg[newOp]--;
					if(indeg[newOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = newOp;
					}
					if(dists[newOp] < newMax) {
						dists[newOp] = newMax;	
						prev[newOp] = curOp;
					}			
				}

				//from MACH
				newOp = mach[curOp];
				if(newOp != 0) {
					assert(indeg[newOp] > 0);
					indeg[newOp]--;
					if(indeg[newOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = newOp;
					}
					//put it in biggest bin that still generates the lowest machine head:
					//if exists at least one bin that is lower the the jobHead of curOp put it in the biggest bin ammong it
					//if all bins are higher the job head put it in the lowest one
					jobHead = dists[_job[curOp]] + inst.P[_job[curOp]];
					assert(jobHead <= dists[curOp]);

					//getting biggest bucket that doesnt delay curOp
					foundNonDelayBucket = false;
					for(unsigned bin=0; bin<inst.K; bin++) {
						if(buckets[inst.operToM[curOp]][bin] < jobHead) {
							if( ! foundNonDelayBucket) {
								chosenBucket = bin;
								foundNonDelayBucket = true;
							} else {
								if(buckets[inst.operToM[curOp]][bin] > buckets[inst.operToM[curOp]][chosenBucket])
									chosenBucket = bin; //bins bucket still smaller the jobHead but bigger then prev chosen
							}
						}
					}

					// all buckets delay curOp? -> get smallest one then
					if( ! foundNonDelayBucket) {
						chosenBucket = 0;
						for(unsigned bin=1; bin<inst.K; bin++) {
							if(buckets[inst.operToM[curOp]][bin] < buckets[inst.operToM[curOp]][chosenBucket])
								chosenBucket = bin;
						}
					}

					//filling chosen bucket
					buckets[inst.operToM[curOp]][chosenBucket] = newMax;


					//setting prev of next oper
					if(dists[newOp] < newMax) {
						//dists[newOp] = newMax; 
						prev[newOp] = curOp;
					}


					//setting heads of next in machine order
					bucketHead = buckets[inst.operToM[curOp]][0];
					for(unsigned bin=1; bin<inst.K; bin++)
						bucketHead = min(bucketHead, buckets[inst.operToM[curOp]][bin]);
					if(dists[newOp] < bucketHead) {
						dists[newOp] = bucketHead; 
						//prev[newOp] = curOp;
					}

				}
			}

			return qAccess<inst.O-1;
		}





		//@return: makes - UINT_MAX if cycle
		//cands: all cand changes in critical path
		unsigned setMetaParallel(vector<unsigned> & dists, unsigned & lastOp, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & bucketLimits, multi_array<unsigned, 2> & bucketTops) const {
			assert(isAlloced());
			assert(dists.size() == inst.O);
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(bucketLimits.shape()[0] == inst.M);
			assert(bucketLimits.shape()[1] == inst.K);
			assert(bucketTops.shape()[0] == inst.M);
			assert(bucketTops.shape()[1] == inst.K);

			unsigned qInsert = 0;
			unsigned qAccess = 0;

			unsigned curOp;
			unsigned curMach;

			unsigned newMax;

			unsigned chosenBucket;
			unsigned bucketHead;
			unsigned jobHead;
			bool foundNonDelayBucket;

			unsigned foundMakes = 0;

			fill(dists.begin(), dists.end(), 0);
			fill(indeg.begin(), indeg.end(), 0);
			fill(bucketLimits.data(),bucketLimits.data()+bucketLimits.num_elements(), 0);
			fill(bucketTops.data(),bucketTops.data()+bucketTops.num_elements(), 0);

			for(unsigned o=1; o<inst.O; o++) {
				if(_job[o] != 0)
					indeg[o]++;
				if(_mach[o] != 0)
					indeg[o]++;
				if(indeg[o] == 0)
					Q[qInsert++] = o;
			}

			while(qAccess < qInsert) {
				//GET NEW CUR OP
				assert(qAccess<Q.size());
				curOp = Q[qAccess++];
				assert(indeg[curOp] == 0);

				//AUX FOR CUROP
				curMach = inst.operToM[curOp];
				jobHead = dists[_job[curOp]] + inst.P[_job[curOp]];

				//RELEASE NEXT OPS
				if(job[curOp] != 0) {
					assert(indeg[job[curOp]] > 0);
					indeg[job[curOp]]--;
					if(indeg[job[curOp]] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = job[curOp];
					}
				}
				if(mach[curOp] != 0) {
					assert(indeg[mach[curOp]] > 0);
					indeg[mach[curOp]]--;
					if(indeg[mach[curOp]] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = mach[curOp];
					}
				}

				//FIND PROPPER BUCKET
				//getting biggest bucket that doesnt delay curOp
				foundNonDelayBucket = false;
				for(unsigned bin=0; bin<inst.K; bin++) {
					if(bucketLimits[curMach][bin] < jobHead) {
						if( ! foundNonDelayBucket) {
							chosenBucket = bin;
							foundNonDelayBucket = true;
						} else {//other non delay bucket previously found
							if(bucketLimits[curMach][bin] > bucketLimits[curMach][chosenBucket])
								chosenBucket = bin; //bin's bucket still smaller the jobHead but bigger then prev chosen
						}
					}
				}
				// all buckets delay curOp? -> get smallest one then
				if( ! foundNonDelayBucket) {
					chosenBucket = 0;
					for(unsigned bin=1; bin<inst.K; bin++) {
						if(bucketLimits[curMach][bin] < bucketLimits[curMach][chosenBucket])
							chosenBucket = bin;
					}
				}
				bucketHead = bucketLimits[curMach][chosenBucket];

				//UPDATE PREV
				if(jobHead >  bucketHead) {
					prev[curOp] = _job[curOp];
					assert( foundNonDelayBucket);
				} else {//jobHead <=  bucketHead
					prev[curOp] = bucketTops[curMach][chosenBucket];
					assert(jobHead==bucketHead   ||   ! foundNonDelayBucket);
				}

				//UPDATE dists  - foundMakes - lastOps
				dists[curOp] = max(jobHead, bucketHead);
				newMax = dists[curOp] + inst.P[curOp];
				if(foundMakes < newMax) {
					foundMakes = newMax;
					lastOp = curOp;
				}

				//UPDATE BUCKET
				bucketLimits[curMach][chosenBucket] = newMax;
				bucketTops[curMach][chosenBucket] = curOp;
			}//while (qAccess < qInsert)

			if(qAccess<inst.O-1)
				return UINT_MAX;
			else 
				return foundMakes;
		}





		//parallel JSP insa
		void insaParallel(vector<unsigned> & indeg, vector<unsigned> & Q, vector<unsigned> & dists, multi_array<unsigned, 2> & buckets) {
			assert(isAlloced());
			assert(isClear());
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(dists.size() == inst.O);

			unsigned curOp;
			unsigned nextOp;

			unsigned curCost;
			unsigned bigJob = UINT_MAX;
			unsigned bigJobCost = 0;

			vector<unsigned> machRoot(inst.M, 0);

			vector<pair<unsigned, unsigned>> costVoper;
			costVoper.reserve(inst.O-inst.M);

			//linearizing jobs - expects JSP
			assert(inst.roots.size() == inst.J);
			for(unsigned r : inst.roots) {
				curOp = r;
				for(unsigned u=0; u<inst.M-1; u++) {
					assert(inst.next[curOp].size() == 1);
					nextOp = inst.next[curOp][0];

					assert(job[curOp] == 0);
					job[curOp] = nextOp;
					assert(_job[nextOp] == 0);
					_job[nextOp] = curOp;

					curOp = nextOp;
				}
			}

			//insert biggest job
		 
			//curCost = 0;
			for(unsigned j=0; j<inst.J; j++) {
				curCost = 0;
				for(unsigned o : inst.jobOpers[j])
					curCost += inst.P[o];
				if(curCost > bigJobCost) {
					bigJobCost = curCost;
					bigJob = j;
				}
			}
			assert(bigJob < inst.J);
		 

			for(unsigned o : inst.jobOpers[bigJob]) {
				assert(machRoot[inst.operToM[o]] == 0);
				machRoot[inst.operToM[o]] = o;
			}


			//ordering the rest of the opers
			for(unsigned o=1; o<inst.O; o++) {
				if(inst.operToJ[o] != bigJob)
					costVoper.push_back(pair<unsigned, unsigned>(inst.P[o], o));
			}
			sort(costVoper.rbegin(), costVoper.rend());

			//inserting missing ops
			for(unsigned pos=0; pos<costVoper.size(); pos++) {
				curOp = costVoper[pos].second;
				insertInsaParallel(curOp, dists, indeg, Q, buckets, machRoot);
			}
		}



		void insertInsaParallel(unsigned curOp, vector<unsigned> & dists, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & buckets, vector<unsigned> & machRoot) {
			unsigned curCost;

			unsigned bestCost;
			unsigned bestPrev;
			unsigned bestNext;
			const unsigned curMach = inst.operToM[curOp];


			assert(mach[curOp] == 0);
			assert(_mach[curOp] == 0);
			assert(machRoot[curMach] != 0);
			assert(_mach[machRoot[curMach]] == 0);

			mach[curOp] = machRoot[curMach];
			_mach[machRoot[curMach]] = curOp;

			bestCost = insaPosCostParallel(curOp, dists, indeg, Q, buckets);
			bestPrev = 0;
			bestNext = machRoot[curMach];

			while(true) {
				curCost = insaPosCostParallel(curOp, dists, indeg, Q, buckets);

				if(curCost < bestCost) {
					bestCost = curCost;
					bestPrev = _mach[curOp];
					bestNext = mach[curOp];		
				}

				if(mach[curOp] == 0)
					break;

				swap(curOp, mach[curOp]);
			}

			assert(_mach[curOp] != 0);
			assert(mach[_mach[curOp]] == curOp);

			//reseting leaf that sould pojnt to cur op at end of run
			mach[_mach[curOp]] = 0;

			assert(bestPrev==0     ||     mach[bestPrev] == bestNext);
			assert(bestNext==0     ||     _mach[bestNext] == bestPrev);

			_mach[curOp] = bestPrev;
			if(bestPrev != 0) mach[bestPrev] = curOp;
			mach[curOp] = bestNext;
			if(bestNext != 0) _mach[bestNext] = curOp;

			if(bestPrev==0)
				machRoot[curMach] = curOp;
		}



		unsigned insaPosCostParallel(unsigned targOp, vector<unsigned> & dists, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & buckets) const {
			unsigned c1;
			unsigned c2;
			c1 = maxDistFromTargParallel(targOp, dists, indeg, Q, buckets);
			if(c1 == UINT_MAX) return UINT_MAX;
			c2 =  maxDistToTargParallel(targOp, dists, indeg, Q, buckets);
			if(c2 == UINT_MAX) return UINT_MAX;
			return c1 + inst.P[targOp] + c2;
		}


		//dist from sink to targOp
		//returns UINT_MAX if has a cycle contaiing targOp
		//expects jobs linearized aready - will ignore precs
		unsigned maxDistFromTargParallel(unsigned targOp, vector<unsigned> & tails, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & buckets) const {
			assert(isAlloced());
			assert(tails.size() == inst.O);
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);

			unsigned qInsert = 0;
			unsigned qAccess = 0;

			unsigned curOp;
			unsigned newOp;

			unsigned newMax;

			unsigned chosenBucket;
			unsigned bucketHead;
			unsigned jobTail;
			unsigned curMach;
			bool foundNonDelayBucket;

			fill(tails.begin(), tails.end(), 0);
			fill(indeg.begin(), indeg.end(), 0);
			fill(buckets.data(),buckets.data()+buckets.num_elements(), 0);

			for(unsigned o=1; o<inst.O; o++) {
				//indeg[o] = inst.next[o].size();
				if(job[o] != 0)
					indeg[o]++;
				if(mach[o] != 0)
					indeg[o]++;
				if(indeg[o] == 0)
					Q[qInsert++] = o;
			}

			//assert(qInsert>0); //if o will skip while and return UINT_MAX as expected

			while(qAccess < qInsert) {
				assert(qAccess<Q.size());
				curOp = Q[qAccess++];
				assert(indeg[curOp] == 0);

				if(curOp == targOp)
					return tails[curOp];

				newMax = tails[curOp] + inst.P[curOp];

				//job order
				newOp = _job[curOp];
				if(newOp != 0) {
					assert(indeg[newOp]>0);
					indeg[newOp]--;
					if(indeg[newOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = newOp;
					}
					tails[newOp] = max(tails[newOp], newMax);
				}
				/*
				//job precs
				for(unsigned n : inst.prev[curOp]) {
				assert(indeg[n]>0);
				indeg[n]--;
				if(indeg[n] == 0) {
				assert(qInsert<Q.size()-1);
				Q[qInsert++] = n;
				}
				tails[n] = max(tails[n], newMax);
				}
				*/

				//from MACH
				newOp = _mach[curOp];
				curMach = inst.operToM[curOp];
				if(newOp != 0) {
					assert(indeg[newOp] > 0);
					indeg[newOp]--;
					if(indeg[newOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = newOp;
					}
					//put it in biggest bin that still generates the lowest machine head:
					//if exists at least one bin that is lower the the jobTail of curOp put it in the biggest bin ammong it
					//if all bins are higher the job head put it in the lowest one
					jobTail = tails[job[curOp]] + inst.P[job[curOp]];
					assert(jobTail <= tails[curOp]);

					//getting biggest bucket that doesnt delay curOp
					foundNonDelayBucket = false;
					for(unsigned bin=0; bin<inst.K; bin++) {
						if(buckets[curMach][bin] < jobTail) {
							if( ! foundNonDelayBucket) {
								chosenBucket = bin;
								foundNonDelayBucket = true;
							} else {
								if(buckets[curMach][bin] > buckets[curMach][chosenBucket])
									chosenBucket = bin; //bins bucket still smaller the jobTail but bigger then prev chosen
							}
						}
					}

					// all buckets delay curOp? -> get smallest one then
					if( ! foundNonDelayBucket) {
						chosenBucket = 0;
						for(unsigned bin=1; bin<inst.K; bin++) {
							if(buckets[curMach][bin] < buckets[curMach][chosenBucket])
								chosenBucket = bin;
						}
					}

					//filling chosen bucket
					buckets[curMach][chosenBucket] = newMax;

					//setting tails of next in machine order
					bucketHead = buckets[curMach][0];
					for(unsigned bin=1; bin<inst.K; bin++)
						bucketHead = min(bucketHead, buckets[curMach][bin]);
					if(tails[newOp] < bucketHead) {
						tails[newOp] = bucketHead;
					}
				}
			}

			return UINT_MAX;
		}



		//dist fro dource to targOp
		//returns UINT_MAX if has a cycle contaiing targOp
		//expects jobs linearized aready - will ignore precs
		unsigned maxDistToTargParallel(unsigned targOp, vector<unsigned> & heads, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & buckets) const {
			assert(isAlloced());
			assert(heads.size() == inst.O);
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);

			unsigned qInsert = 0;
			unsigned qAccess = 0;

			unsigned curOp;
			unsigned newOp;

			unsigned newMax;

			unsigned chosenBucket;
			unsigned bucketHead;
			unsigned jobHead;
			unsigned curMach;
			bool foundNonDelayBucket;

			fill(heads.begin(), heads.end(), 0);
			fill(indeg.begin(), indeg.end(), 0);
			fill(buckets.data(),buckets.data()+buckets.num_elements(), 0);

			for(unsigned o=1; o<inst.O; o++) {
				//indeg[o] = inst.prev[o].size();
				if(_job[o] != 0)
					indeg[o]++;
				if(_mach[o] != 0)
					indeg[o]++;
				if(indeg[o] == 0)
					Q[qInsert++] = o;
			}

			assert(qInsert>0);

			while(qAccess < qInsert) {
				assert(qAccess<Q.size());
				curOp = Q[qAccess++];
				assert(indeg[curOp] == 0);

				if(curOp == targOp)
					return heads[curOp];

				newMax = heads[curOp] + inst.P[curOp];

				//job order
				newOp = job[curOp];
				if(newOp != 0) {
					assert(indeg[newOp]>0);
					indeg[newOp]--;
					if(indeg[newOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = newOp;
					}
					heads[newOp] = max(heads[newOp], newMax);
				}
				/*
				//job precs
				for(unsigned n : inst.next[curOp]) {
				assert(indeg[n]>0);
				indeg[n]--;
				if(indeg[n] == 0) {
				assert(qInsert<Q.size()-1);
				Q[qInsert++] = n;
				}
				heads[n] = max(heads[n], newMax);
				}
				*/

				//from MACH
				newOp = mach[curOp];
				curMach = inst.operToM[curOp];
				if(newOp != 0) {
					assert(indeg[newOp] > 0);
					indeg[newOp]--;
					if(indeg[newOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = newOp;
					}
					//put it in biggest bin that still generates the lowest machine head:
					//if exists at least one bin that is lower the the jobHead of curOp put it in the biggest bin ammong it
					//if all bins are higher the job head put it in the lowest one
					jobHead = heads[_job[curOp]] + inst.P[_job[curOp]];
					assert(jobHead <= heads[curOp]);

					//getting biggest bucket that doesnt delay curOp
					foundNonDelayBucket = false;
					for(unsigned bin=0; bin<inst.K; bin++) {
						if(buckets[curMach][bin] < jobHead) {
							if( ! foundNonDelayBucket) {
								chosenBucket = bin;
								foundNonDelayBucket = true;
							} else {
								if(buckets[curMach][bin] > buckets[curMach][chosenBucket])
									chosenBucket = bin; //bins bucket still smaller the jobHead but bigger then prev chosen
							}
						}
					}

					// all buckets delay curOp? -> get smallest one then
					if( ! foundNonDelayBucket) {
						chosenBucket = 0;
						for(unsigned bin=1; bin<inst.K; bin++) {
							if(buckets[curMach][bin] < buckets[curMach][chosenBucket])
								chosenBucket = bin;
						}
					}

					//filling chosen bucket
					buckets[curMach][chosenBucket] = newMax;

					//setting heads of next in machine order
					bucketHead = buckets[curMach][0];
					for(unsigned bin=1; bin<inst.K; bin++)
						bucketHead = min(bucketHead, buckets[curMach][bin]);
					if(heads[newOp] < bucketHead) {
						heads[newOp] = bucketHead;
					}
				}
			}
	
			return UINT_MAX;
		}


		static void fillCandAllInCritic(vector<pair<unsigned, unsigned>> & cand, const vector<unsigned> & critic) {
			assert( ! critic.empty());
		
			cand.clear();

			for(unsigned pos=0; pos<critic.size()-1; pos++) {
				if(inst.operToJ[critic[pos]] != inst.operToJ[critic[pos+1]]) {
					assert(inst.operToM[critic[pos]] == inst.operToM[critic[pos+1]]);
					cand.push_back(pair<unsigned, unsigned>(critic[pos] , critic[pos+1]));
				}
			}
		}



		void fillCandAllMachs(vector<pair<unsigned, unsigned>> & cand) const {
			cand.clear();
			unsigned root = 0;
			for(unsigned m=0; m<inst.M; m++) {
				for(unsigned o : inst.machOpers[m]) {
					if(_mach[o] == 0) {
						root = o;
						break;
					}
				}
				assert(inst.operToM[root] == m);
				for(unsigned i=0; i<inst.J-1; i++) {
					assert(root != 0);
					assert(mach[root] != 0);
					cand.push_back(pair<unsigned, unsigned>(root, mach[root]));
					root = mach[root];
				}
				assert(mach[root] ==0);
			}
		}



		void bubbleParallel(double p, unsigned perturbSize, vector<vector<unsigned>> & enforce, vector<unsigned> & enfIndeg, vector<unsigned> & indeg, vector<unsigned> & Q, vector<unsigned> & operPos, vector<bool> & reach) {
#ifndef NDEBUG
			assert(enforce.size() == inst.O);
			//for(const vector<unsigned> & v : enforce) assert(v.capacity() == max(inst.J, inst.M));
			assert(enfIndeg.size() == inst.O);
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(operPos.size() == inst.O);
			assert(reach.size() == inst.O);
			assert(p>0.0 && p<1.0);
			assert(perturbSize <= inst.O);
#endif//NDEBUG

			unsigned index;
			unsigned root;
			unsigned curOp;
			unsigned lastInsertOp;
			unsigned oldPrev;
			unsigned oldNext;

			for(unsigned pert=0; pert<perturbSize; pert++) {
				index = rand() % inst.M;

				//finding order that must be enforced so there is no cycle and precs are obeyed
				mustEnforce(enforce, enfIndeg, MACH, index, indeg, Q, operPos, reach);
				//enforce is poset so no cycles and precs are ok, enfIndeg is indeg of such poset

				//finding root
				for(unsigned o : inst.machOpers[index]) {
					if(_mach[o] == 0) {
						root = o;
						break;
					}
				}

				//bubbling
				lastInsertOp = 0;
				while(root != 0) {
					curOp = root;
					assert(curOp==0   ||   inst.operToM[curOp] == index);
					while(curOp != 0) {
						if(random_number<double>(0.0, 1.0) <= p) { //insert curOp in new order
							//buffers
							oldPrev = _mach[curOp];
							oldNext = mach[curOp];

							//inserting in new order
							if(lastInsertOp != 0) 
								mach[lastInsertOp] = curOp;
							_mach[curOp]=lastInsertOp;

							//removing from old order
							mach[oldPrev] = oldNext;
							_mach[oldNext] = oldPrev;

							//updating root
							if(curOp == root)
								root = mach[curOp];

							lastInsertOp = curOp;
						}
						curOp = mach[curOp];
					}
				}//root != 0
				mach[lastInsertOp] = 0;








			}//pert<perturbSize
		}




		//N5 neighbourhood first improvement - randomizes
		void localSearchParallel(vector<unsigned> & dists, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, vector<pair<unsigned, unsigned>> & cands, vector<unsigned> & critic, multi_array<unsigned, 2> & bucketLimits, multi_array<unsigned, 2> & bucketTops) {
#ifndef NDEBUG
		 
			assert(dists.size() == inst.O);
			assert(prev.size() == inst.O);
			assert(prev.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(cands.capacity() == inst.O);
			assert(critic.capacity() == inst.O);
#endif //NDEBUG
			bool cycle;
			unsigned lastOp;
			unsigned thisMakes;

			setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);

			while(true) {
				thisMakes = makes;

				//cout << "\t\tls iter: " << makes << endl;

				computeCritic(critic, lastOp, prev);
				assert( ! critic.empty());
				assert(critic.size() < inst.O);

				fillCandAllMachs(cands);
				//State::computeCritic(critic, lastOp, prev);
				//State::fillCandAllInCritic(cands, critic);

				random_shuffle(cands.begin(), cands.end());

				for(const pair<unsigned, unsigned> & p : cands) {
					swap(p.first, p.second);

					cycle = setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);
					if(cycle)
						makes = UINT_MAX;

					//cout << "\t\t\ttest: " << makes << endl;
					 
					if(thisMakes > makes) {
						break; //first improvement
					} else {
						swap(p.second, p.first); //undoing swap
					}
				}
				if(thisMakes <= makes)
					break;
			}
		}



#endif //PARALLEL_MACHS














































































#ifdef VAR_PARALLEL_MACHS


	void fillExtendedNeigh(vector<pair<unsigned, unsigned>> & cands, vector<unsigned> & critic, const vector<vector<unsigned>> & allDelay) const {
		unsigned x;
		unsigned y;
		unsigned m;
		set<pair<unsigned, unsigned>> cycleCands;
		vector<vector<unsigned>> cylePreds(inst.O);
		vector<unsigned> Q(inst.O);
		vector<unsigned> operPos(inst.O);
		vector<unsigned> indeg(inst.O);
		
		for(unsigned pos=0; pos<critic.size(); pos++) {
			for(unsigned delayOp : allDelay[critic[pos]]) {
				assert(delayOp>0  && delayOp<inst.O);
				x = delayOp;
				y = critic[pos];
				m = _mach[y];
				//cand.push_back(pair<unsigned, unsigned>(delayOp , critic[pos]));//adding actual neighbour - should be making a cycle
				fillCycleCands(x, y, m, cycleCands, cylePreds, Q, operPos, indeg);
			}
		}
		
		cands.clear();
		for(const pair<unsigned, unsigned> & p : cycleCands)
			cands.push_back(p);
		//cands = vector<pair<unsigned, unsigned>>(cycleCands);
	}


	//expects no cycles in state yet - shiftBack(x,y) and sees cycle members
	//fills (adds) cycle cands with the arcs in the cycle if any - there may be others in it before
	//@return would actually form a cycle?
	bool fillCycleCands(unsigned x, unsigned y, unsigned m, set<pair<unsigned, unsigned>> & cycleCands, vector<vector<unsigned>> & cylePreds, vector<unsigned> & Q, vector<unsigned> & operPos, vector<unsigned> & indeg) const {
#ifndef NDEBUG
			assert(x>0 && x<inst.O);
			assert(y>0 && y<inst.O);
			assert(m>0 && m<inst.O);
			assert(cylePreds.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(operPos.size() == inst.O);
			assert(indeg.size() == inst.O);
			unsigned iters = 0;
#endif
			bool foundCycle;
			queue<unsigned> toCheck;
			unsigned curOper;
		
			foundCycle = fillPredecessorForCycle(x, y, m, cylePreds, Q, operPos, indeg);
			
			//retracing the cycles
			
			//setiing starts
			for(unsigned o : cylePreds[y]){
				assert(foundCycle);
				toCheck.push(o);
			}
			
			while( ! toCheck.empty()) {
#ifndef NDEBUG
				iters++;
				assert(iters<100000);//maybe for bigger instances, would not happen with 10x10 or less
#endif
				assert(foundCycle);//or else how got inside while?
				curOper = toCheck.front();
				toCheck.pop();
				assert(cylePreds[curOper].size() <= 2);//can not have more than two predecessors in a cycle other than job and machine
				for(unsigned aPred : cylePreds[curOper]) {
					assert(aPred > 0);
					assert(aPred < inst.O);
					assert(inst.operToM[aPred]==inst.operToM[curOper]   ||   inst.operToJ[aPred]==inst.operToJ[curOper]);//share job or machine
					toCheck.push(aPred);
					if(inst.operToM[aPred]==inst.operToM[curOper]) {//same machine can swap
						cycleCands.insert(pair<unsigned, unsigned>(aPred, curOper));
					}
				}
			}
		
			return foundCycle;
	}






	//is going to exhcange (x,y) arc with (y,x) - it may form a cycle - cyclePreds holds the predecessesors in the cycles if they are fromed
	//if there is a cycle cyclePreds
	//expects no cycle before
	//@return: will it actually make a cycle?
	//before shift the amchine looks like this: . -> x -> . . . . -> m -> y -> . (m == x if x->y)
	//only arc that can add cycle is y->x . check if x can reach y before shift not using (m->y)
	bool fillPredecessorForCycle(unsigned x, unsigned y, unsigned m, vector<vector<unsigned>> & cylePreds, vector<unsigned> & Q, vector<unsigned> & operPos, vector<unsigned> & indeg) const {
#ifndef NDEBUG
		assert(indeg.size() == inst.O);
		assert(Q.size() == inst.O);
		assert(operPos.size() == inst.O);
		assert(cylePreds.size() == inst.O);
		State bufState = *this;
#endif
	
		unsigned qInsert = 0;
		unsigned qAccess = 0;

		unsigned curOp;
		unsigned newOp;
		
		unsigned startPos;
		unsigned endPos;
		
		for(vector<unsigned> & acp : cylePreds)
			acp.clear();

		//Step 1: prepare a topological order in Q, and its inverted table operPos		
		fill(indeg.begin(), indeg.end(), 0);

		for(unsigned o=1; o<inst.O; o++) {
			if(_job[o] != 0)
				indeg[o]++;
			if(_mach[o]!=0)
				indeg[o]++;
			if(indeg[o] == 0)
				Q[qInsert++] = o;
		}
		
		assert(qInsert>0);
		
		//preparing topological order
		while(qAccess < qInsert) {
			assert(qAccess<Q.size());
			curOp = Q[qAccess];
			assert(indeg[curOp] == 0);

			operPos[curOp] = qAccess;

			qAccess++;

			//from JOB
			newOp = job[curOp];
			if(newOp != 0) {
				assert(indeg[newOp] > 0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
			}

			//from MACH
			newOp = mach[curOp];
			if(newOp != 0) {
				assert(indeg[newOp] > 0);
				indeg[newOp]--;
				if(indeg[newOp] == 0) {
					assert(qInsert<Q.size()-1);
					Q[qInsert++] = newOp;
				}
			}
		}
		assert(qAccess == inst.O-1);
		
		//preparing to propagate from x
		assert(Q[operPos[x]] = x);
		assert(Q[operPos[y]] = y);
		endPos = operPos[y];
		curOp = job[x];
		if(curOp != 0) {
			assert(Q[operPos[curOp]] = curOp);
			startPos = operPos[curOp];
			cylePreds[curOp].push_back(x);
		}
		curOp = mach[x];
		assert((x != m    &&   curOp != y)       ||       (x == m    &&   curOp == y));
		if(curOp != 0   &&   curOp != y) {
			assert(x != m);
			assert(Q[operPos[curOp]] = curOp);
			cylePreds[curOp].push_back(x);
			startPos = min(startPos, operPos[curOp]);
		}
		
		//propagating
		for(unsigned aPos=startPos; aPos<=endPos; aPos++) {
			curOp = Q[aPos];
			if( ! cylePreds[curOp].empty()) {//we can reach curOp from x
				if(job[curOp] != 0) {
					assert(operPos[curOp] < operPos[job[curOp]]);
					assert(cylePreds[job[curOp]].size()<=1);
					cylePreds[job[curOp]].push_back(curOp);
				}

				assert(curOp!=m   ||   mach[curOp]==y);//m is by design the oper right before y in the machine
				if(curOp != m) {//ignoring arc m->y
					if(mach[curOp] != 0) {
						assert(operPos[curOp] < operPos[mach[curOp]]);		
						assert(cylePreds[mach[curOp]].size()<=1);	
						cylePreds[mach[curOp]].push_back(curOp);
					}
				}
			}
		}
		
		assert(bufState == *this);
		
		return ! cylePreds[y].empty();//y cyclePreds not empty? -> reachable from x (other than (x,y)) so it will form a cycle
	}

	

	//expects o1 to be before o2
	void shiftForth(unsigned o1, unsigned o2) {
		assert(inst.operToM[o1] == inst.operToM[o2]);

		unsigned next = mach[o1];
		
		while(next != o2) {
			assert(next != 0);
			swap(o1, next);
			assert(next != mach[o1]);
			next = mach[o1];
		}
		
		swap(o1, o2);
	}




	//expects o1 to be before o2
	void shiftBack(unsigned o1, unsigned o2) {
			assert(o1 != o2);
			assert(o1 != 0);
			assert(o2 != 0);
			assert(inst.operToM[o1] == inst.operToM[o2]);
			unsigned prev = _mach[o2];
		
			while(prev != o1) {
				assert(prev != 0);
				swap(prev, o2);
				assert(prev != _mach[o2]);
				prev = _mach[o2];
			}
		
			swap(o1, o2);
		}






		static string printBucket(const vector<vector<pair<unsigned, unsigned>>> & machGaps, unsigned aK) {
			string str = "";

			for(unsigned k=0; k<aK; k++) {
				str.append("\tk: " + unsigStr(k)+" - ");
				for(unsigned x=0; x<machGaps[k].size(); x++) {
					str.append("["+unsigStr(machGaps[k][x].first)+","+unsigStr(machGaps[k][x].second)+"]    ");
				}
				str.append("\n");
			}

			return str;
		}



		static void fillCandParallelN5(vector<pair<unsigned, unsigned>> & cand, const vector<unsigned> & critic, const vector<vector<unsigned>> & allDelay) {
					unsigned lOp;
		unsigned rOp;

		cand.clear();
			vector<unsigned> jobBb;
			vector<unsigned> machBb;
			blockBegins(jobBb, machBb, critic);
		
			
			//MACH BLOCKS
		assert(machBb.size() > 1);//or else in machine lower bound. why still calling neighbourhhod still?
		assert(machBb[0] == 0);
		assert(machBb[1] > 0);
		assert(machBb[1] < critic.size());
		
		//first mach block
		if(machBb[0]+2 <= machBb[1]) {//has 2 or more opers
			lOp = critic[machBb[1]-2];
			rOp = critic[machBb[1]-1];
			assert(inst.operToM[lOp] == inst.operToM[rOp]);
			//cand.push_back(pair<unsigned, unsigned>(lOp, rOp));
			for(unsigned delayOp : allDelay[rOp])
				cand.push_back(pair<unsigned, unsigned>(delayOp, rOp));
		}
		//middle mach blocks
		for(unsigned bbBaseInd=1;  bbBaseInd<machBb.size()-1; bbBaseInd++) {
			assert(machBb[bbBaseInd] < machBb[bbBaseInd+1]);
			assert(machBb[bbBaseInd+1]<critic.size());
			if(machBb[bbBaseInd]+2 <= machBb[bbBaseInd+1]) { //block has size at least 2
				//left edge
				lOp = critic[machBb[bbBaseInd]];
				rOp = critic[machBb[bbBaseInd]+1];
				//cand.push_back(pair<unsigned, unsigned>(lOp, rOp));
				for(unsigned delayOp : allDelay[lOp])
					cand.push_back(pair<unsigned, unsigned>(delayOp, lOp));

				if(machBb[bbBaseInd]+3 <= machBb[bbBaseInd+1]) {//block has size at least 3
					//right edge
					lOp = critic[machBb[bbBaseInd+1]-2];
					rOp = critic[machBb[bbBaseInd+1]-1];
					//cand.push_back(pair<unsigned, unsigned>(lOp, rOp));
					for(unsigned delayOp : allDelay[rOp])
						cand.push_back(pair<unsigned, unsigned>(delayOp, rOp));
				}
			}
		}
		//last mach block
		if(machBb.back()+2<=critic.size()) {//has 2 or more opers
			lOp = critic[machBb.back()];
			rOp = critic[machBb.back()+1];
			//cand.push_back(pair<unsigned, unsigned>(lOp, rOp));//internal edge
			for(unsigned delayOp : allDelay[lOp])
				cand.push_back(pair<unsigned, unsigned>(delayOp, lOp));
		}


		assert( ! cand.empty());
		assert(cand.size() < inst.O);
		
		}




		//if this works consider giving more weight to - garbage
		//	1 - stuff with more workload
		//	2 - stuff that are in a critical path
		static void fillCandFullDelaySweep(vector<pair<unsigned, unsigned>> & cand, const vector<vector<unsigned>> & allDelay) {
			
			for(unsigned o=1; o<inst.O; o++) {
				for(unsigned delayOp : allDelay[o]) {
					assert(delayOp>0  && delayOp<inst.O);
					cand.push_back(pair<unsigned, unsigned>(delayOp , o));
				}
			}
		}
	
	
	
		static void fillCandCritGapFiller(vector<pair<unsigned, unsigned>> & cand, const vector<unsigned> & critic, const vector<vector<unsigned>> & allDelay) {
		
			assert( ! critic.empty());
		
			//cand.clear();
			
			for(unsigned pos=0; pos<critic.size(); pos++) {
				for(unsigned delayOp : allDelay[critic[pos]]) {;
					assert(delayOp>0  && delayOp<inst.O);
					cand.push_back(pair<unsigned, unsigned>(delayOp , critic[pos]));
				}
			}
		}
		
		
		
		
		
		static void fillCandAllInCritic(vector<pair<unsigned, unsigned>> & cand, const vector<unsigned> & critic) {
			assert( ! critic.empty());
		
			//cand.clear();

			for(unsigned pos=0; pos<critic.size()-1; pos++) {
				if(inst.operToJ[critic[pos]] != inst.operToJ[critic[pos+1]]) {
					assert(inst.operToM[critic[pos]] == inst.operToM[critic[pos+1]]);
					cand.push_back(pair<unsigned, unsigned>(critic[pos] , critic[pos+1]));
				}
			}
		}



		void fillCandAllMachs(vector<pair<unsigned, unsigned>> & cand) const {
			cand.clear();
			unsigned root = 0;
			for(unsigned m=0; m<inst.M; m++) {
				for(unsigned o : inst.machOpers[m]) {
					if(_mach[o] == 0) {
						root = o;
						break;
					}
				}
				assert(inst.operToM[root] == m);
				for(unsigned i=0; i<inst.J-1; i++) {
					assert(root != 0);
					assert(mach[root] != 0);
					cand.push_back(pair<unsigned, unsigned>(root, mach[root]));
					root = mach[root];
				}
				assert(mach[root] ==0);
			}
		}




		//@return: makes - UINT_MAX if cycle
		//@allDelay[x]: All operation in same machine as x that held it back when it was scheduled in the amchine. Empty if non delay bucket was available
		unsigned setMetaParallelAllDelay(vector<unsigned> & dists, unsigned & lastOp, vector<unsigned> & prev, vector<vector<unsigned>> & allDelay, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & bucketLimits, multi_array<unsigned, 2> & bucketTops) const {
			assert(isAlloced());
			assert(dists.size() == inst.O);
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(bucketLimits.shape()[0] == inst.M);
			assert(bucketLimits.shape()[1] == inst.maxK);
			assert(bucketTops.shape()[0] == inst.M);
			assert(bucketTops.shape()[1] == inst.maxK);
			
			//cout << toString() << endl;
			
			
			
			assert(allDelay.size() == inst.O);

			unsigned qInsert = 0;
			unsigned qAccess = 0;

			unsigned curOp;
			unsigned curMach;

			unsigned newMax;

			unsigned chosenBucket;
			unsigned bucketHead;
			unsigned jobHead;
			bool foundNonDelayBucket;

			unsigned foundMakes = 0;

			fill(dists.begin(), dists.end(), 0);
			fill(indeg.begin(), indeg.end(), 0);
			fill(bucketLimits.data(),bucketLimits.data()+bucketLimits.num_elements(), 0);
			fill(bucketTops.data(),bucketTops.data()+bucketTops.num_elements(), 0);

			for(unsigned o=1; o<inst.O; o++) {
				if(_job[o] != 0)
					indeg[o]++;
				if(_mach[o] != 0)
					indeg[o]++;
				if(indeg[o] == 0)
					Q[qInsert++] = o;
			}


			while(qAccess < qInsert) {
				//GET NEW CUR OP
				assert(qAccess<Q.size());
				curOp = Q[qAccess++];
				assert(indeg[curOp] == 0);

				//AUX FOR CUR OP
				curMach = inst.operToM[curOp];
				jobHead = dists[_job[curOp]] + inst.P[_job[curOp]];

				//RELEASE NEXT OPS
				if(job[curOp] != 0) {
					assert(indeg[job[curOp]] > 0);
					indeg[job[curOp]]--;
					if(indeg[job[curOp]] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = job[curOp];
					}
				}
				if(mach[curOp] != 0) {
					assert(indeg[mach[curOp]] > 0);
					indeg[mach[curOp]]--;
					if(indeg[mach[curOp]] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = mach[curOp];
					}
				}

				//FIND PROPPER BUCKET
				//getting biggest bucket that doesnt delay curOp
				foundNonDelayBucket = false;
				for(unsigned bin=0; bin<inst.varK[curMach]; bin++) {
					//if(bucketLimits[curMach][bin] <= jobHead) {
					if(bucketLimits[curMach][bin] < jobHead) {
						if( ! foundNonDelayBucket) {//first non-delay bucket found
							chosenBucket = bin;
							foundNonDelayBucket = true;
						} else {//other non delay bucket previously found
							if(bucketLimits[curMach][bin] > bucketLimits[curMach][chosenBucket])
								chosenBucket = bin; //bin's bucket still smaller the jobHead but bigger then prev chosen
						}
					}
				}
				// all buckets delay curOp? -> get smallest one then
				if( ! foundNonDelayBucket) {
					chosenBucket = 0;
					for(unsigned bin=1; bin<inst.varK[curMach]; bin++) {
						if(bucketLimits[curMach][bin] < bucketLimits[curMach][chosenBucket])
							chosenBucket = bin;
					}
				}
				bucketHead = bucketLimits[curMach][chosenBucket];

				//UPDATE PREV
				if(foundNonDelayBucket) {
					prev[curOp] = _job[curOp];
				} else {//jobHead <  bucketHead
					prev[curOp] = bucketTops[curMach][chosenBucket];
					assert(jobHead==bucketHead   ||   ! foundNonDelayBucket);
				}
				/*
				if(jobHead <= bucketHead) {//machine is holding the oper back(job may be as well)
					assert(! foundNonDelayBucket   ||   jobHead == bucketHead);
					prev[curOp] = bucketTops[curMach][chosenBucket];
				} else {//only job is holding oper back
					prev[curOp] = _job[curOp];
				}
				*/
				//UPDATE ALLDELAY
				allDelay[curOp].clear();
				if( ! foundNonDelayBucket) {
					for(unsigned bin=0; bin<inst.varK[curMach]; bin++) {
						//assert(bucketTops[curMach][bin]>0   &&   bucketTops[curMach][bin]<inst.O);
						assert(bucketTops[curMach][bin]<inst.O);
						if(bucketTops[curMach][bin]>0)
							allDelay[curOp].push_back(bucketTops[curMach][bin]);
					}
				} 

				//UPDATE dists  - foundMakes - lastOps
				dists[curOp] = max(jobHead, bucketHead);
				newMax = dists[curOp] + inst.P[curOp];
				if(foundMakes < newMax) {
					foundMakes = newMax;
					lastOp = curOp;
				}

				//UPDATE BUCKET
				bucketLimits[curMach][chosenBucket] = newMax;
				bucketTops[curMach][chosenBucket] = curOp;
			}//while (qAccess < qInsert)


			if(qAccess<inst.O-1)
				return UINT_MAX;
			else 
				return foundMakes;
		}





		//@return: makes - UINT_MAX if cycle
		//cands: all cand changes in critical path
		unsigned setMetaParallel(vector<unsigned> & dists, unsigned & lastOp, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & bucketLimits, multi_array<unsigned, 2> & bucketTops) const {
			assert(isAlloced());
			assert(dists.size() == inst.O);
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(bucketLimits.shape()[0] == inst.M);
			assert(bucketLimits.shape()[1] == inst.maxK);
			assert(bucketTops.shape()[0] == inst.M);
			assert(bucketTops.shape()[1] == inst.maxK);
			
#ifdef ANALYZE_MODE
			observer.setMetaTimes++;
#endif //ANALYZE_MODE

			unsigned qInsert = 0;
			unsigned qAccess = 0;

			unsigned curOp;
			unsigned curMach;

			unsigned newMax;

			unsigned chosenBucket;
			unsigned bucketHead;
			unsigned jobHead;
			bool foundNonDelayBucket;

			unsigned foundMakes = 0;

			fill(dists.begin(), dists.end(), 0);
			fill(indeg.begin(), indeg.end(), 0);
			fill(bucketLimits.data(),bucketLimits.data()+bucketLimits.num_elements(), 0);
			fill(bucketTops.data(),bucketTops.data()+bucketTops.num_elements(), 0);

			for(unsigned o=1; o<inst.O; o++) {
				if(_job[o] != 0)
					indeg[o]++;
				if(_mach[o] != 0)
					indeg[o]++;
				if(indeg[o] == 0)
					Q[qInsert++] = o;
			}

			while(qAccess < qInsert) {
				//GET NEW CUR OP
				assert(qAccess<Q.size());
				curOp = Q[qAccess++];
				assert(indeg[curOp] == 0);

				//AUX FOR CUR OP
				curMach = inst.operToM[curOp];
				jobHead = dists[_job[curOp]] + inst.P[_job[curOp]];

				//RELEASE NEXT OPS
				if(job[curOp] != 0) {
					assert(indeg[job[curOp]] > 0);
					indeg[job[curOp]]--;
					if(indeg[job[curOp]] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = job[curOp];
					}
				}
				if(mach[curOp] != 0) {
					assert(indeg[mach[curOp]] > 0);
					indeg[mach[curOp]]--;
					if(indeg[mach[curOp]] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = mach[curOp];
					}
				}

				//FIND PROPPER BUCKET
				//getting biggest bucket that doesnt delay curOp
				foundNonDelayBucket = false;
				for(unsigned bin=0; bin<inst.varK[curMach]; bin++) {
					if(bucketLimits[curMach][bin] < jobHead) {
						if( ! foundNonDelayBucket) {
							chosenBucket = bin;
							foundNonDelayBucket = true;
						} else {//other non delay bucket previously found
							if(bucketLimits[curMach][bin] > bucketLimits[curMach][chosenBucket])
								chosenBucket = bin; //bin's bucket still smaller the jobHead but bigger then prev chosen
						}
					}
				}
				// all buckets delay curOp? -> get smallest one then
				if( ! foundNonDelayBucket) {
					chosenBucket = 0;
					for(unsigned bin=1; bin<inst.varK[curMach]; bin++) {
						if(bucketLimits[curMach][bin] < bucketLimits[curMach][chosenBucket])
							chosenBucket = bin;
					}
				}
				bucketHead = bucketLimits[curMach][chosenBucket];

				//UPDATE PREV
				if(jobHead >  bucketHead) {
					prev[curOp] = _job[curOp];
					assert( foundNonDelayBucket);
				} else {//jobHead <=  bucketHead
					prev[curOp] = bucketTops[curMach][chosenBucket];
					assert(jobHead==bucketHead   ||   ! foundNonDelayBucket);
				}

				//UPDATE dists  - foundMakes - lastOps
				dists[curOp] = max(jobHead, bucketHead);
				newMax = dists[curOp] + inst.P[curOp];
				if(foundMakes < newMax) {
					foundMakes = newMax;
					lastOp = curOp;
				}

				//UPDATE BUCKET
				bucketLimits[curMach][chosenBucket] = newMax;
				bucketTops[curMach][chosenBucket] = curOp;
			}//while (qAccess < qInsert)

			if(qAccess<inst.O-1)
				return UINT_MAX;
			else 
				return foundMakes;
		}





		//parallel JSP insa
		void insaParallel(vector<unsigned> & indeg, vector<unsigned> & Q, vector<unsigned> & dists, multi_array<unsigned, 2> & buckets) {
			assert(isAlloced());
			assert(isClear());
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(dists.size() == inst.O);

			unsigned curOp;
			unsigned nextOp;

			unsigned curCost;
			unsigned bigJob = UINT_MAX;
			unsigned bigJobCost = 0;

			vector<unsigned> machRoot(inst.M, 0);

			vector<pair<unsigned, unsigned>> costVoper;
			costVoper.reserve(inst.O-inst.M);

			//linearizing jobs - expects JSP
			assert(inst.roots.size() == inst.J);
			for(unsigned r : inst.roots) {
				curOp = r;
				for(unsigned u=0; u<inst.M-1; u++) {
					assert(inst.next[curOp].size() == 1);
					nextOp = inst.next[curOp][0];

					assert(job[curOp] == 0);
					job[curOp] = nextOp;
					assert(_job[nextOp] == 0);
					_job[nextOp] = curOp;

					curOp = nextOp;
				}
			}

			//insert biggest job
		 
			//curCost = 0;
			for(unsigned j=0; j<inst.J; j++) {
				curCost = 0;
				for(unsigned o : inst.jobOpers[j])
					curCost += inst.P[o];
				if(curCost > bigJobCost) {
					bigJobCost = curCost;
					bigJob = j;
				}
			}
			assert(bigJob < inst.J);
		 

			for(unsigned o : inst.jobOpers[bigJob]) {
				assert(machRoot[inst.operToM[o]] == 0);
				machRoot[inst.operToM[o]] = o;
			}


			//ordering the rest of the opers
			for(unsigned o=1; o<inst.O; o++) {
				if(inst.operToJ[o] != bigJob)
					costVoper.push_back(pair<unsigned, unsigned>(inst.P[o], o));
			}
			sort(costVoper.rbegin(), costVoper.rend());

			//inserting missing ops
			for(unsigned pos=0; pos<costVoper.size(); pos++) {
				curOp = costVoper[pos].second;
				insertInsaParallel(curOp, dists, indeg, Q, buckets, machRoot);
			}
		}



		void insertInsaParallel(unsigned curOp, vector<unsigned> & dists, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & buckets, vector<unsigned> & machRoot) {
			unsigned curCost;

			unsigned bestCost;
			unsigned bestPrev;
			unsigned bestNext;
			const unsigned curMach = inst.operToM[curOp];


			assert(mach[curOp] == 0);
			assert(_mach[curOp] == 0);
			assert(machRoot[curMach] != 0);
			assert(_mach[machRoot[curMach]] == 0);

			mach[curOp] = machRoot[curMach];
			_mach[machRoot[curMach]] = curOp;

			bestCost = insaPosCostParallel(curOp, dists, indeg, Q, buckets);
			bestPrev = 0;
			bestNext = machRoot[curMach];

			while(true) {
				curCost = insaPosCostParallel(curOp, dists, indeg, Q, buckets);

				if(curCost < bestCost) {
					bestCost = curCost;
					bestPrev = _mach[curOp];
					bestNext = mach[curOp];		
				}

				if(mach[curOp] == 0)
					break;

				swap(curOp, mach[curOp]);
			}

			assert(_mach[curOp] != 0);
			assert(mach[_mach[curOp]] == curOp);

			//reseting leaf that sould pojnt to cur op at end of run
			mach[_mach[curOp]] = 0;

			assert(bestPrev==0     ||     mach[bestPrev] == bestNext);
			assert(bestNext==0     ||     _mach[bestNext] == bestPrev);

			_mach[curOp] = bestPrev;
			if(bestPrev != 0) mach[bestPrev] = curOp;
			mach[curOp] = bestNext;
			if(bestNext != 0) _mach[bestNext] = curOp;

			if(bestPrev==0)
				machRoot[curMach] = curOp;
		}



		unsigned insaPosCostParallel(unsigned targOp, vector<unsigned> & dists, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & buckets) const {
			unsigned c1;
			unsigned c2;
			c1 = maxDistFromTargParallel(targOp, dists, indeg, Q, buckets);
			if(c1 == UINT_MAX) return UINT_MAX;
			c2 =  maxDistToTargParallel(targOp, dists, indeg, Q, buckets);
			if(c2 == UINT_MAX) return UINT_MAX;
			return c1 + inst.P[targOp] + c2;
		}


		//dist from sink to targOp
		//returns UINT_MAX if has a cycle contaiing targOp
		//expects jobs linearized aready - will ignore precs
		unsigned maxDistFromTargParallel(unsigned targOp, vector<unsigned> & tails, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & buckets) const {
			assert(isAlloced());
			assert(tails.size() == inst.O);
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);

			unsigned qInsert = 0;
			unsigned qAccess = 0;

			unsigned curOp;
			unsigned newOp;

			unsigned newMax;

			unsigned chosenBucket;
			unsigned bucketHead;
			unsigned jobTail;
			unsigned curMach;
			bool foundNonDelayBucket;

			fill(tails.begin(), tails.end(), 0);
			fill(indeg.begin(), indeg.end(), 0);
			fill(buckets.data(),buckets.data()+buckets.num_elements(), 0);

			for(unsigned o=1; o<inst.O; o++) {
				//indeg[o] = inst.next[o].size();
				if(job[o] != 0)
					indeg[o]++;
				if(mach[o] != 0)
					indeg[o]++;
				if(indeg[o] == 0)
					Q[qInsert++] = o;
			}

			//assert(qInsert>0); //if o will skip while and return UINT_MAX as expected

			while(qAccess < qInsert) {
				assert(qAccess<Q.size());
				curOp = Q[qAccess++];
				assert(indeg[curOp] == 0);

				if(curOp == targOp)
					return tails[curOp];

				newMax = tails[curOp] + inst.P[curOp];

				//job order
				newOp = _job[curOp];
				if(newOp != 0) {
					assert(indeg[newOp]>0);
					indeg[newOp]--;
					if(indeg[newOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = newOp;
					}
					tails[newOp] = max(tails[newOp], newMax);
				}


				//from MACH
				newOp = _mach[curOp];
				curMach = inst.operToM[curOp];
				if(newOp != 0) {
					assert(indeg[newOp] > 0);
					indeg[newOp]--;
					if(indeg[newOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = newOp;
					}
					//put it in biggest bin that still generates the lowest machine head:
					//if exists at least one bin that is lower the the jobTail of curOp put it in the biggest bin ammong it
					//if all bins are higher the job head put it in the lowest one
					jobTail = tails[job[curOp]] + inst.P[job[curOp]];
					assert(jobTail <= tails[curOp]);

					//getting biggest bucket that doesnt delay curOp
					foundNonDelayBucket = false;
					for(unsigned bin=0; bin<inst.varK[curMach]; bin++) {
						if(buckets[curMach][bin] < jobTail) {
							if( ! foundNonDelayBucket) {
								chosenBucket = bin;
								foundNonDelayBucket = true;
							} else {
								if(buckets[curMach][bin] > buckets[curMach][chosenBucket])
									chosenBucket = bin; //bins bucket still smaller the jobTail but bigger then prev chosen
							}
						}
					}

					// all buckets delay curOp? -> get smallest one then
					if( ! foundNonDelayBucket) {
						chosenBucket = 0;
						for(unsigned bin=1; bin<inst.varK[curMach]; bin++) {
							if(buckets[curMach][bin] < buckets[curMach][chosenBucket])
								chosenBucket = bin;
						}
					}

					//filling chosen bucket
					buckets[curMach][chosenBucket] = newMax;

					//setting tails of next in machine order
					bucketHead = buckets[curMach][0];
					for(unsigned bin=1; bin<inst.varK[curMach]; bin++)
						bucketHead = min(bucketHead, buckets[curMach][bin]);
					if(tails[newOp] < bucketHead) {
						tails[newOp] = bucketHead;
					}
				}
			}

			return UINT_MAX;
		}



		//dist fro dource to targOp
		//returns UINT_MAX if has a cycle contaiing targOp
		//expects jobs linearized aready - will ignore precs
		unsigned maxDistToTargParallel(unsigned targOp, vector<unsigned> & heads, vector<unsigned> & indeg, vector<unsigned> & Q, multi_array<unsigned, 2> & buckets) const {
			assert(isAlloced());
			assert(heads.size() == inst.O);
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);

			unsigned qInsert = 0;
			unsigned qAccess = 0;

			unsigned curOp;
			unsigned newOp;

			unsigned newMax;

			unsigned chosenBucket;
			unsigned bucketHead;
			unsigned jobHead;
			unsigned curMach;
			bool foundNonDelayBucket;

			fill(heads.begin(), heads.end(), 0);
			fill(indeg.begin(), indeg.end(), 0);
			fill(buckets.data(),buckets.data()+buckets.num_elements(), 0);

			for(unsigned o=1; o<inst.O; o++) {
				//indeg[o] = inst.prev[o].size();
				if(_job[o] != 0)
					indeg[o]++;
				if(_mach[o] != 0)
					indeg[o]++;
				if(indeg[o] == 0)
					Q[qInsert++] = o;
			}

			assert(qInsert>0);

			while(qAccess < qInsert) {
				assert(qAccess<Q.size());
				curOp = Q[qAccess++];
				assert(indeg[curOp] == 0);

				if(curOp == targOp)
					return heads[curOp];

				newMax = heads[curOp] + inst.P[curOp];

				//job order
				newOp = job[curOp];
				if(newOp != 0) {
					assert(indeg[newOp]>0);
					indeg[newOp]--;
					if(indeg[newOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = newOp;
					}
					heads[newOp] = max(heads[newOp], newMax);
				}


				//from MACH
				newOp = mach[curOp];
				curMach = inst.operToM[curOp];
				if(newOp != 0) {
					assert(indeg[newOp] > 0);
					indeg[newOp]--;
					if(indeg[newOp] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = newOp;
					}
					//put it in biggest bin that still generates the lowest machine head:
					//if exists at least one bin that is lower the the jobHead of curOp put it in the biggest bin ammong it
					//if all bins are higher the job head put it in the lowest one
					jobHead = heads[_job[curOp]] + inst.P[_job[curOp]];
					assert(jobHead <= heads[curOp]);

					//getting biggest bucket that doesnt delay curOp
					foundNonDelayBucket = false;
					for(unsigned bin=0; bin<inst.varK[curMach]; bin++) {
						if(buckets[curMach][bin] < jobHead) {
							if( ! foundNonDelayBucket) {
								chosenBucket = bin;
								foundNonDelayBucket = true;
							} else {
								if(buckets[curMach][bin] > buckets[curMach][chosenBucket])
									chosenBucket = bin; //bins bucket still smaller the jobHead but bigger then prev chosen
							}
						}
					}

					// all buckets delay curOp? -> get smallest one then
					if( ! foundNonDelayBucket) {
						chosenBucket = 0;
						for(unsigned bin=1; bin<inst.varK[curMach]; bin++) {
							if(buckets[curMach][bin] < buckets[curMach][chosenBucket])
								chosenBucket = bin;
						}
					}

					//filling chosen bucket
					buckets[curMach][chosenBucket] = newMax;

					//setting heads of next in machine order
					bucketHead = buckets[curMach][0];
					for(unsigned bin=1; bin<inst.varK[curMach]; bin++)
						bucketHead = min(bucketHead, buckets[curMach][bin]);
					if(heads[newOp] < bucketHead) {
						heads[newOp] = bucketHead;
					}
				}
			}
	
			return UINT_MAX;
		}





		void bubbleParallel(double p, unsigned perturbSize, vector<vector<unsigned>> & enforce, vector<unsigned> & enfIndeg, vector<unsigned> & indeg, vector<unsigned> & Q, vector<unsigned> & operPos, vector<bool> & reach) {
#ifndef NDEBUG
			assert(enforce.size() == inst.O);
			//for(const vector<unsigned> & v : enforce) assert(v.capacity() == max(inst.J, inst.M));
			assert(enfIndeg.size() == inst.O);
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(operPos.size() == inst.O);
			assert(reach.size() == inst.O);
			assert(p>0.0 && p<1.0);
			assert(perturbSize <= inst.O);
#endif//NDEBUG

			unsigned index;
			unsigned root;
			unsigned curOp;
			unsigned lastInsertOp;
			unsigned oldPrev;
			unsigned oldNext;

			for(unsigned pert=0; pert<perturbSize; pert++) {
				index = rand() % inst.M;

				//finding order that must be enforced so there is no cycle and precs are obeyed
				mustEnforce(enforce, enfIndeg, MACH, index, indeg, Q, operPos, reach);
				//enforce is poset so no cycles and precs are ok, enfIndeg is indeg of such poset

				//finding root
				for(unsigned o : inst.machOpers[index]) {
					if(_mach[o] == 0) {
						root = o;
						break;
					}
				}


				//bubbling
				lastInsertOp = 0;
				while(root != 0) {
					curOp = root;
					assert(curOp==0   ||   inst.operToM[curOp] == index);
					while(curOp != 0) {
						if(random_number<double>(0.0, 1.0) <= p) { //insert curOp in new order
							//buffers
							oldPrev = _mach[curOp];
							oldNext = mach[curOp];

							//inserting in new order
							if(lastInsertOp != 0) 
								mach[lastInsertOp] = curOp;
							_mach[curOp]=lastInsertOp;

							//removing from old order
							mach[oldPrev] = oldNext;
							_mach[oldNext] = oldPrev;

							//updating root
							if(curOp == root)
								root = mach[curOp];

							lastInsertOp = curOp;
						}
						curOp = mach[curOp];
					}
				}//root != 0
				mach[lastInsertOp] = 0;



			}//pert<perturbSize
		}




		//N5 neighbourhood first improvement - randomizes
		void localSearchParallel(vector<unsigned> & dists, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q, vector<pair<unsigned, unsigned>> & cands, vector<unsigned> & critic, multi_array<unsigned, 2> & bucketLimits, multi_array<unsigned, 2> & bucketTops) {
#ifndef NDEBUG
		 
			assert(dists.size() == inst.O);
			assert(prev.size() == inst.O);
			assert(prev.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(cands.capacity() == inst.O);
			assert(critic.capacity() == inst.O);
#endif //NDEBUG
			bool cycle;
			unsigned lastOp;
			unsigned thisMakes;

			setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);

			while(true) {
				thisMakes = makes;

				//cout << "\t\tls iter: " << makes << endl;

				computeCritic(critic, lastOp, prev);
				assert( ! critic.empty());
				assert(critic.size() < inst.O);

				fillCandAllMachs(cands);
				//State::computeCritic(critic, lastOp, prev);
				//State::fillCandAllInCritic(cands, critic);

				random_shuffle(cands.begin(), cands.end());

				for(const pair<unsigned, unsigned> & p : cands) {
					swap(p.first, p.second);

					cycle = setMetaParallel(dists, lastOp, prev, indeg, Q, bucketLimits, bucketTops);
					if(cycle)
						makes = UINT_MAX;

					//cout << "\t\t\ttest: " << makes << endl;
					 
					if(thisMakes > makes) {
						break; //first improvement
					} else {
						swap(p.second, p.first); //undoing swap
					}
				}
				if(thisMakes <= makes)
					break;
			}
		}
		
		
		
		//Aux for dispatchInOrder - 
		//EXPECTS A FULL ORDER
		void fixJobOrder() {
			unsigned cur;
			
			for(unsigned aRoot : inst.roots) {
				cur = aRoot;
				_job[cur] = 0;
				
				for(unsigned u=0; u<inst.M-1; u++) {
					assert(cur>0 && cur<inst.O);
					assert(inst.next[cur].size() <= 1);
					assert(inst.next[cur].empty()   ||   inst.next[cur][0]>0);
					assert(inst.next[cur].empty()   ||   inst.next[cur][0]<inst.O);
					job[cur] = inst.next[cur][0];
					_job[inst.next[cur][0]] = cur;
					cur = inst.next[cur][0];
				}
				job[cur] = 0;
			}
		
		}
		
		
		void dispatchInOrder(const vector<unsigned> & order) {
#ifndef NDEBUG
			assert(isAlloced());
			assert(order.size() == inst.O-1);
			for(unsigned o : order) {
				assert(o > 0);
				assert(o < inst.O);
			}
#endif //NDEBUG		

			fixJobOrder();


			cout << toString() << endl;
		}
		
		
		
		
		void fillPrioParallel(vector<unsigned> & prio, vector<unsigned> & indeg, vector<unsigned> & Q) const {
			assert(isAlloced());
			assert(indeg.size() == inst.O);
			assert(Q.size() == inst.O);
			assert(prio.size() == inst.O);

			unsigned qInsert = 0;
			unsigned qAccess = 0;

			unsigned curOp;

			fill(prio.begin(), prio.end(), 0);
			fill(indeg.begin(), indeg.end(), 0);

			for(unsigned o=1; o<inst.O; o++) {
				if(_job[o] != 0)
					indeg[o]++;
				if(_mach[o] != 0)
					indeg[o]++;
				if(indeg[o] == 0)
					Q[qInsert++] = o;
			}

			while(qAccess < qInsert) {
				//GET NEW CUR OP
				assert(qAccess<Q.size());
				curOp = Q[qAccess++];
				assert(indeg[curOp] == 0);

				//RELEASE NEXT OPS
				if(job[curOp] != 0) {
					assert(indeg[job[curOp]] > 0);
					indeg[job[curOp]]--;
					if(indeg[job[curOp]] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = job[curOp];
					}
				}
				if(mach[curOp] != 0) {
					assert(indeg[mach[curOp]] > 0);
					indeg[mach[curOp]]--;
					if(indeg[mach[curOp]] == 0) {
						assert(qInsert<Q.size()-1);
						Q[qInsert++] = mach[curOp];
					}
				}


				//UPDATE dists  - foundMakes - lastOps
				prio[curOp] = max(prio[_job[curOp]], prio[_mach[curOp]]) + 1;

			}//while (qAccess < qInsert)

			assert(qAccess==inst.O-1);
		}
		
		
		Schedule reoptFillGaps() const {
			vector<unsigned> prio(inst.O);
			vector<unsigned> indeg(inst.O);
			vector<unsigned> Q(inst.O);	
			
			fillPrioParallel(prio, indeg, Q);
			Schedule s;
			s.dipatchByPiority(prio, SMALLER_GAP);
			
			return s;
		}
		
		
		//receives the order of the machines and sets the state accordingly
		//designed to be used in conjunction with Schedule::getMachOrders()
		//does not adjust the makespan
		void importFromMachOrder(const vector<vector<unsigned>> & allMachOrder) {
			assert(isAlloced());
			assert(allMachOrder.size() == inst.M);
			
			unsigned curOp;
			
			//set jobs
			fixJobOrder();
			
			//set machines
			for(unsigned m=0; m<inst.M; m++) {
				assert(allMachOrder[m].size() == inst.J);
				
				for(unsigned index=0; index<inst.J; index++) {
					curOp = allMachOrder[m][index];
					assert(curOp>0    &&    curOp<inst.O);
					if(index != 0)
						_mach[curOp] = allMachOrder[m][index-1];
					else
						_mach[curOp] = 0;
						
					if(index != inst.J-1)
						mach[curOp] = allMachOrder[m][index+1];
					else
						mach[curOp] = 0;
				}
			}
		}



#endif //VAR_PARALLEL_MACHS




		vector<unsigned> job;
		vector<unsigned> _job;
		vector<unsigned> mach;
		vector<unsigned> _mach;
		unsigned makes;
		unsigned millisecsFound;
		unsigned lPenalty;
		unsigned ePenalty;
		unsigned penalties;
	};
