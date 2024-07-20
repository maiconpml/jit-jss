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

	void testSwap(unsigned o1, unsigned o2) {

		assert(o1 != 0);
		assert(o2 != 0);
		assert(o1 != o2);

		unsigned prevO1 = _mach[o1];
		unsigned postO1 = mach[o1];
		unsigned prevO2 = _mach[o2];
		unsigned postO2 = mach[o2];

		mach[o1] = postO2;
		_mach[o2] = prevO1;

		if (prevO1 != 0) 
			mach[prevO1] = o2;

		if (postO2 != 0) 
			_mach[postO2] = o1;

		if (postO1 == o2) {
			_mach[o1] = o2;

			mach[o2] = o1;

			return;
		}

		_mach[o1] = prevO2;
		mach[o2] = postO1;

		if (postO1 != 0) 
			_mach[postO1] = o2;
		
		if (prevO2 != 0) 
			mach[prevO2] = o1;
	}

	
	void swap(unsigned o1, unsigned o2) {
		assert(o1 != 0);
		assert(o2 != 0);
		assert(o1 != o2);

		unsigned prev;
		unsigned post;

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
			assert(Q[operPos[fromOp]] == fromOp);			
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
		//assert(prevJobOp==0   ||   job[prevJobOp]==postJobOp);
		//assert(postJobOp==0   ||   _job[postJobOp]==prevJobOp);

		job[newOp] = postJobOp;
		if(postJobOp!=0) _job[postJobOp] = newOp;
		_job[newOp] = prevJobOp;
		if(prevJobOp!=0) job[prevJobOp] = newOp;
	}
	void insertMachOper(unsigned newOp, unsigned prevMachOp, unsigned postMachOp) {
		assert(newOp != 0);
		//assert(prevMachOp==0   ||   mach[prevMachOp]==postMachOp);
		//assert(postMachOp==0   ||   _mach[postMachOp]==prevMachOp);

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
		unsigned mach = 0;
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

	static void findCriticOper(vector<vector<unsigned>>& critic, vector<unsigned>& starts, vector<unsigned> _job, vector<unsigned> _mach) {

		for (unsigned i = 0; i < inst.O; ++i) {
			critic[i].clear();
		}

		for (unsigned op = 1; op < inst.O; ++op) {

			if (starts[op] + inst.P[op] > inst.deadlines[op]) {
	
				unsigned auxOp = op;
				vector<unsigned> opCritic;
				
				opCritic.push_back(op);
				
				while (_mach[auxOp] && starts[_mach[auxOp]] + inst.P[_mach[auxOp]] > starts[_job[auxOp]] + inst.P[_job[auxOp]]) {

					opCritic.push_back(_mach[auxOp]);

					auxOp = _mach[auxOp];
				}

				if (opCritic.size() > 1) {
					critic[op] = opCritic;
				}
			}
		}
	}

	static void fillCandidatesTest1(vector<pair<unsigned, unsigned>>& cands, vector<unsigned>& mach) {

		assert(mach.size() == inst.O);
		assert(cands.capacity() == inst.O);

		for (unsigned currentOp = 1; currentOp < inst.O; ++currentOp) {

			if (inst.operToJ[currentOp] != inst.operToJ[mach[currentOp]] && mach[currentOp]) {
				cands.push_back(pair<unsigned, unsigned>(currentOp, mach[currentOp]));
			}
		}
	}

	static void fillCandidatesCritic(vector<pair<unsigned, unsigned>>& cands, vector<vector<unsigned>>& criticOper) {

		assert(cands.capacity() == inst.O);

		for (unsigned currentOp = 1; currentOp < inst.O; ++currentOp) {

			if (criticOper[currentOp].size() > 0) {
				cands.push_back(pair<unsigned, unsigned>(criticOper[currentOp][0], criticOper[currentOp][1]));
				if (criticOper[currentOp].size() > 2) cands.push_back(pair<unsigned, unsigned>(criticOper[currentOp][criticOper[currentOp].size() - 1], criticOper[currentOp][criticOper[currentOp].size() - 2]));
			}
		}
	}

	static void fillCandidatesCritic2(vector<pair<unsigned, unsigned>>& cands, vector<vector<unsigned>>& criticOper) {

		set<pair<unsigned, unsigned>> candsAux;
		// 
		for (unsigned currentOp = 1; currentOp < inst.O; ++currentOp) {
			// verify if there is a critic way
			if (criticOper[currentOp].size() > 0) {

				candsAux.insert(pair<unsigned, unsigned>(criticOper[currentOp][1], criticOper[currentOp][0]));

				if (criticOper[currentOp].size() > 2)candsAux.insert(pair<unsigned, unsigned>(criticOper[currentOp][criticOper[currentOp].size() - 1], criticOper[currentOp][criticOper[currentOp].size() - 2]));
			}
		}

		for (set<pair<unsigned, unsigned>>::iterator it = candsAux.begin(); it != candsAux.end(); ++it) {
			cands.push_back(*it);
		}
	}

	static void fillCandidatesTestA(vector<pair<unsigned, unsigned>>& cands, vector<unsigned>& mach, vector<unsigned>& startTime, vector<unsigned>& operPenalties) {

		assert(mach.size() == inst.O);
		assert(cands.capacity() == inst.O);

		unsigned sumPenaltiesBefore;
		unsigned sumPenaltiesAfter;

		for (unsigned currentOp = 1; currentOp < inst.O; ++currentOp) {

			sumPenaltiesBefore = operPenalties[currentOp] + operPenalties[mach[currentOp]];
			sumPenaltiesAfter = (inst.deadlines[mach[currentOp]] > startTime[currentOp] + inst.P[mach[currentOp]] ? inst.deadlines[mach[currentOp]] - startTime[currentOp] + inst.P[mach[currentOp]] : startTime[currentOp] + inst.P[mach[currentOp]] - inst.deadlines[mach[currentOp]]) + (inst.deadlines[currentOp] > startTime[mach[currentOp]] + inst.P[mach[currentOp]] ? inst.deadlines[currentOp] - startTime[mach[currentOp]] + inst.P[mach[currentOp]] : startTime[mach[currentOp]] + inst.P[mach[currentOp]] - inst.deadlines[currentOp]);

			if (!mach[currentOp]) {
				continue;
			}
			else if ((startTime[currentOp] + inst.P[currentOp]) < inst.deadlines[currentOp] && (startTime[mach[currentOp]] + inst.P[mach[currentOp]]) > inst.deadlines[mach[currentOp]]) {
				cands.push_back(pair<unsigned, unsigned>(currentOp, mach[currentOp]));
			}
			else if ((startTime[currentOp] + inst.P[currentOp]) > inst.deadlines[currentOp] && (startTime[mach[currentOp]] + inst.P[mach[currentOp]]) > inst.deadlines[mach[currentOp]] && sumPenaltiesAfter < sumPenaltiesBefore) {
				cands.push_back(pair<unsigned, unsigned>(currentOp, mach[currentOp]));
			}
			else if ((startTime[currentOp] + inst.P[currentOp]) < inst.deadlines[currentOp] && (startTime[mach[currentOp]] + inst.P[mach[currentOp]]) < inst.deadlines[mach[currentOp]] && sumPenaltiesAfter < sumPenaltiesBefore) {
				cands.push_back(pair<unsigned, unsigned>(currentOp, mach[currentOp]));
			}
		}
	}

	static void fillCandidatesTest2(vector<pair<unsigned, unsigned>> & cands, vector<unsigned> & mach, vector<unsigned> starts) {

		assert(mach.size() == inst.O);
		assert(cands.capacity() == inst.O);

		for (unsigned currentOp = 1; currentOp < inst.O; ++currentOp) {
			assert(inst.operToJ[currentOp] != inst.operToJ[mach[currentOp]]);
			if (!mach[currentOp]) {
				continue;
			}
			else if ((starts[currentOp] + inst.P[currentOp]) < inst.deadlines[currentOp] && (starts[mach[currentOp]] + inst.P[mach[currentOp]]) > inst.deadlines[mach[currentOp]]) {
				cands.push_back(pair<unsigned, unsigned>(currentOp, mach[currentOp]));
			}
		}
	}

	static void fillCandidatesTest3(vector<pair<unsigned, unsigned>>& cands, vector<unsigned>& mach, vector<unsigned> starts) {

		for (unsigned curOp = 1; curOp < inst.O; ++curOp) {
			if (starts[curOp] + inst.P[curOp] >= inst.deadlines[curOp]) continue;

			unsigned nextCurOp = mach[curOp];

			while (nextCurOp != 0) {
				assert(inst.operToJ[curOp] != inst.operToJ[nextCurOp]);

				if (starts[nextCurOp] + inst.P[nextCurOp] > inst.deadlines[nextCurOp]) {
					cands.push_back(pair<unsigned, unsigned>(curOp, nextCurOp));
				}

				nextCurOp = mach[nextCurOp];
			}
		}
	}

	static void fillCandidatesTest4(vector<pair<unsigned, unsigned>>& cands, vector<unsigned>& mach) {

		for (unsigned curOp = 1; curOp < inst.O; ++curOp) {

			unsigned nextCurOp = mach[curOp];

			while (nextCurOp != 0) {
				assert(inst.operToJ[curOp] != inst.operToJ[nextCurOp]);

				cands.push_back(pair<unsigned, unsigned>(curOp, nextCurOp));

				nextCurOp = mach[nextCurOp];
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


	//Luv(0, *) = max {Luv(0,w) + Luv(w, *)  :  w in Q} - Q is {u,...., V}
	//
/*	unsigned guidepost(const vector<unsigned> & heads, const vector<unsigned> & tails) const {


	}
	*/
	
	

	//Balas Vazacopoulos 1998
	//cand: type - direction - bailiff - provost
	

	//Zhang 2007
	//cand: direction - bailiff - provost
	///XXX CHECK WHEN CANDS REPEAT -
	// -- instance block ize two all the same
	// -- other situations???
	//TODO SPECIAL CASES WHERE SIZE IS 2 or 3 WLL REPEAT CANDIDATES (FORTH AND BACK IN ADJ ARE SAME)


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
	static bool detectRepeat(vector<int> & oldValues, unsigned & posCurMakes, unsigned & cycleLastPos, unsigned & cycleL, double newMakes, bool isNewBest, unsigned maxD, unsigned maxLen) {
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

	//shift early operations making its completion time closer of its due date
	void shiftOperations(vector<unsigned>& starts, vector<unsigned>& prev, vector<unsigned>& indeg, vector<unsigned>& Q) {

		vector<unsigned> newStarts(inst.O, 0);
		vector<unsigned> startJobSuccessor(inst.J, UINT_MAX);
		vector<unsigned> startMachSuccessor(inst.M, UINT_MAX);

		unsigned qInsert = Q.size();
		unsigned curOp;
		unsigned newMakes = 0;
		
		while (qInsert > 0) {

			curOp = Q[--qInsert];

			if (!curOp) continue;

			unsigned sMS = startMachSuccessor[inst.operToM[curOp]];
			unsigned sJS = startJobSuccessor[inst.operToJ[curOp]];

			unsigned shiftLimit = min(sJS, sMS);

			unsigned shiftTarget = min(inst.deadlines[curOp], shiftLimit);
			shiftTarget = shiftTarget - inst.P[curOp];

			newStarts[curOp] = max(starts[curOp], shiftTarget);

			newMakes = max(newStarts[curOp]+inst.P[curOp], newMakes);

			startJobSuccessor[inst.operToJ[curOp]] = newStarts[curOp];
			startMachSuccessor[inst.operToM[curOp]] = newStarts[curOp];
 		}

		assert(newMakes >= makes);
		makes = newMakes;

		starts = newStarts;
	}

	//set makes with shift to minimize earliness penalties
	bool setMeta(vector<unsigned>& dists, unsigned& lastOp, vector<unsigned>& prev, vector<unsigned>& indeg, vector<unsigned>& Q) {
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

		for (unsigned o = 1; o < inst.O; o++) {
			if (_job[o] != 0)
				indeg[o]++;
			if (_mach[o] != 0)
				indeg[o]++;
			if (indeg[o] == 0) {
				prev[o] = 0;
				Q[qInsert++] = o;
			}
		}

		//assert(qInsert>0);

		while (qAccess < qInsert) {
			assert(qAccess < Q.size());
			curOp = Q[qAccess++];

			assert(indeg[curOp] == 0);

			newMax = dists[curOp] + inst.P[curOp];
			if (makes < newMax) {
				makes = newMax;
				lastOp = curOp;
			}

			//from JOB
			newOp = job[curOp];
			if (newOp != 0) {
				assert(indeg[newOp] > 0);
				indeg[newOp]--;
				if (indeg[newOp] == 0) {
					assert(qInsert < Q.size() - 1);
					Q[qInsert++] = newOp;
				}
				if (dists[newOp] < newMax) {
					dists[newOp] = newMax;
					prev[newOp] = curOp;
				}
			}

			//from MACH
			newOp = mach[curOp];
			if (newOp != 0) {
				assert(indeg[newOp] > 0);
				indeg[newOp]--;
				if (indeg[newOp] == 0) {
					assert(qInsert < Q.size() - 1);
					Q[qInsert++] = newOp;
				}
				if (dists[newOp] < newMax) {
					dists[newOp] = newMax;
					prev[newOp] = curOp;
				}
			}
		}

#ifdef SHIFT_OPERS
#ifndef NDEBUG
		inst.calcPenalties(dists, ePenalty, lPenalty, operPenalties);
		
		unsigned testPenalties = ePenalty + lPenalty;
#endif //NDEBUG

		shiftOperations(dists, prev, indeg, Q);
#endif //SHIFT_OPERS

		inst.calcPenalties(dists, ePenalty, lPenalty, operPenalties);
		penalties = ePenalty + lPenalty;
		startTime = dists;

#ifdef SHIFT_OPERS
		assert(penalties <= testPenalties);
#endif //SHIFT_OPERS

		return qAccess < inst.O - 1;
	}
	

//	//does not set critic but do set makes
//	//@return: cycle?
//	bool setMeta(vector<unsigned> & dists, unsigned & lastOp, vector<unsigned> & prev, vector<unsigned> & indeg, vector<unsigned> & Q)  {
//
//#ifdef SHIFT_OPERS
//		return setMetaWithShift(dists, lastOp, prev, indeg, Q);
//#endif //SHIFT_OPERS
//
//		assert(isAlloced());
//		assert(dists.size() == inst.O);
//		assert(prev.size() == inst.O);
//		assert(indeg.size() == inst.O);
//		assert(Q.size() == inst.O);
//
//		unsigned qInsert = 0;
//		unsigned qAccess = 0;
//
//		unsigned curOp;
//		unsigned newOp;
//
//		unsigned newMax;
//
//		makes = 0;
//
//		fill(dists.begin(), dists.end(), 0);
//		fill(indeg.begin(), indeg.end(), 0);
//
//		for(unsigned o=1; o<inst.O; o++) {
//			if(_job[o] != 0)
//				indeg[o]++;
//			if(_mach[o] != 0)
//				indeg[o]++;
//			if(indeg[o] == 0) {
//				prev[o] = 0;
//				Q[qInsert++] = o;
//			}
//		}
//
//		//assert(qInsert>0);
//
//		while(qAccess < qInsert) {
//			assert(qAccess<Q.size());
//			curOp = Q[qAccess++];
//
//			assert(indeg[curOp] == 0);
//
//			newMax = dists[curOp] + inst.P[curOp];
//			if(makes < newMax) {
//				makes = newMax;
//				lastOp = curOp;
//			}
//
//			//from JOB
//			newOp = job[curOp];
//			if(newOp != 0) {
//				assert(indeg[newOp] > 0);
//				indeg[newOp]--;
//				if(indeg[newOp] == 0) {
//					assert(qInsert<Q.size()-1);
//					Q[qInsert++] = newOp;
//				}
//				if(dists[newOp] < newMax) {
//					dists[newOp] = newMax;	
//					prev[newOp] = curOp;
//				}			
//			}
//
//			//from MACH
//			newOp = mach[curOp];
//			if(newOp != 0) {
//				assert(indeg[newOp] > 0);
//				indeg[newOp]--;
//				if(indeg[newOp] == 0) {
//					assert(qInsert<Q.size()-1);
//					Q[qInsert++] = newOp;
//				}
//				if(dists[newOp] < newMax) {
//					dists[newOp] = newMax; 
//					prev[newOp] = curOp;
//				}			
//			}
//		}
//		inst.calcPenalties(dists, ePenalty,lPenalty, operPenalties);
//		penalties = ePenalty + lPenalty;
//		startTime = dists;
//		
//		return qAccess<inst.O-1;
//	}

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

	void printPenaltys(){

		vector<unsigned> schedule = genSchedule();
#ifndef PRINT_ONLY_RESULT
		cout << makes << " ";
#endif //PRINT_ONLY_RESULT

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


		vector<unsigned> job;
		vector<unsigned> _job;
		vector<unsigned> mach;
		vector<unsigned> _mach;
		vector<unsigned> startTime;
		vector<double> operPenalties;
		vector<double> tardPenalties;
		unsigned makes;
		unsigned millisecsFound;
		double lPenalty;
		double ePenalty;
		double penalties;
	};
