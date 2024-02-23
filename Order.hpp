/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#pragma once

#include "Settings.hpp"
#include "Utilities.hpp"
#include "Instance.hpp"

class OrderCritic {
public:
	OrderCritic() {  }
	OrderCritic(unsigned _makes, unsigned _entry, unsigned _exit) : cycle(true), makes(_makes), entry(_entry), exit(_exit) {  }
	~OrderCritic() {  }
	
	string toString() const {
		if(cycle)
			return "cycle";
		
		return "makes: " + unsigStr(makes) + " ; entry: " + unsigStr(entry) + " ; exit: " + unsigStr(exit);
	}
	
	bool cycle;
	unsigned makes;
	unsigned entry;
	unsigned exit;
};



class MapValue {
public:
	MapValue() {  }
	MapValue(unsigned _index, unsigned _value) : index(_index), value(_value){  }
	~MapValue() {  }
	
	bool operator < (const MapValue & mv) const { return value < mv.value; }
	bool operator > (const MapValue & mv) const { return value > mv.value; }
	
	bool operator == (const MapValue & mv) const { return mv.index == index; }
	
	string toString() const {
		string str = "";
		
		str.append("<i:"+unsigStr(index));
		str.append("|v:"+unsigStr(value)+">");
		
		return str;
	}
	
	unsigned index;
	unsigned value;
};



class BnbNode {
public:
	BnbNode() {  }
	BnbNode(const vector<unsigned> & _order, const vector<unsigned> & _heads, const vector<unsigned> & _tails, unsigned _lowerBound, unsigned _makes) : order(_order), heads(_heads), tails(_tails), lowerBound(_lowerBound), makes(_makes){  }
	~BnbNode() {  }
	
	bool canImprove() const {  
		assert(lowerBound <= makes);
		assert(cPos <= pPos);
		return lowerBound < makes   &&   cPos < pPos;
	}


	string toString() const {
	
		string str = "lowerBound: "+unsigStr(lowerBound)+" ; makes: "+unsigStr(makes);
		
		str.append(" | order:");
		for(unsigned u : order)
			str.append(","+unsigStr(u));
		
		str.append("  ;  heads:");
		for(unsigned u : heads)
			str.append(","+unsigStr(u));
			
		str.append("  ;  tails:");
		for(unsigned u : tails)
			str.append(","+unsigStr(u));
		
		return str;
	}
	
	bool operator< (const BnbNode & bn) const {
		if(this->lowerBound > bn.lowerBound) {
			return true;
		} else if(this->lowerBound == bn.lowerBound  &&  this->makes > bn.makes) {
			return true;
		} else {
			return false;
		}
	}
	
	vector<unsigned> order;
	vector<unsigned> heads;
	vector<unsigned> tails;
	unsigned lowerBound;
	unsigned makes;
	unsigned cPos;
	unsigned pPos;
};




class BnInfo {
public:
	BnInfo(unsigned _makes, const vector<unsigned> & _order) : defined(true), makes(_makes), order(_order) {  }
	BnInfo() : defined(false) {  }
	~BnInfo() {  }
	
	string toString() const {
		if( ! defined) {
			return "DUMMY";
		}
			
		string str = "";
		
		str.append("  ;  makes "+unsigStr(makes));
		str.append("  ;  order");
		for(unsigned o : order) {
			str.append(" "+unsigStr(o));
		}
		
		return str;
	}
	
	bool operator< (const BnInfo & bn) const { return this->makes > bn.makes; }
	
	bool defined;
	unsigned makes;
	vector<unsigned> order;
};



namespace Carlier {
	//Aux for schrage. Fills queue with all opers that will be ready at current time T
	//@return: new readyLimit of orderOfHeads (delimits the ready opers), readys are up to and including orderOfHeads[readyLimit]
	//Updates: toVisitByTail inserting newly ready elements from old readyLimit until the heads surpasses T
	unsigned fillQueue(const vector<unsigned> & tails, const vector<MapValue> & orderOfHeads, unsigned readyLimit, unsigned T, priority_queue<MapValue> & toVisitByTail) {
		assert(orderOfHeads[readyLimit].value <= T);
		
		unsigned newReadyLimit = readyLimit;
		
		while(newReadyLimit<orderOfHeads.size()-1   &&   orderOfHeads[newReadyLimit+1].value <= T) {
			newReadyLimit++;//before push because readyLimit is up to AND INCLUDING
			toVisitByTail.push(MapValue(orderOfHeads[newReadyLimit].index, tails[orderOfHeads[newReadyLimit].index]));
		}
		
		return newReadyLimit;
	}



	//@return: lower bound for onemach problem
	unsigned preemptSchrage(const vector<unsigned> & heads, const vector<unsigned> & tails, const vector<unsigned> & procTs) {
		const unsigned size = heads.size();
		assert(tails.size() == size);
		assert(procTs.size() == size);
		
		unsigned makes = 0;
		unsigned T;
		unsigned nextRel;
		
		unsigned currOp;
		priority_queue<MapValue> toVisitByTail;
		
		vector<MapValue> orderOfHeads(heads.size());//value is heads, index is index -- all indexes order by head once at the
		unsigned readyLimit;//up to (and including) this index in orderOfHeads all opers are ready and should be already inserted in toVisistByTail
		
		vector<unsigned> timeStillNeeded(size, 0);
		for(unsigned u=0; u<size; u++)
			timeStillNeeded[u] = procTs[u];
		
		//ordering heads from smallest to biggest
		for(unsigned u=0; u<size; u++)
			orderOfHeads[u] = MapValue(u, heads[u]);
		sort(orderOfHeads.begin(), orderOfHeads.end());  //O(n*log(n))
		
		//initializing initial releases
		T = orderOfHeads[0].value; //T is lowest head
		readyLimit=0; //at the very least first element is ready
		toVisitByTail.push(MapValue(orderOfHeads[readyLimit].index, tails[orderOfHeads[readyLimit].index]));
		
		readyLimit = fillQueue(tails, orderOfHeads, readyLimit, T, toVisitByTail);
		
		while( ! toVisitByTail.empty()) {
			currOp = toVisitByTail.top().index;

			assert(timeStillNeeded[currOp] <= procTs[currOp]);
			
			nextRel = (readyLimit+1 < orderOfHeads.size()   ?   orderOfHeads[readyLimit+1].value   :   UINT_MAX);

			if(nextRel < T + timeStillNeeded[currOp]) {
				timeStillNeeded[currOp] -= nextRel-T;
				assert(timeStillNeeded[currOp] > 0);
				T = nextRel;
			} else if(nextRel >= T + timeStillNeeded[currOp]) {
				T += timeStillNeeded[currOp];
				makes = max(makes, T+tails[currOp]);
				timeStillNeeded[currOp] = 0;
				toVisitByTail.pop();
			}
			
			readyLimit = fillQueue(tails, orderOfHeads, readyLimit, T, toVisitByTail);
			
			if(toVisitByTail.empty()   &&   readyLimit != orderOfHeads.size()-1) {//Will have idle time
				T = nextRel;
				readyLimit = fillQueue(tails, orderOfHeads, readyLimit, T, toVisitByTail);
			}
		}
		
		return makes;
	}


		//@return: makes
	unsigned schrage(const vector<unsigned> & heads, const vector<unsigned> & tails, vector<unsigned> & schrageOrder, const vector<unsigned> & procTs) {

		const unsigned size = heads.size();

		assert(size == tails.size());
		assert(size == schrageOrder.size());
		
		unsigned T;
		priority_queue<MapValue> toVisitByTail;//value is tail, index is index -- will go from highest tail to lowest, but cointains indexes ready at current time T
		vector<MapValue> orderOfHeads(heads.size());//value is heads, index is index -- all indexes order by head once at the
		unsigned readyLimit;//up to (and including) this index in orderOfHeads all opers are ready and should be already inserted in toVisistByTail
		
		unsigned currOp;
		
		unsigned makes = 0;
		
		//ordering heads from smallest to biggest
		for(unsigned u=0; u<heads.size(); u++)
			orderOfHeads[u] = MapValue(u, heads[u]);
		sort(orderOfHeads.begin(), orderOfHeads.end());  //O(n*log(n))
		
		T = orderOfHeads[0].value; //T is lowest head
		readyLimit=0; //at the very least first element is ready
		toVisitByTail.push(MapValue(orderOfHeads[readyLimit].index, tails[orderOfHeads[readyLimit].index]));

		readyLimit = fillQueue(tails, orderOfHeads, readyLimit, T, toVisitByTail);
		
		for(unsigned nSched=0; nSched<size; nSched++) {//must schedule exacly size opers, will schedule one each loop
			assert( ! toVisitByTail.empty());

			currOp = toVisitByTail.top().index;
			toVisitByTail.pop();

			schrageOrder[nSched] = currOp;

			T += procTs[currOp];
			makes = max(makes, T+tails[currOp]);
			readyLimit = fillQueue(tails, orderOfHeads, readyLimit, T, toVisitByTail);
			if(toVisitByTail.empty()   &&   nSched<size-1) {//time is not enough to have an unscheduled ready jobs
				assert(readyLimit+1 < orderOfHeads.size());
				T = orderOfHeads[readyLimit+1].value;
				readyLimit = fillQueue(tails, orderOfHeads, readyLimit, T, toVisitByTail);
				assert( ! toVisitByTail.empty());
			}
		}
		assert(toVisitByTail.empty());
		
		return makes;
	}


	//all indexes (params and return) are related to order, not overall op index
	//@return:oper C pos in order, if == critPathEnd it is optimal
	unsigned jSet(unsigned critPathBeg, unsigned critPathEnd, const vector<unsigned> & tails, const vector<unsigned> & order) {
		assert(critPathBeg <= critPathEnd);
		assert(critPathEnd < order.size());
		assert(tails.size() == order.size());
		
		unsigned cPos = critPathEnd; 
		
		for(unsigned pos=critPathBeg; pos<critPathEnd; pos++) {
			if(tails[order[pos]] < tails[order[critPathEnd]])
				cPos = pos;
		}
		
		return cPos;
	}



	/*
	bool cycle;
	unsigned makes;
	unsigned entry;
	unsigned exit;
	*/
	//@return 0:makespan, 1:index in order of entry op, 2:index in order of exit op
	//tuple<unsigned, unsigned, unsigned> critic(const vector<unsigned> heads, const vector<unsigned> tails, const vector<unsigned> & order) const {
	OrderCritic getOrderCritic(const vector<unsigned> & heads, const vector<unsigned> & tails, const vector<unsigned> & order, const vector<unsigned> & procTs)  {

		const unsigned size = heads.size();	

		assert(size == tails.size());
		assert(size == order.size());
		assert(size == procTs.size());
		
		unsigned makes = 0;
		unsigned entry = 0;
		unsigned exit = 0;
		
		unsigned tripaCost = 0;
		unsigned tripaEntry = 0;
		
		for(unsigned u=0; u<size; u++) {
			if(tripaCost   <   heads[order[u]]) {
				tripaEntry = u;
				tripaCost = heads[order[u]];
			}
			tripaCost += procTs[order[u]];

			if(makes < tripaCost+tails[order[u]]     ||     (makes == tripaCost+tails[order[u]]  &&  exit-entry < u-tripaEntry) ) {
				makes = tripaCost+tails[order[u]];
				entry = tripaEntry;
				exit = u;
			}
		}
		
		assert(entry <= exit);
		
		return OrderCritic(makes, entry, exit);
		//return tuple<unsigned, unsigned, unsigned>(makes, entry, exit);
	}



	unsigned lowerBound(const vector<unsigned> & heads, const vector<unsigned> & tails, const vector<unsigned> & order, unsigned cPos, unsigned pPos, const vector<unsigned> & procTs) {

		const unsigned size = procTs.size();
		assert(cPos <= pPos);
		assert(pPos < size);
		assert(heads.size() == size);
		assert(heads.size() == size);
		assert(order.size() == size);
		
		unsigned bestLower;
		
		unsigned minHead = UINT_MAX;
		unsigned minTail = UINT_MAX;
		unsigned cumProcT = 0;
		
		bestLower = preemptSchrage(heads, tails, procTs);
		
		for(unsigned u=cPos+1; u<=pPos; u++) {
			minHead = min(minHead, heads[order[u]]);
			minTail = min(minTail, tails[order[u]]);
			cumProcT += procTs[order[u]];
		}
		
		if(minHead != UINT_MAX   &&   minTail != UINT_MAX) {//cPos==pPos => jSet is empty => solution is optimal
			bestLower = max(bestLower, minHead+cumProcT+minTail); //h(J)
			
			minHead = min(minHead, heads[order[cPos]]);
			minTail = min(minTail, tails[order[cPos]]);
			cumProcT += procTs[order[cPos]];
			
			bestLower = max(bestLower, minHead+cumProcT+minTail); //h(J+{c})
		}
		
		
		return bestLower;
	}


	//will shcrage to find order, and set all fields of bn
	void setBnbNode(const vector<unsigned> & heads, const vector<unsigned> & tails, BnbNode & bn, const vector<unsigned> & procTs) {
		const unsigned size = heads.size();

		assert(bn.order.size() == size);
		assert(tails.size() == size);
		
		//tuple<unsigned, unsigned, unsigned> critInfo;//0:makespan, 1:index in order of entry op, 2:index in order of exit op
		OrderCritic critInfo;
		unsigned cPos;
		
		bn.heads = heads;
		bn.tails = tails;
		
		bn.makes = schrage(heads, tails, bn.order, procTs);
		critInfo = getOrderCritic(heads, tails, bn.order, procTs);
		//assert(get<0>(critInfo) == bn.makes);
		assert(critInfo.makes == bn.makes);
		//cPos = jSet(get<1>(critInfo), get<2>(critInfo), tails, bn.order);
		cPos = jSet(critInfo.entry, critInfo.exit, tails, bn.order);
		//bn.lowerBound = lowerBound(heads, tails, bn.order, cPos, get<2>(critInfo));
		bn.lowerBound = lowerBound(heads, tails, bn.order, cPos, critInfo.exit, procTs);
		bn.cPos = cPos;
		//bn.pPos = get<2>(critInfo);
		bn.pPos = critInfo.exit;
	}



	//will adjust newN heads and tails accorfing to test 1 and 2 from (Carlier 1986)
	//can be done MUCH better  --- is recomuting loads of stuff
	void carlierTrimm(vector<unsigned> & heads, vector<unsigned> & tails, const vector<unsigned> & order, unsigned bestMakes, const vector<unsigned> & procTs) {
		
		const unsigned size = procTs.size();

		assert(size > 0);
		assert(tails.size() == size);
		assert(heads.size() == size);
		assert(order.size() == size);

		OrderCritic critInfo = getOrderCritic(heads, tails, order, procTs);
		unsigned cPos;
		unsigned pPos;

		vector<unsigned> kSetPosV;

		unsigned cumTjSet=0;
		unsigned minHead = UINT_MAX;
		unsigned minTail = UINT_MAX;
		unsigned h_J;
		unsigned kOper;

		cPos = jSet(critInfo.entry, critInfo.exit, tails, order);
		pPos = critInfo.exit;

		assert(cPos <= pPos); 
		if(cPos == pPos)
			return;


		for(unsigned aPos=cPos+1; aPos<=pPos; aPos++) {
			cumTjSet += procTs[order[aPos]];
			minHead = min(minHead, heads[order[aPos]]);
			minTail = min(minTail, tails[order[aPos]]);
		}

		h_J = minHead + cumTjSet + minTail;

		assert(minHead < UINT_MAX);
		assert(minTail < UINT_MAX);
		
		for(unsigned aPos=0; aPos<cPos; aPos++) {
			//order[aPos] is candidate to kSet as it is not in jSet
		    if(procTs[order[aPos]] > bestMakes - h_J)
				kSetPosV.push_back(aPos);
		}

		for(unsigned aPos=pPos+1; aPos<order.size(); aPos++) {
			//order[aPos] is candidate to kSet as it is not in jSet
		    if(procTs[order[aPos]] > bestMakes - h_J)
				kSetPosV.push_back(aPos);
		}

		for(unsigned kPos : kSetPosV) {
			kOper = order[kPos];
			assert(kOper < order.size());
			//test 1
			if( (heads[kOper] + procTs[kOper] + cumTjSet + tails[order[pPos]])    >    bestMakes) {
				heads[kOper] = max(heads[kOper], minHead + cumTjSet);
			}

			//test 2
			if( (minHead + procTs[kOper] + cumTjSet + tails[kOper])  >  bestMakes ) {
				tails[kOper] = max(tails[kOper] , minTail + cumTjSet);
			}
		}
	}



		//will not enforce order o elements with 0 proc time. Not a problem in schedule but can be when fixing cycle
	void forcePrec(const Poset & poset, vector<unsigned> & heads, vector<unsigned> & tails, const vector<unsigned> & procTs) {
		assert(poset.cg.n == heads.size());
		assert(poset.nEl == heads.size());

		queue<unsigned> toVisit;
		vector<unsigned> indeg(poset.nEl);
		unsigned pivot;

		
		//Forward -> update pos
		for(unsigned u=0; u<poset.nEl; u++) {
			indeg[u] = poset.parents[u].size();
		}
		for(unsigned aRoot : poset.roots) {
			assert(poset.parents[aRoot].empty());
			toVisit.push(aRoot);
		}
		assert( ! toVisit.empty());

		while( ! toVisit.empty() ) {
			pivot = toVisit.front();
			toVisit.pop();

			for(unsigned pos : poset.children[pivot]) {
				//pos is ready to schedule only when preds(all pivots leading to it) are. 
				heads[pos] = max(heads[pos], heads[pivot]+procTs[pivot]);
				assert(indeg[pos] > 0);
				indeg[pos]--;
				if(indeg[pos] == 0)
					toVisit.push(pos);
			}
		}


		//Backward -> update pre
		for(unsigned u=0; u<poset.nEl; u++) {
			indeg[u] = poset.children[u].size();
		}
		for(unsigned aLeaf : poset.leafs) {
			assert(poset.children[aLeaf].empty());
			toVisit.push(aLeaf);
		}
		assert( ! toVisit.empty());

		while( ! toVisit.empty() ) {
			pivot = toVisit.front();
			toVisit.pop();

			for(unsigned pre : poset.parents[pivot]) {
			    //if both are ready pre has longer tail so it gets priority
				tails[pre] = max(tails[pre], tails[pivot]+procTs[pivot]);
				assert(indeg[pre] > 0);
				indeg[pre]--;
				if(indeg[pre] == 0)
					toVisit.push(pre);
			}
		}
	}





	bool isPrecInOrder(const Poset & poset, const vector<unsigned> & order, const vector<unsigned> & procTs) {
		assert(poset.cg.n == order.size());
		for(unsigned bailiff=0; bailiff<order.size()-1; bailiff++) {
			for(unsigned provost=bailiff+1; provost<order.size(); provost++) {
				if(poset.cg.path(order[provost], order[bailiff])) {
					if(procTs[order[provost]] != 0   &&   procTs[order[bailiff]] != 0) {
						return false;//the order goes against at least one prec
					}
				}
			}
		}
		return true;
	}




	BnInfo carlierBnb(const Poset & poset, const vector<unsigned> & heads, const vector<unsigned> & tails, bool capIters, const vector<unsigned> & procTs) {

		const unsigned size = procTs.size();

		assert(size == tails.size());
		assert(heads.size() == size);

		priority_queue<pair<BnbNode, unsigned>> toVisit;
	
		BnbNode currN;
		BnbNode bestN;
		BnbNode newN;
		
		//bool precIsOk;
		unsigned nLoops = 0;
		
		unsigned cumT;
		unsigned minHead;
		
		bestN.order.resize(size);
		bestN.heads = heads;
		bestN.tails = tails;
		
		newN.order.resize(size);

		forcePrec(poset, bestN.heads, bestN.tails, procTs); //will adjust heads and tails to force prec

		setBnbNode(bestN.heads, bestN.tails, bestN, procTs); //will shcrage to find order, and set all fields of bn
		
		if(bestN.canImprove()) { //else it is already optimal
			toVisit.push(pair<BnbNode, unsigned>(bestN,rand()));
		}
		
		while( ! toVisit.empty()) {
			nLoops++;
			currN = toVisit.top().first;
			toVisit.pop();

			assert(isPrecInOrder(poset, currN.order, procTs));
			assert(currN.canImprove());//or else shouldn't be in queue
			
			if(currN.lowerBound < bestN.makes) {
				newN.heads = currN.heads;
				newN.tails = currN.tails;
				cumT = 0;
				for(unsigned pos=currN.cPos+1; pos<=currN.pPos; pos++) {
					cumT += procTs[currN.order[pos]];
				}
				//tail(c) = MAX( tail(c) , SUM(procT(i))+tail(p) )
				newN.tails[currN.order[currN.cPos]] = max(newN.tails[currN.order[currN.cPos]]   ,   cumT + newN.tails[currN.order[currN.pPos]]);

				//will adjust newN heads and tails accorfing to test 1 and 2 from (Carlier 1986)
				carlierTrimm(newN.heads, newN.tails, currN.order, bestN.makes, procTs);
					
				setBnbNode(newN.heads, newN.tails, newN, procTs); //will shcrage to find order and set all fields of bn

				if(isPrecInOrder(poset, newN.order, procTs)) {
					if(bestN.makes > newN.makes)
						bestN = newN;
					
					if(newN.canImprove()  &&  newN.lowerBound<bestN.makes) {
						toVisit.push(pair<BnbNode, unsigned>(newN,rand()));
					}
				}
					

				newN.heads = currN.heads;
				newN.tails = currN.tails;
				cumT=0;
				minHead = UINT_MAX;
				for(unsigned pos=currN.cPos+1; pos<=currN.pPos; pos++) {
					cumT += procTs[currN.order[pos]];
					minHead = min(minHead, currN.heads[currN.order[pos]]);
				}
				assert(minHead < UINT_MAX);//J set empty? if so shouldn't be in queue
				//head(c) = MAX( head(c) , MIN(head(i))+SUM(procT(i)) )
				newN.heads[currN.order[currN.cPos]] = max(newN.heads[currN.order[currN.cPos]]   ,   minHead+cumT);

				carlierTrimm(newN.heads, newN.tails, currN.order, bestN.makes, procTs);

				setBnbNode(newN.heads, newN.tails, newN, procTs); //will shcrage to find order and set all fields of bn


				if(isPrecInOrder(poset, newN.order, procTs)) {
					if(bestN.makes > newN.makes)
						bestN = newN;
					
					if(newN.canImprove()  &&  newN.lowerBound<bestN.makes) {
						toVisit.push(pair<BnbNode, unsigned>(newN,rand()));
					}
				}
			}	
		}

		return BnInfo(bestN.makes, bestN.order);
	}

}
