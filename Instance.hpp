/*
  @author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/
#pragma once

#include <vector>
#include <climits>
#include <queue>
#include <algorithm>

#include <boost/multi_array.hpp>

#include "Settings.hpp"
#include "italiano.hpp"
#include "Utilities.hpp"
#include "RandomNumberGenerator.hpp"

using namespace boost;

class Poset {
public:
	Poset(const vector<pair<unsigned, unsigned>> & precs, unsigned _nEl) : nEl(_nEl), children(_nEl), parents(_nEl), cg(_nEl) {
		for(unsigned u=0; u<nEl; u++) {
			leafs.push_back(u);
			roots.push_back(u);
		}
		for(const pair<unsigned, unsigned> & p : precs) {
			addPrec(p.first, p.second);
		}
	}

	Poset(unsigned _nEl) : nEl(_nEl) , children(_nEl), parents(_nEl), cg(_nEl) {
		for(unsigned u=0; u<nEl; u++) {
			leafs.push_back(u);
			roots.push_back(u);
		}
	}
	
	~Poset() {  }

	bool operator!=(const Poset & p) { return !((*this)==p); }
	bool operator==(const Poset & p) {

		if(nEl != p.nEl)
			return false;
		
		assert(children.size() == nEl);
		assert(children.size() == p.children.size());
		for(unsigned a=0; a<nEl; a++) {
			if(children[a].size() != p.children[a].size())
				return false;
			for(unsigned b=0; b<children[a].size(); b++) {
				if(children[a][b] != p.children[a][b])
					return false;
			}
		}

		return true;
	}
	
	//@return: created cycle doing so?
	bool addAll(const vector<pair<unsigned, unsigned>> & precs) {
		bool anyCycle = false;
		for(const pair<unsigned, unsigned> & p : precs){
			anyCycle |= addPrec(p.first, p.second);
		}
		return anyCycle;
	}

	//@return: created cycle doing so?
	bool addPrec(unsigned from, unsigned to) {
		assert( ! cg.hasCycle());

		unsigned u=0;
		if( ! cg.path(from, to)) {

			cg.add(from, to);
		
			children[from].push_back(to);
			parents[to].push_back(from);

			if(children[from].size() == 1) {
				for(u=0; u<leafs.size(); u++) {
					if(leafs[u] == from)
						break;
				}
				assert(leafs[u] == from);//or else u==leafs.size() because from is not in leafs
				leafs.erase(leafs.begin()+u);
				assert( ! leafs.empty());
			}
			if(parents[to].size() == 1) {
				for(u=0; u<roots.size(); u++) {
					if(roots[u] == to)
						break;
				}
				assert(roots[u] == to);//or else u==leafs.size() because from is not in leafs
				roots.erase(roots.begin()+u);
				assert( ! roots.empty());
			}
		}
		return cg.hasCycle();
	}

	string toString() const {		
		string str = (cg.hasCycle() ? "\t\tCycle" : "\t\tNo Cycle");
		str.append("   ;   roots{ ");
		for(unsigned r : roots)
			str.append(unsigStr(r)+" ");
		str.append("}");
			
		str.append(" leafs{ ");
		for(unsigned l : leafs)
			str.append(unsigStr(l)+" ");
		str.append("}");
		
		str.append("\n\tForward: ");
		for(unsigned o=0; o<children.size(); o++) {
			str.append(unsigStr(o)+"->{ ");
			for(unsigned c : children[o]) {
				str.append(unsigStr(c)+" ");
			}
			str.append("} ");
		}
		
		str.append(" \n\tbackward: ");
		for(unsigned o=0; o<parents.size(); o++) {
			str.append("{ ");
			for(unsigned p : parents[o]) {
				str.append(unsigStr(p)+" ");
			}
			str.append("}->"+unsigStr(o)+" ");
		}
		
		return str;
	}
	
	double orderStr() const {
		double os = 0.0;

		for(unsigned from=0; from<nEl; from++) {
			for(unsigned to=0; to<nEl; to++) {
				if(from != to) {
					if(cg.path(from, to)) {
						os++;
					}
				}
			}
		}
		
		unsigned comb = (nEl*(nEl-1))/2;
		
		os /= (double)comb;
		
		return os;
	}
	

	//index here are not global - they are the machine index
	unsigned nEl;
	vector<unsigned> roots;
	vector<unsigned> leafs;
	vector<vector<unsigned>> children; //children[x] == vector of the machine index of the DIRECT children of x
	vector<vector<unsigned>> parents;
	ClosedGraph<MAX_MACHS> cg; //takes into account the transitive closure
};

class Inst {
public:
	Inst() { }
	Inst(const string & filePath) { parse(filePath);}
	~Inst() { }

	string toString() const {
		string str = "";

		str.append("J: " + unsigStr(J));
		str.append(" ; M: " + unsigStr(M));
		str.append(" ; O: " + unsigStr(O));

		str.append("\nP:");	
		for(unsigned u=1 ; u<O; u++)
			str.append(unsigStr(u)+"(" + unsigStr(P[u])+")  ");

		str.append("\njSize:");
		for(unsigned n : jSize)
			str.append(" " + unsigStr(n));
		str.append("\nmSize:");
		for(unsigned n : mSize)
			str.append(" " + unsigStr(n));

		str.append("\njobOpers: ");
		for(const vector<unsigned> & aJob : jobOpers) {
			str.append("\n\t");
			for(unsigned o : aJob)
				str.append(" "+unsigStr(o));
		}
		str.append("\nmachOpers: ");
		for(const vector<unsigned> & aMach : machOpers) {
			str.append("\n\t");
			for(unsigned o : aMach)
				str.append(" "+unsigStr(o));
		}

		str.append("\noperToJ: ");
		for(unsigned pos : operToJ)
			str.append(" " + unsigStr(pos));
		str.append("\noperToM: ");
		for(unsigned pos : operToM)
			str.append(" " + unsigStr(pos));

		str.append("\njmToIndex: ");
		assert(jmToIndex.size() == J);
		for(unsigned j=0; j<J; j++) {
			assert(jmToIndex[j].size() == M);
			str.append("\n\t");
			for(unsigned m=0; m<M; m++) {
				str.append(" "+unsigStr(jmToIndex[j][m]));
			}
		}

		str.append("\nroots: ");
		for(unsigned o : roots)
			str.append(" " + unsigStr(o));
		str.append("\nleafs: ");
		for(unsigned o : leafs)
			str.append(" " + unsigStr(o));

		str.append("\nnext: ");
		assert(next.size() == O);
		for(unsigned o=0; o<O; o++) {
			str.append("    "+unsigStr(o)+":");
			for(unsigned nextO : next[o])
				str.append(","+unsigStr(nextO));
		}

		str.append("\nprev: ");
		assert(prev.size() == O);
		for(unsigned o=0; o<O; o++) {
			str.append("    "+unsigStr(o)+":");
			for(unsigned prevO : prev[o])
				str.append(","+unsigStr(prevO));
		}

		assert(posets.size() == J);
		for(const Poset & aPoset : posets)
			str.append("\n"+aPoset.toString());

		return str;
	}


	void setHeads(vector<unsigned> & heads, vector<unsigned> & indeg, vector<unsigned> & Q) const {
#ifndef NDEBUG
		assert(heads.size() == O);
		assert(indeg.size() == O);
		assert(Q.size() == O);
		bool found;
		fill(Q.begin(), Q.end(), 0);
#endif
		fill(heads.begin(), heads.end(), 0);

		unsigned insertQ = 0;
		unsigned removeQ = 0;
		unsigned curOp;

		for(unsigned o=1; o<O; o++)
			indeg[o] = prev[o].size();

		for(unsigned r : roots) {
			if(indeg[r] == 0) {
				Q[insertQ++] = r;
			}
		}

#ifndef NDEBUG
		for(unsigned o=1; o<O; o++) {
			if(indeg[o] == 0) {
				found = false;
			    for(unsigned r : roots) {
					if(o == r) {
						found = true;
						break;
					}
				}
				assert(found);
			}
		}
#endif

		while(removeQ < insertQ) {
			curOp = Q[removeQ++];

			assert(curOp != 0);
			assert(indeg[curOp] == 0);

			for(unsigned p : next[curOp]) {
				assert(indeg[p] > 0);
				indeg[p]--;
				if(indeg[p] == 0)
					Q[insertQ++] = p;

				heads[p] = max(heads[p], heads[curOp] + P[curOp]);
			}
		}
#ifndef NDEBUG
		assert(insertQ == O-1);//no cycle
		for(unsigned o=1; o<O; o++) {
			if(prev[o].size() == 0)
				assert(heads[o] == 0);//roots
			else
				assert(heads[o] != 0);//not roots
		} 
#endif
	}


	void setTails(vector<unsigned> & tails, vector<unsigned> & indeg, vector<unsigned> & Q) const {
#ifndef NDEBUG
		assert(tails.size() == O);
		assert(indeg.size() == O);
		assert(Q.size() == O);
		bool found;
		fill(Q.begin(), Q.end(), 0);
#endif
		fill(tails.begin(), tails.end(), 0);

		unsigned insertQ = 0;
		unsigned removeQ = 0;
		unsigned curOp;

		for(unsigned o=1; o<O; o++)
			indeg[o] = next[o].size();

		for(unsigned l : leafs) {
			if(indeg[l] == 0) {
				Q[insertQ++] = l;
			}
		}

#ifndef NDEBUG
		for(unsigned o=1; o<O; o++) {
			if(indeg[o] == 0) {
				found = false;
			    for(unsigned l : leafs) {
					if(o == l) {
						found = true;
						break;
					}
				}
				assert(found);
			}
		}
#endif

		while(removeQ < insertQ) {
			curOp = Q[removeQ++];

			assert(curOp != 0);
			assert(indeg[curOp] == 0);

			for(unsigned p : prev[curOp]) {
				assert(indeg[p] > 0);
				indeg[p]--;
				if(indeg[p] == 0)
					Q[insertQ++] = p;

				tails[p] = max(tails[p], tails[curOp] + P[curOp]);
			}
		}

#ifndef NDEBUG
		assert(insertQ == O-1);//no cycle
		for(unsigned o=1; o<O; o++) {
			if(next[o].size() == 0)
				assert(tails[o] == 0);//leafs
			else
				assert(tails[o] != 0);//not leafs
		} 
#endif
	}

	//max for each machine: of min head + min tail + procTs of that machine
	unsigned lowerBoundNasiri(vector<unsigned> & indeg, vector<unsigned> & Q) const {
		assert(O > 0);
		assert(indeg.size() == O);
		assert(Q.size() == O);

	    vector<unsigned> heads(O);
		vector<unsigned> tails(O);
		//vector<unsigned> indeg(O);
		//vector<unsigned> Q(O);

		unsigned minHead;
		unsigned minTail;
		unsigned machT;
		unsigned lb = 0;

		setHeads(heads, indeg, Q);
		setTails(tails,indeg,Q);

		for(const vector<unsigned> & aMach : machOpers) {
			minHead = UINT_MAX;
			minTail = UINT_MAX;
			machT = 0;

			for(unsigned o : aMach) {
				assert(o!=0);
				assert(o<O);
				minHead = min(minHead, heads[o]);
				minTail = min(minTail, tails[o]);
				machT += P[o];
			}
			lb = max(lb, minHead+machT+minTail);
		}
		return lb;
	}

	//max of lowerBoundnasiri and job lower bound
	unsigned lowerBoundTkz(vector<unsigned> & indeg, vector<unsigned> & Q) const {
		unsigned lb = lowerBoundNasiri(indeg, Q);
		unsigned jobSum;

		for(unsigned j=0; j<J; j++) {
			jobSum = 0;
			for(unsigned o : jobOpers[j])
				jobSum += P[o];
			lb = max(lb, jobSum);
		}

		return lb;
	}

	//max of lowerBoundnasiri and job lower bound
	unsigned lowerBoundTkz() const {
		vector<unsigned> indeg(O);
		vector<unsigned> Q(O);
		return lowerBoundTkz(indeg, Q);
	}

	unsigned lowerBoundMss() const {
		unsigned sum;
		unsigned lb = 0;

		for(unsigned j=0; j<J; j++) {
			sum = 0;
			for(unsigned o : jobOpers[j])
				sum += P[o];
			lb = max(lb, sum);
		}
		for(unsigned m=0; m<M; m++) {
			sum = 0;
			for(unsigned o : machOpers[m])
				sum += P[o];
			lb = max(lb, sum);
		}
		
		return lb;
	}

	bool verifySchedule(const vector<unsigned> & starts, unsigned expecMakes) const {
		unsigned foundMakes = 0;

		for(unsigned o1=1; o1<O; o1++) {
			for(unsigned o2=1; o2<O; o2++) {
				if(o1 != o2) {
					//previous were ok?
					if(hasPrec(o1, o2)) {
						if(starts[o1] + P[o1] > starts[o2]) {
							cout << "Prev not OK: " << o1 << "->" << o2 << endl;
							return false;
						}
					}

					//Job lock
					if(operToJ[o1] == operToJ[o2]) {
						if(starts[o1] < starts[o2]) {
							if(starts[o1] + P[o1] > starts[o2]) {
								cout << "job lock not OK: " << o1 << "-" << o2 << endl;
								return false;
							}
						}	else {
							if(starts[o2] + P[o2] > starts[o1]) {
								cout << "job lock not OK: " << o1 << "-" << o2 << endl;
								return false;
							}
						}
					}

					//mach lock
					if(operToM[o1] == operToM[o2]) {
						if(starts[o1] < starts[o2]) {
							if(starts[o1] + P[o1] > starts[o2]) {
								cout << "mach lock not OK: " << o1 << "-" << o2 << endl;
								return false;
							}
						}	else {
							if(starts[o2] + P[o2] > starts[o1]) {
								cout << "mach lock not OK: " << o1 << "-" << o2 << endl;
								return false;
							}
						}
					}
				}
			}
		}

		for(unsigned o=1; o<O; o++) 
			foundMakes = max(foundMakes, starts[o]+P[o]);
		
		if(foundMakes != expecMakes) {
			cout << "makes not OK: " << endl;
			return false;
		}

		return true;
	}

	void iToR2(const vector<unsigned>& starts) {

		cout << "machine task timeStart timeEnd";

		for (unsigned i = 1; i < O; ++i) {
			cout << operToM[i] << " " << operToJ[i] << " " << starts[i] << " " << starts[i] + P[i] << endl;
		}
	}

	void iToR(const vector<unsigned> & starts) const{
    cout << "library(candela)" << endl;
    cout << "data <- list(" << endl;
    for(unsigned m = 0; m < M; m++){
        
			for(unsigned j = 0; j < J;j++){
				int i  = jmToIndex[j][m];
    	  cout<<"  list(label = 'task "<< j <<"',name='machine " << m <<"', level=" <<m<<", start="<< starts[i]<<", end="<< starts[i] +  P[i]<<")";
     	  if( m +1 != M || j+1 != J) cout <<",";
     	  cout << endl;
     	}
   	}
    cout <<")" << endl;
    cout << "candela('GanttChart'," << endl;
    cout << "    data=data, label='name', " << endl;
    cout << "    start='start', end='end', level='level', " << endl;
    cout << "    width=" <<J*100 <<", height="<< M * 100<<")" << endl;
	}

	void calcPenalties(const vector<unsigned> & starts, double&  ePenalty, double& lPenalty, vector<double>& operPenalties){
		ePenalty = 0;
		lPenalty = 0;
		int curDueDate;
		int curStart;
		double curTardiness;
		double curEarliness;
		//this->iToR(starts);

		operPenalties.resize(O, 0);

		for(unsigned o=1; o<O; ++o){

			curStart = starts[o];
			curDueDate = deadlines[o];
			
			curEarliness = max(curDueDate-(curStart+(int)P[o]), 0) * earlPenalties[o];
			curTardiness = max((curStart+(int)P[o])-curDueDate, 0) * tardPenalties[o];

			assert(curEarliness >= 0);
			assert(curTardiness >= 0);
			ePenalty += curEarliness;
			lPenalty += curTardiness;
			operPenalties[o] = curEarliness > 0 ? curEarliness : curTardiness;
		}
	}

	void printPenaltys(const vector<unsigned> & starts, const unsigned & makes){
		
		double sumTardPenaltys = 0;
		double sumEarlPenaltys = 0;
		double sumPenaltys = 0;
		int curDueDate;
		int curStart;
		double curTardiness;
		double curEarliness;

#ifdef RFILE
		iToR2(starts);
		return;
#endif //RFILE

		for(unsigned o=1; o<O; ++o){

			curStart = starts[o];
			curDueDate = deadlines[o];
			
			curEarliness = max(curDueDate-(curStart+(int)P[o]), 0);
			curTardiness = max((curStart+(int)P[o])-curDueDate, 0);

			assert(curEarliness >= 0);
			assert(curTardiness >= 0);

			sumEarlPenaltys += curEarliness*earlPenalties[o];
			sumTardPenaltys += curTardiness*tardPenalties[o];
		}

		sumPenaltys = sumEarlPenaltys + sumTardPenaltys;

#ifdef PRINT_SCHEDULE

		for (unsigned j = 0; j < J; ++j) {

			printf("%02d - ", j);

			for (unsigned m = 0; m < M; ++m) {
				int i = jmToIndex[j][m];
				curStart = starts[i];
				curDueDate = deadlines[i];
				curEarliness = max(curDueDate - (curStart + (int)P[i]), 0);
				curTardiness = max((curStart + (int)P[i]) - curDueDate, 0);
				printf("[%03u|%03u|%03u|%2.0f|%2.0f] ", starts[i], P[i], deadlines[i], curEarliness, curTardiness);
			}
			cout << endl;
		}

		return;
#endif //PRINT_SCHEDULE

#ifndef PRINT_ONLY_RESULT
		cout << sumPenaltys << " " << sumEarlPenaltys << " " << sumTardPenaltys << " ";
#endif //PRINT_ONLY_RESULT
#ifdef PRINT_ONLY_RESULT
		cout << sumPenaltys;
#endif // PRINT_ONLY_RESULT


#ifdef NEIGHBOURS_NB
		double mean = 0;
		for (double n : neigh) {
			mean += n;
		}
		if (neigh.size()) {
			mean = mean / neigh.size();
		}
		cout << mean << " ";
#endif // NEIGHBOURS_NB

		cout << endl;
	}


	//o1 preceedes o2?
	bool hasPrec(unsigned o1, unsigned o2) const {
		assert(o1 < O   &&   o1!=0);
		assert(o2 < O   &&   o2!=0);
		assert(o1!=0   &&   o2!=0);

		if(operToJ[o1] != operToJ[o2])
			return false;
		return posets[operToJ[o1]].cg.path(operToM[o1], operToM[o2]);
	}

	void parse(const string & filePath) {
		P.clear();
		operToJ.clear(); 
		operToM.clear();
		jSize.clear();
		mSize.clear();
		jobOpers.clear();
		machOpers.clear();
		deadlines.clear();
		earlPenalties.clear();
		tardPenalties.clear();
		
		string bufferStr;
		double bufferT;
		vector<pair<unsigned, unsigned>> order;
		unsigned fstM, sndM;
		unsigned fstO, sndO;
		vector<unsigned> fromDummy;
		vector<unsigned> toDummy;
	
		ifstream streamFromFile;
		streamFromFile.open(filePath);
		if( ! streamFromFile.is_open())
			throw errorText("Could not open instance file: "+filePath, "Instance.hpp", "Inst::parse()");

		streamFromFile >> J;
		
		streamFromFile >> M;

		P.push_back(0);
		deadlines.push_back(0);
		earlPenalties.push_back(0);
		tardPenalties.push_back(0); 
		operToJ.push_back(UINT_MAX);
		operToM.push_back(UINT_MAX);

		//initializing
		O = 1;
		jSize.resize(J,0);
		mSize.resize(M,0);

		jmToIndex.resize(extents[J][M]);
		jobOpers.resize(J);
		machOpers.resize(M);

		unsigned auxM;

		for(unsigned j=0; j<J; j++) {
			for(unsigned m=0; m<M; m++) {
				streamFromFile >> auxM;
				streamFromFile >> bufferT;

				assert(bufferT != 0);

				if(bufferT != 0) {
					P.push_back(bufferT);
					assert(P[O] == bufferT);
					jmToIndex[j][auxM] = O;
					jobOpers[j].push_back(O);
					machOpers[auxM].push_back(O);

					streamFromFile >> bufferT;
					deadlines.push_back(bufferT);
					assert(deadlines[O] == bufferT);
					streamFromFile >> bufferT;
					earlPenaltys.push_back(bufferT);
					assert(earlPenaltys[O] == bufferT);
					streamFromFile >> bufferT;
					tardPenaltys.push_back(bufferT);
					assert(tardPenaltys[O] == bufferT);

					O++;
					operToJ.push_back(j);
					operToM.push_back(auxM);
					jSize[j]++;
					mSize[auxM]++;
				} else {
					jmToIndex[j][auxM] = 0;
				}

			}
		}
		assert(O <= J*M+1);
		assert(operToJ.size() == O);
		assert(operToM.size() == O);

		streamFromFile.close();

		assert(!streamFromFile.is_open());

		streamFromFile.open(filePath);

		assert(streamFromFile.is_open());

		//Meta that needs O and P
		next.resize(O);
		prev.resize(O);

		streamFromFile >> bufferT;
		streamFromFile >> bufferT;

		//Reading precedences
		for(unsigned j=0; j<J; j++) {

			order.clear();

			fromDummy.clear();
			toDummy.clear();
			
			streamFromFile >> fstM;
			streamFromFile >> bufferT;
			streamFromFile >> bufferT;
			streamFromFile >> bufferT;
			streamFromFile >> bufferT;
			for(unsigned u=1; u<M; u++) {
				assert(fstM < M);
				fstO = jmToIndex[j][fstM];
				assert(fstO < O);
				
				streamFromFile >> sndM;
				assert(sndM < M);
				sndO = jmToIndex[j][sndM];
				assert(sndO < O);

				if(fstO!=0   &&   sndO!=0) {
					next[fstO].push_back(sndO);
					prev[sndO].push_back(fstO);
				}

				order.push_back(pair<unsigned, unsigned>(fstM, sndM));
				fstM = sndM;
				streamFromFile >> bufferT;
				streamFromFile >> bufferT;
				streamFromFile >> bufferT;
				streamFromFile >> bufferT;
			}

			posets.push_back(Poset(order, M));
		}

		assert(M <= MAX_MACHS);
		for(unsigned j=0; j<J; j++) {
			assert(! posets[j].roots.empty());
			for(unsigned m : posets[j].roots) {
				if(jmToIndex[j][m] != 0)
					roots.push_back(jmToIndex[j][m]);
			}
		}
		for(unsigned j=0; j<J; j++) {
			assert(! posets[j].leafs.empty());
			for(unsigned m : posets[j].leafs) {
				if(jmToIndex[j][m] != 0)
					leafs.push_back(jmToIndex[j][m]);
			}
		}

		streamFromFile.close();
	}

	unsigned J;
	unsigned M;
	//includes dummy
	unsigned O; //<sum(jSize)+sum(mSize)+1> if no reentrant and no 0 time -- actual opers start  1. 0 is dummy
	vector<unsigned> P; //<O> processing times

	vector<unsigned> jSize; //<J> --  number of opers of job j is jSize[j]
	vector<unsigned> mSize; //<M> 

	vector<vector<unsigned>> jobOpers; //<J, ?> jobs[x] is vector of opers in job x
	vector<vector<unsigned>> machOpers; //<M, ?>

	vector<unsigned> operToJ; // <O> -- oper o is in operToJ[o] job
	vector<unsigned> operToM; // <O> 

	
	vector<unsigned> roots;//<?>
	vector<unsigned> leafs;

	//XXX TODO MAYBE CHANGE FOR EXTENDED_JSS to vec instead of matrix
	vector<vector<unsigned>> next; //<O,?>
	vector<vector<unsigned>> prev;

	//JIT-JSS informations
	vector<unsigned> deadlines; //<O> deadlines[o] is due-date of an operation o
	vector<double> earlPenalties; //<O> earlPenalties[o] is earliness penalty of operation o
	vector<double> tardPenalties;//<O>	tardPenalties[o] is tardiness penalty of operation o


	//precedences

#ifdef VANILLA_PSS
	vector<Poset> posets; // italianos of job orders careful useing directly - use hasPrec - O(1) but expensive
#endif //VANILLA_PSS


#ifdef VANILLA_PSS
	multi_array<unsigned,2> jmToIndex; //jmToIndex[j][m] = o -- may be dummy if doesnt exist
#endif //VANILLA_PSS
	

};
Inst inst;
