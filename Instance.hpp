/*
  @author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/
#pragma once

#include <boost/multi_array.hpp>

#include "Settings.hpp"
#include "italiano.hpp"
#include "Utilities.hpp"

using namespace boost;

class Poset {
public:
	Poset(const vector<pair<unsigned, unsigned>>& precs, unsigned _nEl);

	Poset(unsigned _nEl);
	
	~Poset();

	bool operator!=(const Poset& p);
	bool operator==(const Poset& p);
	
	//@return: created cycle doing so?
	bool addAll(const vector<pair<unsigned, unsigned>>& precs);

	//@return: created cycle doing so?
	bool addPrec(unsigned from, unsigned to);

	string toString() const;
	
	double orderStr() const;
	

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

	string toString() const;

	void setHeads(vector<unsigned>& heads, vector<unsigned>& indeg, vector<unsigned>& Q) const;

	void setTails(vector<unsigned>& tails, vector<unsigned>& indeg, vector<unsigned>& Q) const;

	//max for each machine: of min head + min tail + procTs of that machine
	unsigned lowerBoundNasiri(vector<unsigned>& indeg, vector<unsigned>& Q) const;

	//max of lowerBoundnasiri and job lower bound
	unsigned lowerBoundTkz(vector<unsigned>& indeg, vector<unsigned>& Q) const;

	//max of lowerBoundnasiri and job lower bound
	unsigned lowerBoundTkz() const;

	unsigned lowerBoundMss() const;

	bool verifySchedule(const vector<unsigned>& starts, unsigned expecMakes) const;

	void generateGantt(const vector<unsigned>& starts);

	void calcPenalties(const vector<unsigned>& starts, double& ePenalty, double& lPenalty, vector<double>& operPenalties);

	void printOutForTest(const vector<unsigned>& starts);

	void printPenalties(const vector<unsigned>& starts, const unsigned& makes);

	bool hasPrec(unsigned o1, unsigned o2) const;

	void parse(const string& filePath);

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
	vector<Poset> posets; // italianos of job orders careful useing directly - use hasPrec - O(1) but expensive

	multi_array<unsigned,2> jmToIndex; //jmToIndex[j][m] = o -- may be dummy if doesnt exist
};

extern Inst inst;
