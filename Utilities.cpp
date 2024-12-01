#include <math.h>

#include <sys/stat.h>

#include "Utilities.hpp"

bool lexiOrder(const string& s1, const string& s2) {

	char cs1[100];
	char cs2[100];

	strcpy(cs1, s1.c_str());
	strcpy(cs2, s2.c_str());

	return lexicographical_compare(cs1, cs1 + s1.size(), cs2, cs2 + s2.size());
}

bool doesFileExist(const string& filename) {
	struct stat buf;
	if (stat(filename.c_str(), &buf) != -1) {
		return true;
	}
	return false;
}


bool metropolis(unsigned currMakes, unsigned childMakes, double T) {
	assert(T > 0.0);
	if (childMakes <= currMakes)
		return true;

	if (drand48() <= exp(((double)currMakes - (double)childMakes) / T)) {
		//cerr << "I accepted " << childMakes << " compared to " << currMakes << " with temp. " << endl; 
		return true;
	}
	else {
		return false;
	}
}