/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#pragma once

#include "Settings.hpp"

#include <stdio.h>
#include <string>

#include <iostream>
#include <sstream>
#include <fstream>

#include <vector>

#include <dirent.h> 
#include <algorithm> 

#include <cstring>

#include <sys/stat.h>

inline string doubleStr(const long double aDouble) {
	ostringstream strs;
	strs << aDouble;
	string str = strs.str();
	return str;
}

inline string intStr(const long long aInt) {
    stringstream mySs;
    mySs << aInt;
    return mySs.str();
}

inline string unsigStr(const unsigned aInt) {
	if(aInt == UINFTY)
		return "I";
    stringstream mySs;
    mySs << aInt;
    return mySs.str();
}

inline const string errorText(const string & message, const string & fileName, const string & methodName) {
	return "\nERROR. "+message+" in file: "+fileName+" ; method: "+methodName;
}

bool hasTermination(const string & name, const string & term) {
	if(name.size() < term.size())
		return false;
		
	if(name.compare(".") == 0   &&   name.compare(".."))
		return false;
			
	string s;
	s.resize(term.size());
		
	for(unsigned u=0; u<s.size(); u++) {
		if(term[term.size()-1-u] != name[name.size()-1-u]) {
			return false;
		}
	}
		
	return true;
}

bool lexiOrder(const string & s1, const string & s2) { 
	
	char cs1[100];
	char cs2[100];
	
	strcpy (cs1, s1.c_str());
	strcpy (cs2, s2.c_str());
	
	return lexicographical_compare(cs1,cs1+s1.size(),cs2,cs2+s2.size());
}
	
vector<string> filesInDir(const string & dirPath, const string & termination) {
	vector<string> instFileNames;
		
	DIR *d;
	struct dirent *dir;
	d = opendir(dirPath.c_str());
	string s;
		
	if( ! d) {
		throw errorText("Can't open inst dir: " + dirPath, "Batch.hpp", "batch()");
	}
		
	dir = readdir(d);
	while (dir != NULL) {
		s = dir->d_name;
		if(hasTermination(s, termination)) {
			instFileNames.push_back(s);
		}
		dir = readdir(d);
	}
	
	closedir(d);
	
	sort (instFileNames.begin(), instFileNames.end(), lexiOrder); 
	return instFileNames;
}

bool doesFileExist(const string & filename) {
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1) {
        return true;
    }
    return false;
}


bool metropolis(unsigned currMakes, unsigned childMakes, double T) {
	assert(T > 0.0);
	if(childMakes<=currMakes)
		return true;

	if (drand48() <= exp(((double)currMakes-(double)childMakes)/T)) { 
		//cerr << "I accepted " << childMakes << " compared to " << currMakes << " with temp. " << endl; 
		return true; 
	} else { 
		return false; 
	}
}







