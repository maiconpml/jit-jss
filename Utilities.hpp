/*
@author: Tadeu Knewitz Zubaran tkzubaran@gmail.com
*/

#pragma once

#include <sstream>
#include <cstring>

#include "Settings.hpp"

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

bool lexiOrder(const string& s1, const string& s2);

bool doesFileExist(const string& filename);


bool metropolis(unsigned currMakes, unsigned childMakes, double T);







