#pragma once

#include <string>

class Parameters {
public:

	string instPath;
	string name;
	double maxMillisecs;
	unsigned seed;
	unsigned tenure;
	unsigned initialjumpLimit;
	unsigned decreaseDivisor;
	unsigned bjSize;
	unsigned maxD;
	unsigned maxC;
	bool timeLog;
	double scaleTime;
	string searchTypeStr;
	bool onlyMakesLowerBound;
	unsigned schedulerType;
};