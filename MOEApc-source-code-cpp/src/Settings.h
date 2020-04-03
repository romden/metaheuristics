#ifndef Settings_H
#define Settings_H

#include <fstream>
#include <string>

class CSettings
{
public:
	CSettings(void);
	~CSettings(void);

	void initialization(void);

	double CR;
	double F;
    double pm;
    double etam;
	double delta;
	int T;
    int mu;
	int maxEval;
};

#endif