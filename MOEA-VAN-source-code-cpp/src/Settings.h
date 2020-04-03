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

	double CR; // DE crossover probability
	double F; //  scaling factor 
    double pm; // mutation probability
    double etam; // mutation distribution index
	double delta; // probability for mating pool
	int T; // neighborhood size
    int p; // parameter to compute Lp norm
	int K; // number of closest to check dominance
	int mu; // population size
	int maxEval; // stopping criterion
};

#endif