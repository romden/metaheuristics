#ifndef MOEA_H
#define MOEA_H

#include <fstream>

#include "Problem.h"
#include "Settings.h"
#include "Individual.h"
#include "Util.h"

class CMOEA
{
public:

	// constructor
	CMOEA(void);

	// destructor
	~CMOEA(void);
		
	int poolSize;

	int * parents;
	int * indexes;
	double * refPoint;
	double * tempArray;
	double ** mxAnlges;
	CIndividual * population;

	CProblem problem;
	CSettings settings;

	CUtil util;
	
	// INITIALIZATION
	void initialization(void);

	void memoryAllocation(void);
	void initIndividuals(void);
	bool updateReferencePoint(void);
	void objectiveTranslation(int head, int tail);
	void angleComputation(int head, int tail);	

	// Evolution
	void run(void);

	void matingSelection(void);
	void variation(void);
	void polynomialMutation(void);
	void DEoperator(void);
	void update(void);
	void environmentalSelection(void);
		
	void findClosest(int *);
	bool dominates(double * a, double * b, int m);	
	void replace(int);	
	bool dominateClosest(int k);

	// finalization
	void finalization(void);	

	void saveOutcome(void);
	void memoryRelease(void);

};


#endif