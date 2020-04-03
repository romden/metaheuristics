#ifndef MOEA_H
#define MOEA_H

#include <fstream>
#include <math.h>

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

	double ** grid; // grid based on polar coordinates
	int numberOfGrids;
	double deltaAngle;

	int poolSize;
	int ** neighborhood;
	int * whole;
	int * matingPool;
	double * refPoint;
	bool update;

	CIndividual * population;
	CIndividual offspring;	

	CProblem problem;
	CSettings settings;

	CUtil util;

	// INITIALIZATION

	void initGrid(void);

	void initNeighborhoods(void);

	void initIndividuals(void);

	void initReferencePoint(void);

	void initPolarCoordinates(void);

	void calculateAngles(double * f, int m, double rho, double * theta);


	// EVOLUTIONARY PROCESS

	void run(void);

	// MATING SELECTION

	void matingSelection(int);

	// VARIATION

	void variation(int);

	void polynomialMutation(void);

	void DEoperator(int);

	// UPDATE

	void updateReferencePoint(void);

	void updatePolarCoordinates(void);
    
	// SELECTION

	void environmentalSelection(void);

	bool dominates(double * a, double * b, int m);

	bool belongs2grid(double *, int);

	void replace(int);

	// FINALIZATION

	void saveOutcome(void);

	void freeMemory(void);
};


#endif
