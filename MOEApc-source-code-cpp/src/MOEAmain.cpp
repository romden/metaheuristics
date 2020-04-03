#include "MOEA.h"
#include <iostream>


CMOEA::CMOEA(void)
{
}


CMOEA::~CMOEA(void)
{
}


void CMOEA::run(void)
{
	// INITIALIZATION

	// initialize grid
	initGrid();

	// initialize neighborhood
	initNeighborhoods();

	// initialize population
	initIndividuals();

	// initialize reference point
	initReferencePoint();

	// initialize polar coordinates
	initPolarCoordinates();


	// EVOLUTIONARY PROCESS
	int g, idx;
	for(g = 0; g < settings.maxEval; g++)
	{
		// select individual producing offspring
		idx = util.randi(numberOfGrids);

		// select mating pool based on probability
		matingSelection(idx);

		// variation to produce offspring
		variation(idx);

		// update reference point
		updateReferencePoint();

		// update projections
		updatePolarCoordinates();

		// update population
		environmentalSelection();

	} // for


	// FINALIZATION

	// save results
	saveOutcome();

	// release memory
	freeMemory();
}



