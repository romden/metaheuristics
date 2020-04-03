#include "MOEA.h"


CMOEA::CMOEA(void)
{
}


CMOEA::~CMOEA(void)
{
}


void CMOEA::run(void)
{
	// INITIALIZATION
	initialization();

	// EVOLUTIONARY PROCESS
	for(int g = 0; g < settings.maxEval; g++)
	{
		// select mating pool based on probability
		matingSelection();

		// variation to produce offspring
		variation();

		// update ref point and scalar products
		update();

		// select new population
		environmentalSelection();

	} // for


	// FINALIZATION
	finalization();
}
