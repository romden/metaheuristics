#include "MOEA.h"


void CMOEA::matingSelection(int current)
{
	// mating selection based on probability
	if (util.rand() <= settings.delta)
	{
		matingPool = neighborhood[current]; // neighborhood
		poolSize = settings.T;
	}
	else
	{
		matingPool = whole; // whole population
		poolSize = numberOfGrids;
	}

 }