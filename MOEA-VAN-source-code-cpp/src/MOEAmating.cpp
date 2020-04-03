#include "MOEA.h"

void CMOEA::matingSelection()
{
	// select individual producing offspring
	parents[0] = util.randi(settings.mu);

	for(int i = 0; i < settings.mu; i++)
	{
		indexes[i] = i; 	
		tempArray[i] = -mxAnlges[parents[0]][i];
	}

	// mating selection based on probability
	if (util.rand() <= settings.delta)
	{
		// neighborhood
		poolSize = settings.T; 
		// find nearest neighboring subproblems
		util.minFastSort(tempArray, indexes, settings.mu, settings.T);
	}
	else
	{
		// whole population
		poolSize = settings.mu;
	}


	// select two parents from mating pool
	parents[1] = indexes[util.randi(poolSize)];
	while(parents[1] == parents[0]){
		parents[1] = indexes[util.randi(poolSize)];
	}
	parents[2] = indexes[util.randi(poolSize)];
	while((parents[2] == parents[1]) || (parents[2] == parents[0])){
		parents[2] = indexes[util.randi(poolSize)];
	}

}