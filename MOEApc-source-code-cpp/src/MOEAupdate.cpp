#include "MOEA.h"



void CMOEA::updateReferencePoint(void)
{
	// ref point is not updated
	update = false;

	for(int j = 0; j < problem.m; j++)
	{
		if(offspring.f[j] < refPoint[j])
		{
			refPoint[j] = offspring.f[j];
			// ref point is updated
			update = true;
		}
	}
}



void CMOEA::updatePolarCoordinates(void)
{
	int i, j;

	if(update)
	{		
		for(i = 0; i < numberOfGrids; i++)
		{
			// calc dist to ref
			population[i].rho = util.distance(population[i].f, refPoint, problem.m);

			// calc angles
			calculateAngles(population[i].f, problem.m, population[i].rho, population[i].theta);
		}		
	}

	// calc dist to ref for offspring
	offspring.rho = util.distance(offspring.f, refPoint, problem.m);

	// calculate polar coordinate for offspring
	calculateAngles(offspring.f, problem.m, offspring.rho, offspring.theta);
}

