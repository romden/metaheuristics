#include "MOEA.h"
#include <iostream>


// function to perform environmental selection
void CMOEA::environmentalSelection(void)
{
	int i, j;
	update = true;

	while(update)
	{
		update = false;

		for(i = 0; i < numberOfGrids; i++)
		{
			if(belongs2grid(offspring.theta, i))
			{

				if(!belongs2grid(population[i].theta, i) && !dominates(population[i].f, offspring.f, problem.m))
				{
					replace(i);
					update = true;
				}
				else if(belongs2grid(population[i].theta, i) && offspring.rho < population[i].rho)
				{
					replace(i);
				}

				break;
			}
		}
	}

}


bool CMOEA::belongs2grid(double * theta, int gridNumber)
{
	bool belongs = true;
	double lbA, ubA;

	for(int i = 0; i < problem.m-1; i++)
	{
		lbA = grid[gridNumber][i];
		ubA = grid[gridNumber][i] + deltaAngle;

		if(!(lbA <= theta[i] && theta[i] <= ubA))
		{
			belongs = false;
			break;
		}
	}

	return belongs;
}




// function to check dominance
bool CMOEA::dominates(double * a, double * b, int m)
{
	bool dom = true;

	// check whether a dominates b
	for(int j = 0; j < problem.m; j++)
	{
		if(b[j] < a[j])
		{
			dom = false;
			break;
		}
	}

	return dom;
}


// function to replace idx-th individual by offspring
void CMOEA::replace(int idx)
{
	int i;
	double temp;

	// vars
	for(i = 0; i < problem.n; i++)
	{
		temp = population[idx].x[i];
		population[idx].x[i] = offspring.x[i];
		offspring.x[i] = temp;
	}

	// objs
	for(i = 0; i < problem.m; i++)
	{
		// original 
		temp = population[idx].f[i];
		population[idx].f[i] = offspring.f[i];
		offspring.f[i] = temp;
	}

	// theta
	for(i = 0; i < problem.m-1; i++)
	{
		temp = population[idx].theta[i];
		population[idx].theta[i] = offspring.theta[i];
		offspring.theta[i] = temp;
	}

	// distance to ref point
	temp = population[idx].rho;
	population[idx].rho = offspring.rho;
	offspring.rho = temp;

}