#include "MOEA.h"


void CMOEA::variation(int current)
{
	DEoperator(current);

	polynomialMutation();

	problem.evaluate(offspring);
}


// apply polynomial mutation on offspring
void CMOEA::DEoperator(int current)
{
	// select two parents from mating pool
	int p1 = matingPool[util.randi(poolSize)];
	while(p1 == current){
		p1 = matingPool[util.randi(poolSize)];
	}
	int p2 = matingPool[util.randi(poolSize)];
	while((p2 == p1) || (p2 == current)){
		p2 = matingPool[util.randi(poolSize)];
	}

	// apply DE operator
	int jr = util.randi(problem.n);
	for(int j = 0; j < problem.n; j++)
	{
		// appy DE operator
		if((util.rand() <= settings.CR) || (j == jr)){
			offspring.x[j] = population[current].x[j] + settings.F*(population[p1].x[j] - population[p2].x[j]);
		}
		else{
			offspring.x[j] = population[current].x[j];
		}

		// bounds
		if(offspring.x[j] < problem.lb[j]){
			offspring.x[j] = problem.lb[j] + util.rand()*(population[current].x[j] - problem.lb[j]);
		}
		if(offspring.x[j] > problem.ub[j]){
			offspring.x[j] = problem.ub[j] - util.rand()*(problem.ub[j] - population[current].x[j]);
		}

	} 
}


// apply polynomial mutation on offspring
void CMOEA::polynomialMutation()
{
	double u, delta;

	for(int j = 0; j < problem.n; j++)
	{
		// apply polinomial mutation
		if(util.rand() <= settings.pm)
		{
			u = util.rand();

			if(u <= 0.5){
				delta = pow(2.0*u, 1.0/(settings.etam + 1.0)) - 1.0;
			}
			else{
				delta = 1.0 - pow(2.0*(1.0 - u), 1.0/(settings.etam + 1.0));
			}

			offspring.x[j] += (problem.ub[j] - problem.lb[j])*delta;
		}

		// bounds
		if(offspring.x[j] < problem.lb[j]){
			offspring.x[j] = problem.lb[j];
		}
		if(offspring.x[j] > problem.ub[j]){
			offspring.x[j] = problem.ub[j];
		}
	}

}