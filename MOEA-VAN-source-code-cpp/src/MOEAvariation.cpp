#include "MOEA.h"


void CMOEA::variation()
{
	DEoperator();

	polynomialMutation();

	problem.evaluate(population[settings.mu]);
}


// apply polynomial mutation on offspring
void CMOEA::DEoperator()
{
	// select parents from mating pool
	int p0 = parents[0];
	int p1 = parents[1];
	int p2 = parents[2];

	// apply DE operator
	int jr = util.randi(problem.n);
	for(int j = 0; j < problem.n; j++)
	{
		// appy DE operator
		if((util.rand() <= settings.CR) || (j == jr)){
			population[settings.mu].x[j] = population[p0].x[j] + settings.F*(population[p1].x[j] - population[p2].x[j]);
		}
		else{
			population[settings.mu].x[j] = population[p0].x[j];
		}

		// bounds
		if(population[settings.mu].x[j] < problem.lb[j]){
			population[settings.mu].x[j] = problem.lb[j] + util.rand()*(population[p0].x[j] - problem.lb[j]);
		}
		if(population[settings.mu].x[j] > problem.ub[j]){
			population[settings.mu].x[j] = problem.ub[j] - util.rand()*(problem.ub[j] - population[p0].x[j]);
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

			population[settings.mu].x[j] += (problem.ub[j] - problem.lb[j])*delta;
		}

		// bounds
		if(population[settings.mu].x[j] < problem.lb[j]){
			population[settings.mu].x[j] = problem.lb[j];
		}
		if(population[settings.mu].x[j] > problem.ub[j]){
			population[settings.mu].x[j] = problem.ub[j];
		}
	}

}