#include "MOEA.h"


void CMOEA::initialization(void)
{
	memoryAllocation();

	initIndividuals();

	objectiveTranslation(0, settings.mu);

	angleComputation(0, settings.mu);
}


// function to allocate memory for variables
void CMOEA::memoryAllocation(void)
{
	// allocate memory for population
	population = new CIndividual[settings.mu+1];

	// allocate mamory for distance mx
	mxAnlges = new double*[settings.mu+1];

	for(int i = 0; i < settings.mu+1; i++)
	{  
		// allocate memory for individuals variables
		population[i].x = new double[problem.n];
		population[i].f = new double[problem.m];
		population[i].nf = new double[problem.m];

		// allocate memory for kernel matrix
		mxAnlges[i] = new double[settings.mu];
	}

	// allocate memory for temporary array
	tempArray = new double[settings.mu];

	// allocate memory for reference point
	refPoint = new double[problem.m];

	// allocate memory for mating pool arrays
	indexes = new int[settings.mu];

	// allocate memory for parents
	parents = new int[3];
}


void CMOEA::initIndividuals(void)
{
	int i, j;

	// init ref point
	for(j = 0; j < problem.m; j++){
		refPoint[j] = std::numeric_limits<double>::infinity();
	}

	// generate initial population	
	for(i = 0; i < settings.mu; i++)
	{
		// generate idividual in the decision space
		for(j = 0; j < problem.n; j++){
			population[i].x[j] = util.rand(problem.lb[j], problem.ub[j]);
		}

		// evaluate individual except offspring 
		problem.evaluate(population[i]);

		// update ref point
		for(j = 0; j < problem.m; j++)
		{
			if(population[i].f[j] < refPoint[j]){
				refPoint[j] = population[i].f[j];
			}
		}
	}
}



