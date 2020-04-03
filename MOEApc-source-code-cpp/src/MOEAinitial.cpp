#include "MOEA.h"
#include <iostream>


void CMOEA::initGrid(void)
{
	int i, j;
	int m = problem.m;

	// calculate number of divisions
	int numberOfDivisions = (int)ceil( pow( settings.mu, 1.0/(m-1.0) ) );

	// number of grids
	numberOfGrids = (int)pow(numberOfDivisions, (double)(m-1));

	// calculate angle range
	deltaAngle = 3.1415926535897932384626433832795/(2.0*numberOfDivisions);

	// angles
	double * angles = new double[numberOfDivisions];
	for(i = 0; i < numberOfDivisions; i++){
		angles[i] = deltaAngle*i;		
	}

	// allocate mem for grid
	grid = new double*[numberOfGrids];
	for(i = 0; i < numberOfGrids; i++){
		grid[i] = new double[m-1];		
	}

	// generate grid
	int idx;
	for(j = 0; j < m-1; j++)
	{
		idx = 0;
		for(i = 0; i < numberOfGrids; i++)
		{
			grid[i][j] =  angles[idx];       

			if(fmod((double)(i+1), pow(numberOfDivisions, (double)(j+1))) == 0){
				idx = 0;
			}
			else if(fmod((double)(i+1), pow(numberOfDivisions, (double)j)) == 0){
				idx++;
			}
		}
	}

}


void CMOEA::initNeighborhoods(void)
{
	neighborhood = new int*[numberOfGrids];

	whole = new int[numberOfGrids];

	double * distVector = new double[numberOfGrids];

	int * idx = new int[numberOfGrids];

	int i, j;
	for(i = 0; i < numberOfGrids; i++)
	{
		whole[i] = i;

		neighborhood[i] = new int[settings.T];

		// calculate the distances based on weight vectors
		for(j = 0; j < numberOfGrids; j++)
		{
			distVector[j] = util.distance(grid[i], grid[j], problem.m-1);
			idx[j] = j;
		} 

		// find nearest neighboring subproblems
		util.minFastSort(distVector, idx, numberOfGrids, settings.T);

		for(j = 0; j < settings.T; j++){
			neighborhood[i][j] = idx[j];
		}
	}

	// free memory
	delete [] distVector;
	delete [] idx;

}



void CMOEA::initIndividuals(void)
{

	// allocate memory for population
	population = new CIndividual[numberOfGrids];

	int i, j;
	// generate initial population	
	for(i = 0; i < numberOfGrids; i++)
	{            
		population[i].x = new double[problem.n];

		population[i].f = new double[problem.m];

		population[i].theta = new double[problem.m-1];

		// generate idividual in the decision space
		for(j = 0; j < problem.n; j++){
			population[i].x[j] = util.rand(problem.lb[j], problem.ub[j]);
		}

		// evaluate individual
		problem.evaluate(population[i]);
	}

	// initialize offspring
	offspring.x = new double[problem.n];

	offspring.f = new double[problem.m]; 

	offspring.theta = new double[problem.m-1];

}



void CMOEA::initReferencePoint(void)
{
	int i, j;

	// initialize reference point
	refPoint = new double[problem.m];

	for(j = 0; j < problem.m; j++){
		refPoint[j] = population[0].f[j]; 			
	}

	for(i = 1; i < numberOfGrids; i++)
	{
		// update reference point
		for(j = 0; j < problem.m; j++)
		{
			if(population[i].f[j] < refPoint[j]){
				refPoint[j] = population[i].f[j];
			}
		}
	}
}


void CMOEA::initPolarCoordinates(void)
{
	for(int i = 0; i < numberOfGrids; i++)
	{
		// dist to ref point
		population[i].rho = util.distance(population[i].f, refPoint, problem.m);

		// calc angles
		calculateAngles(population[i].f, problem.m, population[i].rho, population[i].theta);
	}

}



void CMOEA::calculateAngles(double * y, int m, double rho, double * theta)
{
	int i, j;
	double u;
	for(i = 0; i < m-1; i++)
	{
		u = (y[m-1-i] - refPoint[m-1-i])/rho;

		for(j = 0; j < i; j++){
			u /= cos(theta[j]);
		}

		theta[i] = asin(u);
	}

}

