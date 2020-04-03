#include "MOEA.h"


void CMOEA::saveOutcome(void)
{
	// open files to output
	std::ofstream fout1, fout2;
	fout1.open("VAR.dat");
	fout2.open("FUN.dat");

	// declare
	int i, j;

	// write to file
	for(i = 0; i < numberOfGrids; i++)
	{
		// save variables
		for(j = 0; j < problem.n; j++){
			fout1  << " " << population[i].x[j];
		}
		fout1 << "\n";

		// save function values
		for(j = 0; j < problem.m; j++){
			fout2 << " " << population[i].f[j];
		}
		fout2 << "\n";		
	}

	// close file
	fout1.close();
	fout2.close();
}


void CMOEA::freeMemory(void)
{
	delete [] refPoint;
	delete [] whole;


	for (int i = 0; i < numberOfGrids; i++)
	{
		delete [] population[i].x;
		delete [] population[i].f;
		delete [] population[i].theta;

		delete [] grid[i];
		delete [] neighborhood[i];
	}


	delete [] population;

	delete [] grid;
	delete [] neighborhood;

	delete [] offspring.x;
	delete [] offspring.f;
	delete [] offspring.theta;
}