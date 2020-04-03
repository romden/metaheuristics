#include "MOEA.h"


void CMOEA::finalization(void)
{
	// save results
	saveOutcome();

	// release memory
	memoryRelease();
}


void CMOEA::saveOutcome(void)
{
	// open files to output
	std::ofstream fout1, fout2;
	fout1.open("VAR.dat");
	fout2.open("FUN.dat");

	// declare
	int i, j;

	// write to file
	for(i = 0; i < settings.mu; i++)
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


void CMOEA::memoryRelease(void)
{
	for(int i = 0; i < settings.mu+1; i++)
	{
		delete [] population[i].x;
		delete [] population[i].f;
		delete [] population[i].nf;
		delete [] mxAnlges[i];
	}

	delete [] population;
	delete [] mxAnlges;

	delete [] tempArray;
	delete [] refPoint;
	delete [] indexes;	
	delete [] parents;
}