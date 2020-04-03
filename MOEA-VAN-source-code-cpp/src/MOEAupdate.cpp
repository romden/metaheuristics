#include "MOEA.h"



void CMOEA::update(void)
{
	int head = settings.mu;
	int tail = settings.mu+1;

	if(updateReferencePoint()){
		head = 0;
	}

	// translate objectives and compute norm
	objectiveTranslation(head, tail);

	// compute angles between objective vectors
	angleComputation(head, tail);
}


bool CMOEA::updateReferencePoint(void)
{
	bool updating = false; // ref point is not updated

	for(int j = 0; j < problem.m; j++)
	{
		if(population[settings.mu].f[j] < refPoint[j])
		{
			refPoint[j] = population[settings.mu].f[j];

			updating = true;// ref point is updated
		}
	}

	return updating;
}


void CMOEA::objectiveTranslation(int head, int tail)
{
	int i, j;

	// project to the hyperplane 
	for(i = head; i < tail; i++)
	{
		// normalize objective vector
		for(j = 0; j < problem.m; j++){		
			population[i].nf[j] = population[i].f[j] - refPoint[j];
		}

		// calculate distance to ref point
		population[i].dist2ref = util.Lpnorm(population[i].nf, problem.m, settings.p);
	}
}


void CMOEA::angleComputation(int head, int tail)
{
	int i, j;

	for(i = head; i < tail; i++)
	{
		// dist between individuals
		for(j = 0; j < settings.mu; j++){
			mxAnlges[i][j] = util.cosAngle(population[i].nf, population[j].nf, problem.m);
		}
	}
}
