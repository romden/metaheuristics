#include "MOEA.h"
#include <iostream>


// select new population
void CMOEA::environmentalSelection(void)
{

	if(!dominateClosest(settings.K))
	{	
		findClosest(indexes);

		int toReplace = indexes[util.randi(2)];

		if(population[indexes[0]].dist2ref < population[indexes[1]].dist2ref){
			toReplace = indexes[1];
		}
		else if(population[indexes[1]].dist2ref < population[indexes[0]].dist2ref){
			toReplace = indexes[0];
		}

		if(toReplace < settings.mu){
			replace(toReplace);
		}
	}

} // environmentalSelection


// check whether k closest dominate offspring
bool CMOEA::dominateClosest(int k)
{
	int i;
	for(i = 0; i < settings.mu; i++){	
		indexes[i] = i;
		tempArray[i] = -mxAnlges[settings.mu][i];
	}
	// sort
	util.minFastSort(tempArray, indexes, settings.mu, k);

	for(i = 0; i < k; i++)
	{
		if(dominates(population[indexes[i]].nf, population[settings.mu].nf, problem.m)){
			return true;
		}
	}

	return false;
}


// dominance relaton
bool CMOEA::dominates(double * a, double * b, int m)
{
	// check whether a dominates b
	for(int j = 0; j < problem.m; j++)
	{
		if(b[j] < a[j]){
			return false;
		}
	}

	return true;
}


// find 2 closest vectors
void CMOEA::findClosest(int * idx)
{	
	double maxValue = -std::numeric_limits<double>::infinity();

	int i, j;
	for(i = 0; i < settings.mu+1; i++)
	{
		for(j = 0; j < settings.mu; j++)
		{
			if(j == i){
				continue;
			}
			if(mxAnlges[i][j] > maxValue)
			{
				maxValue = mxAnlges[i][j];
				idx[0] = i;
				idx[1] = j;
			}
		}
	}

}


// replace the worst individual
void CMOEA::replace(int idx)
{
	int i;

	// vars
	for(i = 0; i < problem.n; i++)
	{
		tempArray[i] = population[idx].x[i];
		population[idx].x[i] = population[settings.mu].x[i];
		population[settings.mu].x[i] = tempArray[i];
	}

	// obj and norm obj
	for(i = 0; i < problem.m; i++)
	{
		tempArray[i] = population[idx].f[i];
		population[idx].f[i] = population[settings.mu].f[i];
		population[settings.mu].f[i] = tempArray[i];

		tempArray[i] = population[idx].nf[i];
		population[idx].nf[i] = population[settings.mu].nf[i];
		population[settings.mu].nf[i] = tempArray[i];
	}

	// dists to ref
	tempArray[0] = population[idx].dist2ref;
	population[idx].dist2ref = population[settings.mu].dist2ref;
	population[settings.mu].dist2ref = tempArray[0];


	// update matrix
	for(i = 0; i < settings.mu; i++){
		tempArray[i] = mxAnlges[settings.mu][i];
		mxAnlges[settings.mu][i] = mxAnlges[idx][i];
		mxAnlges[idx][i] = tempArray[i];
	}
	mxAnlges[settings.mu][idx] = mxAnlges[idx][idx];

	for(i = 0; i < settings.mu; i++){
		mxAnlges[i][idx] = tempArray[i];
	}
	mxAnlges[idx][idx] = util.cosAngle(population[idx].nf, population[idx].nf, problem.m);
}