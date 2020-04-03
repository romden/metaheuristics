#ifndef UTIL_H
#define UTIL_H

#include <cstdlib>
#include <ctime> 
#include <math.h>


class CUtil
{    
public:
	CUtil::CUtil(void);
	CUtil::~CUtil(void);

	int randi(int);
	double rand(void);
	double rand(double, double);

	void minFastSort(double * x, int * idx, int n, int m);
	void randomPermutation(int * perm, int size);

	double distance(double * vector1, double * vector2, int dim);
	double innerProduct(double * vector1, double * vector2, int dim);
	double cosAngle(double * vector1, double * vector2, int dim);
	double Lpnorm(double * vector, int dim, int p);
};


#endif