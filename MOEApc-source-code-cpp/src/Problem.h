#ifndef Problem_H
#define Problem_H

#include <fstream>
#include <string>

#include "Individual.h"

#include "problems\LZ09\LZ09.h"
#include "problems\LGZ\LGZ.h"



class CProblem
{
public:
	CProblem(void);
	~CProblem(void);

	std::string name;

	int m;
	int n;
	double * lb;
	double * ub;

	void (*feval)(const double *, int, double *, int);

	void evaluate(CIndividual);

	void initialization(void);
};

#endif