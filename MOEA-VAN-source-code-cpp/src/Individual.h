#ifndef Individual_H
#define Individual_H

class CIndividual
{
public:
	double * x;
	double * f;

	double * nf; // normalized objectives
	double dist2ref; // distance to reference point

public:
	CIndividual(void);
	~CIndividual(void);
};


#endif